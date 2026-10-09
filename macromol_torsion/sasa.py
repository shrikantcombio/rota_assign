"""
SasaAnalyzer Module
-------------------
Solvent Accessible Surface Area (SASA) calculation using NACCESS.
Classifies residues into Interface ('I'), Non-Interface ('N'), Surface, and Core.
"""

from pathlib import Path
import os
import shutil
import subprocess
from macromol_torsion.structure_reader import StructureReader


class SasaAnalyzer:
    """
    Interface to NACCESS executable for SASA calculation and interface classification.
    """

    def __init__(self, naccess_binary=None):
        if naccess_binary:
            self.naccess_binary = str(naccess_binary)
        else:
            default_bin = Path("/home/labuser/naccess2.1.1/naccess")
            if default_bin.exists():
                self.naccess_binary = str(default_bin)
            else:
                self.naccess_binary = shutil.which("naccess") or "naccess"

    def run_naccess(self, structure_file, output_dir=None):
        file_p = Path(structure_file)
        if file_p.suffix.lower() == '.rsa':
            return file_p, file_p.with_suffix('.asa')

        pdb_file = StructureReader.write_temp_pdb(file_p)
        work_dir = Path(output_dir) if output_dir else pdb_file.parent

        cwd = os.getcwd()
        os.chdir(work_dir)
        try:
            cmd = f"{self.naccess_binary} {pdb_file.resolve()}"
            subprocess.run(cmd, shell=True, capture_output=True, text=True)
        finally:
            os.chdir(cwd)

        rsa_p = work_dir / f"{pdb_file.stem}.rsa"
        asa_p = work_dir / f"{pdb_file.stem}.asa"
        return rsa_p, asa_p

    def parse_rsa(self, rsa_file_path):
        rsa_p = Path(rsa_file_path)
        if not rsa_p.exists():
            raise FileNotFoundError(f"RSA file not found: {rsa_p}")

        with open(rsa_p, 'r') as f:
            lines = [line for line in f if line.startswith("RES")]
        return lines

    def get_surface(self, unbound_rsa_lines):
        """
        Following Guharoy et al. (2010) and Miller et al. (1987), surface residues
        are defined if the relative accessibility of residue 'X' in NACCESS is > 8.0%
        from its accessibility in the extended Ala-X-Ala tripeptide conformation.
        """
        surface = []
        for residue in unbound_rsa_lines:
            if residue.startswith("RES"):
                try:
                    res_label = residue[4:14].strip()
                    chain = residue[8]
                    res_seq = residue[9:14].strip()
                    rel_sasa = float(residue[22:28].strip())
                    if rel_sasa > 8.0:
                        surface.append(res_label)
                        surface.append((chain, res_seq))
                except (ValueError, IndexError):
                    pass
        return surface

    def get_core(self, unbound_rsa_lines):
        core = []
        for residue in unbound_rsa_lines:
            if residue.startswith("RES"):
                try:
                    rel_sasa = float(residue[22:28].strip())
                    if rel_sasa <= 8.0:
                        res_label = residue[4:14].strip()
                        chain = residue[8]
                        res_seq = residue[9:14].strip()
                        core.append(res_label)
                        core.append((chain, res_seq))
                except (ValueError, IndexError):
                    pass
        return core

    def get_interface(self, bound_rsa_lines, bound_chains, unbound_rsa_lines):
        """
        Following Barik et al. (2012 PRince) and Belapure et al. (2023 IRAA),
        residues are classified as interface residues if their SASA decreases
        by at least 1.0 Angstrom^2 (delta SASA >= 1.0 A^2) upon complex formation.
        """
        if isinstance(bound_chains, (list, tuple, set)):
            chains_set = set(bound_chains)
        else:
            s_chains = str(bound_chains).replace(',', '').replace('+', '').replace(' ', '')
            chains_set = set(s_chains)

        bound_dict = {}
        for line in bound_rsa_lines:
            if line.startswith("RES"):
                res_label = line[4:14].strip()
                chain = line[8]
                res_seq = line[9:14].strip()
                try:
                    abs_sasa = float(line[14:22].strip())
                    bound_dict[(chain, res_seq)] = (res_label, abs_sasa)
                except ValueError:
                    pass

        unbound_dict = {}
        for line in unbound_rsa_lines:
            if line.startswith("RES"):
                chain = line[8]
                res_seq = line[9:14].strip()
                try:
                    abs_sasa = float(line[14:22].strip())
                    unbound_dict[(chain, res_seq)] = abs_sasa
                except ValueError:
                    pass

        interface = []
        for (chain, res_seq), (res_label, b_sasa) in bound_dict.items():
            if (not chains_set or chain in chains_set) and (chain, res_seq) in unbound_dict:
                u_sasa = unbound_dict[(chain, res_seq)]
                delta_sasa = u_sasa - b_sasa
                if delta_sasa >= 1.0:
                    interface.append(res_label)
                    interface.append((chain, res_seq))

        return interface

    def get_non_interface_surface(self, bound_rsa_lines, bound_chains, unbound_rsa_lines):
        """
        Surface residues, other than the interface residues, are considered as non-interface residues.
        """
        surface = self.get_surface(unbound_rsa_lines)
        interface = set(self.get_interface(bound_rsa_lines, bound_chains, unbound_rsa_lines))

        return [res for res in surface if res not in interface]

    def parse_int_file(self, int_file_path):
        """
        Parses an .int interface file (PRince/NACCESS interface format) and returns
        a list of res_label strings and (chain, res_seq) tuples for interface residues.
        """
        int_p = Path(int_file_path)
        if not int_p.exists():
            return []

        interface = []
        with open(int_p, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')) and len(line) >= 26:
                    resname = line[17:20].strip()
                    chain = line[21].strip()
                    res_seq = line[22:26].strip()
                    if chain and res_seq:
                        res_label = f"{resname:<3s} {chain:1s} {res_seq:>4s}".strip()
                        interface.append(res_label)
                        interface.append((chain, res_seq))
                        if res_seq.isdigit():
                            interface.append((chain, str(int(res_seq))))
                            interface.append((chain, int(res_seq)))
        return list(set(interface))

    def classify_residues(self, bound_structure, unbound_structure, bound_chains, output_dir=None, int_file=None, int_type='protein'):
        """
        Classifies residues into interface, non-interface, surface, and core.
        Searches for interface files in 263 subdirectories ({PDB}_p{protein}_r{rna}).
        int_type: 'protein' ({PDB}P.int), 'rna' ({PDB}R.int), or 'complex' ({PDB}.int).
        """
        b_p = Path(bound_structure)
        pdb_tag = b_p.stem.split('_')[0].upper()

        int_candidates = []
        if int_file:
            int_candidates.append(Path(int_file))

        # 1. Search in 263 subdirectories matching {PDB}_p* (e.g. 1C0A_pA_rB)
        parent_dir = b_p.parent
        subdirs = [d for d in parent_dir.glob(f"{pdb_tag}_p*") if d.is_dir()]
        for sd in subdirs:
            if int_type == 'protein':
                int_candidates.extend([sd / f"{pdb_tag}P.int", sd / f"{pdb_tag}.int"])
            elif int_type == 'rna':
                int_candidates.extend([sd / f"{pdb_tag}R.int", sd / f"{pdb_tag}.int"])
            else:
                int_candidates.extend([sd / f"{pdb_tag}.int", sd / f"{pdb_tag}P.int", sd / f"{pdb_tag}R.int"])

        # 2. Check structure directory directly
        if int_type == 'protein':
            int_candidates.extend([
                parent_dir / f"{pdb_tag}P.int",
                parent_dir / f"{pdb_tag}.int",
                Path("/home/labuser/Projects/PhD_projects/rota_assign/dataset") / f"{pdb_tag}P.int",
                Path("/home/labuser/Projects/PhD_projects/rota_assign/dataset") / f"{pdb_tag}.int",
            ])
        elif int_type == 'rna':
            int_candidates.extend([
                parent_dir / f"{pdb_tag}R.int",
                parent_dir / f"{pdb_tag}.int",
                Path("/home/labuser/Projects/PhD_projects/rota_assign/dataset") / f"{pdb_tag}R.int",
                Path("/home/labuser/Projects/PhD_projects/rota_assign/dataset") / f"{pdb_tag}.int",
            ])
        else:
            int_candidates.extend([
                parent_dir / f"{pdb_tag}.int",
                parent_dir / f"{pdb_tag}P.int",
                parent_dir / f"{pdb_tag}R.int",
            ])

        for cand in int_candidates:
            if cand.exists():
                interface = self.parse_int_file(cand)
                return interface, [], [], []

        # 3. Fallback to running NACCESS SASA calculation
        try:
            b_rsa, _ = self.run_naccess(bound_structure, output_dir=output_dir)
            u_rsa, _ = self.run_naccess(unbound_structure, output_dir=output_dir)

            b_lines = self.parse_rsa(b_rsa)
            u_lines = self.parse_rsa(u_rsa)

            interface = self.get_interface(b_lines, bound_chains, u_lines)
            non_int_surf = self.get_non_interface_surface(b_lines, bound_chains, u_lines)
            surface = self.get_surface(u_lines)
            core = self.get_core(u_lines)

            return interface, non_int_surf, surface, core
        except Exception:
            return [], [], [], []

