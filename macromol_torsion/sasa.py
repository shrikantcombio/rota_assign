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
        rel_sasa = lambda a, r: 100.0 * float(a - r) / a if float(a) > 0 else 0.0
        surface = []
        for residue in unbound_rsa_lines:
            if residue.startswith("RES"):
                try:
                    abs_sasa = float(residue[15:22].strip())
                    total_sasa = float(residue[23:28].strip())
                    if rel_sasa(abs_sasa, total_sasa) > 8.0:
                        surface.append(residue[4:14].strip())
                except (ValueError, ZeroDivisionError):
                    pass
        return surface

    def get_core(self, unbound_rsa_lines):
        core = []
        for residue in unbound_rsa_lines:
            if residue.startswith("RES"):
                try:
                    if float(residue[15:22].strip()) == 0.0:
                        core.append(residue[4:14].strip())
                except ValueError:
                    pass
        return core

    def get_interface(self, bound_rsa_lines, bound_chains, unbound_rsa_lines):
        """
        Calculates interface residues by matching (chain, resSeq) between bound complex and unbound structure.
        """
        chains_set = set(bound_chains) if isinstance(bound_chains, (list, tuple, set)) else set(str(bound_chains))

        bound_dict = {}
        for line in bound_rsa_lines:
            if line.startswith("RES"):
                res_label = line[4:14].strip()
                chain = line[8]
                res_seq = line[9:14].strip()
                try:
                    abs_sasa = float(line[15:22].strip())
                    bound_dict[(chain, res_seq)] = (res_label, abs_sasa)
                except ValueError:
                    pass

        unbound_dict = {}
        for line in unbound_rsa_lines:
            if line.startswith("RES"):
                chain = line[8]
                res_seq = line[9:14].strip()
                try:
                    abs_sasa = float(line[15:22].strip())
                    unbound_dict[(chain, res_seq)] = abs_sasa
                except ValueError:
                    pass

        interface = []
        for (chain, res_seq), (res_label, b_sasa) in bound_dict.items():
            if (not chains_set or chain in chains_set) and (chain, res_seq) in unbound_dict:
                u_sasa = unbound_dict[(chain, res_seq)]
                if abs(b_sasa - u_sasa) > 0.1:
                    interface.append(res_label)
                    interface.append((chain, res_seq))

        return interface

    def get_non_interface_surface(self, bound_rsa_lines, bound_chains, unbound_rsa_lines):
        surface = self.get_surface(unbound_rsa_lines)
        interface = self.get_interface(bound_rsa_lines, bound_chains, unbound_rsa_lines)

        interface_set = set(interface)
        return [res for res in surface if res not in interface_set]

    def classify_residues(self, bound_structure, unbound_structure, bound_chains, output_dir=None):
        b_rsa, _ = self.run_naccess(bound_structure, output_dir=output_dir)
        u_rsa, _ = self.run_naccess(unbound_structure, output_dir=output_dir)

        b_lines = self.parse_rsa(b_rsa)
        u_lines = self.parse_rsa(u_rsa)

        interface = self.get_interface(b_lines, bound_chains, u_lines)
        non_int_surf = self.get_non_interface_surface(b_lines, bound_chains, u_lines)
        surface = self.get_surface(u_lines)
        core = self.get_core(u_lines)

        return interface, non_int_surf, surface, core
