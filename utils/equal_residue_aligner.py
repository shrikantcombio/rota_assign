"""
equal_residue_aligner.py
-------------------------
Class EqualResidueAligner for aligning bound (complex) and unbound structure pairs using ClustalW2.
Extracts equal residue positions for protein chains and preserves full-length RNA chains for downstream rotamer and torsion angle calculations.

Features:
- Extracts fasta sequence from PDB and CIF files
- Runs ClustalW2 sequence alignment
- Maps identical (*) and similar (:, .) residue positions
- Filters common atom records for matching protein residue pairs
- Retains full-length RNA chains in bound equal structure output files
- Writes curated equal-residue dataset files with structured naming:
  bound_unbound_C.pdb (complex) and bound_unbound_U.pdb (unbound) under dataset/bound_unbound/
"""

from pathlib import Path
import sys
import os
import shutil
import subprocess
import re
from collections import OrderedDict

# Import StructureReader for PDB & CIF parsing
ROTA_ASSIGN_DIR = Path(__file__).resolve().parent.parent
if str(ROTA_ASSIGN_DIR) not in sys.path:
    sys.path.insert(0, str(ROTA_ASSIGN_DIR))

try:
    from macromol_torsion.structure_reader import StructureReader
except ImportError:
    StructureReader = None

AA_3TO1 = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
    'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
    'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
    'MSE': 'M', ' A': 'A', ' C': 'C', ' G': 'G', ' U': 'U',
    'ADE': 'A', 'CYT': 'C', 'GUA': 'G', 'URA': 'U'
}

PROTEIN_BACKBONE = ['N', 'CA', 'C', 'O']
RNA_BACKBONE = ["O3'", "P", "O5'", "C5'", "C4'", "C3'"]
RNA_RESNAMES = {'A', 'G', 'C', 'U', 'RA', 'RG', 'RC', 'RU', 'ADE', 'GUA', 'CYT', 'URA'}


class EqualResidueAligner:
    """
    Aligns bound and unbound structures using ClustalW2 and extracts matching equal-residue datasets.
    """

    def __init__(self, clustalw2_binary=None, default_output_dir="dataset"):
        if clustalw2_binary:
            self.clustalw2_binary = str(clustalw2_binary)
        else:
            self.clustalw2_binary = (os.environ.get("CLUSTALW2_BIN") or shutil.which("clustalw2")
                                     or shutil.which("clustalw") or "clustalw2")

        self.default_output_dir = Path(default_output_dir)

    def extract_sequence(self, structure_file, chains):
        """
        Extracts 1-letter sequence dictionary per chain from PDB or CIF structure file.
        :returns: (stem_name, chain_seq_dict, chain_list)
        """
        file_p = Path(structure_file)
        stem_name = file_p.stem.replace('_clean', '')

        if StructureReader is not None:
            atoms = StructureReader.parse_structure(file_p)
        else:
            atoms = []

        chains_set = set(chains) if isinstance(chains, (list, tuple, set)) else set(str(chains))
        chain_seq_dict = OrderedDict()
        chain_list = []

        seen_res = set()

        for a in atoms:
            if a.get('record') in ['ATOM', 'HETATM']:
                chain = a['chain']
                if chain in chains_set:
                    aname = a['atom_name']
                    if aname in ['CA', "C4'", 'P']:
                        res_key = (chain, a['resSeq'], a['icode'])
                        if res_key not in seen_res:
                            seen_res.add(res_key)
                            resname_3 = a['resname'].strip().upper()
                            aa_1 = AA_3TO1.get(resname_3, 'X')
                            if chain not in chain_seq_dict:
                                chain_seq_dict[chain] = []
                                chain_list.append(chain)
                            chain_seq_dict[chain].append(aa_1)

        for ch in chain_seq_dict:
            chain_seq_dict[ch] = "".join(chain_seq_dict[ch])

        return stem_name, chain_seq_dict, chain_list

    def run_clustalw2(self, fasta_file):
        """
        Executes ClustalW2 multiple sequence alignment on input fasta file.
        Uses -INFILE= and -OUTFILE= parameters required by ClustalW2.
        :returns: Output .aln alignment file path.
        """
        fa_p = Path(fasta_file)
        aln_p = fa_p.with_suffix('.aln')

        cmd = f"{self.clustalw2_binary} -INFILE={fa_p.resolve()} -OUTFILE={aln_p.resolve()}"
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if not aln_p.exists() or aln_p.stat().st_size == 0 or res.returncode != 0:
            sys.stderr.write(f"Warning: ClustalW2 binary '{self.clustalw2_binary}' returned code {res.returncode}.\n")

        return aln_p

    def parse_alignment_file(self, aln_file, bound_tag, unbound_tag):
        """
        Parses ClustalW2 .aln file to identify identical (*) and similar (:, .) residue positions.
        """
        aln_p = Path(aln_file)
        if not aln_p.exists():
            raise FileNotFoundError(f"Alignment file not found: {aln_p}")

        status = []
        b_chars = []
        ub_chars = []

        with open(aln_p, 'r') as f:
            for line in f:
                if line.startswith('                '):
                    status.extend(line[16:76])
                elif line.startswith(bound_tag):
                    b_chars.extend(line[16:76])
                elif line.startswith(unbound_tag):
                    ub_chars.extend(line[16:76])

        b_stat_index = []
        ub_stat_index = []

        for i in range(min(len(status), len(b_chars), len(ub_chars))):
            if b_chars[i].isalpha():
                b_stat_index.append(status[i])
            if ub_chars[i].isalpha():
                ub_stat_index.append(status[i])

        b_identical = [idx for idx, s in enumerate(b_stat_index, start=1) if s == '*']
        ub_identical = [idx for idx, s in enumerate(ub_stat_index, start=1) if s == '*']

        b_similar = [idx for idx, s in enumerate(b_stat_index, start=1) if s in [':', '.']]
        ub_similar = [idx for idx, s in enumerate(ub_stat_index, start=1) if s in [':', '.']]

        return b_identical, ub_identical, b_similar, ub_similar

    def align_and_extract_equal_residues(self, bound_file, bound_chains, unbound_file, unbound_chains, rna_chains=None, output_dir=None):
        """
        Extracts sequences, runs ClustalW2, identifies equal residues for protein,
        and retains full-length RNA chains in bound equal structure.
        Saves output files with structured naming:
        dataset/bound_unbound/bound_unbound_C.pdb and dataset/bound_unbound/bound_unbound_U.pdb
        """
        bound_p = Path(bound_file)
        unbound_p = Path(unbound_file)

        b_name, b_dict, b_list = self.extract_sequence(bound_p, bound_chains)
        u_name, u_dict, u_list = self.extract_sequence(unbound_p, unbound_chains)

        pair_tag = f"{b_name}_{u_name}"
        out_dir = Path(output_dir) if output_dir else self.default_output_dir / pair_tag
        out_dir.mkdir(parents=True, exist_ok=True)

        # 2. Write FASTA file
        fasta_file = out_dir / f"{pair_tag}.fa"
        with open(fasta_file, 'w') as fa:
            for c in b_list:
                fa.write(f">{b_name}:{c}\n{b_dict[c]}\n")
            for c in u_list:
                fa.write(f">{u_name}:{c}\n{u_dict[c]}\n")

        # 3. Run ClustalW2
        aln_file = self.run_clustalw2(fasta_file)

        # Fallback to python alignment if ClustalW2 failed or file doesn't exist
        if not aln_file.exists() or aln_file.stat().st_size == 0:
            aln_file = self._fallback_python_alignment(b_dict, u_dict, b_name, u_name, fasta_file)

        # 4. Parse alignment
        try:
            b_ident, ub_ident, b_sim, ub_sim = self.parse_alignment_file(aln_file, b_name, u_name)
        except Exception:
            b_ident, ub_ident, b_sim, ub_sim = [], [], [], []

        # 5. Extract atoms and match residues
        b_atoms = StructureReader.parse_structure(bound_p) if StructureReader else []
        ub_atoms = StructureReader.parse_structure(unbound_p) if StructureReader else []

        # Filter equal protein residue atoms
        b_curated_atoms = self._filter_equal_residue_atoms(b_atoms, bound_chains, b_ident, b_sim, ub_atoms, unbound_chains)
        ub_curated_atoms = self._filter_equal_residue_atoms(ub_atoms, unbound_chains, ub_ident, ub_sim, b_atoms, bound_chains)

        # 6. Preserve full-length RNA chain(s) from bound_file in b_curated_atoms
        rna_atoms = self._extract_full_length_rna_atoms(b_atoms, rna_chains=rna_chains, protein_chains=bound_chains)
        if rna_atoms:
            b_curated_atoms.extend(rna_atoms)

        # 7. Write structured equal-residue output files: bound_unbound_C.pdb and bound_unbound_U.pdb
        ext = bound_p.suffix.lower()
        out_bound = out_dir / f"{pair_tag}_C{ext}"
        out_unbound = out_dir / f"{pair_tag}_U{ext}"

        from clean_structure import StructureCleaner
        sc = StructureCleaner()
        if ext in ['.cif', '.mmcif']:
            sc.write_cif(b_curated_atoms, out_bound, structure_id=f"{pair_tag}_C")
            sc.write_cif(ub_curated_atoms, out_unbound, structure_id=f"{pair_tag}_U")
        else:
            sc.write_pdb(b_curated_atoms, out_bound)
            sc.write_pdb(ub_curated_atoms, out_unbound)

        return {
            'BOUND_EQUAL': out_bound,
            'UNBOUND_EQUAL': out_unbound,
            'FASTA': fasta_file,
            'ALIGNMENT': aln_file
        }

    def _extract_full_length_rna_atoms(self, atoms, rna_chains=None, protein_chains=None):
        """
        Extracts all full-length RNA/nucleic acid atom records from structure.
        """
        rna_set = set(rna_chains) if isinstance(rna_chains, (list, tuple, set)) else (set(str(rna_chains)) if rna_chains else set())
        prot_set = set(protein_chains) if isinstance(protein_chains, (list, tuple, set)) else (set(str(protein_chains)) if protein_chains else set())

        rna_atoms = []
        for a in atoms:
            if a.get('record') in ['ATOM', 'HETATM']:
                chain = a['chain']
                resname = a['resname'].strip().upper()

                if rna_set and chain in rna_set:
                    rna_atoms.append(a)
                elif not rna_set and chain not in prot_set and resname in RNA_RESNAMES:
                    rna_atoms.append(a)

        return rna_atoms

    def _fallback_python_alignment(self, b_dict, u_dict, b_name, u_name, fasta_file):
        aln_file = fasta_file.with_suffix('.aln')
        b_seq = "".join(b_dict.values())
        u_seq = "".join(u_dict.values())

        with open(aln_file, 'w') as f:
            f.write(f"CLUSTAL 2.1 multiple sequence alignment\n\n\n")
            n = max(len(b_seq), len(u_seq))
            for i in range(0, n, 60):
                b_sub = b_seq[i:i+60]
                u_sub = u_seq[i:i+60]

                stat = []
                for b_c, u_c in zip(b_sub, u_sub):
                    if b_c == u_c:
                        stat.append('*')
                    else:
                        stat.append(':')
                stat_str = "".join(stat)

                f.write(f"{b_name:<16s}{b_sub}\n")
                f.write(f"{u_name:<16s}{u_sub}\n")
                f.write(f"{' '*16}{stat_str}\n\n")

        return aln_file

    def _filter_equal_residue_atoms(self, atoms, chains, ident_indices, sim_indices, ref_atoms, ref_chains):
        chains_set = set(chains) if isinstance(chains, (list, tuple, set)) else set(str(chains))

        res_order = []
        res_atoms = OrderedDict()

        for a in atoms:
            if a.get('record') in ['ATOM', 'HETATM'] and a['chain'] in chains_set:
                res_key = (a['resname'], a['chain'], a['resSeq'], a['icode'])
                if res_key not in res_atoms:
                    res_order.append(res_key)
                    res_atoms[res_key] = OrderedDict()
                res_atoms[res_key][a['atom_name']] = a

        ident_set = set(ident_indices)
        sim_set = set(sim_indices)

        backbone_set = set(PROTEIN_BACKBONE + RNA_BACKBONE)

        curated = []
        for idx, res_key in enumerate(res_order, start=1):
            if idx in ident_set:
                for a_dict in res_atoms[res_key].values():
                    curated.append(a_dict)
            elif idx in sim_set:
                for aname, a_dict in res_atoms[res_key].items():
                    if aname in backbone_set:
                        curated.append(a_dict)

        return curated


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract and align equal residue datasets between bound and unbound structure pairs using ClustalW2.")
    parser.add_argument("--bound", required=True, help="Bound structure file (.pdb or .cif)")
    parser.add_argument("--bound-chains", required=True, help="Chain ID(s) for bound structure")
    parser.add_argument("--unbound", required=True, help="Unbound structure file (.pdb or .cif)")
    parser.add_argument("--unbound-chains", required=True, help="Chain ID(s) for unbound structure")
    parser.add_argument("--rna-chains", default=None, help="Chain ID(s) for RNA in bound structure")
    parser.add_argument("-o", "--output-dir", default="dataset", help="Output directory for equal-residue dataset files")
    parser.add_argument("--clustalw2", default=None, help="Path to clustalw2 executable (default: $CLUSTALW2_BIN or PATH)")

    args = parser.parse_args()

    aligner = EqualResidueAligner(clustalw2_binary=args.clustalw2, default_output_dir=args.output_dir)
    results = aligner.align_and_extract_equal_residues(
        args.bound, args.bound_chains,
        args.unbound, args.unbound_chains,
        rna_chains=args.rna_chains,
        output_dir=args.output_dir
    )

    print("Equal residue alignment completed:")
    print("Bound Equal Structure:  ", results['BOUND_EQUAL'])
    print("Unbound Equal Structure:", results['UNBOUND_EQUAL'])
    print("FASTA File:             ", results['FASTA'])
    print("Alignment File:         ", results['ALIGNMENT'])


if __name__ == "__main__":
    main()
