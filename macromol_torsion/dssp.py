"""
DsspAnalyzer Module
-------------------
Secondary structure assignment using DSSP and 3-state Rost & Sander correction mapping.
"""

from pathlib import Path
import os
import shutil
import subprocess
from macromol_torsion.structure_reader import StructureReader

AA_1TO3 = {
    "R": "ARG", "K": "LYS", "D": "ASP", "E": "GLU", "N": "ASN", "Q": "GLN", "H": "HIS",
    "S": "SER", "T": "THR", "Y": "TYR", "W": "TRP", "F": "PHE", "C": "CYS", "M": "MET",
    "L": "LEU", "I": "ILE", "P": "PRO", "V": "VAL", "A": "ALA", "G": "GLY"
}


class DsspAnalyzer:
    """
    Interface to DSSP executable for secondary structure assignment.
    """

    def __init__(self, dssp_binary=None):
        if dssp_binary:
            self.dssp_binary = str(dssp_binary)
        else:
            default_bin = Path("/home/labuser/bin/dssp")
            if default_bin.exists():
                self.dssp_binary = str(default_bin)
            else:
                self.dssp_binary = shutil.which("dssp") or shutil.which("mkdssp") or "dssp"

    def run_dssp(self, structure_file, output_dssp=None):
        file_p = Path(structure_file)
        if file_p.suffix.lower() == '.dssp':
            return file_p

        stem = file_p.stem.split('.')[0].split('_')[0]
        dssp_cache_dir = Path("/home/labuser/Projects/results/PRDBv35_clean/proRNA_data/dssp")
        if output_dssp:
            out_p = Path(output_dssp)
        elif dssp_cache_dir.exists():
            out_p = dssp_cache_dir / f"{stem}.dssp"
        else:
            out_p = file_p.with_suffix('.dssp')

        if out_p.exists() and out_p.stat().st_size > 0:
            return out_p

        dssp_same = file_p.with_suffix('.dssp')
        if dssp_same.exists() and dssp_same.stat().st_size > 0:
            return dssp_same

        out_p.parent.mkdir(parents=True, exist_ok=True)
        pdb_file = StructureReader.write_temp_pdb(file_p)

        cmd = f"{self.dssp_binary} -i {pdb_file} -o {out_p}"
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if not out_p.exists() or out_p.stat().st_size == 0:
            cmd_fallback = f"{self.dssp_binary} {pdb_file} {out_p}"
            subprocess.run(cmd_fallback, shell=True, capture_output=True, text=True)

        return out_p

    def parse_dssp(self, dssp_file_path):
        dssp_p = Path(dssp_file_path)
        if not dssp_p.exists():
            raise FileNotFoundError(f"DSSP file not found: {dssp_p}")

        with open(dssp_p, 'r') as f:
            lines = f.readlines()

        res_details = []
        flag = False
        for line in lines:
            if line.startswith("  #"):
                flag = True
                continue
            if flag:
                if len(line) > 16 and line[:5].strip().isdigit() and line[13] != "!":
                    res_details.append(line)
        return res_details

    def corrected_assignment(self, dssp_file_path):
        """
        Secondary structure 3-state classification and Rost & Sander (1993) correction rules:
        - Helices (H): 'H' (alpha-helix), 'G' (3_10-helix), 'I' (pi-helix)
        - Strands (S): 'E' (extended strand)
        - Loops (L): 'B' (isolated beta-bridge), 'T' (turn), 'S' (bend), ' ', '-', 'C' (coils)

        Rost & Sander (1993) correction rules for isolated beta-bridges ('B'):
        1. An isolated beta-bridge ('B') preceded or followed by a coil ('_') is considered as two strands: 'B_' -> 'SS' or '_B' -> 'SS'.
        2. Two isolated beta-bridges interrupted by a coil ('B_B') is considered as three loops: 'B_B' -> 'LLL'.
        """
        dssp_lines = self.parse_dssp(dssp_file_path)
        n = len(dssp_lines)
        if n == 0:
            return []

        raw_codes = [line[16] for line in dssp_lines]

        # Step 1: Initial 3-state mapping
        mapped = []
        for c in raw_codes:
            if c in ['H', 'G', 'I']:
                mapped.append('H')
            elif c == 'E':
                mapped.append('S')
            else:
                mapped.append('L')

        final_ss = list(mapped)

        # Step 2: Apply Rost & Sander (1993) rules
        i = 0
        while i < n:
            c = raw_codes[i]
            if c == 'B':
                # Rule 2: Two isolated beta-bridges interrupted by coil (B _ B) -> 3 loops (LLL)
                if i + 2 < n and raw_codes[i + 2] == 'B' and raw_codes[i + 1] in ['-', ' ', 'S', 'T']:
                    final_ss[i] = 'L'
                    final_ss[i + 1] = 'L'
                    final_ss[i + 2] = 'L'
                    i += 3
                    continue

                # Rule 1: Isolated beta-bridge preceded or followed by coil (B _ = SS or _ B = SS)
                next_is_coil = (i + 1 < n and raw_codes[i + 1] in ['-', ' ', 'S', 'T'])
                prev_is_coil = (i - 1 >= 0 and raw_codes[i - 1] in ['-', ' ', 'S', 'T'])
                if next_is_coil or prev_is_coil:
                    final_ss[i] = 'S'
                    if next_is_coil and i + 1 < n:
                        final_ss[i + 1] = 'S'
                    if prev_is_coil and i - 1 >= 0:
                        final_ss[i - 1] = 'S'
            i += 1

        return [(dssp_lines[k], final_ss[k]) for k in range(n)]

    def get_secondary_structure(self, structure_file, output_dssp=None):
        dssp_file = self.run_dssp(structure_file, output_dssp=output_dssp)
        assigned_pairs = self.corrected_assignment(dssp_file)

        ss_dict = {}
        for raw_line, ss_code in assigned_pairs:
            try:
                res_seq = raw_line[5:10].strip()
                chain = raw_line[11]
                aa_1 = raw_line[13]

                aa_3 = AA_1TO3.get(aa_1, "UNK")
                res_label = f"{aa_3:<3s} {chain:1s} {res_seq:>4s}".strip()

                # Store by full res_label AND by (chain, res_seq) tuple
                ss_dict[res_label] = ss_code
                ss_dict[(chain, res_seq)] = ss_code
            except Exception:
                continue

        return ss_dict
