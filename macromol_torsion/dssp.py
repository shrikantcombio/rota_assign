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

        pdb_file = StructureReader.write_temp_pdb(file_p)
        out_p = Path(output_dssp) if output_dssp else file_p.with_suffix('.dssp')

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
        helix_tag = ["H", "G"]
        strand_tag = ["E"]
        coils_tag = ["I", "-", "S", "T", " "]

        dssp_lines = self.parse_dssp(dssp_file_path)
        result = []

        i = 0
        n = len(dssp_lines)
        while i < n:
            ss_char = dssp_lines[i][16]
            line_raw = dssp_lines[i]

            if ss_char not in helix_tag + strand_tag + coils_tag:
                if ss_char == 'B' and (i < n - 1 and dssp_lines[i + 1][16] == '-'):
                    if i < n - 2 and dssp_lines[i + 2][16] == 'B':
                        result.append((line_raw, "L"))
                        result.append((dssp_lines[i + 1], "L"))
                        result.append((dssp_lines[i + 2], "L"))
                        i += 3
                    else:
                        result.append((line_raw, "S"))
                        result.append((dssp_lines[i + 1], "S"))
                        i += 2
                else:
                    result.append((line_raw, "L"))
                    i += 1
            elif ss_char in helix_tag:
                result.append((line_raw, "H"))
                i += 1
            elif ss_char in strand_tag:
                result.append((line_raw, "S"))
                i += 1
            else:
                result.append((line_raw, "L"))
                i += 1

        return result

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
