#!/usr/bin/python3
"""
Refactored pdb_module.py wrapping EqualResidueAligner.
Maintains full backward compatibility for get_res_details, pdbTofasta, clustalw2, mapping_residue, and main().
New implementation: equal_residue_aligner.py
"""

from pathlib import Path
import sys
import os

UTILS_DIR = Path(__file__).resolve().parent
if str(UTILS_DIR) not in sys.path:
    sys.path.insert(0, str(UTILS_DIR))

from equal_residue_aligner import EqualResidueAligner, AA_3TO1

aa3to1 = AA_3TO1


class get_res_details:
    def __init__(self, data=None, pdb=None, chains=None):
        aligner = EqualResidueAligner()
        self.result = []
        if data is not None:
            self.result = data
        elif pdb is not None:
            pdb_p = Path(pdb) if pdb.endswith('.pdb') else Path(f"{pdb}.pdb")
            if pdb_p.exists():
                with open(pdb_p, 'r') as f:
                    for line in f:
                        if line.startswith("ATOM") and (chains is None or line[21] in chains):
                            self.result.append(line)

        self.resTag_index = []
        self.res_dict = {}
        self._parse_res_details()

    def _parse_res_details(self):
        for line in self.result:
            if line.startswith('ATOM'):
                resNo = line[21:27]
                resTag = line[17:27].strip()
                if resTag not in self.resTag_index:
                    self.resTag_index.append(resTag)
                    self._get_atoms(resTag, resNo)

    def _get_atoms(self, resTag, resno):
        atomList = []
        atomDict = {}
        for line in self.result:
            if line.startswith("ATOM") and line[21:27] == resno:
                atomNm = line[13:17]
                if atomNm not in atomList:
                    atomList.append(atomNm)
                    atomDict[atomNm] = line
        self.res_dict[resTag] = atomDict


class clustalw2:
    def __init__(self, filename):
        self.filename = filename
        self.do_clustalw2()

    def do_clustalw2(self):
        aligner = EqualResidueAligner()
        aligner.run_clustalw2(self.filename)


def pdbTofasta(pdb_file, chains):
    aligner = EqualResidueAligner()
    b_name, b_dict, b_list = aligner.extract_sequence(pdb_file, chains)
    
    fasta_file = f"{b_name}.fa"
    with open(fasta_file, 'w') as fa:
        for c in b_list:
            if c in chains:
                fa.write(f">{b_name}:{c}\n{b_dict[c]}\n")

    return b_name, b_dict, b_list


def mapping_residue(alignFile, bound, bChain, unbound, ubChain):
    aligner = EqualResidueAligner()
    res = aligner.align_and_extract_equal_residues(
        bound_file=f"{bound}.pdb" if not bound.endswith('.pdb') else bound,
        bound_chains=bChain,
        unbound_file=f"{unbound}.pdb" if not unbound.endswith('.pdb') else unbound,
        unbound_chains=ubChain,
        output_dir=bound
    )
    return res


def main():
    list_file = Path('List_1.txt')
    if not list_file.exists():
        print("List_1.txt not found.")
        return

    aligner = EqualResidueAligner()

    with open(list_file, 'r') as fileList:
        for f in fileList:
            parts = f.split()
            if len(parts) < 3:
                continue

            complx = parts[0]
            p_chain = parts[1].split(':')[0]
            unbound = parts[2].split(':')[0]
            ub_chain = parts[2].split(':')[1] if ':' in parts[2] else 'A'

            aligner.align_and_extract_equal_residues(
                bound_file=f"{complx}.pdb",
                bound_chains=p_chain,
                unbound_file=f"{unbound}.pdb",
                unbound_chains=ub_chain,
                output_dir=complx
            )


if __name__ == "__main__":
    main()
