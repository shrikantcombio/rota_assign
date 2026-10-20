#!/usr/bin/python
"""
Refactored pdb_HETATM.py wrapping clean_structure.StructureCleaner.
Provides full backward compatibility for cleaning HETATM modified residues in PDB files.
"""

from pathlib import Path
import sys
import os
import glob

CLEAN_PDB_DIR = Path(__file__).resolve().parent  # bundled structure_cleaner.py
if CLEAN_PDB_DIR.exists() and str(CLEAN_PDB_DIR) not in sys.path:
    sys.path.insert(0, str(CLEAN_PDB_DIR))

from clean_structure import StructureCleaner, ModResMap

modres = ModResMap.MODRES_MAP


def hetatm(pdb_path):
    """Refactored hetatm function converting modified residues to ATOM records."""
    cleaner = StructureCleaner(keep_waters=True, keep_ions=True)
    in_p = Path(pdb_path)
    out_mod = in_p.with_name(f"{in_p.stem}_mod.pdb")
    out_log = in_p.with_name(f"{in_p.stem}.log")

    out_p, logs = cleaner.clean_structure(in_p, out_mod)

    with open(out_log, 'w') as f_log:
        f_log.write(f"{pdb_path}:\t The following HETATMs were modified for this structure:\n")
        for l in logs:
            f_log.write(f"{l}\n")
        f_log.write("\n")


def main():
    names = glob.glob(r'????.pdb')
    names.sort()

    het_cleaned_dir = Path("HETCleaned")
    het_cleaned_dir.mkdir(exist_ok=True)

    cleaner = StructureCleaner()

    for i, file_name in enumerate(names, start=1):
        print(f"{i} {file_name}")
        in_p = Path(file_name)
        out_p = het_cleaned_dir / f"{in_p.stem}_cln.pdb"
        cleaner.clean_structure(in_p, out_p)


if __name__ == "__main__":
    main()
