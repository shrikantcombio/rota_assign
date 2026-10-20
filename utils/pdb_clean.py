#!/usr/bin/python
"""
Refactored pdb_clean.py wrapping clean_structure.StructureCleaner.
Provides full backward compatibility for cleaning modified residues, waters, ions, and alternate occupancies.
"""

from pathlib import Path
import sys
import os
import glob

# Ensure clean_pdb directory is in path
CLEAN_PDB_DIR = Path(__file__).resolve().parent  # bundled structure_cleaner.py
if CLEAN_PDB_DIR.exists() and str(CLEAN_PDB_DIR) not in sys.path:
    sys.path.insert(0, str(CLEAN_PDB_DIR))

from clean_structure import StructureCleaner, ModResMap

modres = ModResMap.MODRES_MAP


def verify_extension(name):
    """Verifies file extension for .pdb files."""
    if name.endswith(".pdb"):
        name = name[:-4]
    return name


def occupancy(pdb_lines):
    """Legacy wrapper for occupancy resolution."""
    cleaner = StructureCleaner()
    atoms = cleaner.parse_file(pdb_lines if isinstance(pdb_lines, (str, Path)) else pdb_lines)
    cleaned_atoms = cleaner._resolve_occupancies(atoms)
    
    # Format back to PDB lines
    lines = []
    for a in cleaned_atoms:
        if a.get('record') == 'TER':
            lines.append("TER\n")
            continue
        aname = a['atom_name']
        aname_fmt = f" {aname:<3s}" if len(aname) < 4 and not aname.startswith(('1', '2', '3', '4')) else f"{aname:<4s}"
        resname = f"{a['resname']:<3s}"[:3]
        chain = f"{a['chain']:1s}"[:1]
        try:
            res_seq = int(a['resSeq'])
        except Exception:
            res_seq = 1
        icode = f"{a['icode']:1s}"[:1]
        x, y, z = a['xyz']
        occ = float(a.get('occupancy', 1.0))
        elem = f"{a.get('elem', 'C'):>2s}"[:2]

        line = (
            f"{a.get('record', 'ATOM'):<6s}{1:>5d} {aname_fmt}{a.get('alt_loc', ' '):1s}{resname} "
            f"{chain}{res_seq:>4d}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f} 0.00          {elem}\n"
        )
        lines.append(line)
    return lines


def hetatm(pdb_lines):
    """Legacy wrapper for modified residue (HETATM) conversion."""
    cleaner = StructureCleaner(keep_waters=True, keep_ions=True)
    atoms = cleaner.parse_file(pdb_lines if isinstance(pdb_lines, (str, Path)) else pdb_lines)
    cleaned_atoms, _ = cleaner.clean_atoms(atoms)
    
    lines = []
    for a in cleaned_atoms:
        if a.get('record') == 'TER':
            lines.append("TER\n")
            continue
        aname = a['atom_name']
        aname_fmt = f" {aname:<3s}" if len(aname) < 4 and not aname.startswith(('1', '2', '3', '4')) else f"{aname:<4s}"
        resname = f"{a['resname']:<3s}"[:3]
        chain = f"{a['chain']:1s}"[:1]
        try:
            res_seq = int(a['resSeq'])
        except Exception:
            res_seq = 1
        icode = f"{a['icode']:1s}"[:1]
        x, y, z = a['xyz']
        occ = float(a.get('occupancy', 1.0))
        elem = f"{a.get('elem', 'C'):>2s}"[:2]

        line = (
            f"{a.get('record', 'ATOM'):<6s}{1:>5d} {aname_fmt}{a.get('alt_loc', ' '):1s}{resname} "
            f"{chain}{res_seq:>4d}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f} 0.00          {elem}\n"
        )
        lines.append(line)
    return lines


def main():
    cleaner = StructureCleaner()
    names = glob.glob(r'????.pdb')
    names.sort()
    for name in names:
        base_name = verify_extension(name)
        out_p, logs = cleaner.clean_structure(f"{base_name}.pdb", f"{base_name}.pdb")
        print(f"Cleaned {name} -> {out_p} ({len(logs)} modified residue conversions)")


if __name__ == "__main__":
    main()
