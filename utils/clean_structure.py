"""
utils/clean_structure.py
-------------------------
Imports StructureCleaner from /home/labuser/clean_pdb/clean_structure.py using importlib
to avoid circular import issues while maintaining full backward compatibility.
"""

from pathlib import Path
import sys
import importlib.util

CLEAN_PDB_DIR = Path("/home/labuser/clean_pdb")
CLEAN_PDB_FILE = CLEAN_PDB_DIR / "clean_structure.py"

if CLEAN_PDB_FILE.exists():
    spec = importlib.util.spec_from_file_location("clean_pdb_module", CLEAN_PDB_FILE)
    clean_pdb_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(clean_pdb_module)

    StructureCleaner = clean_pdb_module.StructureCleaner
    ModResMap = clean_pdb_module.MODRES_MAP
else:
    raise ImportError(f"Cannot locate clean_structure.py at {CLEAN_PDB_FILE}")

__all__ = ["StructureCleaner", "ModResMap"]
