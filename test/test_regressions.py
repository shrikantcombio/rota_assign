#!/usr/bin/env python3
"""
test_regressions.py
-------------------
Unit tests for defects fixed in October 2026. Each test fails on the earlier code:

  - bound/unbound residues were joined by label, so relabelled or renumbered unbound chains
    lost or mis-paired residues (macromol_torsion.pairing)
  - phi/psi/omega were computed across chain breaks and chain boundaries
  - residues with insertion codes shared a label
  - rotamer assignment ignored the last chi (LEU tp vs tt were indistinguishable)
  - missing chis shifted later chis into the wrong position
  - proline was named with 6-state letters, so both ring puckers were the off-rotamer 'T'
  - symmetric terminal chis (PHE/TYR chi2, ...) were folded and then compared on a 360-degree circle

Run:  python test/test_regressions.py
"""

import math
import sys
import tempfile
import unittest
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_DIR))
sys.path.insert(0, str(BASE_DIR / "script"))

from macromol_torsion.pairing import pair_equal_residues, pair_residues_across_entries
from macromol_torsion.protein_backbone import ProteinBackboneCalculator
import get_rotastate


def _pdb_line(serial, name, resname, chain, resseq, xyz, icode=" "):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {name:<4s} {resname:3s} {chain}{resseq:4d}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {name[0]}\n")


def _write_chain(lines, chain, start_serial, residues, offset=(0.0, 0.0, 0.0)):
    """Idealised extended backbone: residue i at x = 3.8 * i (C(i)-N(i+1) = 1.33 A)."""
    serial = start_serial
    for k, (resname, resseq, x0, icode) in enumerate(residues):
        ox, oy, oz = offset
        for name, dx, dy in (("N", 0.0, 0.0), ("CA", 1.2, 0.8), ("C", 2.47, 0.0)):
            lines.append(_pdb_line(serial, name, resname, chain, resseq, (x0 + dx + ox, dy + oy, oz), icode))
            serial += 1
    return serial


class PairingTests(unittest.TestCase):
    def test_relabelled_unbound_chains_are_paired(self):
        bound = ["ILE I   18", "ASN I   19", "LYS J    5", "SER J    6"]
        unbound = ["ILE A   18", "ASN A   19", "LYS B    5", "SER B    6"]
        mapping, stats = pair_equal_residues(bound, unbound)
        self.assertEqual(stats["paired"], 4)
        self.assertEqual(mapping["LYS J    5"], "LYS B    5")

    def test_renumbered_unbound_chain_is_paired_by_sequence(self):
        bound = ["PRO A    1", "GLU A    2", "LEU A    3"]
        unbound = ["PRO A    2", "GLU A    3", "LEU A    4"]
        mapping, _ = pair_equal_residues(bound, unbound)
        self.assertEqual(mapping["PRO A    1"], "PRO A    2")

    def test_different_residue_types_are_not_paired(self):
        mapping, stats = pair_equal_residues(["ALA A    1", "GLY A    2"], ["ALA A    1", "SER A    2"])
        self.assertNotIn("GLY A    2", mapping)
        self.assertEqual(stats["unpaired_bound"], 1)

    def test_cross_entry_pairing_picks_matching_chain(self):
        seq = ["MET", "LYS", "THR", "ALA", "TYR", "ILE"]
        bound = [f"{r} D {i + 1:4d}" for i, r in enumerate(seq)]
        unbound = ([f"{r} H {i + 101:4d}" for i, r in enumerate(["GLY"] * 6)]
                   + [f"{r} L {i + 10:4d}" for i, r in enumerate(seq)])
        mapping, stats = pair_residues_across_entries(bound, unbound, "D", "HL")
        self.assertEqual(stats["chain_map"], "D>L")
        self.assertEqual(mapping["MET D    1"], "MET L   10")


    def test_homo_oligomer_copies_share_the_unbound_monomer(self):
        seq = ["SER", "GLN", "PHE", "ARG", "VAL"]
        bound = [f"{r} D {i + 1:4d}" for i, r in enumerate(seq)] + [f"{r} E {i + 1:4d}" for i, r in enumerate(seq)]
        unbound = [f"{r} A {i + 1:4d}" for i, r in enumerate(seq)]
        mapping, stats = pair_residues_across_entries(bound, unbound, "DE", "A")
        self.assertEqual(stats["paired"], 10)
        self.assertEqual(stats["chain_map"], "D>AE>A")


class BackboneTests(unittest.TestCase):
    def _calc(self, lines):
        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as f:
            f.writelines(lines)
            path = f.name
        return ProteinBackboneCalculator().calculate_dihedrals(path)

    def test_no_dihedrals_across_chain_boundary_or_break(self):
        lines = []
        serial = _write_chain(lines, "A", 1, [("ALA", 1, 0.0, " "), ("ALA", 2, 3.8, " "), ("ALA", 3, 7.6, " ")])
        # chain B starts right after chain A in the file; residue 3 of B follows a 20 A gap
        _write_chain(lines, "B", serial, [("GLY", 1, 11.4, " "), ("GLY", 2, 15.2, " "), ("GLY", 3, 39.0, " ")])
        phi, psi, omega, _, _, labels = self._calc(lines)
        idx = {lbl: i for i, lbl in enumerate(labels)}
        self.assertTrue(math.isnan(phi[idx["GLY B    1"]]), "phi of a chain's first residue must be undefined")
        self.assertTrue(math.isnan(psi[idx["ALA A    3"]]), "psi of a chain's last residue must be undefined")
        self.assertTrue(math.isnan(psi[idx["GLY B    2"]]), "psi across a chain break must be undefined")
        self.assertTrue(math.isnan(phi[idx["GLY B    3"]]), "phi across a chain break must be undefined")
        self.assertFalse(math.isnan(psi[idx["ALA A    2"]]))

    def test_insertion_codes_give_distinct_labels(self):
        lines = []
        _write_chain(lines, "A", 1, [("ALA", 52, 0.0, " "), ("SER", 52, 3.8, "A"), ("ALA", 53, 7.6, " ")])
        *_, labels = self._calc(lines)
        self.assertEqual(len(set(labels)), 3)
        self.assertIn("SER A   52A", labels)


class RotamerTests(unittest.TestCase):
    def test_classify_uses_every_rotameric_chi(self):
        self.assertEqual(get_rotastate.classify_rotamer("LEU", [-65.0, 175.0]), "mt")
        self.assertEqual(get_rotastate.classify_rotamer("LEU", [-65.0, 65.0]), "mp")
        self.assertEqual(get_rotastate.classify_rotamer("LYS", [168.7, 58.7, 55.1, -78.4]), "tppm")

    def test_sp2_terminal_chi_is_not_a_letter(self):
        self.assertEqual(get_rotastate.classify_rotamer("ASP", [-60.0, -20.0]), "m")
        self.assertEqual(get_rotastate.classify_rotamer("ARG", [-65.0, 180.0, 180.0, 85.0, 0.0]), "mtt")

    def test_missing_chi_is_not_shifted(self):
        # chi1 missing: must be NA, not "chi2 treated as chi1"
        self.assertEqual(get_rotastate.classify_rotamer("LEU", [float("nan"), 175.0]), "NA")
        self.assertEqual(get_rotastate.get_aa_rotamer("LEU", "L", "bound", "N", None, 175.0,
                                                      phi=-65.0, psi=-40.0), "NA")

    def test_library_assignment_distinguishes_last_chi(self):
        a = get_rotastate.get_aa_rotamer("LEU", "L", "bound", "N", -65.0, 175.0, phi=-65.0, psi=-40.0)
        b = get_rotastate.get_aa_rotamer("LEU", "L", "bound", "N", -177.0, 65.0, phi=-65.0, psi=-40.0)
        self.assertNotEqual(a, b)
        self.assertTrue(a.isalpha() and b.isalpha(), "LEU names are letters only (Lovell)")

    def test_symmetric_chi_compared_with_180_period(self):
        # PHE chi2 = 91 is the ring orientation of -89: the same library row must match both
        self.assertAlmostEqual(get_rotastate._chi_diff("PHE", 1, 91.0, -89.0), 0.0)
        self.assertAlmostEqual(get_rotastate._chi_diff("PHE", 1, 89.0, -89.0), -2.0)
        self.assertAlmostEqual(get_rotastate._chi_diff("PHE", 0, 89.0, -89.0), 178.0)
        lib = {"m": [[-65.0, 95.0, 0.6], [-65.0, 10.0, 0.4]]}
        self.assertEqual(get_rotastate._best_rotamer("PHE", [-65.0, -86.0], lib),
                         get_rotastate._best_rotamer("PHE", [-65.0, 94.0], lib))

    def test_proline_named_by_ring_pucker(self):
        self.assertEqual(get_rotastate.classify_rotamer("PRO", [28.0, -35.0]), "endo")
        self.assertEqual(get_rotastate.classify_rotamer("PRO", [-35.0, 40.0]), "exo")
        self.assertEqual(get_rotastate.get_aa_rotamer("PRO", "L", "bound", "N", 25.0, -35.0,
                                                      phi=-65.0, psi=-35.0), "endo")
        self.assertEqual(get_rotastate.get_aa_rotamer("PRO", "L", "bound", "N", -28.0, 38.0,
                                                      phi=-65.0, psi=-35.0), "exo")

    def test_v2_libraries_are_the_default(self):
        import build_rotalib
        self.assertEqual(build_rotalib.LIBRARY_VERSION, "v2")
        self.assertTrue(build_rotalib.LIBRARY_SPECS["BBD"][1].endswith("_v2"))


if __name__ == "__main__":
    unittest.main()
