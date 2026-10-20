"""
ProteinBackboneCalculator Module
--------------------------------
Calculates protein mainchain backbone torsion angles (Phi, Psi, Omega, Tau, Xi).
"""

from pathlib import Path
from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader

# Longest C(i)-N(i+1) distance still treated as a peptide bond (ideal 1.33 A)
MAX_PEPTIDE_BOND = 2.0

STANDARD_AA = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
]


class ProteinBackboneCalculator:
    """
    Calculator for protein backbone torsion angles.
    """

    def __init__(self, allowed_residues=None):
        self.allowed_residues = allowed_residues or STANDARD_AA

    def calculate_dihedrals(self, structure_input):
        if isinstance(structure_input, (str, Path)):
            atoms = StructureReader.parse_structure(structure_input)
        elif isinstance(structure_input, list) and len(structure_input) > 0 and isinstance(structure_input[0], dict):
            atoms = structure_input
        else:
            atoms = StructureReader.parse_structure(structure_input)

        main_chain = ["N", "CA", "C"]

        residue_atoms = {}
        residue_order = []

        for a in atoms:
            if a['atom_name'] in main_chain and a['record'] == 'ATOM':
                res_key = (a['resname'], a['chain'], a['resSeq'], a['icode'])
                if res_key not in residue_atoms:
                    residue_order.append(res_key)
                    residue_atoms[res_key] = {}
                if a['atom_name'] not in residue_atoms[res_key]:
                    residue_atoms[res_key][a['atom_name']] = a['xyz']

        N_coords = []
        CA_coords = []
        C_coords = []
        chains = []
        valid_residues = []

        for res_key in residue_order:
            r_dict = residue_atoms[res_key]
            if "N" in r_dict and "CA" in r_dict and "C" in r_dict:
                N_coords.append(r_dict["N"])
                CA_coords.append(r_dict["CA"])
                C_coords.append(r_dict["C"])
                resname, chain, seq, icode = res_key
                chains.append(chain)
                formatted_label = f"{resname:<3s} {chain:1s} {seq:>4s}{(icode or '').strip()}".strip()
                valid_residues.append(formatted_label)

        n_res = len(valid_residues)

        # Residue i is peptide-bonded to i + 1 only within a chain and with C(i)-N(i+1) <= 2.0 A;
        # dihedrals spanning a chain break or a chain boundary are undefined (NaN)
        bonded = [
            chains[i] == chains[i + 1] and PDBGeometry.distance(C_coords[i], N_coords[i + 1]) <= MAX_PEPTIDE_BOND
            for i in range(n_res - 1)
        ]

        def _dihedral(*pts):
            try:
                return PDBGeometry.get_dihedral(*pts)
            except Exception:
                return float("nan")

        phi, psi, omega, tau, xi, labels = [], [], [], [], [], []
        for i in range(n_res):
            prev_ok = i > 0 and bonded[i - 1]
            next_ok = i < n_res - 1 and bonded[i]
            phi.append(_dihedral(C_coords[i - 1], N_coords[i], CA_coords[i], C_coords[i]) if prev_ok else float("nan"))
            psi.append(_dihedral(N_coords[i], CA_coords[i], C_coords[i], N_coords[i + 1]) if next_ok else float("nan"))
            omega.append(_dihedral(CA_coords[i], C_coords[i], N_coords[i + 1], CA_coords[i + 1]) if next_ok else float("nan"))
            tau.append(PDBGeometry.tau_angle(N_coords[i], CA_coords[i], C_coords[i]))
            xi_ok = prev_ok and next_ok and i < n_res - 2 and bonded[i + 1]
            xi.append(_dihedral(CA_coords[i - 1], CA_coords[i], CA_coords[i + 1], CA_coords[i + 2]) if xi_ok else float("nan"))
            labels.append(valid_residues[i])

        return phi, psi, omega, tau, xi, labels
