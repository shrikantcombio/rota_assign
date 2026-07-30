"""
ProteinBackboneCalculator Module
--------------------------------
Calculates protein mainchain backbone torsion angles (Phi, Psi, Omega, Tau, Xi).
"""

from pathlib import Path
from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader

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
        valid_residues = []

        for res_key in residue_order:
            r_dict = residue_atoms[res_key]
            if "N" in r_dict and "CA" in r_dict and "C" in r_dict:
                N_coords.append(r_dict["N"])
                CA_coords.append(r_dict["CA"])
                C_coords.append(r_dict["C"])
                resname, chain, seq, icode = res_key
                formatted_label = f"{resname:<3s} {chain:1s} {seq:>4s}".strip()
                valid_residues.append(formatted_label)

        phi, psi, omega, tau, xi, labels = [], [], [], [], [], []
        n_res = len(valid_residues)

        for i in range(1, n_res - 1):
            try:
                p_angle = PDBGeometry.get_dihedral(C_coords[i - 1], N_coords[i], CA_coords[i], C_coords[i])
                s_angle = PDBGeometry.get_dihedral(N_coords[i], CA_coords[i], C_coords[i], N_coords[i + 1])
                o_angle = PDBGeometry.get_dihedral(CA_coords[i], C_coords[i], N_coords[i + 1], CA_coords[i + 1])
                t_angle = PDBGeometry.tau_angle(N_coords[i], CA_coords[i], C_coords[i])

                phi.append(p_angle)
                psi.append(s_angle)
                omega.append(o_angle)
                tau.append(t_angle)
                labels.append(valid_residues[i])
            except Exception:
                continue

        for i in range(1, n_res - 2):
            try:
                x_angle = PDBGeometry.get_dihedral(CA_coords[i - 1], CA_coords[i], CA_coords[i + 1], CA_coords[i + 2])
                xi.append(x_angle)
            except Exception:
                xi.append(0.0)

        return phi, psi, omega, tau, xi, labels
