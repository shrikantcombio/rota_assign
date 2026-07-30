"""
RNAPseudoTorsionCalculator Module
----------------------------------
Calculates standard (eta, theta) and modified (eta', theta') pseudo-torsion angles for RNA residues.
"""

from pathlib import Path
from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader


class RNAPseudoTorsionCalculator:
    """
    Calculator for RNA pseudo-torsion angles.
    """

    @classmethod
    def calculate_pseudo_torsions(cls, structure_input):
        """
        Calculates pseudo-torsion angles: eta, theta, eta_prime (etad), theta_prime (thetad).
        :returns: (labels, eta, theta, etad, thetad)
        """
        if isinstance(structure_input, (str, Path)):
            atoms = StructureReader.parse_structure(structure_input)
        elif isinstance(structure_input, list) and len(structure_input) > 0 and isinstance(structure_input[0], dict):
            atoms = structure_input
        else:
            atoms = StructureReader.parse_structure(structure_input)

        residue_atoms = {}
        residue_order = []

        for a in atoms:
            aname = a['atom_name']
            if aname in ["P", "C4'", "C1'"]:
                res_key = (a['resname'], a['chain'], a['resSeq'], a['icode'])
                if res_key not in residue_atoms:
                    residue_order.append(res_key)
                    residue_atoms[res_key] = {}
                if aname not in residue_atoms[res_key]:
                    residue_atoms[res_key][aname] = a['xyz']

        P_coords, C4_coords, C1_coords = [], [], []
        valid_residues = []

        for res_key in residue_order:
            r = residue_atoms[res_key]
            p_xyz = r.get("P")
            c4_xyz = r.get("C4'")
            c1_xyz = r.get("C1'")

            if p_xyz and c4_xyz and c1_xyz:
                P_coords.append(p_xyz)
                C4_coords.append(c4_xyz)
                C1_coords.append(c1_xyz)

                resname, chain, seq, icode = res_key
                formatted_label = f"{resname:<3s} {chain:1s} {seq:>4s}".strip()
                valid_residues.append(formatted_label)

        labels, eta, theta, etad, thetad = [], [], [], [], []
        n_res = len(valid_residues)

        for i in range(1, n_res - 1):
            try:
                e_val = PDBGeometry.get_dihedral(C4_coords[i - 1], P_coords[i], C4_coords[i], P_coords[i + 1])
                t_val = PDBGeometry.get_dihedral(P_coords[i], C4_coords[i], P_coords[i + 1], C4_coords[i + 1])

                ed_val = PDBGeometry.get_dihedral(C1_coords[i - 1], P_coords[i], C1_coords[i], P_coords[i + 1])
                td_val = PDBGeometry.get_dihedral(P_coords[i], C4_coords[i], P_coords[i + 1], C4_coords[i + 1])

                eta.append(e_val)
                theta.append(t_val)
                etad.append(ed_val)
                thetad.append(td_val)
                labels.append(valid_residues[i])
            except Exception:
                continue

        return labels, eta, theta, etad, thetad
