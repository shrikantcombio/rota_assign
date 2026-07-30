"""
RNABackboneCalculator Module
----------------------------
Calculates RNA backbone torsion angles (alpha, beta, gamma, delta, epsilon, zeta)
and 7D dinucleotide suite torsion angles (delta[i-1], epsilon[i-1], zeta[i-1], alpha[i], beta[i], gamma[i], delta[i]).
"""

from pathlib import Path
from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader


class RNABackboneCalculator:
    """
    Calculator for RNA backbone 6-bond and 7D dinucleotide suite torsion angles.
    """

    BACKBONE_ATOMS = ["O3'", "P  ", "O5'", "C5'", "C4'", "C3'"]

    @classmethod
    def parse_backbone_coords(cls, structure_input):
        if isinstance(structure_input, (str, Path)):
            atoms = StructureReader.parse_structure(structure_input)
        elif isinstance(structure_input, list) and len(structure_input) > 0 and isinstance(structure_input[0], dict):
            atoms = structure_input
        else:
            atoms = StructureReader.parse_structure(structure_input)

        residue_atoms = {}
        residue_order = []

        for a in atoms:
            # Check for RNA backbone atoms
            aname = a['atom_name']
            if aname in ["O3'", "P", "O5'", "C5'", "C4'", "C3'"]:
                if aname == "P":
                    aname_std = "P"
                else:
                    aname_std = aname

                res_key = (a['resname'], a['chain'], a['resSeq'], a['icode'])
                if res_key not in residue_atoms:
                    residue_order.append(res_key)
                    residue_atoms[res_key] = {}
                if aname_std not in residue_atoms[res_key]:
                    residue_atoms[res_key][aname_std] = a['xyz']

        O3_coords, P_coords, O5_coords, C5_coords, C4_coords, C3_coords = [], [], [], [], [], []
        valid_residues = []

        for res_key in residue_order:
            r = residue_atoms[res_key]
            # Match atom names flexibility ('P' or 'P  ')
            p_xyz = r.get("P") or r.get("P  ")
            o3_xyz = r.get("O3'") or r.get("O3*")
            o5_xyz = r.get("O5'") or r.get("O5*")
            c5_xyz = r.get("C5'") or r.get("C5*")
            c4_xyz = r.get("C4'") or r.get("C4*")
            c3_xyz = r.get("C3'") or r.get("C3*")

            if p_xyz and o3_xyz and o5_xyz and c5_xyz and c4_xyz and c3_xyz:
                O3_coords.append(o3_xyz)
                P_coords.append(p_xyz)
                O5_coords.append(o5_xyz)
                C5_coords.append(c5_xyz)
                C4_coords.append(c4_xyz)
                C3_coords.append(c3_xyz)

                resname, chain, seq, icode = res_key
                formatted_label = f"{resname:<3s} {chain:1s} {seq:>4s}".strip()
                valid_residues.append(formatted_label)

        return O3_coords, P_coords, O5_coords, C5_coords, C4_coords, C3_coords, valid_residues

    def calculate_backbone(self, structure_input):
        """
        Calculates standard 6-bond RNA backbone torsion angles:
        alpha, beta, gamma, delta, epsilon, zeta.
        :returns: (labels, alpha, beta, gamma, delta, epsilon, zeta)
        """
        O3, P, O5, C5, C4, C3, valid_residues = self.parse_backbone_coords(structure_input)

        labels, alpha, beta, gamma, delta, epsilon, zeta = [], [], [], [], [], [], []
        n_res = len(valid_residues)

        for i in range(1, n_res - 1):
            try:
                alp = PDBGeometry.get_dihedral(O3[i - 1], P[i], O5[i], C5[i])
                bet = PDBGeometry.get_dihedral(P[i], O5[i], C5[i], C4[i])
                gam = PDBGeometry.get_dihedral(O5[i], C5[i], C4[i], C3[i])
                del_a = PDBGeometry.get_dihedral(C5[i], C4[i], C3[i], O3[i])
                eps = PDBGeometry.get_dihedral(C4[i], C3[i], O3[i], P[i + 1])
                zet = PDBGeometry.get_dihedral(C3[i], O3[i], P[i + 1], O5[i + 1])

                alpha.append(alp)
                beta.append(bet)
                gamma.append(gam)
                delta.append(del_a)
                epsilon.append(eps)
                zeta.append(zet)
                labels.append(valid_residues[i])
            except Exception:
                continue

        return labels, alpha, beta, gamma, delta, epsilon, zeta

    def calculate_7d_suite(self, structure_input):
        """
        Calculates 7D dinucleotide suite torsion angles:
        delm (delta[i-1]), epsm (epsilon[i-1]), zetm (zeta[i-1]), alpha[i], beta[i], gamma[i], delta[i].
        :returns: (labels, alpha, beta, gamma, delta, epsilon, zeta, delm, epsm, zetm)
        """
        O3, P, O5, C5, C4, C3, valid_residues = self.parse_backbone_coords(structure_input)

        labels = []
        alpha, beta, gamma, delta, epsilon, zeta = [], [], [], [], [], []
        delm, epsm, zetm = [], [], []

        n_res = len(valid_residues)

        for i in range(2, n_res - 1):
            try:
                alp = PDBGeometry.get_dihedral(O3[i - 1], P[i], O5[i], C5[i])
                bet = PDBGeometry.get_dihedral(P[i], O5[i], C5[i], C4[i])
                gam = PDBGeometry.get_dihedral(O5[i], C5[i], C4[i], C3[i])
                del_a = PDBGeometry.get_dihedral(C5[i], C4[i], C3[i], O3[i])
                eps = PDBGeometry.get_dihedral(C4[i], C3[i], O3[i], P[i + 1])
                zet = PDBGeometry.get_dihedral(C3[i], O3[i], P[i + 1], O5[i + 1])

                dm = PDBGeometry.get_dihedral(C5[i - 1], C4[i - 1], C3[i - 1], O3[i - 1])
                em = PDBGeometry.get_dihedral(C4[i - 1], C3[i - 1], O3[i - 1], P[i])
                zm = PDBGeometry.get_dihedral(C3[i - 1], O3[i - 1], P[i], O5[i])

                alpha.append(alp)
                beta.append(bet)
                gamma.append(gam)
                delta.append(del_a)
                epsilon.append(eps)
                zeta.append(zet)
                delm.append(dm)
                epsm.append(em)
                zetm.append(zm)

                labels.append(valid_residues[i])
            except Exception:
                continue

        return labels, alpha, beta, gamma, delta, epsilon, zeta, delm, epsm, zetm
