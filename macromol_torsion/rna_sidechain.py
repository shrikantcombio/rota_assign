"""
RNASidechainCalculator Module
-----------------------------
Calculates RNA glycosidic side-chain Chi angle (chi1/chi2) and inter-base geometric parameters
(C1'-C1' distance, N1/N9-N1/N9 distance, mu pseudo-torsion angle).
"""

from pathlib import Path
from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader


class RNASidechainCalculator:
    """
    Calculator for RNA glycosidic Chi angles and inter-nucleotide distances.
    """

    PURINE_RES = ["A", "G", "ADE", "GUA"]
    PYRIMIDINE_RES = ["C", "U", "CYT", "URA"]

    @classmethod
    def parse_rna_sidechain_atoms(cls, structure_input):
        if isinstance(structure_input, (str, Path)):
            atoms = StructureReader.parse_structure(structure_input)
        elif isinstance(structure_input, list) and len(structure_input) > 0 and isinstance(structure_input[0], dict):
            atoms = structure_input
        else:
            atoms = StructureReader.parse_structure(structure_input)

        residue_atoms = {}
        residue_order = []

        needed_atoms = ["O4'", "C1'", "C2'", "N9", "C4", "N1", "C2"]

        for a in atoms:
            aname = a['atom_name']
            if aname in needed_atoms:
                res_key = (a['resname'], a['chain'], a['resSeq'], a['icode'])
                if res_key not in residue_atoms:
                    residue_order.append(res_key)
                    residue_atoms[res_key] = {}
                if aname not in residue_atoms[res_key]:
                    residue_atoms[res_key][aname] = a['xyz']

        return residue_atoms, residue_order

    def calculate_sidechain_and_distances(self, structure_input):
        """
        Calculates RNA glycosidic Chi angles, inter-nucleotide distances, and mu pseudo-torsions.
        :returns: (labels, chi1, chi2, nn_dist, cc_dist, mu)
        """
        residue_atoms, residue_order = self.parse_rna_sidechain_atoms(structure_input)

        O4_coords, C1_coords, C2_coords = [], [], []
        N_coords, C_base_coords = [], []
        valid_labels = []

        for res_key in residue_order:
            r = residue_atoms[res_key]
            resname, chain, seq, icode = res_key
            resname_upper = resname.upper()

            o4 = r.get("O4'")
            c1 = r.get("C1'")
            c2 = r.get("C2'")

            is_purine = resname_upper in self.PURINE_RES or resname_upper.endswith("A") or resname_upper.endswith("G")
            n_base = r.get("N9") if is_purine else r.get("N1")
            c_base = r.get("C4") if is_purine else r.get("C2")

            if not n_base:
                n_base = r.get("N9") or r.get("N1")
            if not c_base:
                c_base = r.get("C4") or r.get("C2")

            if o4 and c1 and c2 and n_base and c_base:
                O4_coords.append(o4)
                C1_coords.append(c1)
                C2_coords.append(c2)
                N_coords.append(n_base)
                C_base_coords.append(c_base)

                formatted_label = f"{resname:<3s} {chain:1s} {seq:>4s}".strip()
                valid_labels.append(formatted_label)

        labels, chi1, chi2 = [], [], []
        nn_dist, cc_dist, mu = [], [], []

        n_res = len(valid_labels)

        # 1. Glycosidic Chi angles & mu pseudo-torsion for adjacent pairs
        for i in range(1, n_res):
            try:
                c1_prev = PDBGeometry.get_dihedral(O4_coords[i - 1], C1_coords[i - 1], N_coords[i - 1], C_base_coords[i - 1])
                c2_curr = PDBGeometry.get_dihedral(O4_coords[i], C1_coords[i], N_coords[i], C_base_coords[i])
                mu_val = PDBGeometry.get_dihedral(N_coords[i - 1], C1_coords[i - 1], C1_coords[i], N_coords[i])

                nn = PDBGeometry.distance(N_coords[i - 1], N_coords[i])
                cc = PDBGeometry.distance(C1_coords[i - 1], C1_coords[i])

                chi1.append(c1_prev)
                chi2.append(c2_curr)
                mu.append(mu_val)
                nn_dist.append(round(float(nn), 2))
                cc_dist.append(round(float(cc), 2))
                labels.append(valid_labels[i - 1])
            except Exception:
                continue

        return labels, chi1, chi2, nn_dist, cc_dist, mu
