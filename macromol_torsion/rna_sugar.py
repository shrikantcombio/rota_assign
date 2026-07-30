"""
RNASugarPuckerCalculator Module
-------------------------------
Calculates RNA sugar ring torsion angles (nu0, nu1, nu2, nu3, nu4), Altona & Sundaralingam
pseudorotation phase angle P, maximum amplitude tau_m, and 10-state sugar pucker classification.
"""

import math
from math import sin, cos, radians, atan2, degrees
from pathlib import Path
from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader


class RNASugarPuckerCalculator:
    """
    Calculator for RNA sugar ring pucker torsions and pseudorotation parameters.
    """

    SUGAR_ATOMS = ["O4'", "C1'", "C2'", "C3'", "C4'"]

    @classmethod
    def calculate_ring_torsions(cls, structure_input):
        """
        Calculates ring torsions V0, V1, V2, V3, V4.
        :returns: (labels, V0, V1, V2, V3, V4)
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
            if aname in cls.SUGAR_ATOMS:
                res_key = (a['resname'], a['chain'], a['resSeq'], a['icode'])
                if res_key not in residue_atoms:
                    residue_order.append(res_key)
                    residue_atoms[res_key] = {}
                if aname not in residue_atoms[res_key]:
                    residue_atoms[res_key][aname] = a['xyz']

        V0, V1, V2, V3, V4, labels = [], [], [], [], [], []

        for res_key in residue_order:
            r = residue_atoms[res_key]
            o4 = r.get("O4'")
            c1 = r.get("C1'")
            c2 = r.get("C2'")
            c3 = r.get("C3'")
            c4 = r.get("C4'")

            if o4 and c1 and c2 and c3 and c4:
                try:
                    v0_val = PDBGeometry.get_dihedral(c4, o4, c1, c2)
                    v1_val = PDBGeometry.get_dihedral(o4, c1, c2, c3)
                    v2_val = PDBGeometry.get_dihedral(c1, c2, c3, c4)
                    v3_val = PDBGeometry.get_dihedral(c2, c3, c4, o4)
                    v4_val = PDBGeometry.get_dihedral(c3, c4, o4, c1)

                    V0.append(v0_val)
                    V1.append(v1_val)
                    V2.append(v2_val)
                    V3.append(v3_val)
                    V4.append(v4_val)

                    resname, chain, seq, icode = res_key
                    formatted_label = f"{resname:<3s} {chain:1s} {seq:>4s}".strip()
                    labels.append(formatted_label)
                except Exception:
                    continue

        return labels, V0, V1, V2, V3, V4

    @classmethod
    def calculate_pseudorotation(cls, structure_input):
        """
        Calculates Altona & Sundaralingam pseudorotation parameters:
        phase angle P, maximum amplitude tau_m, and sugar pucker conformation string.
        :returns: List of dicts per residue: [{'label': str, 'phase': float, 'amplitude': float, 'pucker': str}, ...]
        """
        labels, V0, V1, V2, V3, V4 = cls.calculate_ring_torsions(structure_input)

        constant = 2 * (math.sin(math.radians(36.0)) + math.sin(math.radians(72.0)))
        results = []

        for i in range(len(labels)):
            v0_f = float(V0[i])
            v1_f = float(V1[i])
            v2_f = float(V2[i])
            v3_f = float(V3[i])
            v4_f = float(V4[i])

            phase_rad = math.atan2((v4_f + v1_f) - (v3_f + v0_f), v2_f * constant)
            amplitude = v2_f / math.cos(phase_rad) if math.cos(phase_rad) != 0 else v2_f

            phase_deg = (math.degrees(phase_rad) + 72.0) % 360.0
            if phase_deg < 0:
                phase_deg += 360.0

            if phase_deg < 36.0:
                pucker = "C3'-endo"
            elif phase_deg < 72.0:
                pucker = "C4'-exo"
            elif phase_deg < 108.0:
                pucker = "O4'-endo"
            elif phase_deg < 144.0:
                pucker = "C1'-exo"
            elif phase_deg < 180.0:
                pucker = "C2'-endo"
            elif phase_deg < 216.0:
                pucker = "C3'-exo"
            elif phase_deg < 252.0:
                pucker = "C4'-endo"
            elif phase_deg < 288.0:
                pucker = "O4'-exo"
            elif phase_deg < 324.0:
                pucker = "C1'-endo"
            else:
                pucker = "C2'-exo"

            results.append({
                'label': labels[i],
                'phase': round(float(phase_deg), 2),
                'amplitude': round(float(amplitude), 2),
                'pucker': pucker,
                'v0': v0_f, 'v1': v1_f, 'v2': v2_f, 'v3': v3_f, 'v4': v4_f
            })

        return results
