"""
ProteinSidechainCalculator Module
----------------------------------
Calculates side-chain torsion angles (Chi1, Chi2, Chi3, Chi4, Chi5) for protein residues.
"""

from pathlib import Path
from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader


class ProteinSidechainCalculator:
    """
    Calculator for protein side-chain Chi torsion angles.
    """

    CHI1_ATOMS = ["N", "CA", "CB", "CG", "SG", "CG1", "OG", "OG1"]
    CHI2_ATOMS = ["CA", "CB", "CG", "CG1", "CD", "OD1", "ND1", "CD1", "SD"]
    CHI3_ATOMS = ["CB", "CG", "CD", "SD", "NE", "OE1", "CE"]
    CHI4_ATOMS = ["CG", "CD", "NE", "CZ", "CE", "NZ"]
    CHI5_ATOMS = ["CD", "NE", "CZ", "NH1"]

    @classmethod
    def _parse_atoms_by_residue(cls, structure_input, needed_atoms):
        if isinstance(structure_input, (str, Path)):
            atoms = StructureReader.parse_structure(structure_input)
        elif isinstance(structure_input, list) and len(structure_input) > 0 and isinstance(structure_input[0], dict):
            atoms = structure_input
        else:
            atoms = StructureReader.parse_structure(structure_input)

        res_dict = {}
        res_list = []

        for a in atoms:
            if a['record'] in ['ATOM', 'HETATM'] and a['resname'] not in ['ALA', 'GLY']:
                res_key = a['res_label']
                if res_key not in res_dict:
                    res_list.append(res_key)
                    res_dict[res_key] = {}
                if a['atom_name'] in needed_atoms and a['atom_name'] not in res_dict[res_key]:
                    res_dict[res_key][a['atom_name']] = a['xyz']

        return res_dict, res_list

    @classmethod
    def calculate_chi1(cls, structure_input):
        res_dict, res_list = cls._parse_atoms_by_residue(structure_input, cls.CHI1_ATOMS)
        key_list, chi1 = [], []

        for key in res_list:
            r = res_dict[key]
            n_coord = r.get("N")
            ca_coord = r.get("CA")
            cb_coord = r.get("CB")
            cg_coord = next((r[k] for k in ["CG", "SG", "CG1", "OG", "OG1"] if k in r), None)

            if n_coord and ca_coord and cb_coord and cg_coord:
                try:
                    val = PDBGeometry.get_dihedral(n_coord, ca_coord, cb_coord, cg_coord)
                    key_list.append(key)
                    chi1.append(val)
                except Exception:
                    continue

        return key_list, chi1

    @classmethod
    def calculate_chi2(cls, structure_input):
        res_dict, res_list = cls._parse_atoms_by_residue(structure_input, cls.CHI2_ATOMS)
        key_list, chi2 = [], []

        for key in res_list:
            r = res_dict[key]
            ca_coord = r.get("CA")
            cb_coord = r.get("CB")
            cg_coord = next((r[k] for k in ["CG", "CG1"] if k in r), None)
            cd_coord = next((r[k] for k in ["CD", "OD1", "ND1", "CD1", "SD"] if k in r), None)

            if ca_coord and cb_coord and cg_coord and cd_coord:
                try:
                    val = PDBGeometry.get_dihedral(ca_coord, cb_coord, cg_coord, cd_coord)
                    key_list.append(key)
                    chi2.append(val)
                except Exception:
                    continue

        return key_list, chi2

    @classmethod
    def calculate_chi3(cls, structure_input):
        res_dict, res_list = cls._parse_atoms_by_residue(structure_input, cls.CHI3_ATOMS)
        key_list, chi3 = [], []

        for key in res_list:
            r = res_dict[key]
            cb_coord = r.get("CB")
            cg_coord = r.get("CG")
            cd_coord = next((r[k] for k in ["CD", "SD"] if k in r), None)
            ne_coord = next((r[k] for k in ["NE", "OE1", "CE"] if k in r), None)

            if cb_coord and cg_coord and cd_coord and ne_coord:
                try:
                    val = PDBGeometry.get_dihedral(cb_coord, cg_coord, cd_coord, ne_coord)
                    key_list.append(key)
                    chi3.append(val)
                except Exception:
                    continue

        return key_list, chi3

    @classmethod
    def calculate_chi4(cls, structure_input):
        res_dict, res_list = cls._parse_atoms_by_residue(structure_input, cls.CHI4_ATOMS)
        key_list, chi4 = [], []

        for key in res_list:
            r = res_dict[key]
            cg_coord = r.get("CG")
            cd_coord = r.get("CD")
            ne_coord = next((r[k] for k in ["NE", "CE"] if k in r), None)
            cz_coord = next((r[k] for k in ["CZ", "NZ"] if k in r), None)

            if cg_coord and cd_coord and ne_coord and cz_coord:
                try:
                    val = PDBGeometry.get_dihedral(cg_coord, cd_coord, ne_coord, cz_coord)
                    key_list.append(key)
                    chi4.append(val)
                except Exception:
                    continue

        return key_list, chi4

    @classmethod
    def calculate_chi5(cls, structure_input):
        res_dict, res_list = cls._parse_atoms_by_residue(structure_input, cls.CHI5_ATOMS)
        key_list, chi5 = [], []

        for key in res_list:
            r = res_dict[key]
            cd_coord = r.get("CD")
            ne_coord = r.get("NE")
            cz_coord = r.get("CZ")
            nh1_coord = r.get("NH1")

            if cd_coord and ne_coord and cz_coord and nh1_coord:
                try:
                    val = PDBGeometry.get_dihedral(cd_coord, ne_coord, cz_coord, nh1_coord)
                    key_list.append(key)
                    chi5.append(val)
                except Exception:
                    continue

        return key_list, chi5

    def calculate_all_chis(self, structure_input):
        k1, v1 = self.calculate_chi1(structure_input)
        k2, v2 = self.calculate_chi2(structure_input)
        k3, v3 = self.calculate_chi3(structure_input)
        k4, v4 = self.calculate_chi4(structure_input)
        k5, v5 = self.calculate_chi5(structure_input)

        return {
            'CHI1': dict(zip(k1, v1)),
            'CHI2': dict(zip(k2, v2)),
            'CHI3': dict(zip(k3, v3)),
            'CHI4': dict(zip(k4, v4)),
            'CHI5': dict(zip(k5, v5))
        }
