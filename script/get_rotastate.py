"""
get_rotastate.py
----------------
Module for calculating and assigning amino acid rotamer states (e.g. t-17, M145, m92, no_chi, MA)
based on chi angles, phi/psi dihedrals (BBD library), secondary structure (SS), bound/unbound state, and interface status.
Uses building blocks from build_rotalib.py.
"""

from pathlib import Path
import sys
import numpy as np
import json

# Path setup for script directory
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

ROOT_DIR = SCRIPT_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

import build_rotalib

aa_chis_dict = build_rotalib.aa_chis_dict

# Load SSD rotamer libraries
(
    full_rota_lib_dict_b_I,
    full_rota_lib_dict_b_N,
    full_rota_lib_dict_ub_I,
    full_rota_lib_dict_ub_N,
    rota_lib_dict_b_I,
    rota_lib_dict_b_N,
    rota_lib_dict_ub_I,
    rota_lib_dict_ub_N
) = build_rotalib.run_all_files()

# Load BBD rotamer libraries
(
    bbd_b_I,
    bbd_b_N,
    bbd_ub_I,
    bbd_ub_N
) = build_rotalib.run_bbd_files()


def wrap_180(angle):
    """
    Wraps angle in degrees to [-180, 180].
    """
    return angle - np.ceil(angle / 360.0 - 0.5) * 360.0


def angles_in_valid_range(aa, chis, rotamer, rotastate, rotastate_values, width=30):
    """
    Checks if given chi angles fall within width window of library rotamer state values.
    """
    temp_possible_states_values = []

    for each_rotastate_value in rotastate_values:
        if all(
            np.abs(wrap_180(float(chis[each_chi_index]) - float(each_rotastate_value[each_chi_index]))) < width
            for each_chi_index in range(len(each_rotastate_value) - 1)
        ):
            temp_possible_states_values.append(each_rotastate_value)

    if len(temp_possible_states_values) == 0:
        return [rotamer, 0]
    else:
        sorted_temp_possible_states_values = sorted(temp_possible_states_values, key=lambda x: x[-1])
        assigned_val = round(sorted_temp_possible_states_values[0][-2])
        prob_val = sorted_temp_possible_states_values[0][-1]
        return [rotastate + str(assigned_val), prob_val]


def find_nearest_valid_conf(aa, target_phi, target_psi, lib_dict):
    """
    Finds nearest (PHI, PSI) conformation grid key for amino acid in BBD dictionary.
    """
    if aa not in lib_dict:
        return None
    valid_confs = [c for c, val in lib_dict[aa].items() if len(val) > 0]
    if not valid_confs:
        return None

    best_conf = None
    min_dist = float('inf')

    for conf in valid_confs:
        try:
            phi_str, psi_str = conf.split(',')
            c_phi, c_psi = float(phi_str), float(psi_str)
            d_phi = wrap_180(target_phi - c_phi)
            d_psi = wrap_180(target_psi - c_psi)
            dist = d_phi**2 + d_psi**2
            if dist < min_dist:
                min_dist = dist
                best_conf = conf
        except ValueError:
            continue

    return best_conf


class RotamerStateAssigner:
    """
    Object-Oriented Rotamer State Assigner for protein amino acid residues.
    """

    def __init__(self):
        self.aa_chis_dict = aa_chis_dict
        self.full_b_I = full_rota_lib_dict_b_I
        self.full_b_N = full_rota_lib_dict_b_N
        self.full_ub_I = full_rota_lib_dict_ub_I
        self.full_ub_N = full_rota_lib_dict_ub_N

        self.bbd_b_I = bbd_b_I
        self.bbd_b_N = bbd_b_N
        self.bbd_ub_I = bbd_ub_I
        self.bbd_ub_N = bbd_ub_N

    def get_rotamer(self, aa, ss, b_ub, interface, chi_1=None, chi_2=None, chi_3=None, chi_4=None, chi_5=None, phi=None, psi=None):
        return get_aa_rotamer(aa, ss, b_ub, interface, chi_1, chi_2, chi_3, chi_4, chi_5, phi=phi, psi=psi)


def get_aa_rotamer(aa, ss, b_ub, interface, chi_1=None, chi_2=None, chi_3=None, chi_4=None, chi_5=None, phi=None, psi=None):
    """
    Assigns rotamer state for an amino acid residue.
    Uses Backbone Dependent (BBD) library when phi and psi angles are provided,
    and falls back to Secondary Structure Dependent (SSD) library otherwise.
    """
    chis_temp = [chi_1, chi_2, chi_3, chi_4, chi_5]
    chis = [val for val in chis_temp if val is not None and not (isinstance(val, float) and np.isnan(val))]

    aa_chis = aa_chis_dict.get(aa, 0)

    if len(chis) == 0 and aa in ['ALA', 'GLY']:
        return 'no_chi'
    elif len(chis) < aa_chis:
        return 'MA'

    b_ub_clean = str(b_ub).lower().strip()
    interface_clean = str(interface).upper().strip()

    # 1. Try BBD rotamer assignment if phi and psi are present
    if phi is not None and psi is not None and not (np.isnan(phi) or np.isnan(psi)):
        if b_ub_clean in ['bound', 'b', 'complex']:
            bbd_lib = bbd_b_I if interface_clean == 'I' else bbd_b_N
        else:
            bbd_lib = bbd_ub_I if interface_clean == 'I' else bbd_ub_N

        if bbd_lib:
            nearest_conf = find_nearest_valid_conf(aa, float(phi), float(psi), bbd_lib)
            if nearest_conf:
                aa_rota_library = bbd_lib[aa][nearest_conf]
                final_assign_rotastates = []
                for rotastate, rotastate_value in aa_rota_library.items():
                    rotamer = '!' * aa_chis
                    assigned = angles_in_valid_range(aa, chis, rotamer, rotastate, rotastate_value)
                    final_assign_rotastates.append(assigned)

                if final_assign_rotastates:
                    sorted_states = sorted(final_assign_rotastates, key=lambda x: x[-1], reverse=True)
                    if sorted_states[0][1] > 0:
                        return sorted_states[0][0]

    # 2. SSD rotamer assignment fallback
    if b_ub_clean in ['bound', 'b', 'complex']:
        lib_dict = full_rota_lib_dict_b_I if interface_clean == 'I' else full_rota_lib_dict_b_N
    else:
        lib_dict = full_rota_lib_dict_ub_I if interface_clean == 'I' else full_rota_lib_dict_ub_N

    if not lib_dict or aa not in lib_dict:
        return '!' * aa_chis

    aa_rota_library = lib_dict[aa].get(ss, {})
    if not aa_rota_library:
        return '!' * aa_chis

    final_assign_rotastates = []
    for rotastate, rotastate_value in aa_rota_library.items():
        rotamer = '!' * aa_chis
        assigned = angles_in_valid_range(aa, chis, rotamer, rotastate, rotastate_value)
        final_assign_rotastates.append(assigned)

    if not final_assign_rotastates:
        return '!' * aa_chis

    sorted_states = sorted(final_assign_rotastates, key=lambda x: x[-1], reverse=True)
    return sorted_states[0][0]


if __name__ == "__main__":
    test_state = get_aa_rotamer('ASN', 'L', 'bound', 'I', -74.42, 101.18, phi=74.42, psi=-101.18)
    print("Test BBD Rotamer Assignment (ASN, L, bound, I, chi1=-74.42, chi2=101.18, phi=74.42, psi=-101.18):", test_state)
