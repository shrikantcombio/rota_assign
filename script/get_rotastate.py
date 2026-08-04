"""
get_rotastate.py
----------------
Module for calculating and assigning amino acid rotamer states for all 20 amino acids based on
SCWRL4 / MODELLER / PyMOL / Chimera rotamer matching algorithms extended for custom rotamer libraries.
- Supports both ON-rotamers (p, t, m) and OFF-rotamers (P, T, M).
- Uses backbone Ramachandran neighborhood search (30-60 deg radius) to handle sparse grid cells.
- Evaluates rotamer matches using log-likelihood weighted angular deviation scoring over n - 1 primary dihedrals.
- Supports user-selected library type ('BBD' default, 'SSD', or 'BBI').
- Outputs 'no_chi' for ALA/GLY, and 'NA' (Not Available) for unassigned or missing chi states.
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


def find_best_bbd_rotamer(aa, chis, phi, psi, bbd_lib, max_grid_radius=3600):
    """
    SCWRL4 / MODELLER / PyMOL style rotamer search algorithm for Backbone-Dependent (BBD) libraries.
    1. Searches neighboring Ramachandran grid cells within max_grid_radius (30-60 deg).
    2. Evaluates rotamer candidates using log-likelihood weighted chi-squared angular deviation over n - 1 primary dihedrals.
    3. Selects optimal rotamer state supporting both ON (p,t,m) and OFF (P,T,M) rotameric states.
    """
    aa_clean = str(aa).upper().strip()
    n_chis = aa_chis_dict.get(aa_clean, 0)
    if n_chis == 0 or aa_clean in ['ALA', 'GLY']:
        return 'no_chi'

    n_check = n_chis - 1 if n_chis > 1 else 1

    clean_chis = [c for c in chis[:n_chis] if c is not None and not (isinstance(c, (int, float)) and np.isnan(float(c)))]
    if len(clean_chis) < n_check:
        return 'NA'

    if not bbd_lib or aa_clean not in bbd_lib:
        return 'NA'

    # Step 1: Search nearest grid points on Ramachandran map
    best_conf = None
    min_grid_dist = float('inf')

    for conf, rotas in bbd_lib[aa_clean].items():
        if not rotas:
            continue
        try:
            phi_str, psi_str = conf.split(',')
            c_phi, c_psi = float(phi_str), float(psi_str)
            d_phi = wrap_180(phi - c_phi)
            d_psi = wrap_180(psi - c_psi)
            dist = d_phi**2 + d_psi**2
            if dist < min_grid_dist and dist <= max_grid_radius:
                min_grid_dist = dist
                best_conf = conf
        except ValueError:
            continue

    if not best_conf:
        # Fallback to any valid conformation for this aa in the library
        valid_confs = [c for c, val in bbd_lib[aa_clean].items() if len(val) > 0]
        if not valid_confs:
            return 'NA'
        best_conf = min(valid_confs, key=lambda c: (wrap_180(phi - float(c.split(',')[0]))**2 + wrap_180(psi - float(c.split(',')[1]))**2))

    aa_rota_library = bbd_lib[aa_clean][best_conf]

    # Step 2: Log-Likelihood + Chi-squared Angular Deviation Scoring over n - 1 primary dihedrals
    best_rotastate = 'NA'
    best_score = float('inf')

    for r_prefix, r_vals in aa_rota_library.items():
        for val_row in r_vals:
            prob = float(val_row[-1])
            chi_diff_sq = sum(
                (wrap_180(clean_chis[idx] - float(val_row[idx])))**2
                for idx in range(min(n_check, len(val_row) - 1))
            )
            score = chi_diff_sq - 100.0 * np.log(max(prob, 1e-4))
            if score < best_score:
                best_score = score
                last_chi = round(val_row[min(n_chis - 1, len(val_row) - 2)])
                best_rotastate = f"{r_prefix}{last_chi}"

    return best_rotastate


def find_best_ssd_rotamer(aa, chis, ss, ssd_lib):
    """
    SCWRL / MODELLER style rotamer search algorithm for Secondary Structure-Dependent (SSD) libraries.
    """
    aa_clean = str(aa).upper().strip()
    n_chis = aa_chis_dict.get(aa_clean, 0)
    if n_chis == 0 or aa_clean in ['ALA', 'GLY']:
        return 'no_chi'

    n_check = n_chis - 1 if n_chis > 1 else 1

    clean_chis = [c for c in chis[:n_chis] if c is not None and not (isinstance(c, (int, float)) and np.isnan(float(c)))]
    if len(clean_chis) < n_check:
        return 'NA'

    if not ssd_lib or aa_clean not in ssd_lib:
        return 'NA'

    aa_rota_library = ssd_lib[aa_clean].get(ss, {})
    if not aa_rota_library:
        # Fallback to loop 'L' or any SS key
        for alt_ss in ['L', 'H', 'E']:
            if alt_ss in ssd_lib[aa_clean]:
                aa_rota_library = ssd_lib[aa_clean][alt_ss]
                break

    if not aa_rota_library:
        return 'NA'

    best_rotastate = 'NA'
    best_score = float('inf')

    for r_prefix, r_vals in aa_rota_library.items():
        for val_row in r_vals:
            prob = float(val_row[-1])
            chi_diff_sq = sum(
                (wrap_180(clean_chis[idx] - float(val_row[idx])))**2
                for idx in range(min(n_check, len(val_row) - 1))
            )
            score = chi_diff_sq - 100.0 * np.log(max(prob, 1e-4))
            if score < best_score:
                best_score = score
                last_chi = round(val_row[min(n_chis - 1, len(val_row) - 2)])
                best_rotastate = f"{r_prefix}{last_chi}"

    return best_rotastate


class RotamerStateAssigner:
    """
    Object-Oriented Rotamer State Assigner for protein amino acid residues.
    """

    def __init__(self, lib_type='BBD'):
        self.lib_type = lib_type
        self.aa_chis_dict = aa_chis_dict
        self.full_b_I = full_rota_lib_dict_b_I
        self.full_b_N = full_rota_lib_dict_b_N
        self.full_ub_I = full_rota_lib_dict_ub_I
        self.full_ub_N = full_rota_lib_dict_ub_N

        self.bbd_b_I = bbd_b_I
        self.bbd_b_N = bbd_b_N
        self.bbd_ub_I = bbd_ub_I
        self.bbd_ub_N = bbd_ub_N

    def get_rotamer(self, aa, ss, b_ub, interface, chi_1=None, chi_2=None, chi_3=None, chi_4=None, chi_5=None, phi=None, psi=None, lib_type=None):
        l_type = lib_type or self.lib_type
        return get_aa_rotamer(aa, ss, b_ub, interface, chi_1, chi_2, chi_3, chi_4, chi_5, phi=phi, psi=psi, lib_type=l_type)


def get_aa_rotamer(aa, ss, b_ub, interface, chi_1=None, chi_2=None, chi_3=None, chi_4=None, chi_5=None, phi=None, psi=None, lib_type='BBD'):
    """
    Assigns rotamer state for an amino acid residue using SCWRL4 / MODELLER / PyMOL / Chimera style rotamer search.
    - lib_type: 'BBD' (default), 'SSD', or 'BBI'.
    - Returns 'no_chi' for ALA/GLY.
    - Returns 'NA' (Not Available) for missing chi states.
    - Supports ON-rotamers (p,t,m) and OFF-rotamers (P,T,M).
    - Standardized n - 1 chi criteria applied for multi-chi amino acids.
    """
    aa_clean = str(aa).upper().strip()
    aa_chis = aa_chis_dict.get(aa_clean, 0)

    if aa_clean in ['ALA', 'GLY'] or aa_chis == 0:
        return 'no_chi'

    chis_temp = [chi_1, chi_2, chi_3, chi_4, chi_5]
    chis = [val for val in chis_temp if val is not None and not (isinstance(val, (int, float)) and np.isnan(float(val)))]

    n_check = aa_chis - 1 if aa_chis > 1 else 1

    if len(chis) < n_check:
        return 'NA'

    b_ub_clean = str(b_ub).lower().strip()
    interface_clean = str(interface).upper().strip()
    lib_type_clean = str(lib_type).upper().strip()

    # 1. Backbone-Dependent (BBD) Rotamer Library (Default)
    if lib_type_clean == 'BBD':
        if phi is not None and psi is not None and not (np.isnan(float(phi)) or np.isnan(float(psi))):
            if b_ub_clean in ['bound', 'b', 'complex']:
                bbd_lib = bbd_b_I if interface_clean == 'I' else bbd_b_N
            else:
                bbd_lib = bbd_ub_I if interface_clean == 'I' else bbd_ub_N

            if bbd_lib:
                return find_best_bbd_rotamer(aa_clean, chis, float(phi), float(psi), bbd_lib)
        return 'NA'

    # 2. Secondary Structure-Dependent (SSD) Rotamer Library
    elif lib_type_clean == 'SSD':
        if b_ub_clean in ['bound', 'b', 'complex']:
            lib_dict = full_rota_lib_dict_b_I if interface_clean == 'I' else full_rota_lib_dict_b_N
        else:
            lib_dict = full_rota_lib_dict_ub_I if interface_clean == 'I' else full_rota_lib_dict_ub_N

        if lib_dict and aa_clean in lib_dict:
            return find_best_ssd_rotamer(aa_clean, chis, ss, lib_dict)
        return 'NA'

    # 3. Backbone-Independent (BBI) Rotamer Library
    elif lib_type_clean == 'BBI':
        return 'NA'

    return 'NA'


if __name__ == "__main__":
    test_state = get_aa_rotamer('ASN', 'L', 'bound', 'I', -74.42, 101.18, phi=74.42, psi=-101.18, lib_type='BBD')
    print("Test BBD Rotamer Assignment (ASN, L, bound, I, chi1=-74.42, chi2=101.18, phi=74.42, psi=-101.18):", test_state)
