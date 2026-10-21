"""
get_rotastate.py
----------------
Module for calculating and assigning amino acid rotamer states for all 20 amino acids based on
SCWRL4 / MODELLER / PyMOL / Chimera rotamer matching algorithms extended for custom rotamer libraries.
- Supports both ON-rotamers (p, t, m) and OFF-rotamers (P, T, M).
- Uses backbone Ramachandran neighborhood search (30-60 deg radius) to handle sparse grid cells.
- Evaluates rotamer matches using log-likelihood weighted angular deviation scoring over all chis (up to chi4).
- Names rotamers Lovell/Richardson style: letters for rotameric chis, the value for a terminal sp2 chi.
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

# Load BBI rotamer libraries
(
    bbi_b_I,
    bbi_b_N,
    bbi_ub_I,
    bbi_ub_N
) = build_rotalib.run_bbi_files()


# Residues whose last chi is a terminal sp2 torsion: named with that chi's value (Lovell et al. 2000,
# e.g. ASP 'm-20', GLN 'mt-30', ARG 'mtt180'). All other chis are rotameric and named by letter only.
SP2_TERMINAL_AA = {'ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'PHE', 'TYR', 'TRP', 'ARG'}
# Terminal chis with two-fold symmetry (1-indexed): compared with a 180-degree period and named
# after folding into [-90, 90)
SYMMETRIC_CHI_IDX = build_rotalib.SYMMETRIC_CHI_IDX
MAX_COMPARED_CHIS = 4  # ARG chi5 (planar guanidinium) is not part of the rotamer definition


def wrap_180(angle):
    """
    Wraps angle in degrees to [-180, 180].
    """
    return angle - np.ceil(angle / 360.0 - 0.5) * 360.0


def _is_missing(value):
    return value is None or (isinstance(value, (int, float, np.floating)) and np.isnan(float(value)))


def _n_compared_chis(aa):
    return min(aa_chis_dict.get(aa, 0), MAX_COMPARED_CHIS)


def _comparable(aa, idx, angle):
    """Angle used in rotamer names: symmetric terminal chis are folded into [-90, 90)."""
    if SYMMETRIC_CHI_IDX.get(aa) == idx + 1:
        return build_rotalib.fold_symmetric_chi(angle)
    return float(angle)


def _chi_diff(aa, idx, a, b):
    """
    Circular difference a - b of chi `idx` (0-based). Two-fold symmetric terminal chis have a
    180-degree period: PHE chi2 = 89 and -89 differ by 2 degrees, not 178.
    """
    if SYMMETRIC_CHI_IDX.get(aa) == idx + 1:
        return (float(a) - float(b) + 90.0) % 180.0 - 90.0
    return wrap_180(float(a) - float(b))


def _best_rotamer(aa, chis, aa_rota_library):
    """
    Scores every library rotamer against the residue's chis and names the best one.

    score = sum of squared circular chi deviations - 100 * ln(P(rotamer)); all chis of the
    residue (up to chi4) are compared. Library rows are [chi1, ..., chiN, prob].
    """
    n_cmp = _n_compared_chis(aa)
    obs = [float(chis[i]) for i in range(n_cmp)]

    best_rotastate, best_score = 'NA', float('inf')
    for r_prefix, r_vals in aa_rota_library.items():
        for val_row in r_vals:
            if len(val_row) < n_cmp + 1:
                continue
            prob = float(val_row[-1])
            chi_diff_sq = sum(
                _chi_diff(aa, i, obs[i], val_row[i]) ** 2 for i in range(n_cmp)
            )
            score = chi_diff_sq - 100.0 * np.log(max(prob, 1e-4))
            if score < best_score:
                best_score = score
                best_rotastate = _rotamer_name(aa, r_prefix, val_row, n_cmp)
    return best_rotastate


def _rotamer_name(aa, r_prefix, val_row, n_cmp):
    """Lovell/Richardson style name from the library prefix (letters of chi1..chi(n-1)) and row."""
    if n_cmp == 1:
        return r_prefix
    last = val_row[n_cmp - 1]
    if aa in SP2_TERMINAL_AA:
        return f"{r_prefix}{round(_comparable(aa, n_cmp - 1, last))}"
    return f"{r_prefix}{build_rotalib.assign_ptm_state(last)}"


def _residue_chis(aa, chis):
    """chi1..chiN for comparison, or None if any of them is missing (positions are preserved)."""
    n_cmp = _n_compared_chis(aa)
    vals = list(chis[:n_cmp])
    if len(vals) < n_cmp or any(_is_missing(v) for v in vals):
        return None
    return [float(v) for v in vals]


def classify_rotamer(aa, chis):
    """
    Rotamer class read directly from a residue's own chi angles (no library): one letter per
    rotameric chi from build_rotalib.assign_ptm_state (p/t/m on-rotamer, P/T/M off-rotamer).
    For residues with a terminal sp2 chi (SP2_TERMINAL_AA) that chi is not part of the class,
    matching the letter part of the Lovell names produced by the library assignment.
    Returns 'no_chi' for ALA/GLY and 'NA' when a required chi is missing.
    """
    aa_clean = str(aa).upper().strip()
    if aa_chis_dict.get(aa_clean, 0) == 0:
        return 'no_chi'
    res_chis = _residue_chis(aa_clean, chis)
    if res_chis is None:
        return 'NA'
    n_letters = len(res_chis) - 1 if (aa_clean in SP2_TERMINAL_AA and len(res_chis) > 1) else len(res_chis)
    return "".join(build_rotalib.assign_ptm_state(c) for c in res_chis[:n_letters])


def find_best_bbd_rotamer(aa, chis, phi, psi, bbd_lib, max_grid_radius=3600):
    """
    Rotamer search for Backbone-Dependent (BBD) libraries.
    1. Uses the nearest populated Ramachandran grid cell within sqrt(max_grid_radius) degrees.
    2. Scores library rotamers by squared chi deviation and log-probability (see _best_rotamer).
    """
    aa_clean = str(aa).upper().strip()
    if aa_chis_dict.get(aa_clean, 0) == 0:
        return 'no_chi'

    res_chis = _residue_chis(aa_clean, chis)
    if res_chis is None or not bbd_lib or aa_clean not in bbd_lib:
        return 'NA'

    best_conf, min_grid_dist = None, float('inf')
    for conf, rotas in bbd_lib[aa_clean].items():
        if not rotas:
            continue
        try:
            c_phi, c_psi = (float(x) for x in conf.split(','))
        except ValueError:
            continue
        dist = wrap_180(phi - c_phi) ** 2 + wrap_180(psi - c_psi) ** 2
        if dist < min_grid_dist and dist <= max_grid_radius:
            min_grid_dist, best_conf = dist, conf

    if best_conf is None:
        return 'NA'
    return _best_rotamer(aa_clean, res_chis, bbd_lib[aa_clean][best_conf])


def find_best_ssd_rotamer(aa, chis, ss, ssd_lib):
    """
    Rotamer search for Secondary Structure-Dependent (SSD) and Backbone-Independent libraries.
    """
    aa_clean = str(aa).upper().strip()
    if aa_chis_dict.get(aa_clean, 0) == 0:
        return 'no_chi'

    res_chis = _residue_chis(aa_clean, chis)
    if res_chis is None or not ssd_lib or aa_clean not in ssd_lib:
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
    return _best_rotamer(aa_clean, res_chis, aa_rota_library)


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
    - All chis up to chi4 are compared; names follow Lovell/Richardson (see _rotamer_name).
    """
    aa_clean = str(aa).upper().strip()
    aa_chis = aa_chis_dict.get(aa_clean, 0)

    if aa_clean in ['ALA', 'GLY'] or aa_chis == 0:
        return 'no_chi'

    chis = [chi_1, chi_2, chi_3, chi_4, chi_5]
    if _residue_chis(aa_clean, chis) is None:
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
        if b_ub_clean in ['bound', 'b', 'complex']:
            lib_dict = bbi_b_I if interface_clean == 'I' else bbi_b_N
        else:
            lib_dict = bbi_ub_I if interface_clean == 'I' else bbi_ub_N

        if lib_dict and aa_clean in lib_dict:
            return find_best_ssd_rotamer(aa_clean, chis, 'ALL', lib_dict)
        return 'NA'

    return 'NA'


if __name__ == "__main__":
    test_state = get_aa_rotamer('ASN', 'L', 'bound', 'I', -74.42, 101.18, phi=74.42, psi=-101.18, lib_type='BBD')
    print("Test BBD Rotamer Assignment (ASN, L, bound, I, chi1=-74.42, chi2=101.18, phi=74.42, psi=-101.18):", test_state)
