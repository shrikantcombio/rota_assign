"""
Rotamer Library Builder Module
Object-Oriented implementation for building rotamer library dictionaries
from Backbone Dependent (BBD), Backbone Independent (BBI), and Secondary Structure Dependent (SSD) libraries.
Loads existing json dictionaries from rota_lib_dict if available.
"""

from pathlib import Path
import os
import sys
import json
import copy
import pandas as pd
import numpy as np


AA_CHIS_DICT = {
    'ARG': 5, 'LYS': 4, 'GLU': 3, 'GLN': 3, 'MET': 3,
    'ASP': 2, 'ASN': 2, 'HIS': 2, 'LEU': 2, 'ILE': 2,
    'PHE': 2, 'PRO': 2, 'TRP': 2, 'TYR': 2, 'CYS': 1,
    'SER': 1, 'THR': 1, 'VAL': 1, 'ALA': 0, 'GLY': 0
}

aa_chis_dict = AA_CHIS_DICT

AA_LIST = [
    "ARG", "ASN", "ASP", "GLU", "GLN", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "TRP", "TYR", "CYS", "SER", "THR", "VAL", "GLY", "ALA"
]

aaList = AA_LIST


def assign_ptm_state(val):
    val = float(val)
    if -30 < val <= 30:
        return 'T'
    elif 30 < val <= 90:
        return 'p'
    elif 90 < val <= 150:
        return 'P'
    elif val > 150 or val <= -150:
        return 't'
    elif -150 < val <= -90:
        return 'M'
    elif -90 < val <= -30:
        return 'm'
    else:
        return str(int(val))


SYMMETRIC_CHI_IDX = {'ASP': 2, 'GLU': 3, 'PHE': 2, 'TYR': 2}  # 1-indexed, two-fold symmetric (period 180)

def fold_symmetric_chi(angle):
    angle = float(angle)
    if angle >= 90.0:
        return angle - 180.0
    if angle < -90.0:
        return angle + 180.0
    return angle

def proline_pucker(chi1):
    """
    Proline ring pucker from chi1: 'endo' (C-gamma endo, chi1 > 0, chi2 < 0) or 'exo' (C-gamma exo,
    chi1 < 0, chi2 > 0). The pyrrolidine ring has no staggered chi minima, so the 6-state letters
    (which would call both puckers the off-rotamer 'T') do not apply.
    """
    return 'endo' if float(chi1) > 0 else 'exo'


def chi_state(row):
    aa = row['AA']
    if aa == 'PRO':
        return proline_pucker(row['Chi1'])
    n_chis = AA_CHIS_DICT.get(aa, 0)
    sym_idx = SYMMETRIC_CHI_IDX.get(aa, None)

    def get_chi(i):  # i is 1-indexed
        val = float(row[f'Chi{i}'])
        if sym_idx and i == sym_idx:
            val = fold_symmetric_chi(val)
        return assign_ptm_state(val)

    if n_chis == 0:
        return ""
    elif n_chis == 1:
        return get_chi(1)
    elif n_chis == 2:
        return get_chi(1) + get_chi(2)
    elif n_chis == 3:
        return get_chi(1) + get_chi(2) + get_chi(3)
    else:
        return get_chi(1) + get_chi(2) + get_chi(3) + get_chi(4)


def is_on_rotamer(aa, rota_state):
    if not rota_state:
        return True
    n_chis = AA_CHIS_DICT.get(aa, 0)
    if n_chis != 1:
        prefix = rota_state[:2]
        return not any(c in prefix for c in ['P', 'T', 'M'])
    else:
        return not any(c in rota_state for c in ['P', 'T', 'M'])


def _read_clean_lines(file_path):
    with open(file_path, 'r') as file:
        lines = [line.strip() for line in file if line.strip() and not line.strip().strip('\"\'').startswith('#')]
    if not lines:
        raise ValueError(f"No valid data lines found in {file_path}")

    first_line = lines[0]
    sep = ',' if ',' in first_line else ('\t' if '\t' in first_line else None)
    data = [line.split(sep) for line in lines]

    df_raw = pd.DataFrame(data)
    df_raw = df_raw.map(lambda s: s.strip() if isinstance(s, str) else s)

    if str(df_raw.iloc[0, 0]).upper() in ['RES', 'AA', 'RESIDUE']:
        df_raw = df_raw.iloc[1:].reset_index(drop=True)

    return df_raw


class BaseRotamerLibBuilder:
    def __init__(self, aa_chis_dict=None):
        self.aa_chis_dict = aa_chis_dict or copy.deepcopy(AA_CHIS_DICT)

    def load_data(self, file_path):
        raise NotImplementedError("Subclasses must implement load_data()")

    def process_df(self, df):
        df = df.copy()
        df['rota_state'] = df.apply(chi_state, axis=1)

        is_not_one_chi = df['AA'].map(lambda aa: self.aa_chis_dict.get(aa, 0) != 1)

        mask_on_multi = ~df['rota_state'].str[:2].str.contains(r'[A-Z]', na=False)
        mask_on_single = ~df['rota_state'].str.contains(r'[A-Z]', na=False)

        mask_on = np.where(is_not_one_chi, mask_on_multi, mask_on_single)

        on_rota_state_df = df[mask_on].copy()
        off_rota_state_df = df[~mask_on].copy()

        return df, on_rota_state_df, off_rota_state_df, []

    def build_lib_dict(self, df):
        if df.empty:
            return {}

        final_rota_lib = {}
        all_aa = list(df['AA'].unique())

        for each_aa in all_aa:
            final_rota_lib[each_aa] = {}

        for idx in range(len(df)):
            row = df.iloc[idx]
            aa = row['AA']
            conf = row['Conformation']
            r_state = row['rota_state']
            prob = float(row['prob'])

            n_chis = min(self.aa_chis_dict.get(aa, 0), 4)

            if aa != 'PRO' and self.aa_chis_dict.get(aa, 0) != 1 and len(r_state) > 1:
                key_prefix = r_state[:-1]
            else:
                key_prefix = r_state

            chi_vals = []
            for chi_idx in range(1, n_chis + 1):
                chi_vals.append(round(float(row[f'Chi{chi_idx}']), 3))
            chi_vals.append(round(prob, 3))

            if conf not in final_rota_lib[aa]:
                final_rota_lib[aa][conf] = {}

            if key_prefix not in final_rota_lib[aa][conf]:
                final_rota_lib[aa][conf][key_prefix] = []

            final_rota_lib[aa][conf][key_prefix].append(chi_vals)

        return final_rota_lib

    def export_json(self, data_dict, output_path):
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, 'w') as file:
            json.dump(data_dict, file, indent=4)

    def build_library(self, file_path, output_dir=None, name_prefix=None):
        file_p = Path(file_path)
        df = self.load_data(file_p)
        full_df, on_df, off_df, prob_arr = self.process_df(df)

        full_dict = self.build_lib_dict(full_df)
        on_dict = self.build_lib_dict(on_df)
        off_dict = self.build_lib_dict(off_df)

        if output_dir:
            out_d = Path(output_dir)
            prefix = name_prefix or file_p.stem
            self.export_json(full_dict, out_d / f"full_{prefix}.json")
            self.export_json(on_dict, out_d / f"{prefix}.json")

        return full_dict, on_dict, off_dict, prob_arr


class SSDRotamerLibBuilder(BaseRotamerLibBuilder):
    def load_data(self, file_path):
        file_p = Path(file_path)
        if file_p.suffix.lower() == '.json':
            df_raw = pd.read_json(file_p)
            df = pd.DataFrame()
            df['AA'] = df_raw['AA'].astype(str).str.strip()
            df['Conformation'] = df_raw['Conformation'] if 'Conformation' in df_raw.columns else (df_raw['SS'] if 'SS' in df_raw.columns else 'L')
            df['prob'] = df_raw['PROB'] if 'PROB' in df_raw.columns else df_raw['prob']
            for i in range(1, 5):
                chi_col = next((c for c in [f'B_CHI{i}', f'U_CHI{i}', f'CHI{i}', f'Chi{i}'] if c in df_raw.columns), None)
                df[f'Chi{i}'] = df_raw[chi_col] if chi_col else 0.0
            return df
        else:
            df_raw = _read_clean_lines(file_p)
            df = pd.DataFrame()
            df['AA'] = df_raw.iloc[:, 0].astype(str).str.strip()
            df['Conformation'] = df_raw.iloc[:, 1].astype(str).str.strip()
            df['Chi1'] = df_raw.iloc[:, 2]
            df['Chi2'] = df_raw.iloc[:, 3]
            df['Chi3'] = df_raw.iloc[:, 4]
            df['Chi4'] = df_raw.iloc[:, 5]
            df['prob'] = df_raw.iloc[:, -1]
            return df


class BBDRotamerLibBuilder(BaseRotamerLibBuilder):
    def load_data(self, file_path):
        file_p = Path(file_path)
        if file_p.suffix.lower() == '.json':
            df_raw = pd.read_json(file_p)
            df = pd.DataFrame()
            df['AA'] = df_raw['AA'].astype(str).str.strip()
            if 'PHI' in df_raw.columns and 'PSI' in df_raw.columns:
                phi = df_raw['PHI'].apply(lambda x: str(int(round(float(x)))))
                psi = df_raw['PSI'].apply(lambda x: str(int(round(float(x)))))
                df['Conformation'] = phi + ',' + psi
            elif 'Conformation' in df_raw.columns:
                df['Conformation'] = df_raw['Conformation'].astype(str).str.strip()
            else:
                df['Conformation'] = '0,0'

            df['prob'] = df_raw['PROB'] if 'PROB' in df_raw.columns else df_raw['prob']
            for i in range(1, 5):
                chi_col = next((c for c in [f'B_CHI{i}', f'U_CHI{i}', f'CHI{i}', f'Chi{i}'] if c in df_raw.columns), None)
                df[f'Chi{i}'] = df_raw[chi_col] if chi_col else 0.0
            return df
        else:
            df_raw = _read_clean_lines(file_p)
            df = pd.DataFrame()
            df['AA'] = df_raw.iloc[:, 0].astype(str).str.strip()
            phi = df_raw.iloc[:, 1].apply(lambda x: str(int(round(float(x)))))
            psi = df_raw.iloc[:, 2].apply(lambda x: str(int(round(float(x)))))
            df['Conformation'] = phi + ',' + psi
            df['prob'] = df_raw.iloc[:, 5]
            df['Chi1'] = df_raw.iloc[:, 6]
            df['Chi2'] = df_raw.iloc[:, 7]
            df['Chi3'] = df_raw.iloc[:, 8]
            df['Chi4'] = df_raw.iloc[:, 9]
            return df


class BBIRotamerLibBuilder(BaseRotamerLibBuilder):
    def load_data(self, file_path):
        file_p = Path(file_path)
        if file_p.suffix.lower() == '.json':
            df_raw = pd.read_json(file_p)
            df = pd.DataFrame()
            df['AA'] = df_raw['AA'].astype(str).str.strip()
            df['Conformation'] = 'ALL'
            df['prob'] = df_raw['PROB'] if 'PROB' in df_raw.columns else df_raw['prob']
            for i in range(1, 5):
                chi_col = next((c for c in [f'B_CHI{i}', f'U_CHI{i}', f'CHI{i}', f'Chi{i}'] if c in df_raw.columns), None)
                df[f'Chi{i}'] = df_raw[chi_col] if chi_col else 0.0
            return df
        else:
            df_raw = _read_clean_lines(file_p)
            df = pd.DataFrame()
            df['AA'] = df_raw.iloc[:, 0].astype(str).str.strip()
            df['Conformation'] = 'ALL'
            df['prob'] = df_raw.iloc[:, 2]
            df['Chi1'] = df_raw.iloc[:, 3]
            df['Chi2'] = df_raw.iloc[:, 4]
            df['Chi3'] = df_raw.iloc[:, 5]
            df['Chi4'] = df_raw.iloc[:, 6]
            return df


class RotamerLibraryPipeline:
    BUILDERS = {
        'SSD': SSDRotamerLibBuilder,
        'BBD': BBDRotamerLibBuilder,
        'BBI': BBIRotamerLibBuilder,
    }

    def __init__(self, output_dir=None):
        self.output_dir = Path(output_dir) if output_dir else None

    def get_builder(self, lib_type):
        lib_type_upper = lib_type.upper()
        if lib_type_upper not in self.BUILDERS:
            raise ValueError(f"Unknown library type: {lib_type}. Must be one of ['SSD', 'BBD', 'BBI']")
        return self.BUILDERS[lib_type_upper]()

    def process_file(self, file_path, lib_type, output_dir=None, name_prefix=None):
        builder = self.get_builder(lib_type)
        out_dir = output_dir or self.output_dir
        return builder.build_library(file_path, output_dir=out_dir, name_prefix=name_prefix)

    def process_directory(self, dir_path, lib_type, output_dir=None):
        dir_p = Path(dir_path)
        out_dir = output_dir or self.output_dir
        builder = self.get_builder(lib_type)

        results = {}
        for file_p in dir_p.glob("*"):
            if file_p.suffix.lower() in ['.csv', '.tsv', '.txt', '.lib']:
                prefix = file_p.stem
                results[prefix] = builder.build_library(file_p, output_dir=out_dir, name_prefix=prefix)
        return results


# Source rotamer libraries shipped with rota_assign (copies of the RBPs_rotamer_lib releases).
# Set ROTA_ASSIGN_LIBRARY_DIR to use another directory with the same layout.
LIBRARY_DIR = Path(os.environ.get("ROTA_ASSIGN_LIBRARY_DIR", Path(__file__).resolve().parent.parent / "library"))

# BBD/BBI library release: 'v2' (default; terminal sp2 chis binned by residue-specific ranges) or
# 'v1' (every chi in three 120-degree bins; used for the original manuscript analyses).
# The SSD library has a single release.
LIBRARY_VERSIONS = ("v1", "v2")
LIBRARY_VERSION = os.environ.get("ROTA_ASSIGN_LIBRARY_VERSION", "v2")
if LIBRARY_VERSION not in LIBRARY_VERSIONS:
    raise ValueError(f"ROTA_ASSIGN_LIBRARY_VERSION must be one of {LIBRARY_VERSIONS}, got {LIBRARY_VERSION!r}")
_VERSION_SUFFIX = "" if LIBRARY_VERSION == "v1" else f"_{LIBRARY_VERSION}"

# v1 caches stay in rota_lib_dict/<type>/, later releases in rota_lib_dict/<version>/<type>/
_CACHE_ROOT = Path(__file__).resolve().parent.parent / "rota_lib_dict"
CACHE_DIR = _CACHE_ROOT if LIBRARY_VERSION == "v1" else _CACHE_ROOT / LIBRARY_VERSION

# (cache subdir, source subdir, builder, {cache name: source file prefix})
LIBRARY_SPECS = {
    "BBD": ("bbd", f"RBPs_BBD_rotamer_lib{_VERSION_SUFFIX}", "BBD", {
        "bbd_b_I": "RBPs_I_bbd_B_rotamer_lib", "bbd_b_N": "RBPs_N_bbd_B_rotamer_lib",
        "bbd_ub_I": "RBPs_I_bbd_U_rotamer_lib", "bbd_ub_N": "RBPs_N_bbd_U_rotamer_lib"}),
    "BBI": ("bbi", f"RBPs_BBI_rotamer_lib{_VERSION_SUFFIX}", "BBI", {
        "bbi_b_I": "RBPs_I_bbi_B_rotamer_lib", "bbi_b_N": "RBPs_N_bbi_B_rotamer_lib",
        "bbi_ub_I": "RBPs_I_bbi_U_rotamer_lib", "bbi_ub_N": "RBPs_N_bbi_U_rotamer_lib"}),
    "SSD": ("ssd", "RBPs_SSD_rotamer_lib", "SSD", {
        "b_I": "InterfaceResidueBound", "b_N": "NonInterfaceResidueBound",
        "ub_I": "InterfaceResidueUnbound", "ub_N": "NonInterfaceResidueUnbound"}),
}


def _source_file(src_dir, prefix):
    for ext in (".json", ".csv", ".tsv"):
        p = src_dir / f"{prefix}{ext}"
        if p.exists():
            return p
    raise FileNotFoundError(f"No rotamer library file for {prefix} in {src_dir}")


# Bump when the dictionary layout built from the sources changes, so caches are rebuilt
# (2: PRO keyed by ring pucker)
CACHE_FORMAT = "2"


def _sources_digest(files):
    import hashlib
    h = hashlib.sha256(CACHE_FORMAT.encode())
    for f in files:
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def load_libraries(lib_type, force_rebuild=False, save_to_cache=True):
    """
    Rotamer library dictionaries for lib_type ('BBD', 'BBI' or 'SSD'), as
    {name: (full_dict, on_rotamer_dict)}.

    The JSON cache in rota_lib_dict/<type>/ is used only while it was built from the current
    source files (SHA-256 recorded in source.sha256); otherwise it is rebuilt.
    """
    cache_sub, src_sub, builder_type, names = LIBRARY_SPECS[lib_type]
    src_dir = LIBRARY_DIR / src_sub
    cache_dir = CACHE_DIR / cache_sub
    sources = {name: _source_file(src_dir, prefix) for name, prefix in names.items()}
    digest = _sources_digest(sources.values())
    manifest = cache_dir / "source.sha256"

    cache_files = {n: (cache_dir / f"full_rota_lib_{n}.json", cache_dir / f"rota_lib_{n}.json") for n in names}
    cached = (not force_rebuild and manifest.exists() and manifest.read_text().strip() == digest
              and all(f.exists() for pair in cache_files.values() for f in pair))
    if cached:
        return {n: tuple(json.loads(f.read_text()) for f in pair) for n, pair in cache_files.items()}

    builder = RotamerLibraryPipeline().get_builder(builder_type)
    libs = {}
    for name, src in sources.items():
        full_dict, on_dict, _, _ = builder.build_library(src)
        libs[name] = (full_dict, on_dict)
    if save_to_cache:
        cache_dir.mkdir(parents=True, exist_ok=True)
        for name, (full_dict, on_dict) in libs.items():
            builder.export_json(full_dict, cache_files[name][0])
            builder.export_json(on_dict, cache_files[name][1])
        manifest.write_text(digest + "\n")
    return libs


def run_all_files(force_rebuild=False, prompt_user=False):
    """SSD libraries: (full b_I, full b_N, full ub_I, full ub_N, on b_I, on b_N, on ub_I, on ub_N)."""
    libs = load_libraries("SSD", force_rebuild=force_rebuild)
    order = ("b_I", "b_N", "ub_I", "ub_N")
    return tuple(libs[n][0] for n in order) + tuple(libs[n][1] for n in order)


def run_bbd_files(force_rebuild=False, prompt_user=False, save_to_cache=True):
    """BBD libraries (full dictionaries): (b_I, b_N, ub_I, ub_N)."""
    libs = load_libraries("BBD", force_rebuild=force_rebuild, save_to_cache=save_to_cache)
    return tuple(libs[n][0] for n in ("bbd_b_I", "bbd_b_N", "bbd_ub_I", "bbd_ub_N"))


def run_bbi_files(force_rebuild=False, prompt_user=False, save_to_cache=True):
    """BBI libraries (full dictionaries): (b_I, b_N, ub_I, ub_N)."""
    libs = load_libraries("BBI", force_rebuild=force_rebuild, save_to_cache=save_to_cache)
    return tuple(libs[n][0] for n in ("bbi_b_I", "bbi_b_N", "bbi_ub_I", "bbi_ub_N"))


if __name__ == '__main__':
    run_bbd_files(force_rebuild=True)
    run_bbi_files(force_rebuild=True)
