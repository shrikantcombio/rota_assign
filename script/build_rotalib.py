"""
Rotamer Library Builder Module
Object-Oriented implementation for building rotamer library dictionaries
from Backbone Dependent (BBD), Backbone Independent (BBI), and Secondary Structure Dependent (SSD) libraries.
Loads existing json dictionaries from rota_lib_dict if available.
"""

from pathlib import Path
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


def chi_state(row):
    aa = row['AA']
    n_chis = AA_CHIS_DICT.get(aa, 0)
    if n_chis == 0:
        return ""
    elif n_chis == 1:
        return assign_ptm_state(row['Chi1'])
    elif n_chis == 2:
        return assign_ptm_state(row['Chi1']) + assign_ptm_state(row['Chi2'])
    elif n_chis == 3:
        return assign_ptm_state(row['Chi1']) + assign_ptm_state(row['Chi2']) + assign_ptm_state(row['Chi3'])
    else:
        return (assign_ptm_state(row['Chi1']) +
                assign_ptm_state(row['Chi2']) +
                assign_ptm_state(row['Chi3']) +
                assign_ptm_state(row['Chi4']))


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

            if self.aa_chis_dict.get(aa, 0) != 1 and len(r_state) > 1:
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
            self.export_json(full_dict, out_d / f"full_{prefix}.txt")
            self.export_json(on_dict, out_d / f"{prefix}.txt")

        return full_dict, on_dict, off_dict, prob_arr


class SSDRotamerLibBuilder(BaseRotamerLibBuilder):
    def load_data(self, file_path):
        file_p = Path(file_path)
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


def run_all_files():
    current_dir = Path(__file__).resolve().parent
    root_dir = current_dir.parent

    dump_dir = root_dir / 'rota_lib_dict'
    if not dump_dir.exists():
        dump_dir = current_dir / 'rota_lib_dict'

    required_json_files = [
        'full_rota_lib_b_I.txt', 'full_rota_lib_b_N.txt',
        'full_rota_lib_ub_I.txt', 'full_rota_lib_ub_N.txt',
        'rota_lib_b_I.txt', 'rota_lib_b_N.txt',
        'rota_lib_ub_I.txt', 'rota_lib_ub_N.txt'
    ]

    if dump_dir.exists() and all((dump_dir / fname).exists() for fname in required_json_files):
        dicts = []
        for fname in required_json_files:
            with open(dump_dir / fname, 'r') as f:
                dicts.append(json.load(f))
        return tuple(dicts)

    return ({}, {}, {}, {}, {}, {}, {}, {})


def run_bbd_files():
    """
    Loads BBD rotamer library dictionaries from json files in rota_lib_dict/bbd directory,
    or builds them directly from library/RBPs_BBD_rotamer_lib/*.csv if json files do not exist.
    """
    current_dir = Path(__file__).resolve().parent
    root_dir = current_dir.parent

    dump_dir = root_dir / 'rota_lib_dict' / 'bbd'

    bbd_files = [
        'full_rota_lib_bbd_b_I.txt',
        'full_rota_lib_bbd_b_N.txt',
        'full_rota_lib_bbd_ub_I.txt',
        'full_rota_lib_bbd_ub_N.txt'
    ]

    if dump_dir.exists() and all((dump_dir / fname).exists() for fname in bbd_files):
        dicts = []
        for fname in bbd_files:
            with open(dump_dir / fname, 'r') as f:
                dicts.append(json.load(f))
        return tuple(dicts)

    # Fallback to building directly from library/RBPs_BBD_rotamer_lib CSVs
    csv_dir = root_dir / 'library' / 'RBPs_BBD_rotamer_lib'
    if csv_dir.exists():
        builder = BBDRotamerLibBuilder()
        b_I, _, _, _ = builder.build_library(csv_dir / 'RBPs_I_bbd_B_rotamer_lib.csv')
        b_N, _, _, _ = builder.build_library(csv_dir / 'RBPs_N_bbd_B_rotamer_lib.csv')
        ub_I, _, _, _ = builder.build_library(csv_dir / 'RBPs_I_bbd_U_rotamer_lib.csv')
        ub_N, _, _, _ = builder.build_library(csv_dir / 'RBPs_N_bbd_U_rotamer_lib.csv')
        return (b_I, b_N, ub_I, ub_N)

    return ({}, {}, {}, {})


if __name__ == '__main__':
    run_all_files()
