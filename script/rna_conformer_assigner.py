"""
rna_conformer_assigner.py
-------------------------
RNA Dinucleotide Conformer Assigner module using StructAnnotate.
Calculates NtC conformer cluster ID, CANA cluster string, and confalA score
from RNA 7D suite torsion angle dataframes or CSV files.
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np

# Path setup for StructAnnotate
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
SA_DIR = ROOT_DIR / "StructAnnotate"

if str(SA_DIR) not in sys.path:
    sys.path.insert(0, str(SA_DIR))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    from rna_annotator.confal import assign_conformer, prepare_conformer_library
except ImportError:
    raise ImportError(f"Cannot import StructAnnotate.rna_annotator.confal. Ensure StructAnnotate is located at {SA_DIR}")

# Mapping for NtC cluster IDs to standard NtC codes and CANA strings
NTC_MAP = {
    0: ("AA00", "AAA"),
    1: ("AA01", "AAw"),
    2: ("AA02", "AAA"),
    3: ("AA03", "AAA"),
    4: ("AA04", "AAA"),
    5: ("AA05", "AAA"),
    6: ("AA06", "AAw"),
    7: ("AA07", "AAA"),
    8: ("AA08", "AAA"),
    9: ("AA09", "AAA"),
    10: ("AA10", "AAw"),
    11: ("AA11", "AAw"),
    12: ("AA12", "AAu"),
    13: ("AB01", "A-B"),
    14: ("AB02", "A-B"),
    15: ("AB03", "A-B"),
    16: ("AB04", "A-B"),
    17: ("AB05", "A-B"),
    18: ("BA01", "B-A"),
    19: ("BA02", "B-A"),
    20: ("BB01", "BBB"),
    21: ("OP01", "OPN"),
    22: ("OP02", "OPN"),
    23: ("OP03", "OPN"),
    24: ("OP04", "OPN"),
    25: ("OP05", "OPN"),
    26: ("OP06", "OPN"),
    27: ("OP07", "OPN")
}


class RNAConformerAssigner:
    """
    Assigner class for computing RNA dinucleotide conformer states (NtC, CANA, confalA).
    """

    def __init__(self, library_csv=None):
        if library_csv:
            lib_p = Path(library_csv)
        else:
            lib_p = SA_DIR / "library" / "conformer_library.csv"

        if not lib_p.exists():
            raise FileNotFoundError(f"StructAnnotate conformer library CSV not found: {lib_p}")

        self.library_path = lib_p
        self.library = prepare_conformer_library(lib_p)

    def assign_suite(self, observed_angles):
        """
        Assigns NtC, CANA, and confalA score for a single dictionary of observed torsion angles.
        :param observed_angles: dict with keys matching ['delta1', 'epsilon1', 'zeta1', 'alpha2', 'beta2', 'gamma2', 'delta2']
        :returns: (ntc_code, cana_str, confal_score)
        """
        best_cluster, best_score, best_meta = assign_conformer(observed_angles, self.library)

        if best_cluster is not None and isinstance(best_cluster, (int, str)):
            try:
                cluster_int = int(best_cluster)
                ntc_code, cana_str = NTC_MAP.get(cluster_int, (f"AA{cluster_int:02d}", "AAA"))
            except ValueError:
                ntc_code = str(best_cluster)
                cana_str = best_meta.get("cana", "AAA")
        else:
            ntc_code = "NANT"
            cana_str = "NAN"
            best_score = 0.0

        return ntc_code, cana_str, round(float(best_score), 1)

    def assign_df(self, df):
        """
        Assigns NtC, CANA, and confalA columns to an RNA 7D suite DataFrame.
        """
        df = df.copy()
        cols = {c.upper(): c for c in df.columns}

        delm_c = cols.get('DELM', 'DELM')
        epsm_c = cols.get('EPSM', 'EPSM')
        zetm_c = cols.get('ZETM', 'ZETM')
        alpha_c = cols.get('ALPHA', 'ALPHA')
        beta_c = cols.get('BETA', 'BETA')
        gamma_c = cols.get('GAMMA', 'GAMMA')
        delta_c = cols.get('DELTA', 'DELTA')

        ntc_list = []
        cana_list = []
        confal_list = []

        for idx, row in df.iterrows():
            obs = {
                'delta1': float(row[delm_c]) % 360.0 if delm_c in row and not pd.isna(row[delm_c]) else None,
                'epsilon1': float(row[epsm_c]) % 360.0 if epsm_c in row and not pd.isna(row[epsm_c]) else None,
                'zeta1': float(row[zetm_c]) % 360.0 if zetm_c in row and not pd.isna(row[zetm_c]) else None,
                'alpha2': float(row[alpha_c]) % 360.0 if alpha_c in row and not pd.isna(row[alpha_c]) else None,
                'beta2': float(row[beta_c]) % 360.0 if beta_c in row and not pd.isna(row[beta_c]) else None,
                'gamma2': float(row[gamma_c]) % 360.0 if gamma_c in row and not pd.isna(row[gamma_c]) else None,
                'delta2': float(row[delta_c]) % 360.0 if delta_c in row and not pd.isna(row[delta_c]) else None,
            }

            ntc, cana, confal = self.assign_suite(obs)
            ntc_list.append(ntc)
            cana_list.append(cana)
            confal_list.append(confal)

        df['NtC'] = ntc_list
        df['CANA'] = cana_list
        df['confalA'] = confal_list

        return df

    def assign_file(self, input_file, output_file=None):
        """
        Processes an input CSV/TSV RNA 7D suite file and appends NtC, CANA, confalA columns.
        """
        in_p = Path(input_file)
        if not in_p.exists():
            raise FileNotFoundError(f"Input RNA file not found: {in_p}")

        sep = ',' if in_p.suffix.lower() == '.csv' else '\t'
        df = pd.read_csv(in_p, sep=sep, engine='python')

        df_out = self.assign_df(df)

        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            out_sep = ',' if out_p.suffix.lower() == '.csv' else '\t'
            df_out.to_csv(out_p, sep=out_sep, index=False)

        return df_out


def assign_rna_conformers(input_target, output_file=None, library_csv=None):
    """
    Helper function for assigning RNA dinucleotide conformers to DataFrame or file.
    """
    assigner = RNAConformerAssigner(library_csv=library_csv)
    if isinstance(input_target, (str, Path)):
        return assigner.assign_file(input_target, output_file=output_file)
    elif isinstance(input_target, pd.DataFrame):
        df_out = assigner.assign_df(input_target)
        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            df_out.to_csv(out_p, index=False)
        return df_out
    else:
        raise TypeError("input_target must be a file path or pandas DataFrame")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Assign RNA dinucleotide conformers (NtC, CANA, confalA) using StructAnnotate.")
    parser.add_argument("input", help="Input RNA 7D suite CSV/TSV file")
    parser.add_argument("-o", "--output", default=None, help="Output file path")

    args = parser.parse_args()
    df_res = assign_rna_conformers(args.input, output_file=args.output)
    print(f"Assigned RNA conformers shape: {df_res.shape}")
