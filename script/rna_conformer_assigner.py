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

# Mapping for NtC cluster IDs to standard NtC codes and CANA strings (from RNP_NtC_mean.csv)
NTC_MAP = {
    0: ("AA00", "AAA"),
    1: ("AB05", "A-B"),
    2: ("AA01", "AAw"),
    3: ("OP11", "OPN"),
    4: ("OP22", "OPN"),
    5: ("BB16", "BBw"),
    6: ("OP33", "OPN"),
    7: ("BB04", "B12"),
    8: ("OP15", "OPN"),
    9: ("AA10", "AAw"),
    10: ("AA07", "AAu"),
    11: ("OP21", "OPN"),
    12: ("ZZ1S", "ZZZ"),
    13: ("AB06", "A-B"),
    14: ("OP25", "OPN"),
    15: ("BB21", "miB"),
    16: ("IC02", "ICL"),
    17: ("OP32", "OPN"),
    18: ("OP17", "OPN"),
    19: ("IC08", "ICL"),
    20: ("OP29", "OPN"),
    21: ("AA11", "AAw"),
    22: ("OP26", "OPN"),
    23: ("BB22", "miB"),
    24: ("BA16", "B-A"),
    25: ("IC09", "ICL"),
    26: ("AB07", "A-B"),
    27: ("ZZ01", "ZZZ")
}


def normalize_torsion_angle(val, angle_name=""):
    """
    Normalizes observed torsion angles from various conventions into [0, 360) format.
    Handles signed/negative angles in [-180, 180] for delta, beta, and gamma.
    """
    if pd.isna(val) or val is None:
        return None
    v = float(val)
    # Delta (sugar pucker) is always positive (typically 80-150 deg)
    if "delta" in angle_name.lower() or angle_name.lower() in ["delm", "delta"]:
        return abs(v)
    # Beta and Gamma signed angles (-180..0)
    if "beta" in angle_name.lower() or "gamma" in angle_name.lower():
        if v < 0:
            return abs(v)
    if v < 0:
        return (v % 360.0)
    return v


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
                default_ntc, default_cana = NTC_MAP.get(cluster_int, (f"AA{cluster_int:02d}", "AAA"))
                ntc_code = str(best_meta.get("ntc", default_ntc))
                cana_str = str(best_meta.get("cana", default_cana))
            except ValueError:
                ntc_code = str(best_meta.get("ntc", best_cluster))
                cana_str = str(best_meta.get("cana", "AAA"))
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
                'delta1': normalize_torsion_angle(row[delm_c], 'delta1') if delm_c in row else None,
                'epsilon1': normalize_torsion_angle(row[epsm_c], 'epsilon1') if epsm_c in row else None,
                'zeta1': normalize_torsion_angle(row[zetm_c], 'zeta1') if zetm_c in row else None,
                'alpha2': normalize_torsion_angle(row[alpha_c], 'alpha2') if alpha_c in row else None,
                'beta2': normalize_torsion_angle(row[beta_c], 'beta2') if beta_c in row else None,
                'gamma2': normalize_torsion_angle(row[gamma_c], 'gamma2') if gamma_c in row else None,
                'delta2': normalize_torsion_angle(row[delta_c], 'delta2') if delta_c in row else None,
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
