"""
analyze_rotastates.py
---------------------
Analysis module for Arginine forks, protein rotamer state transitions upon RNA binding (bound vs unbound),
amino acid chemical class preferences, and cross-pairing preferences between protein rotamers and RNA dinucleotide conformers.
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np
from collections import Counter

# Path setup
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
TEST_DIR = ROOT_DIR / "test"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# Chemical Class Definitions
AMINO_ACID_GROUPS = {
    'CHARGED': ['ARG', 'LYS', 'HIS', 'ASP', 'GLU'],
    'POLAR': ['SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'],
    'HYDROPHOBIC': ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP'],
    'AROMATIC': ['PHE', 'TYR', 'TRP', 'HIS']
}


def get_chemical_group(aa):
    """
    Returns primary chemical class for a 3-letter amino acid code.
    """
    aa_upper = str(aa).upper().strip()
    if aa_upper in AMINO_ACID_GROUPS['CHARGED']:
        return 'CHARGED'
    elif aa_upper in AMINO_ACID_GROUPS['POLAR']:
        return 'POLAR'
    elif aa_upper in AMINO_ACID_GROUPS['HYDROPHOBIC']:
        return 'HYDROPHOBIC'
    else:
        return 'OTHER'


def classify_rna_moiety(nuc_atom):
    """
    Classifies RNA atom into Phosphate, Ribose Sugar, or Nitrogenous Base.
    """
    atom = str(nuc_atom).upper().strip()
    if atom in ['OP1', 'OP2', 'OP3', 'P', "O5'", "C5'"]:
        return 'Phosphate'
    elif atom in ["O2'", "O3'", "O4'", "C1'", "C2'", "C3'", "C4'"]:
        return 'Ribose Sugar'
    elif atom in ['N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'O2', 'O4', 'O6', 'C2', 'C4', 'C5', 'C6', 'C8']:
        return 'Nitrogenous Base'
    else:
        return 'Other'


def _extract_aa(row):
    if 'AA' in row and not pd.isna(row['AA']):
        return str(row['AA']).strip()
    elif 'LABEL' in row and not pd.isna(row['LABEL']):
        parts = str(row['LABEL']).split()
        return parts[0] if len(parts) > 0 else 'GLY'
    else:
        return 'GLY'


class ArginineForkAnalyzer:
    """
    Analyzer for Arginine forks (bidentate / multi-contact hydrogen bonding between Arginine and RNA).
    """

    def analyze(self, df):
        """
        Extracts Arginine fork residues where count(interactions per (PDB_ID, AA_chain, AA_pos)) > 1.
        :param df: Final merged DataFrame or interactions DataFrame
        :returns: (df_arg_forks, summary_dict)
        """
        df_copy = df.copy()
        df_copy['AA_clean'] = df_copy.apply(_extract_aa, axis=1)

        # Filter Arginine residues
        arg_df = df_copy[df_copy['AA_clean'].str.upper() == 'ARG'].reset_index(drop=True)

        if arg_df.empty:
            return pd.DataFrame(), {'total_arg_interactions': 0, 'arg_fork_residues_count': 0}

        # Identify bidentate / fork interactions (count > 1 per residue)
        group_cols = ['PDB_ID', 'AA_chain', 'AA_pos']
        valid_cols = [c for c in group_cols if c in arg_df.columns]

        if not valid_cols:
            valid_cols = ['AA_chain', 'AA_pos']

        counts = arg_df.groupby(valid_cols)[valid_cols[0]].transform('count')
        arg_forks_df = arg_df[counts > 1].copy()

        # Classify RNA target moieties
        if 'Nuc_atom' in arg_forks_df.columns:
            arg_forks_df['RNA_Moiety'] = arg_forks_df['Nuc_atom'].apply(classify_rna_moiety)

        # Summary calculations
        unique_fork_residues = arg_forks_df[valid_cols].drop_duplicates()
        atom_dist = Counter(arg_forks_df['AA_atom']) if 'AA_atom' in arg_forks_df.columns else {}
        rotastate_dist = Counter(arg_forks_df['b_rotastate']) if 'b_rotastate' in arg_forks_df.columns else {}
        moiety_dist = Counter(arg_forks_df['RNA_Moiety']) if 'RNA_Moiety' in arg_forks_df.columns else {}

        summary = {
            'total_arg_interactions': len(arg_df),
            'arg_fork_interactions': len(arg_forks_df),
            'arg_fork_residues_count': len(unique_fork_residues),
            'arg_atom_distribution': dict(atom_dist),
            'rotastate_distribution': dict(rotastate_dist),
            'rna_moiety_distribution': dict(moiety_dist)
        }

        return arg_forks_df, summary


class RotamerTransitionAnalyzer:
    """
    Analyzer for rotamer state shifts upon RNA binding (unbound ub_rotastate -> bound b_rotastate).
    """

    def analyze(self, df):
        """
        Analyzes rotamer transitions and groups by amino acid chemical class.
        :param df: DataFrame containing 'AA', 'b_rotastate', and 'ub_rotastate'
        """
        df_copy = df.copy()

        if 'b_rotastate' not in df_copy.columns or 'ub_rotastate' not in df_copy.columns:
            return pd.DataFrame(), pd.DataFrame()

        df_copy['AA_clean'] = df_copy.apply(_extract_aa, axis=1)

        df_valid = df_copy[
            (~df_copy['b_rotastate'].astype(str).str.contains('!|MA|NA|no_chi', na=False)) &
            (~df_copy['ub_rotastate'].astype(str).str.contains('!|MA|NA|no_chi', na=False))
        ].copy()

        if df_valid.empty:
            return pd.DataFrame(), pd.DataFrame()

        df_valid['Chemical_Group'] = df_valid['AA_clean'].apply(get_chemical_group)
        df_valid['Rotamer_Shift'] = df_valid['b_rotastate'] != df_valid['ub_rotastate']
        df_valid['Transition'] = df_valid['ub_rotastate'] + " -> " + df_valid['b_rotastate']

        group_summary = []
        for grp, grp_df in df_valid.groupby('Chemical_Group'):
            n_tot = len(grp_df)
            n_shifted = grp_df['Rotamer_Shift'].sum()
            n_preserved = n_tot - n_shifted
            shift_rate = (n_shifted / n_tot * 100) if n_tot > 0 else 0.0

            group_summary.append({
                'Chemical_Group': grp,
                'Total_Residues': n_tot,
                'Preserved_Rotamer_Count': n_preserved,
                'Shifted_Rotamer_Count': n_shifted,
                'Induced_Fit_Rate_Pct': round(shift_rate, 2),
                'Pre_Organized_Rate_Pct': round(100.0 - shift_rate, 2)
            })

        df_summary = pd.DataFrame(group_summary)
        return df_valid, df_summary


class RotamerConformerPreferenceAnalyzer:
    """
    Analyzer for cross-preferences between protein rotamer states and RNA dinucleotide conformers (NtC / CANA).
    """

    def analyze(self, df):
        """
        Cross-tabulates protein rotastates with RNA NtC / CANA conformers.
        """
        df_copy = df.copy()

        if 'b_rotastate' not in df_copy.columns or 'NtC' not in df_copy.columns:
            return pd.DataFrame()

        df_copy['AA_clean'] = df_copy.apply(_extract_aa, axis=1)

        df_valid = df_copy[
            (~df_copy['b_rotastate'].astype(str).str.contains('!|MA|NA|no_chi', na=False)) &
            (df_copy['NtC'] != 'NANT')
        ].copy()

        if df_valid.empty:
            return pd.DataFrame()

        ct_ntc = pd.crosstab(df_valid['AA_clean'] + '_' + df_valid['b_rotastate'], df_valid['NtC'], margins=True)
        return ct_ntc


def run_full_rotastate_analysis(input_csv, output_dir=None):
    """
    Executes full rotamer preference and Arginine fork analysis, exporting output CSVs.
    """
    in_p = Path(input_csv)
    out_d = Path(output_dir) if output_dir else in_p.parent
    out_d.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_p)

    # 1. Arginine Forks Analysis
    fork_analyzer = ArginineForkAnalyzer()
    df_forks, fork_summary = fork_analyzer.analyze(df)

    out_fork_file = out_d / f"{in_p.stem.split('_')[0]}_arg_forks_analysis.csv"
    if not df_forks.empty:
        df_forks.to_csv(out_fork_file, index=False)

    # 2. Rotamer Transitions
    trans_analyzer = RotamerTransitionAnalyzer()
    df_trans, df_group_summary = trans_analyzer.analyze(df)

    out_trans_file = out_d / f"{in_p.stem.split('_')[0]}_rotamer_transitions.csv"
    if isinstance(df_group_summary, pd.DataFrame) and not df_group_summary.empty:
        df_group_summary.to_csv(out_trans_file, index=False)

    # 3. Cross Preference Analysis
    pref_analyzer = RotamerConformerPreferenceAnalyzer()
    df_pref = pref_analyzer.analyze(df)

    out_pref_file = out_d / f"{in_p.stem.split('_')[0]}_rotamer_conformer_analysis.csv"
    if isinstance(df_pref, pd.DataFrame) and not df_pref.empty:
        df_pref.to_csv(out_pref_file)

    return {
        'FORK_SUMMARY': fork_summary,
        'FORK_FILE': out_fork_file,
        'TRANSITION_FILE': out_trans_file,
        'PREFERENCE_FILE': out_pref_file
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Analyze Arginine Forks and Rotamer Preference Shifts.")
    parser.add_argument("--input", required=True, help="Input merged CSV file (e.g. test/1ASY_final_output.csv)")
    parser.add_argument("-o", "--output-dir", default=None, help="Output directory for analysis reports")

    args = parser.parse_args()

    results = run_full_rotastate_analysis(args.input, output_dir=args.output_dir)

    print("=== Arginine Forks Analysis Summary ===")
    for k, v in results['FORK_SUMMARY'].items():
        print(f"  {k}: {v}")

    print("\nAnalysis CSV Files Generated:")
    print("  Arginine Forks File:         ", results['FORK_FILE'])
    print("  Rotamer Transitions File:    ", results['TRANSITION_FILE'])
    print("  Rotamer-Conformer Cross File:", results['PREFERENCE_FILE'])


if __name__ == "__main__":
    main()
