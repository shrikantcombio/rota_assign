"""
rotaassign.py
-------------
Module for batch assigning rotameric states to protein residues in complex, bound, unbound, or both states.
Provides class ProteinRotamerAssigner with flexible DataFrame/CSV file input handling and options for state assignment.
Supports Backbone Dependent (BBD) rotamer assignment when phi and psi dihedrals are present.
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np

# Path setup for script directory
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

ROOT_DIR = SCRIPT_DIR.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from script.get_rotastate import get_aa_rotamer, RotamerStateAssigner


class ProteinRotamerAssigner:
    """
    Class for assigning rotamer states to protein dataset DataFrames or CSV files.
    """

    def __init__(self):
        self.state_assigner = RotamerStateAssigner()

    def assign_df(self, df, mode='both'):
        """
        Assigns rotamer states to pandas DataFrame based on mode ('bound', 'unbound', or 'both').
        """
        df = df.copy()
        mode_clean = str(mode).lower().strip()

        cols = {c.upper(): c for c in df.columns}

        label_col = cols.get('LABEL', 'LABEL')
        interface_col = cols.get('SASA', cols.get('INTERFACE', 'SASA'))

        if mode_clean in ['bound', 'complex', 'both']:
            b_ss_col = cols.get('B_SS', cols.get('SS_B', 'B_SS'))
            b_phi_col = cols.get('B_PHI', cols.get('PHI_B', 'B_PHI'))
            b_psi_col = cols.get('B_PSI', cols.get('PSI_B', 'B_PSI'))

            b_chi1 = cols.get('B_CHI1', 'B_CHI1')
            b_chi2 = cols.get('B_CHI2', 'B_CHI2')
            b_chi3 = cols.get('B_CHI3', 'B_CHI3')
            b_chi4 = cols.get('B_CHI4', 'B_CHI4')
            b_chi5 = cols.get('B_CHI5', 'B_CHI5')

            b_states = []
            for idx, row in df.iterrows():
                if 'AA' in row and not pd.isna(row['AA']):
                    aa = str(row['AA']).strip()
                elif label_col in row and not pd.isna(row[label_col]):
                    parts = str(row[label_col]).split()
                    aa = parts[0] if len(parts) > 0 else 'GLY'
                else:
                    aa = 'GLY'

                ss = row.get(b_ss_col, 'L')
                inter = row.get(interface_col, 'N')
                c1 = row.get(b_chi1, None)
                c2 = row.get(b_chi2, None)
                c3 = row.get(b_chi3, None)
                c4 = row.get(b_chi4, None)
                c5 = row.get(b_chi5, None)

                phi = row.get(b_phi_col, None)
                psi = row.get(b_psi_col, None)

                chis = [c1, c2, c3, c4, c5]
                chis_clean = [None if (x is None or (isinstance(x, (int, float)) and (float(x) == 0.0 or np.isnan(float(x))))) else float(x) for x in chis]

                r_state = get_aa_rotamer(aa, ss, 'bound', inter, chis_clean[0], chis_clean[1], chis_clean[2], chis_clean[3], chis_clean[4], phi=phi, psi=psi)
                b_states.append(r_state)

            df['b_rotastate'] = b_states

        if mode_clean in ['unbound', 'both']:
            u_ss_col = cols.get('U_SS', cols.get('SS_UB', 'U_SS'))
            u_phi_col = cols.get('U_PHI', cols.get('PHI_UB', 'U_PHI'))
            u_psi_col = cols.get('U_PSI', cols.get('PSI_UB', 'U_PSI'))

            u_chi1 = cols.get('U_CHI1', 'U_CHI1')
            u_chi2 = cols.get('U_CHI2', 'U_CHI2')
            u_chi3 = cols.get('U_CHI3', 'U_CHI3')
            u_chi4 = cols.get('U_CHI4', 'U_CHI4')
            u_chi5 = cols.get('U_CHI5', 'U_CHI5')

            u_states = []
            for idx, row in df.iterrows():
                if 'AA' in row and not pd.isna(row['AA']):
                    aa = str(row['AA']).strip()
                elif label_col in row and not pd.isna(row[label_col]):
                    parts = str(row[label_col]).split()
                    aa = parts[0] if len(parts) > 0 else 'GLY'
                else:
                    aa = 'GLY'

                ss = row.get(u_ss_col, 'L')
                inter = row.get(interface_col, 'N')
                c1 = row.get(u_chi1, None)
                c2 = row.get(u_chi2, None)
                c3 = row.get(u_chi3, None)
                c4 = row.get(u_chi4, None)
                c5 = row.get(u_chi5, None)

                phi = row.get(u_phi_col, None)
                psi = row.get(u_psi_col, None)

                chis = [c1, c2, c3, c4, c5]
                chis_clean = [None if (x is None or (isinstance(x, (int, float)) and (float(x) == 0.0 or np.isnan(float(x))))) else float(x) for x in chis]

                r_state = get_aa_rotamer(aa, ss, 'unbound', inter, chis_clean[0], chis_clean[1], chis_clean[2], chis_clean[3], chis_clean[4], phi=phi, psi=psi)
                u_states.append(r_state)

            df['ub_rotastate'] = u_states

        return df

    def assign_file(self, input_file, output_file=None, mode='both'):
        in_p = Path(input_file)
        if not in_p.exists():
            raise FileNotFoundError(f"Input file not found: {in_p}")

        sep = ',' if in_p.suffix.lower() == '.csv' else '\t'
        df = pd.read_csv(in_p, sep=sep, engine='python')

        df_out = self.assign_df(df, mode=mode)

        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            out_sep = ',' if out_p.suffix.lower() == '.csv' else '\t'
            df_out.to_csv(out_p, sep=out_sep, index=False)

        return df_out


def assign_rotamer_states(input_target, mode='both', output_file=None):
    assigner = ProteinRotamerAssigner()
    if isinstance(input_target, (str, Path)):
        return assigner.assign_file(input_target, output_file=output_file, mode=mode)
    elif isinstance(input_target, pd.DataFrame):
        df_out = assigner.assign_df(input_target, mode=mode)
        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            df_out.to_csv(out_p, index=False)
        return df_out
    else:
        raise TypeError("input_target must be a file path or pandas DataFrame")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Assign rotameric states to protein dataset files.")
    parser.add_argument("input", help="Input CSV/TSV file containing protein dihedrals")
    parser.add_argument("-o", "--output", default=None, help="Output file path")
    parser.add_argument("-m", "--mode", choices=['bound', 'unbound', 'both'], default='both', help="Assignment mode")

    args = parser.parse_args()
    df_res = assign_rotamer_states(args.input, mode=args.mode, output_file=args.output)
    print(f"Assigned rotamer states shape: {df_res.shape}")
