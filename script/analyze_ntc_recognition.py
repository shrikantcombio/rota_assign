"""
analyze_ntc_recognition.py
--------------------------
Comprehensive statistical analysis engine for RNA Dinucleotide Conformers (NtCs)
in Protein-RNA recognition across the PRDBv35 benchmark (263 complexes).

Answers:
1. Interface vs. Non-Interface NtC distribution & enrichment (Fisher's exact test, Odds Ratios).
2. Prevalent interface NtCs, Amino Acid preferences, and Reciprocal Recognition matrices.
2a. Unified interaction grading across all 14 interaction types (H-bonds, Stacking, Water, Pi, Salt-bridges).
2b. RNA Secondary Structure associations with NtCs and Amino Acid recognition.
3. Interaction Energetics (E_int in kcal/mol) distributions and stabilizing potentials per NtC class.
"""

import os
import sys
import json
import math
import argparse
from pathlib import Path
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
from scipy import stats

# Protein-RNA results directory; set ROTA_ASSIGN_PRORNA_DIR
PRO_RNA_RESULTS_DIR = Path(os.environ.get("ROTA_ASSIGN_PRORNA_DIR", Path(__file__).resolve().parent.parent / "dataset")) / "results"


# Path setup
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# Standard 20 Amino Acids
STANDARD_AAS = ['ARG', 'LYS', 'HIS', 'ASP', 'GLU', 'ASN', 'GLN', 'SER', 'THR', 'TYR',
                'PHE', 'TRP', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'GLY', 'CYS']

# Chemical Group Mapping
AA_GROUP_MAP = {
    'ARG': 'BASIC_POSITIVE', 'LYS': 'BASIC_POSITIVE', 'HIS': 'BASIC_POSITIVE',
    'ASP': 'ACIDIC_NEGATIVE', 'GLU': 'ACIDIC_NEGATIVE',
    'ASN': 'POLAR_UNCHARGED', 'GLN': 'POLAR_UNCHARGED', 'SER': 'POLAR_UNCHARGED', 'THR': 'POLAR_UNCHARGED',
    'TYR': 'AROMATIC', 'PHE': 'AROMATIC', 'TRP': 'AROMATIC',
    'ALA': 'HYDROPHOBIC', 'VAL': 'HYDROPHOBIC', 'LEU': 'HYDROPHOBIC', 'ILE': 'HYDROPHOBIC', 'MET': 'HYDROPHOBIC',
    'PRO': 'SPECIAL', 'GLY': 'SPECIAL', 'CYS': 'SPECIAL'
}

# Interaction Type Hierarchical Weights
INTERACTION_WEIGHTS = {
    'H-BOND': 1.0,
    'SALT-BRIDGE': 1.2,
    'STACK': 1.0,
    'CATION-PI': 1.0,
    'WEAK_H-BOND': 0.5,
    'WATER': 0.6,
    'C-H-PI_ACYCLIC': 0.6,
    'C-H-PI_CYCLIC': 0.7,
    'O-H-PI_ACYCLIC': 0.6,
    'O-H-PI_CYCLIC': 0.7,
    'N-H-PI_ACYCLIC': 0.6,
    'N-H-PI_CYCLIC': 0.7,
    'SUGAR-PI': 0.6,
    'ION': 0.8
}


def compute_shannon_entropy(prob_list):
    """Computes Shannon entropy in bits for a probability distribution."""
    ent = 0.0
    for p in prob_list:
        if p > 0:
            ent -= p * math.log2(p)
    return round(ent, 3)


def benjamini_hochberg(p_values):
    """Computes Benjamini-Hochberg FDR corrected q-values."""
    p_vals = np.array(p_values)
    n = len(p_vals)
    sorted_indices = np.argsort(p_vals)
    sorted_p = p_vals[sorted_indices]
    q_vals = np.zeros(n)
    min_q = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        q = (sorted_p[i] * n) / rank
        min_q = min(min_q, q)
        q_vals[sorted_indices[i]] = min(1.0, min_q)
    return q_vals.tolist()


class NtcRecognitionAnalyzer:
    """Statistical Analyzer for RNA dinucleotide conformer recognition."""

    def __init__(self, benchmark_json_path):
        self.json_path = Path(benchmark_json_path)
        with open(self.json_path, 'r', encoding='utf-8') as f:
            self.data = json.load(f)

        self.complexes = self.data.get('complexes', {})
        self.metadata = self.data.get('metadata', {})
        self._extract_flat_data()

    def _extract_flat_data(self):
        """Extracts flat records for steps and interactions."""
        self.step_records = []
        self.interaction_records = []

        for pdb, comp in self.complexes.items():
            for step in comp.get('steps', []):
                ntc = step.get('ntc', 'NANT')
                if not ntc or ntc in ['NANT', 'NAN', '']:
                    continue

                is_int = step.get('is_interface', False)
                ss_info = step.get('secondary_structure', {})
                ss_type = ss_info.get('ss_type', 'unpaired')
                if not ss_type or ss_type in ['.', '-', '']:
                    ss_type = 'unpaired'

                step_rec = {
                    'pdb_id': pdb,
                    'step_id': step.get('step_id', ''),
                    'ntc': ntc,
                    'cana': step.get('cana', ''),
                    'confal_a': step.get('confal_a', 0.0),
                    'is_interface': is_int,
                    'sasa': 'I' if is_int else 'N',
                    'ss_type': ss_type,
                    'bp_type': ss_info.get('bp_type', '.'),
                    'num_interactions': step.get('num_interactions', 0)
                }
                self.step_records.append(step_rec)

                for interaction in step.get('interactions', []):
                    prot = interaction.get('protein_residue', {})
                    aa = str(prot.get('aa', '')).upper().strip()
                    if not aa:
                        continue

                    int_type = interaction.get('interaction_type', 'H-BOND')
                    energy = interaction.get('energy_kcalmol', 0.0)
                    score = interaction.get('score', 0.0)
                    dist = interaction.get('dist', '-')

                    int_rec = {
                        'pdb_id': pdb,
                        'step_id': step.get('step_id', ''),
                        'ntc': ntc,
                        'cana': step.get('cana', ''),
                        'confal_a': step.get('confal_a', 0.0),
                        'ss_type': ss_type,
                        'aa': aa,
                        'aa_class': AA_GROUP_MAP.get(aa, 'OTHER'),
                        'aa_rotastate': prot.get('rotastate', 'NA'),
                        'interaction_type': int_type,
                        'edge': interaction.get('edge', '-'),
                        'score': score,
                        'energy_kcalmol': energy,
                        'dist': dist
                    }
                    self.interaction_records.append(int_rec)

        self.df_steps = pd.DataFrame(self.step_records)
        self.df_ints = pd.DataFrame(self.interaction_records)
        print(f"[INFO] Extracted {len(self.df_steps)} valid RNA dinucleotide steps and {len(self.df_ints)} interactions.")

    def analyze_interface_vs_noninterface(self):
        """
        Question 1: Analyzes distribution and statistical enrichment of 28 NtCs
        in interface ('I') vs non-interface ('N') regions.
        """
        ct = pd.crosstab(self.df_steps['ntc'], self.df_steps['sasa'])
        if 'I' not in ct.columns:
            ct['I'] = 0
        if 'N' not in ct.columns:
            ct['N'] = 0

        tot_I = ct['I'].sum()
        tot_N = ct['N'].sum()

        results = []
        p_values = []

        for ntc, row in ct.iterrows():
            a = int(row['I'])
            b = int(row['N'])
            c = tot_I - a
            d = tot_N - b

            # Odds Ratio
            or_val = (a * d) / (b * c) if (b * c) > 0 else (999.0 if a > 0 else 0.0)
            log2_or = round(math.log2(or_val), 3) if or_val > 0 else -99.0

            # Fisher's exact test
            odds_ratio_stat, p_val = stats.fisher_exact([[a, b], [c, d]])
            p_values.append(p_val)

            pct_I = round((a / tot_I * 100), 2) if tot_I > 0 else 0.0
            pct_N = round((b / tot_N * 100), 2) if tot_N > 0 else 0.0

            results.append({
                'NtC': ntc,
                'Interface_Count': a,
                'NonInterface_Count': b,
                'Total_Count': a + b,
                'Interface_Pct': pct_I,
                'NonInterface_Pct': pct_N,
                'Odds_Ratio': round(or_val, 3),
                'Log2_Odds_Ratio': log2_or,
                'P_Value': p_val
            })

        q_values = benjamini_hochberg(p_values)
        for i, res in enumerate(results):
            res['FDR_Q_Value'] = round(q_values[i], 5)
            # Classification
            if res['FDR_Q_Value'] < 0.05 and res['Log2_Odds_Ratio'] > 0.3:
                res['Interface_Preference'] = 'Significantly Enriched'
            elif res['FDR_Q_Value'] < 0.05 and res['Log2_Odds_Ratio'] < -0.3:
                res['Interface_Preference'] = 'Significantly Depleted'
            else:
                res['Interface_Preference'] = 'Neutral / Proportional'

        df_res = pd.DataFrame(results).sort_values(by='Interface_Count', ascending=False).reset_index(drop=True)
        return df_res

    def analyze_amino_acid_preferences(self):
        """
        Questions 2 & 2b: Cross-tabulates NtCs with 20 Amino Acids and Chemical Classes.
        Computes Log-Odds propensity scores and Shannon specificity entropy.
        """
        # Filter interactions to standard amino acids
        df_std = self.df_ints[self.df_ints['aa'].isin(STANDARD_AAS)].copy()

        # 1. NtC x AA Contingency Table
        ct_aa = pd.crosstab(df_std['ntc'], df_std['aa'])
        # Reindex to all standard AAs
        for aa in STANDARD_AAS:
            if aa not in ct_aa.columns:
                ct_aa[aa] = 0
        ct_aa = ct_aa[STANDARD_AAS]

        tot_int = len(df_std)
        p_aa_marginal = df_std['aa'].value_counts() / tot_int
        p_ntc_marginal = df_std['ntc'].value_counts() / tot_int

        # Compute Propensity Matrix: log2( P(AA | NtC) / P(AA) )
        propensity_matrix = pd.DataFrame(index=ct_aa.index, columns=STANDARD_AAS, dtype=float)
        specificity_list = []

        for ntc in ct_aa.index:
            row = ct_aa.loc[ntc]
            n_ntc_ints = row.sum()
            if n_ntc_ints == 0:
                propensity_matrix.loc[ntc] = 0.0
                specificity_list.append({'NtC': ntc, 'Total_Interactions': 0, 'Shannon_Entropy': 0.0, 'Specificity': 'None'})
                continue

            p_aa_given_ntc = row / n_ntc_ints
            entropy = compute_shannon_entropy(p_aa_given_ntc.values)

            # Top preferred AA
            top_aa = p_aa_given_ntc.idxmax()
            top_aa_pct = round(p_aa_given_ntc.max() * 100, 1)

            specificity_desc = 'High Specificity' if entropy < 3.2 else ('Moderate Specificity' if entropy < 3.8 else 'Broad / Promiscuous')

            specificity_list.append({
                'NtC': ntc,
                'Total_Interactions': int(n_ntc_ints),
                'Shannon_Entropy_Bits': entropy,
                'Max_Theoretical_Entropy': round(math.log2(20), 3),
                'Top_Preferred_AA': f"{top_aa} ({top_aa_pct}%)",
                'Specificity_Grade': specificity_desc
            })

            for aa in STANDARD_AAS:
                obs = row[aa]
                exp = n_ntc_ints * p_aa_marginal.get(aa, 1e-6)
                if obs > 0 and exp > 0:
                    prop = math.log2(obs / exp)
                elif obs == 0:
                    prop = -3.5  # Penalize unobserved contacts
                else:
                    prop = 0.0
                propensity_matrix.loc[ntc, aa] = round(prop, 2)

        df_spec = pd.DataFrame(specificity_list).sort_values(by='Total_Interactions', ascending=False).reset_index(drop=True)

        # 2. NtC x Chemical Class Contingency & Propensity
        ct_class = pd.crosstab(df_std['ntc'], df_std['aa_class'])
        p_class_marginal = df_std['aa_class'].value_counts() / tot_int

        prop_class_matrix = pd.DataFrame(index=ct_class.index, columns=ct_class.columns, dtype=float)
        for ntc in ct_class.index:
            row = ct_class.loc[ntc]
            n_ntc = row.sum()
            for grp in ct_class.columns:
                obs = row[grp]
                exp = n_ntc * p_class_marginal.get(grp, 1e-6)
                if obs > 0 and exp > 0:
                    prop = math.log2(obs / exp)
                elif obs == 0:
                    prop = -3.5
                else:
                    prop = 0.0
                prop_class_matrix.loc[ntc, grp] = round(prop, 2)

        return ct_aa, propensity_matrix, df_spec, ct_class, prop_class_matrix

    def analyze_interaction_types(self):
        """
        Question 2a: Unifies and grades all 14 interaction types per NtC.
        Computes interaction composition, primary direct contacts vs secondary vs water,
        and weighted composite contact index.
        """
        ct_types = pd.crosstab(self.df_ints['ntc'], self.df_ints['interaction_type'])

        breakdown_rows = []
        for ntc, row in ct_types.iterrows():
            tot = int(row.sum())
            if tot == 0:
                continue

            hb = int(row.get('H-BOND', 0))
            whb = int(row.get('WEAK_H-BOND', 0))
            wat = int(row.get('WATER', 0))
            sb = int(row.get('SALT-BRIDGE', 0))
            stk = int(row.get('STACK', 0))
            cat_pi = int(row.get('CATION-PI', 0))

            # All Pi interactions (C-H-PI, O-H-PI, N-H-PI, SUGAR-PI)
            pi_sum = sum(int(row.get(k, 0)) for k in row.index if 'PI' in k and k != 'CATION-PI')

            # Direct strong contacts = H-BOND + SALT-BRIDGE + STACK + CATION-PI
            direct_strong = hb + sb + stk + cat_pi

            # Calculate weighted composite contact score for this NtC
            weighted_score = 0.0
            for itype, count in row.items():
                w = INTERACTION_WEIGHTS.get(itype, 0.5)
                weighted_score += count * w

            avg_weighted_score_per_contact = round(weighted_score / tot, 3)

            breakdown_rows.append({
                'NtC': ntc,
                'Total_Contacts': tot,
                'H_Bond_Count': hb,
                'H_Bond_Pct': round(hb / tot * 100, 1),
                'Weak_H_Bond_Count': whb,
                'Weak_H_Bond_Pct': round(whb / tot * 100, 1),
                'Water_Mediated_Count': wat,
                'Water_Mediated_Pct': round(wat / tot * 100, 1),
                'Stacking_Count': stk,
                'Stacking_Pct': round(stk / tot * 100, 1),
                'Salt_Bridge_Count': sb,
                'Salt_Bridge_Pct': round(sb / tot * 100, 1),
                'Pi_Interactions_Count': pi_sum,
                'Pi_Interactions_Pct': round(pi_sum / tot * 100, 1),
                'Direct_Strong_Ratio': round(direct_strong / tot, 3),
                'Avg_Weighted_Contact_Quality': avg_weighted_score_per_contact
            })

        df_breakdown = pd.DataFrame(breakdown_rows).sort_values(by='Total_Contacts', ascending=False).reset_index(drop=True)
        return ct_types, df_breakdown

    def analyze_secondary_structure_associations(self):
        """
        User Addition (2b): Analyzes how RNA secondary structure elements
        (stem, hairpin, junction, bulge, single-stranded) associate with different NtCs
        and modulate protein recognition.
        """
        # NtC x Secondary Structure contingency across all steps
        ct_ss = pd.crosstab(self.df_steps['ntc'], self.df_steps['ss_type'])
        
        # Interaction counts by secondary structure type
        ct_ss_ints = pd.crosstab(self.df_ints['ss_type'], self.df_ints['aa_class'])

        ss_summary = []
        for ntc, row in ct_ss.iterrows():
            tot_steps = int(row.sum())
            if tot_steps == 0:
                continue

            # Dominant SS
            top_ss = row.idxmax()
            top_ss_pct = round(row.max() / tot_steps * 100, 1)

            ss_summary.append({
                'NtC': ntc,
                'Total_Steps': tot_steps,
                'Dominant_Secondary_Structure': f"{top_ss} ({top_ss_pct}%)",
                'Stem_Count': int(row.get('stem', 0)),
                'Hairpin_Count': int(row.get('hairpin', 0)),
                'Junction_Count': int(row.get('junction', 0)),
                'Bulge_Count': int(row.get('bulge', 0)),
                'SingleStranded_Count': int(row.get('single-stranded', 0)),
                'InternalLoop_Count': int(row.get('internal-loop', 0))
            })

        df_ss_summary = pd.DataFrame(ss_summary).sort_values(by='Total_Steps', ascending=False).reset_index(drop=True)
        return ct_ss, df_ss_summary, ct_ss_ints

    def analyze_interaction_energetics(self):
        """
        Question 3: Evaluates interaction energies (E_int in kcal/mol) per NtC class.
        Calculates energy statistics, net stabilizing capacity, and CONFAL score correlation.
        """
        energy_valid = self.df_ints[
            (self.df_ints['energy_kcalmol'] != 0.0) &
            (~pd.isna(self.df_ints['energy_kcalmol']))
        ].copy()

        energy_rows = []
        for ntc, grp in energy_valid.groupby('ntc'):
            e_vals = grp['energy_kcalmol'].values
            n_contacts = len(e_vals)
            if n_contacts < 3:
                continue

            mean_e = round(float(np.mean(e_vals)), 3)
            median_e = round(float(np.median(e_vals)), 3)
            std_e = round(float(np.std(e_vals)), 3)
            min_e = round(float(np.min(e_vals)), 3)  # Most stabilizing
            q25 = round(float(np.percentile(e_vals, 25)), 3)
            q75 = round(float(np.percentile(e_vals, 75)), 3)

            # Fraction of strongly favorable contacts (E < -10 kcal/mol)
            favorable_count = int(np.sum(e_vals < -5.0))
            favorable_pct = round(favorable_count / n_contacts * 100, 1)

            # Average CONFAL score of contacted steps
            avg_confal = round(float(grp['confal_a'].mean()), 1)

            energy_rows.append({
                'NtC': ntc,
                'Energy_Contacts_Count': n_contacts,
                'Mean_E_int_kcalmol': mean_e,
                'Median_E_int_kcalmol': median_e,
                'Std_E_int': std_e,
                'Min_E_int_kcalmol': min_e,
                'IQR_E_int': round(q75 - q25, 3),
                'Favorable_Contacts_Pct': favorable_pct,
                'Mean_CONFAL_Score': avg_confal
            })

        df_energy = pd.DataFrame(energy_rows).sort_values(by='Energy_Contacts_Count', ascending=False).reset_index(drop=True)
        return df_energy


def run_full_ntc_analysis(benchmark_json_path, output_dir=None):
    """Executes full NtC statistical analysis pipeline and exports tables and report."""
    analyzer = NtcRecognitionAnalyzer(benchmark_json_path)
    out_d = Path(output_dir) if output_dir else Path(benchmark_json_path).parent
    out_d.mkdir(parents=True, exist_ok=True)

    print("\n--- Running Statistical Analyses ---")

    # 1. Interface vs Non-Interface
    print("1. Analyzing Interface vs. Non-Interface NtC Distribution...")
    df_int_vs_non = analyzer.analyze_interface_vs_noninterface()
    out_q1 = out_d / "ntc_interface_vs_noninterface_distribution.csv"
    df_int_vs_non.to_csv(out_q1, index=False)

    # 2. Amino Acid Preferences & Reciprocal Recognition
    print("2. Analyzing Amino Acid Preferences & Recognition Specificity...")
    ct_aa, prop_aa, df_spec, ct_class, prop_class = analyzer.analyze_amino_acid_preferences()
    out_q2_aa = out_d / "ntc_amino_acid_contingency.csv"
    ct_aa.to_csv(out_q2_aa)
    out_q2_prop = out_d / "ntc_amino_acid_propensity_matrix.csv"
    prop_aa.to_csv(out_q2_prop)
    out_q2_spec = out_d / "ntc_recognition_specificity.csv"
    df_spec.to_csv(out_q2_spec, index=False)
    out_q2_class = out_d / "ntc_chemical_class_propensity.csv"
    prop_class.to_csv(out_q2_class)

    # 2a. Interaction Types
    print("3. Analyzing Interaction Types and Grading...")
    ct_types, df_types = analyzer.analyze_interaction_types()
    out_q2a = out_d / "ntc_interaction_type_breakdown.csv"
    df_types.to_csv(out_q2a, index=False)

    # 2b. Secondary Structure
    print("4. Analyzing Secondary Structure Associations...")
    ct_ss, df_ss, ct_ss_ints = analyzer.analyze_secondary_structure_associations()
    out_q2b_ss = out_d / "ntc_secondary_structure_distribution.csv"
    df_ss.to_csv(out_q2b_ss, index=False)
    out_q2b_ints = out_d / "secondary_structure_aa_class_matrix.csv"
    ct_ss_ints.to_csv(out_q2b_ints)

    # 3. Energetics
    print("5. Analyzing Interaction Energetics...")
    df_energy = analyzer.analyze_interaction_energetics()
    out_q3 = out_d / "ntc_energetics_summary.csv"
    df_energy.to_csv(out_q3, index=False)

    # Compile Comprehensive JSON Report
    report = {
        "analysis_metadata": {
            "benchmark_complexes": analyzer.metadata.get('total_complexes', 263),
            "total_rna_steps": len(analyzer.df_steps),
            "total_interface_steps": int(analyzer.df_steps['is_interface'].sum()),
            "total_non_interface_steps": int((~analyzer.df_steps['is_interface']).sum()),
            "total_interactions_analyzed": len(analyzer.df_ints)
        },
        "q1_interface_distribution_top5": df_int_vs_non.head(5).to_dict(orient='records'),
        "q1_significantly_enriched_ntcs": df_int_vs_non[df_int_vs_non['Interface_Preference'] == 'Significantly Enriched'].to_dict(orient='records'),
        "q1_significantly_depleted_ntcs": df_int_vs_non[df_int_vs_non['Interface_Preference'] == 'Significantly Depleted'].to_dict(orient='records'),
        "q2_recognition_specificity_top5": df_spec.head(5).to_dict(orient='records'),
        "q2a_interaction_breakdown_top5": df_types.head(5).to_dict(orient='records'),
        "q2b_secondary_structure_top5": df_ss.head(5).to_dict(orient='records'),
        "q3_energetics_top5": df_energy.head(5).to_dict(orient='records')
    }

    out_report = out_d / "ntc_recognition_analysis_summary.json"
    with open(out_report, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2)

    print(f"\n[SUCCESS] Statistical analysis complete! Output reports generated in: {out_d}")
    return report


def main():
    parser = argparse.ArgumentParser(description="Statistical analysis of RNA dinucleotide conformer (NtC) recognition by proteins.")
    parser.add_argument("--json", default=str(PRO_RNA_RESULTS_DIR / "PRDBv35_unified_ntc_interactions.json"), help="Path to unified benchmark JSON")
    parser.add_argument("-o", "--output-dir", default=str(PRO_RNA_RESULTS_DIR), help="Directory to save output CSVs and reports")
    args = parser.parse_args()

    run_full_ntc_analysis(args.json, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
