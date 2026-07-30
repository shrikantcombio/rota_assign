"""
merge_pipeline.py
-----------------
Master end-to-end processing pipeline module.
1. Calculates protein & RNA dihedrals, DSSP, and NACCESS SASA -> dataset/results/1ASY_protein_torsions.csv
2. Assigns protein rotamer states using BBD rotamer library (library/RBPs_BBD_rotamer_lib) -> dataset/results/1ASY_protein_rotastates.csv
3. Assigns RNA dinucleotide conformers (NtC, CANA, confalA) using StructAnnotate -> dataset/results/1ASY_rna_conformers.csv
4. Uses interaction details (including NAxplorer dataset/1ASY_interactions_NAxplorer.csv) -> dataset/results/1ASY_interactions.csv
5. Exports final merged dataset into dataset/results/1ASY_final_output.csv matching schema:
   PDB_ID,AA,AA_atom,AA_chain,AA_pos,Nuc,Nuc_atom,Nuc_chain,Nuc_pos,hb_dis,NtC,CANA,confalA,b_phi,b_psi,b_ss,b_rotastate
6. Runs Arginine forks & rotamer transition preference analysis.
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np

# Path setup
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
RESULTS_DIR = ROOT_DIR / "dataset" / "results"
TEST_DIR = ROOT_DIR / "test"
DATASET_DIR = ROOT_DIR / "dataset"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from macromol_torsion import ComplexTorsionCalculator
from script.rotaassign import ProteinRotamerAssigner
from script.rna_conformer_assigner import RNAConformerAssigner
from script.analyze_rotastates import run_full_rotastate_analysis


class ComplexMergePipeline:
    """
    Master pipeline for running torsion calculations, rotamer assignment, RNA conformer annotation,
    and merging interaction data into final output datasets under dataset/results/.
    """

    def __init__(self, output_dir=None):
        self.output_dir = Path(output_dir) if output_dir else RESULTS_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.rotamer_assigner = ProteinRotamerAssigner()
        self.rna_assigner = RNAConformerAssigner()
        self.calculator = ComplexTorsionCalculator()

    def process_interactions_file(self, interactions_source, target_pdb="1ASY"):
        """
        Loads interaction details from NAxplorer dataset/1ASY_interactions_NAxplorer.csv or standard CSV
        and writes intermediate dataset/results/1ASY_interactions.csv.
        """
        src_p = Path(interactions_source)
        if not src_p.exists():
            raise FileNotFoundError(f"Interactions source file not found: {src_p}")

        try:
            df_raw = pd.read_csv(src_p, sep=',')
            if len(df_raw.columns) == 1 and '\t' in str(df_raw.columns[0]):
                df_raw = pd.read_csv(src_p, sep='\t')
        except Exception:
            df_raw = pd.read_csv(src_p, sep='\t')

        req_cols = ['PDB_ID', 'AA', 'AA_atom', 'AA_chain', 'AA_pos', 'Nuc', 'Nuc_atom', 'Nuc_chain', 'Nuc_pos', 'hb_dis']

        if all(col in df_raw.columns for col in req_cols):
            df_int = df_raw[req_cols].copy()
        elif 'rna_nt' in df_raw.columns and 'prot_aa' in df_raw.columns:
            rows = []
            for idx, r in df_raw.iterrows():
                pdb_id = str(r.get('pdb_id', target_pdb)).strip().upper()
                rna_parts = str(r['rna_nt']).split(':')
                prot_parts = str(r['prot_aa']).split(':')
                atom_parts = str(r.get('atoms', ':')).split(':')

                nuc_chain = rna_parts[0] if len(rna_parts) > 0 else 'R'
                nuc = rna_parts[1] if len(rna_parts) > 1 else 'U'
                nuc_pos = rna_parts[2] if len(rna_parts) > 2 else '0'

                aa_chain = prot_parts[0] if len(prot_parts) > 0 else 'A'
                aa = prot_parts[1] if len(prot_parts) > 1 else 'GLY'
                aa_pos = prot_parts[2] if len(prot_parts) > 2 else '0'

                nuc_atom = atom_parts[0] if len(atom_parts) > 0 else ''
                aa_atom = atom_parts[1] if len(atom_parts) > 1 else ''

                hb_dis = r.get('dist', r.get('hb_dis', 0.0))

                rows.append({
                    'PDB_ID': pdb_id,
                    'AA': aa,
                    'AA_atom': aa_atom,
                    'AA_chain': aa_chain,
                    'AA_pos': aa_pos,
                    'Nuc': nuc,
                    'Nuc_atom': nuc_atom,
                    'Nuc_chain': nuc_chain,
                    'Nuc_pos': nuc_pos,
                    'hb_dis': hb_dis
                })
            df_int = pd.DataFrame(rows)
        else:
            df_int = df_raw.copy()

        out_int_file = self.output_dir / f"{target_pdb}_interactions.csv"
        df_int.to_csv(out_int_file, index=False)
        return out_int_file, df_int

    def run_pipeline(self, bound_file, unbound_file, protein_chains="A", rna_chains="R", interactions_file=None, pdb_id="1ASY"):
        """
        Executes end-to-end pipeline:
        1. Analyzes complex dihedrals, DSSP, and SASA -> 1ASY_protein_torsions.csv
        2. Assigns protein rotamer states using BBD library -> intermediate 1ASY_protein_rotastates.csv
        3. Assigns RNA dinucleotide conformers -> intermediate 1ASY_rna_conformers.csv
        4. Processes interaction details -> intermediate 1ASY_interactions.csv
        5. Saves final merged 1ASY_final_output.csv
        6. Runs Arginine forks & rotamer transition analysis reports
        """
        bound_p = Path(bound_file)
        unbound_p = Path(unbound_file)

        # Step 1: Analyze complex dihedrals & SASA
        pair_tag = f"{bound_p.stem.split('_')[0]}_{unbound_p.stem.split('_')[0]}"
        output_pair_dir = self.output_dir / "pairs" / pair_tag
        output_pair_dir.mkdir(parents=True, exist_ok=True)

        calc_results = self.calculator.analyze_complex(
            bound_file=bound_p,
            unbound_file=unbound_p,
            protein_chains=protein_chains,
            rna_chains=rna_chains,
            orig_bound_file=bound_p,
            orig_unbound_file=unbound_p,
            output_file=output_pair_dir / f"{pair_tag}.csv"
        )

        df_protein = calc_results['PROTEIN']
        df_rna = calc_results['RNA']

        # Save intermediate raw protein torsions file
        prot_torsion_file = self.output_dir / f"{pdb_id}_protein_torsions.csv"
        df_protein.to_csv(prot_torsion_file, index=False)

        # Step 2: Assign Protein rotamer states (intermediate file: 1ASY_protein_rotastates.csv)
        df_protein_assigned = self.rotamer_assigner.assign_df(df_protein, mode='both')
        prot_rota_file = self.output_dir / f"{pdb_id}_protein_rotastates.csv"
        df_protein_assigned.to_csv(prot_rota_file, index=False)

        # Step 3: Assign RNA dinucleotide conformers (intermediate file: 1ASY_rna_conformers.csv)
        df_rna_assigned = self.rna_assigner.assign_df(df_rna)
        rna_conf_file = self.output_dir / f"{pdb_id}_rna_conformers.csv"
        df_rna_assigned.to_csv(rna_conf_file, index=False)

        # Step 4: Extract/process interactions (intermediate file: 1ASY_interactions.csv)
        if interactions_file:
            int_file_path, df_interactions = self.process_interactions_file(interactions_file, target_pdb=pdb_id)
        else:
            default_naxplorer = DATASET_DIR / f"{pdb_id}_interactions_NAxplorer.csv"
            default_test_int = TEST_DIR / f"{pdb_id}_interactions.csv"
            default_hb = ROOT_DIR / "hb_pro_rna_merge_files" / f"{pdb_id}.csv"

            if default_naxplorer.exists():
                int_file_path, df_interactions = self.process_interactions_file(default_naxplorer, target_pdb=pdb_id)
            elif default_test_int.exists():
                int_file_path, df_interactions = self.process_interactions_file(default_test_int, target_pdb=pdb_id)
            elif default_hb.exists():
                int_file_path, df_interactions = self.process_interactions_file(default_hb, target_pdb=pdb_id)
            else:
                df_interactions = None
                int_file_path = None

        # Step 5: Merge into final output DataFrame
        if df_interactions is not None:
            df_final = df_interactions.copy()

            prot_map = {}
            for idx, row in df_protein_assigned.iterrows():
                parts = str(row['LABEL']).split()
                chain = parts[1] if len(parts) > 1 else 'A'
                pos = parts[2] if len(parts) > 2 else '0'
                prot_map[(chain, str(pos))] = {
                    'b_phi': row.get('B_PHI', 0.0),
                    'b_psi': row.get('B_PSI', 0.0),
                    'b_ss': row.get('B_SS', 'L'),
                    'b_rotastate': row.get('b_rotastate', 'MA')
                }

            rna_map = {}
            for idx, row in df_rna_assigned.iterrows():
                parts = str(row['LABEL']).split()
                chain = parts[1] if len(parts) > 1 else 'R'
                pos = parts[2] if len(parts) > 2 else '0'
                rna_map[(chain, str(pos))] = {
                    'NtC': row.get('NtC', 'NANT'),
                    'CANA': row.get('CANA', 'NAN'),
                    'confalA': row.get('confalA', 0.0)
                }

            ntc_col, cana_col, confal_col = [], [], []
            phi_col, psi_col, ss_col, rota_col = [], [], [], []

            for idx, row in df_final.iterrows():
                p_chain = str(row.get('AA_chain', 'A'))
                p_pos = str(row.get('AA_pos', '0'))

                r_chain = str(row.get('Nuc_chain', 'R'))
                r_pos = str(row.get('Nuc_pos', '0'))

                p_data = prot_map.get((p_chain, p_pos), {
                    'b_phi': row.get('b_phi', 0.0),
                    'b_psi': row.get('b_psi', 0.0),
                    'b_ss': row.get('b_ss', 'L'),
                    'b_rotastate': row.get('b_rotastate', 'MA')
                })

                r_data = rna_map.get((r_chain, r_pos), {
                    'NtC': row.get('NtC', 'NANT'),
                    'CANA': row.get('CANA', 'NAN'),
                    'confalA': row.get('confalA', 0.0)
                })

                ntc_col.append(r_data['NtC'])
                cana_col.append(r_data['CANA'])
                confal_col.append(r_data['confalA'])

                phi_col.append(p_data['b_phi'])
                psi_col.append(p_data['b_psi'])
                ss_col.append(p_data['b_ss'])
                rota_col.append(p_data['b_rotastate'])

            df_final['NtC'] = ntc_col
            df_final['CANA'] = cana_col
            df_final['confalA'] = confal_col

            df_final['b_phi'] = phi_col
            df_final['b_psi'] = psi_col
            df_final['b_ss'] = ss_col
            df_final['b_rotastate'] = rota_col

            final_out_file = self.output_dir / f"{pdb_id}_final_output.csv"
            df_final.to_csv(final_out_file, index=False)

            # Step 6: Run Arginine forks & rotamer transition analysis reports
            analysis_results = run_full_rotastate_analysis(final_out_file, output_dir=self.output_dir)
        else:
            final_out_file = None
            df_final = None
            analysis_results = None

        return {
            'PROTEIN_TORSIONS_FILE': prot_torsion_file,
            'PROTEIN_ROTA_FILE': prot_rota_file,
            'RNA_CONF_FILE': rna_conf_file,
            'INTERACTIONS_FILE': int_file_path,
            'FINAL_OUTPUT_FILE': final_out_file,
            'FINAL_DF': df_final,
            'ANALYSIS': analysis_results
        }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="End-to-End Complex Merge Pipeline")
    parser.add_argument("--pdb", default="1ASY", help="PDB ID to process (default: 1ASY)")
    parser.add_argument("--bound", default=None, help="Bound PDB file path")
    parser.add_argument("--unbound", default=None, help="Unbound PDB file path")
    parser.add_argument("--interactions", default=None, help="Interactions CSV file path")
    parser.add_argument("-o", "--output-dir", default=None, help="Output directory (default: dataset/results)")

    args = parser.parse_args()

    out_dir = Path(args.output_dir) if args.output_dir else RESULTS_DIR

    bound_file = Path(args.bound) if args.bound else DATASET_DIR / "1ASY_1EOV" / "1ASY_1EOV_C.pdb"
    if not bound_file.exists():
        bound_file = DATASET_DIR / "1ASY.pdb"
    if not bound_file.exists():
        bound_file = TEST_DIR / "dataset" / "1ASY.pdb"

    unbound_file = Path(args.unbound) if args.unbound else DATASET_DIR / "1ASY_1EOV" / "1ASY_1EOV_U.pdb"
    if not unbound_file.exists():
        unbound_file = DATASET_DIR / "1EOV.pdb"
    if not unbound_file.exists():
        unbound_file = TEST_DIR / "dataset" / "1EOV.pdb"

    int_file = Path(args.interactions) if args.interactions else DATASET_DIR / f"{args.pdb}_interactions_NAxplorer.csv"
    if not int_file.exists():
        int_file = TEST_DIR / f"{args.pdb}_interactions.csv"

    pipeline = ComplexMergePipeline(output_dir=out_dir)
    res = pipeline.run_pipeline(bound_file, unbound_file, interactions_file=int_file if int_file.exists() else None, pdb_id=args.pdb)

    print("Pipeline Execution Completed Successfully!")
    print("Intermediate Protein Torsions File:", res['PROTEIN_TORSIONS_FILE'])
    print("Intermediate Protein Rotastates File:", res['PROTEIN_ROTA_FILE'])
    print("Intermediate RNA Conformers File:  ", res['RNA_CONF_FILE'])
    if res['INTERACTIONS_FILE']:
        print("Intermediate Interactions File:    ", res['INTERACTIONS_FILE'])
    if res['FINAL_OUTPUT_FILE']:
        print("Final Output File:                 ", res['FINAL_OUTPUT_FILE'])
    if res['ANALYSIS']:
        print("\nAnalysis Reports Output:")
        print("  Arginine Forks File:             ", res['ANALYSIS']['FORK_FILE'])
        print("  Rotamer Transitions File:        ", res['ANALYSIS']['TRANSITION_FILE'])
        print("  Rotamer-Conformer Cross File:    ", res['ANALYSIS']['PREFERENCE_FILE'])


if __name__ == "__main__":
    main()
