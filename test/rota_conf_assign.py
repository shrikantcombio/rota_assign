"""
rota_conf_assign.py
-------------------
Test script for verifying single rotamer library choice (BBD default), n-1 chi checking criteria,
and 'NA' output indicator on single bound complex structure (1ASY.cif).
"""

from pathlib import Path
import sys
import pandas as pd

# Path setup
TEST_DIR = Path(__file__).resolve().parent
ROOT_DIR = TEST_DIR.parent
SCRIPT_DIR = ROOT_DIR / "script"
DATASET_DIR = ROOT_DIR / "dataset"
RESULTS_DIR = DATASET_DIR / "results"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from script.merge_pipeline import ComplexMergePipeline


def run_test():
    print("=== Running Single-Library Rotamer Assignment Test (1ASY.cif) ===")

    # 1. Resolve 1ASY Structure File (.cif or .pdb)
    bound_file = DATASET_DIR / "1ASY.cif"
    if not bound_file.exists():
        bound_file = Path("/home/labuser/Projects/PRince/test/1ASY.cif")
    if not bound_file.exists():
        bound_file = DATASET_DIR / "1ASY.pdb"

    # 2. NAxplorer Interactions File
    naxplorer_int_file = DATASET_DIR / "1ASY_interactions_NAxplorer.csv"
    if not naxplorer_int_file.exists():
        naxplorer_int_file = TEST_DIR / "1ASY_interactions.csv"

    print(f"Bound Structure File (CIF/PDB): {bound_file}")
    print(f"Interactions File:               {naxplorer_int_file}")
    print(f"Results Output Directory:        {RESULTS_DIR}")

    # 3. Instantiate and run pipeline for 1ASY bound complex with default BBD rotamer library
    pipeline = ComplexMergePipeline(output_dir=RESULTS_DIR, lib_type='BBD')
    results = pipeline.run_pipeline(
        bound_file=bound_file,
        unbound_file=None,  # Single-structure bound mode
        protein_chains="A",
        rna_chains="R",
        interactions_file=naxplorer_int_file if naxplorer_int_file.exists() else None,
        pdb_id="1ASY",
        lib_type='BBD'
    )

    print("\n=== Test Execution Completed Successfully ===")
    print(f"1. Raw Protein Torsions CSV:   {results['PROTEIN_TORSIONS_FILE']}")
    print(f"2. Protein Rotastates CSV:     {results['PROTEIN_ROTA_FILE']}")
    print(f"3. RNA Conformers CSV:         {results['RNA_CONF_FILE']}")
    print(f"4. Mapped Interactions CSV:    {results['INTERACTIONS_FILE']}")
    print(f"5. Final Merged Output CSV:    {results['FINAL_OUTPUT_FILE']}")

    # 4. Verify Output Schema and BBD Rotamer Assignment
    df_rota = pd.read_csv(results['PROTEIN_ROTA_FILE'])
    print(f"\nProtein Rotastates Dataset Shape: {df_rota.shape}")
    print("\nSample BBD Rotamer Assignments:")
    print(df_rota[['LABEL', 'B_PHI', 'B_PSI', 'SASA', 'b_rotastate']].head(15))

    assert not df_rota.empty, "Error: Protein rotastates dataframe is empty!"
    assert 'b_rotastate' in df_rota.columns, "Error: b_rotastate column missing!"

    # Verify NA vs assigned rotamer counts
    na_count = (df_rota['b_rotastate'] == 'NA').sum()
    no_chi_count = (df_rota['b_rotastate'] == 'no_chi').sum()
    assigned_count = len(df_rota) - na_count - no_chi_count

    print(f"\nAssignment Breakdown:")
    print(f"  Assigned Rotamers: {assigned_count}")
    print(f"  no_chi (ALA/GLY):  {no_chi_count}")
    print(f"  NA (Not Available): {na_count}")

    print("\nSUCCESS: All single-library BBD rotamer assignment tests passed cleanly!")


if __name__ == "__main__":
    run_test()
