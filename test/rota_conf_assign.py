"""
rota_conf_assign.py
-------------------
Test script for verifying full end-to-end protein & RNA torsion calculation,
rotamer & RNA conformer assignment, NAxplorer interaction mapping, and export under dataset/results/.
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
    print("=== Running Pipeline Test: rota_conf_assign.py ===")

    # 1. Resolve Bound and Unbound Structure Files
    bound_file = DATASET_DIR / "1ASY_1EOV_C.pdb"
    if not bound_file.exists():
        bound_file = DATASET_DIR / "1ASY.pdb"
    if not bound_file.exists():
        bound_file = TEST_DIR / "dataset" / "1ASY_1EOV_C.pdb"

    unbound_file = DATASET_DIR / "1ASY_1EOV_U.pdb"
    if not unbound_file.exists():
        unbound_file = DATASET_DIR / "1EOV.pdb"
    if not unbound_file.exists():
        unbound_file = TEST_DIR / "dataset" / "1ASY_1EOV_U.pdb"

    # 2. NAxplorer Interactions File
    naxplorer_int_file = DATASET_DIR / "1ASY_interactions_NAxplorer.csv"
    if not naxplorer_int_file.exists():
        naxplorer_int_file = TEST_DIR / "1ASY_interactions.csv"

    print(f"Bound Structure File:      {bound_file}")
    print(f"Unbound Structure File:    {unbound_file}")
    print(f"NAxplorer Interactions File:{naxplorer_int_file}")
    print(f"Results Output Directory:  {RESULTS_DIR}")

    # 3. Instantiate and run pipeline
    pipeline = ComplexMergePipeline(output_dir=RESULTS_DIR)
    results = pipeline.run_pipeline(
        bound_file=bound_file,
        unbound_file=unbound_file,
        protein_chains="A",
        rna_chains="R",
        interactions_file=naxplorer_int_file,
        pdb_id="1ASY"
    )

    print("\n=== Test Execution Completed Successfully ===")
    print(f"1. Raw Protein Torsions CSV:   {results['PROTEIN_TORSIONS_FILE']}")
    print(f"2. Protein Rotastates CSV:     {results['PROTEIN_ROTA_FILE']}")
    print(f"3. RNA Conformers CSV:         {results['RNA_CONF_FILE']}")
    print(f"4. Mapped Interactions CSV:    {results['INTERACTIONS_FILE']}")
    print(f"5. Final Merged Output CSV:    {results['FINAL_OUTPUT_FILE']}")

    # 4. Verify Final Output Schema and Data
    df_final = pd.read_csv(results['FINAL_OUTPUT_FILE'])
    print(f"\nFinal Merged Output Dataset Shape: {df_final.shape}")
    print(f"Target Columns: {list(df_final.columns)}")
    print("\nFirst 5 Rows of Final Output Dataset:")
    print(df_final.head())

    assert not df_final.empty, "Error: Final merged dataframe is empty!"
    assert 'b_rotastate' in df_final.columns, "Error: b_rotastate column missing!"
    assert 'NtC' in df_final.columns, "Error: NtC column missing!"

    print("\nSUCCESS: All tests passed cleanly!")


if __name__ == "__main__":
    run_test()
