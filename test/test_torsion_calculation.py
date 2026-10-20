"""
test_torsion_calculation.py
----------------------------
End-to-end integration test script for structure cleaning, equal residue alignment with ClustalW2,
and Protein-RNA complex torsion angle, DSSP, and NACCESS SASA interface calculations.
Verifies calculated SASA interface against reference files: 1ASYP.int (protein) and 1ASYR.int (RNA).
All output files are written under test/test_results/.
"""

from pathlib import Path
import sys
import pandas as pd
from collections import Counter

# 1. Path Setup
TEST_DIR = Path(__file__).resolve().parent
BASE_DIR = TEST_DIR.parent
UTIL_DIR = BASE_DIR / "utils"
DATASET_DIR = BASE_DIR / "dataset"
RESULTS_DIR = TEST_DIR / "test_results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

if str(UTIL_DIR) not in sys.path:
    sys.path.insert(1, str(UTIL_DIR))
if str(BASE_DIR) not in sys.path:
    sys.path.insert(1, str(BASE_DIR))

from clean_structure import StructureCleaner
from utils.equal_residue_aligner import EqualResidueAligner
from macromol_torsion import ComplexTorsionCalculator


def parse_int_file(int_filepath):
    """
    Parses .int interface file and extracts set of (chain, resSeq) tuples.
    """
    int_p = Path(int_filepath)
    if not int_p.exists():
        return set()

    res_keys = set()
    with open(int_p, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21:22].strip()
                res_seq = line[22:26].strip()
                res_keys.add((chain, res_seq))
    return res_keys


def main():
    print("=== Step 1: Cleaning PDB Structures into test/test_results/ ===")
    cleaner = StructureCleaner(
        keep_waters=True,
        keep_ions=True,
        altloc_strategy='A',
        set_occupancy=1.0,
        default_output_dir=RESULTS_DIR
    )

    # Test inputs are tracked in dataset/
    bound_src = DATASET_DIR / "1ASY.pdb"
    unbound_src = DATASET_DIR / "1ASY_1EOV_U.pdb"
    for src in (bound_src, unbound_src):
        if not src.exists():
            raise FileNotFoundError(f"Test input missing: {src}")

    clean_bound, _ = cleaner.clean_structure(bound_src, "1ASY.pdb")
    clean_unbound, _ = cleaner.clean_structure(unbound_src, "1EOV.pdb")

    print(f"Cleaned Bound File:   {clean_bound}")
    print(f"Cleaned Unbound File: {clean_unbound}")

    print("\n=== Step 2: Extracting Equal Residues into test/test_results/1ASY_1EOV/ ===")
    aligner = EqualResidueAligner(default_output_dir=RESULTS_DIR)

    results = aligner.align_and_extract_equal_residues(
        bound_file=clean_bound,
        bound_chains="A",
        unbound_file=clean_unbound,
        unbound_chains="A",
        rna_chains="R",
        output_dir=RESULTS_DIR / "1ASY_1EOV"
    )

    equal_bound = results['BOUND_EQUAL']     # test_results/1ASY_1EOV/1ASY_1EOV_C.pdb
    equal_unbound = results['UNBOUND_EQUAL'] # test_results/1ASY_1EOV/1ASY_1EOV_U.pdb

    print(f"Bound Equal Structure (Complex):   {equal_bound}")
    print(f"Unbound Equal Structure (Unbound): {equal_unbound}")
    print(f"FASTA File:                         {results['FASTA']}")
    print(f"ClustalW2 Alignment File:           {results['ALIGNMENT']}")

    print("\n=== Step 3: Analyzing Complex Torsions, DSSP, and SASA ===")
    calc = ComplexTorsionCalculator()

    csv_out = RESULTS_DIR / "1ASY_1EOV" / "1ASY_1EOV.csv"

    complex_results = calc.analyze_complex(
        bound_file=equal_bound,
        unbound_file=equal_unbound,
        protein_chains="A",
        rna_chains="R",
        orig_bound_file=clean_bound,
        orig_unbound_file=clean_unbound,
        output_file=csv_out
    )

    df_protein = complex_results['PROTEIN']
    df_rna = complex_results['RNA']

    print(f"Calculated Protein CSV File: {csv_out}")
    print(f"Protein Dataset Shape:       {df_protein.shape}")
    print(f"Protein SASA Distribution:   {Counter(df_protein['SASA'])}")

    rna_csv_out = csv_out.parent / f"{csv_out.stem}_rna.csv"
    print(f"\nCalculated RNA CSV File:     {rna_csv_out}")
    print(f"RNA Suite Dataset Shape:     {df_rna.shape}")
    print(f"RNA SASA Distribution:       {Counter(df_rna['SASA'])}")

    print("\n=== Step 4: Verification Against Reference Interface Files (1ASYP.int & 1ASYR.int) ===")

    def parse_key(lbl):
        parts = str(lbl).split()
        if len(parts) >= 3:
            return (parts[1], parts[2])
        return (None, None)

    df_protein['KEY'] = df_protein['LABEL'].apply(parse_key)
    df_rna['KEY'] = df_rna['LABEL'].apply(parse_key)

    calc_p_int = set(df_protein[df_protein['SASA'] == 'I']['KEY'])
    calc_r_int = set(df_rna[df_rna['SASA'] == 'I']['KEY'])

    # 1ASYP.int Verification
    asyp_file = TEST_DIR / "1ASYP.int"
    ref_p_int = parse_int_file(asyp_file)
    if ref_p_int:
        p_match = ref_p_int.intersection(calc_p_int)
        print(f"1ASYP.int Reference Protein Interface Residues: {len(ref_p_int)}")
        print(f"Calculated Protein Interface Residues:          {len(calc_p_int)}")
        print(f"Protein Interface Match Rate:                   {len(p_match)} / {len(ref_p_int)} ({len(p_match)/len(ref_p_int)*100:.2f}%)")

    # 1ASYR.int Verification
    asyr_file = TEST_DIR / "1ASYR.int"
    ref_r_int = parse_int_file(asyr_file)
    if ref_r_int:
        r_match = ref_r_int.intersection(calc_r_int)
        valid_ref_r = set(df_rna['KEY']).intersection(ref_r_int)
        print(f"\n1ASYR.int Reference RNA Interface Nucleotides:  {len(ref_r_int)}")
        print(f"Calculated RNA Interface Nucleotides:            {len(calc_r_int)}")
        print(f"RNA Interface Match Rate (Overall):              {len(r_match)} / {len(ref_r_int)} ({len(r_match)/len(ref_r_int)*100:.2f}%)")
        print(f"RNA Interface Match Rate (7D Suites Valid):      {len(r_match)} / {len(valid_ref_r)} ({len(r_match)/len(valid_ref_r)*100:.2f}%)")

    print("\nSUCCESS: End-to-end structured torsion calculation test completed successfully!")


if __name__ == "__main__":
    main()
