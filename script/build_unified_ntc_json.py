"""
build_unified_ntc_json.py
-------------------------
Constructs a unified, standardized JSON dataset uniting all RNA dinucleotide conformer (NtC) steps
(both interface 'I' and non-interface 'N') with their structural, torsion, DSSR secondary structure,
and protein interaction details.

Generates:
1. Per-complex JSON: {results_dir}/{PDB}_unified_ntc_interactions.json
2. Dataset-wide Benchmark JSON: {results_dir}/PRDBv35_unified_ntc_interactions.json
"""

import os
import sys
import json
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

# Protein-RNA results directory; set ROTA_ASSIGN_PRORNA_DIR
PRO_RNA_RESULTS_DIR = Path(os.environ.get("ROTA_ASSIGN_PRORNA_DIR", Path(__file__).resolve().parent.parent / "dataset")) / "results"


# Path setup
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# Chemical Class Definitions for Amino Acids
AMINO_ACID_CLASSES = {
    'ARG': 'BASIC_POSITIVE',
    'LYS': 'BASIC_POSITIVE',
    'HIS': 'BASIC_POSITIVE',
    'ASP': 'ACIDIC_NEGATIVE',
    'GLU': 'ACIDIC_NEGATIVE',
    'ASN': 'POLAR',
    'GLN': 'POLAR',
    'SER': 'POLAR',
    'THR': 'POLAR',
    'TYR': 'AROMATIC_POLAR',
    'PHE': 'AROMATIC_HYDROPHOBIC',
    'TRP': 'AROMATIC_HYDROPHOBIC',
    'ALA': 'HYDROPHOBIC_ALIPHATIC',
    'VAL': 'HYDROPHOBIC_ALIPHATIC',
    'LEU': 'HYDROPHOBIC_ALIPHATIC',
    'ILE': 'HYDROPHOBIC_ALIPHATIC',
    'MET': 'HYDROPHOBIC_ALIPHATIC',
    'PRO': 'SPECIAL_PROLINE',
    'GLY': 'SPECIAL_GLYCINE',
    'CYS': 'SPECIAL_CYSTEINE'
}


def get_aa_chemical_class(aa_code):
    """Returns chemical classification for amino acid 3-letter code."""
    aa = str(aa_code).upper().strip()
    return AMINO_ACID_CLASSES.get(aa, 'OTHER')


def safe_float(val, default=0.0):
    """Safely converts value to float or returns default."""
    try:
        if pd.isna(val) or str(val).strip() in ['-', '', 'nan', 'NaN', 'None', '?']:
            return default
        return round(float(val), 4)
    except (ValueError, TypeError):
        return default


def safe_str(val, default=''):
    """Safely converts value to clean string."""
    if pd.isna(val) or str(val).strip() in ['nan', 'NaN', 'None']:
        return default
    return str(val).strip()


def build_complex_unified_json(pdb_id, results_dir):
    """
    Builds a unified dictionary for a single complex combining all RNA dinucleotide conformers
    (both interface 'I' and non-interface 'N') with any associated protein interactions.
    """
    results_path = Path(results_dir)
    clean_pdb = pdb_id.strip().upper()
    rna_conf_file = results_path / f"{clean_pdb}_rna_conformers.csv"
    final_output_file = results_path / f"{clean_pdb}_final_output.csv"

    if not rna_conf_file.exists():
        return None

    df_rna = pd.read_csv(rna_conf_file)
    df_final = pd.read_csv(final_output_file) if final_output_file.exists() else pd.DataFrame()

    # Index interactions by (Nuc_chain, str(Nuc_pos))
    interactions_by_pos = {}
    if not df_final.empty:
        for _, row in df_final.iterrows():
            nuc_chain = safe_str(row.get('Nuc_chain', 'R'))
            nuc_pos = safe_str(row.get('Nuc_pos', '0'))
            key = (nuc_chain, nuc_pos)
            if key not in interactions_by_pos:
                interactions_by_pos[key] = []

            aa = safe_str(row.get('AA', ''))
            aa_class = get_aa_chemical_class(aa)

            int_obj = {
                "protein_residue": {
                    "aa": aa,
                    "chain": safe_str(row.get('AA_chain', 'A')),
                    "pos": safe_str(row.get('AA_pos', '0')),
                    "atom": safe_str(row.get('AA_atom', '')),
                    "phi": safe_float(row.get('b_phi', 0.0)),
                    "psi": safe_float(row.get('b_psi', 0.0)),
                    "ss": safe_str(row.get('b_ss', 'L')),
                    "rotastate": safe_str(row.get('b_rotastate', 'NA')),
                    "chemical_class": aa_class
                },
                "rna_residue": {
                    "nuc": safe_str(row.get('Nuc', '')),
                    "chain": nuc_chain,
                    "pos": nuc_pos,
                    "atom": safe_str(row.get('Nuc_atom', ''))
                },
                "interaction_type": safe_str(row.get('interaction_type', 'H-BOND')),
                "edge": safe_str(row.get('edge', '-')),
                "dist": safe_str(row.get('dist', '-')),
                "score": safe_float(row.get('score', 0.0)),
                "energy_kcalmol": safe_float(row.get('E_int_kcalmol', 0.0)),
                "mediator_info": safe_str(row.get('mediator_info', '-')),
                "status": safe_str(row.get('status', 'ok')),
                "source": safe_str(row.get('source', 'NAxplorer_energy'))
            }
            interactions_by_pos[key].append(int_obj)

    steps = []
    num_interface_steps = 0
    num_non_interface_steps = 0
    total_assigned_interactions = 0

    for idx, row in df_rna.iterrows():
        label_raw = safe_str(row.get('LABEL', ''))
        parts = label_raw.split()
        res_name = parts[0] if len(parts) > 0 else 'N'
        chain = parts[1] if len(parts) > 1 else 'R'
        pos = parts[2] if len(parts) > 2 else str(idx)

        sasa_status = safe_str(row.get('SASA', 'N')).upper()
        if sasa_status not in ['I', 'N']:
            sasa_status = 'N'

        is_interface = (sasa_status == 'I')
        if is_interface:
            num_interface_steps += 1
        else:
            num_non_interface_steps += 1

        # Look up interactions for this nucleotide suite
        step_interactions = interactions_by_pos.get((chain, pos), [])
        total_assigned_interactions += len(step_interactions)

        step_obj = {
            "step_id": f"{clean_pdb}_{chain}_{pos}",
            "step_index": idx,
            "chain": chain,
            "res_pos": pos,
            "res_name": res_name,
            "suite_label": label_raw,
            "sasa": sasa_status,
            "is_interface": is_interface,
            "ntc": safe_str(row.get('NtC', 'NANT')),
            "cana": safe_str(row.get('CANA', 'NAN')),
            "confal_a": safe_float(row.get('confalA', 0.0)),
            "torsions": {
                "delm": safe_float(row.get('DELM', 0.0)),
                "epsm": safe_float(row.get('EPSM', 0.0)),
                "zetm": safe_float(row.get('ZETM', 0.0)),
                "alpha": safe_float(row.get('ALPHA', 0.0)),
                "beta": safe_float(row.get('BETA', 0.0)),
                "gamma": safe_float(row.get('GAMMA', 0.0)),
                "delta": safe_float(row.get('DELTA', 0.0)),
                "chi1": safe_float(row.get('CHI1', 0.0)),
                "chi2": safe_float(row.get('CHI2', 0.0)),
                "nn_dist": safe_float(row.get('NN_DIST', 0.0)),
                "cc_dist": safe_float(row.get('CC_DIST', 0.0)),
                "mu": safe_float(row.get('MU', 0.0))
            },
            "secondary_structure": {
                "ss_dbn": safe_str(row.get('SS_DBN', '.')),
                "ss_type": safe_str(row.get('SS_TYPE', 'unpaired')),
                "paired_with": safe_str(row.get('PAIRED_WITH', '.')),
                "bp_type": safe_str(row.get('BP_TYPE', '.')),
                "dssr_summary": safe_str(row.get('DSSR_SUMMARY', ''))
            },
            "num_interactions": len(step_interactions),
            "interactions": step_interactions
        }
        steps.append(step_obj)

    complex_data = {
        "pdb_id": clean_pdb,
        "num_total_steps": len(steps),
        "num_interface_steps": num_interface_steps,
        "num_non_interface_steps": num_non_interface_steps,
        "num_interactions": total_assigned_interactions,
        "steps": steps
    }

    return complex_data


def run_build_unified_json(results_dir=str(PRO_RNA_RESULTS_DIR), pdb=None):
    """
    Builds per-complex and dataset-wide unified JSON files.
    """
    res_path = Path(results_dir)
    if not res_path.exists():
        raise FileNotFoundError(f"Results directory does not exist: {results_dir}")

    if pdb:
        pdbs = [pdb.strip().upper()]
    else:
        # Find all complexes with _rna_conformers.csv
        pdbs = sorted([f.stem.replace('_rna_conformers', '') for f in res_path.glob('*_rna_conformers.csv')])

    print(f"[INFO] Building unified NtC JSON for {len(pdbs)} complexes in {results_dir}...")

    all_complexes = {}
    tot_steps = 0
    tot_int_steps = 0
    tot_non_int_steps = 0
    tot_interactions = 0

    for i, p in enumerate(pdbs, 1):
        comp_data = build_complex_unified_json(p, res_path)
        if comp_data is None:
            continue

        # Save per-complex JSON
        out_single = res_path / f"{p}_unified_ntc_interactions.json"
        with open(out_single, 'w', encoding='utf-8') as f:
            json.dump(comp_data, f, indent=2)

        all_complexes[p] = comp_data
        tot_steps += comp_data["num_total_steps"]
        tot_int_steps += comp_data["num_interface_steps"]
        tot_non_int_steps += comp_data["num_non_interface_steps"]
        tot_interactions += comp_data["num_interactions"]

        if i % 25 == 0 or i == len(pdbs):
            print(f"  Processed [{i}/{len(pdbs)}] complexes...")

    # Build dataset-wide benchmark JSON if in batch mode
    if not pdb or len(pdbs) > 1:
        benchmark_json_path = res_path / "PRDBv35_unified_ntc_interactions.json"
        benchmark_data = {
            "metadata": {
                "benchmark_name": "PRDBv35_clean",
                "total_complexes": len(all_complexes),
                "total_rna_dinucleotide_steps": tot_steps,
                "total_interface_steps": tot_int_steps,
                "total_non_interface_steps": tot_non_int_steps,
                "total_mapped_interactions": tot_interactions
            },
            "complexes": all_complexes
        }
        with open(benchmark_json_path, 'w', encoding='utf-8') as f:
            json.dump(benchmark_data, f, indent=2)

        print(f"\n[SUCCESS] Dataset-wide benchmark JSON saved: {benchmark_json_path}")
        print(f"  Total Complexes: {len(all_complexes)}")
        print(f"  Total RNA Steps: {tot_steps} (Interface: {tot_int_steps}, Non-Interface: {tot_non_int_steps})")
        print(f"  Total Mapped Interactions: {tot_interactions}")

    return all_complexes


def main():
    parser = argparse.ArgumentParser(description="Construct unified JSON dataset combining interface and non-interface RNA conformers.")
    parser.add_argument("--pdb", default=None, help="Process single PDB ID")
    parser.add_argument("-d", "--dir", default=str(PRO_RNA_RESULTS_DIR), help="Directory containing results CSVs")
    args = parser.parse_args()

    run_build_unified_json(results_dir=args.dir, pdb=args.pdb)


if __name__ == "__main__":
    main()
