# Rotamer & RNA Conformer Assignment Scripts Package

This directory contains the object-oriented Python modules for building rotamer libraries, assigning amino acid rotamer states, annotating RNA dinucleotide conformers, analyzing Arginine forks & rotamer preference shifts, mapping NAxplorer interactions, and merging Protein-RNA interactions into structured output datasets under `dataset/results/`.

---

## Package Architecture

```text
script/
├── __init__.py                # Package initialization
├── build_rotalib.py           # Object-Oriented Rotamer Library Builder (SSD, BBD, BBI)
├── get_rotastate.py           # Core Rotamer State Assigner (get_aa_rotamer & BBD lookup)
├── rotaassign.py              # Batch Protein Rotamer State Assigner (bound/unbound/both)
├── rna_conformer_assigner.py  # RNA Dinucleotide Conformer Assigner (StructAnnotate integration)
├── analyze_rotastates.py      # Arginine Forks & Rotamer Transition Preference Analyzer
├── merge_pipeline.py          # End-to-end master merging pipeline
└── README.md                  # Package documentation
```

---

## Core Modules & API Usage

### 1. `build_rotalib.py`
Builds rotamer library dictionaries from Secondary Structure Dependent (SSD), Backbone Dependent (BBD), and Backbone Independent (BBI) library files.

**Usage:**
```python
from script.build_rotalib import SSDRotamerLibBuilder, RotamerLibraryPipeline

# Process a single SSD library file
builder = SSDRotamerLibBuilder()
full_df, on_df, off_df, prob_arr = builder.get_rota_stats("path/to/lib.txt")
lib_dict = builder.build_lib_dict(on_df)
builder.export_json(lib_dict, "path/to/output.json")
```

---

### 2. `get_rotastate.py`
Determines the rotamer state for an amino acid residue based on its chi angles, backbone $(\phi, \psi)$ dihedrals (BBD library lookup), secondary structure (`b_ss`), bound/unbound state (`b_ub`), and interface status (`I` vs `N`).

**Usage:**
```python
from script.get_rotastate import get_aa_rotamer

rotastate = get_aa_rotamer(
    aa='ASN',
    ss='L',
    b_ub='bound',
    interface='I',
    chi_1=-74.42,
    chi_2=101.18,
    phi=74.42,
    psi=-101.18
)
print("Assigned Rotastate:", rotastate)  # Output: 't17'
```

---

### 3. `rotaassign.py`
Batch assigns rotamer states to protein dataset DataFrames or CSV files using Backbone Dependent (BBD) rotamer libraries (`library/RBPs_BBD_rotamer_lib`). Supports modes: `'bound'`, `'unbound'`, or `'both'`.

**Usage:**
```python
from script.rotaassign import assign_rotamer_states

# Assign rotamer states to a CSV dataset file
df_assigned = assign_rotamer_states(
    input_target="dataset/1ASY_1EOV/1ASY_1EOV.csv",
    mode="both",
    output_file="dataset/results/1ASY_protein_rotastates.csv"
)
```

---

### 4. `rna_conformer_assigner.py`
Integrates with **StructAnnotate** to calculate `NtC` (conformer cluster code), `CANA` (conformer category), and `confalA` (CONFAL score) for RNA 7D suite dinucleotide steps.

**Usage:**
```python
from script.rna_conformer_assigner import assign_rna_conformers

df_rna_assigned = assign_rna_conformers(
    input_target="dataset/1ASY_1EOV/1ASY_1EOV_rna.csv",
    output_file="dataset/results/1ASY_rna_conformers.csv"
)
```

---

### 5. `analyze_rotastates.py`
Analyzes Arginine forks (bidentate/multi-contact interactions with RNA), amino acid rotamer shifts upon RNA binding (unbound to bound), chemical class group preferences (Charged, Polar, Hydrophobic, Aromatic), and protein rotamer vs RNA conformer cross-tabulations.

**Usage:**
```python
from script.analyze_rotastates import run_full_rotastate_analysis

results = run_full_rotastate_analysis(
    input_csv="dataset/results/1ASY_final_output.csv",
    output_dir="dataset/results/"
)
```

---

### 6. `merge_pipeline.py` & Test Script (`test/rota_conf_assign.py`)
Master pipeline executing end-to-end processing for a Protein-RNA complex pair:
1. Calculates protein & RNA dihedrals, DSSP, and NACCESS SASA ($\rightarrow$ `dataset/results/1ASY_protein_torsions.csv`).
2. Assigns protein rotamer states using BBD library ($\rightarrow$ `dataset/results/1ASY_protein_rotastates.csv`).
3. Assigns RNA dinucleotide conformers using StructAnnotate ($\rightarrow$ `dataset/results/1ASY_rna_conformers.csv`).
4. Maps interaction entries from `dataset/1ASY_interactions_NAxplorer.csv` ($\rightarrow$ `dataset/results/1ASY_interactions.csv`).
5. Merges all fields into **`dataset/results/1ASY_final_output.csv`**.
6. Runs Arginine forks & rotamer transition analysis reports.

**Usage via Command Line:**
```bash
python test/rota_conf_assign.py
```
