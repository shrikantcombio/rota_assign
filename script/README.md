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
├── prince_subzone.py          # Core / Rim Interface Sub-zoning Assigner (PRince dASA)
├── merge_pipeline.py          # End-to-end master merging pipeline
└── README.md                  # Package documentation
```

---

## Core Modules & API Usage

### 1. `build_rotalib.py`
Builds rotamer library dictionaries from Secondary Structure Dependent (`SSD`), Backbone Dependent (`BBD`), and Backbone Independent (`BBI`) library files (supporting `.json`, `.csv`, `.tsv`).
- Caches pre-built JSON dictionary files under `rota_lib_dict/{bbd,bbi,ssd}/` (v1 libraries) or `rota_lib_dict/v2/{bbd,bbi,ssd}/` (v2, default; `ROTA_ASSIGN_LIBRARY_VERSION`).
- Includes `should_rebuild_library()` check & user override prompt logic when existing JSON dictionary files are found.

**Usage:**
```python
from script.build_rotalib import BBDRotamerLibBuilder, run_bbd_files, run_bbi_files

# 1. Load cached JSON rotamer dictionaries (or build if missing)
bbd_b_I, bbd_b_N, bbd_ub_I, bbd_ub_N = run_bbd_files(force_rebuild=False)
bbi_b_I, bbi_b_N, bbi_ub_I, bbi_ub_N = run_bbi_files(force_rebuild=False)

# 2. Build library directly from a JSON input file
builder = BBDRotamerLibBuilder()
full_dict, on_dict, off_dict, _ = builder.build_library("RBPs_rotamer_lib/RBPs_BBD_rotamer_lib/RBPs_I_bbd_B_rotamer_lib.json")
```

---

### 2. `get_rotastate.py`
Determines the rotamer state for an amino acid residue using SCWRL4 / MODELLER style rotamer matching over $n - 1$ primary dihedrals.
- Supports library types: `'BBD'` (default), `'BBI'`, or `'SSD'`.
- Applies symmetric terminal $\chi$ folding (ASP $\chi_2$, GLU $\chi_3$, PHE $\chi_2$, TYR $\chi_2 \rightarrow [-90^\circ, 90^\circ)$).
- Enforces a $60^\circ$ grid radius distance cap ($\text{dist}^2 \le 3600$) to prevent distant Ramachandran grid cell fallbacks.

**Usage:**
```python
from script.get_rotastate import get_aa_rotamer

# Backbone-Dependent (BBD) Assignment
bbd_state = get_aa_rotamer(
    aa='ARG', ss='L', b_ub='bound', interface='I',
    chi_1=-60.0, chi_2=175.0, chi_3=180.0, chi_4=-80.0,
    phi=-60.0, psi=-30.0, lib_type='BBD'
)
print("BBD Rotamer State:", bbd_state)  # Output: 'mtt-79'

# Backbone-Independent (BBI) Assignment
bbi_state = get_aa_rotamer(
    aa='ARG', ss='L', b_ub='bound', interface='I',
    chi_1=-60.0, chi_2=175.0, chi_3=180.0, chi_4=-80.0,
    lib_type='BBI'
)
print("BBI Rotamer State:", bbi_state)  # Output: 'mtt-83'
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

### 6. `prince_subzone.py`
Assigns **Core Interface** ($\Delta\text{ASA} \ge 15.0\ \text{\AA}^2$) vs **Rim Interface** ($0 < \Delta\text{ASA} < 15.0\ \text{\AA}^2$) vs **Non-Interface** status using per-atom buried surface area ($\Delta\text{ASA}$) extracted from PRince `*.int` / `*.intnew` files.
- Literature basis: Chakrabarti & Janin (2002), Levy (2010), Jones & Thornton (1996).
- Supports threshold sensitivity analysis across 10.0, 15.0, 20.0, 25.0, and 30.0 $\text{\AA}^2$ cutoffs.

**Usage:**
```python
from script.prince_subzone import load_prince_subzones, get_subzone_label, compute_sensitivity_analysis

# Load PRince subzone lookup map
subzone_map = load_prince_subzones("path/to/1ASYP.int", cutoff=15.0)

# Classify residue
label = get_subzone_label(subzone_map, chain="A", resnum="119", sasa_flag="I")

# Run cutoff sensitivity analysis
df_sensitivity = compute_sensitivity_analysis("path/to/1ASYP.int")
```

---

### 7. `merge_pipeline.py` & Test Script (`test/rota_conf_assign.py`)
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
