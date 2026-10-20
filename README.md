# Protein-RNA Complex Rotamer & RNA Conformer Assignment Framework (`rota_assign`)

A comprehensive, modular Python framework for structural analysis of Protein-RNA complexes. The framework automates end-to-end processing—from raw PDB structure cleaning and sequence alignment to protein backbone/sidechain dihedral calculation, Backbone Dependent (BBD) rotamer assignment, StructAnnotate RNA dinucleotide conformer annotation, NAxplorer interaction mapping, and Arginine fork / rotamer shift analysis.

Structure file (CIF/PDB)
        ↓
ComplexTorsionCalculator (macromol_torsion)
        ├── protein_torsions.csv  (φ,ψ,ω, χ1–χ5, SS, SASA)
        └── rna_7d_suite           (α,β,γ,δ,ε,ζ,χ, η,θ,η',θ', pseudorotation)
        ↓                               ↓
ProteinRotamerAssigner          RNAConformerAssigner
(get_rotastate → find_best_bbd)  (StructAnnotate confal → NtC/CANA/confalA)
        ↓                               ↓
        └──────────── merge_pipeline ───┘
                              ↓
              interactions file (NAxplorer CSV)
                              ↓
         map_complex_interactions_and_conformers()
                              ↓
              29-column final output CSV
                              ↓
              run_full_rotastate_analysis()
               ├── ArginineForkAnalyzer
               ├── RotamerTransitionAnalyzer
               └── RotamerConformerPreferenceAnalyzer
---

## 🛠️ System Requirements & External Dependencies

Running `rota_assign` requires Python 3.9+, core scientific Python libraries, and three external bioinformatic system binaries for SASA calculation, secondary structure assignment, and sequence alignment.

### 1. External System Binaries (Must be installed & accessible in `$PATH` or configured path)

| Tool / Binary | Version | Purpose in Pipeline | Detection / Default Paths |
| :--- | :--- | :--- | :--- |
| **NACCESS** (`naccess`) | `2.1.1` | Calculates Solvent Accessible Surface Area (SASA) and classifies residues into Interface (`I`) vs. Non-Interface (`N`) states based on $\Delta\text{SASA} \ge 1.0\ \text{Å}^2$. | `$NACCESS_BIN`, else `naccess` on `$PATH`. |
| **DSSP** (`mkdssp` / `dssp`) | `3.x` / `4.x` | Assigns protein 3-state secondary structure (`H`: Helix, `S`: Sheet, `L`: Loop). | `$DSSP_BIN`, else `mkdssp` / `dssp` on `$PATH`. |
| **ClustalW2** (`clustalw2` / `clustalw`) | `2.1` | Performs pairwise sequence alignment between bound complex and unbound free protein structures to extract equal residue length datasets. | `$CLUSTALW2_BIN`, else `clustalw2` / `clustalw` on `$PATH`. |

---

### 2. Core Python Dependencies

| Package | Minimum Version | Purpose |
| :--- | :--- | :--- |
| `python` | `>= 3.9` | Execution environment (tested on Python 3.11). |
| `biopython` | `>= 1.79` | Parsing PDB and mmCIF structure files (`Bio.PDB`). |
| `pandas` | `>= 1.3.0` | Data manipulation, CSV processing, dataframe merging. |
| `numpy` | `>= 1.21.0` | Vectorized angular math, dihedral angle calculations, distance metrics. |
| `scipy` | `>= 1.7.0` | Statistical testing (Mann-Whitney U, $\chi^2$ contingency, matrix ops). |
| `matplotlib` | `>= 3.4.0` | Visualization and figure generation. |
| `seaborn` | `>= 0.11.0` | Statistical plot formatting and distribution plots. |

---

### 3. Bundled Package Modules & Reference Libraries

The repository includes pre-packaged reference libraries and third-party tools:
- **StructAnnotate** (`StructAnnotate/`, git submodule of https://github.com/Eredtywe/StructAnnotate): RNA 7D suite dinucleotide conformer assignment (`NtC`, `CANA`, `confalA`).
- **Rotamer libraries** (`library/RBPs_{BBD,BBI,SSD}_rotamer_lib/`): copies of the RBPs_rotamer_lib releases. `script/build_rotalib.py` caches them as dictionaries in `rota_lib_dict/` and rebuilds the cache whenever the source files change (SHA-256 in `rota_lib_dict/<type>/source.sha256`). Set `ROTA_ASSIGN_LIBRARY_DIR` to use another library directory with the same layout.
- **Structure Cleaner** (`utils/structure_cleaner.py`, imported via `utils/clean_structure.py`): PDB/mmCIF cleaner (modified residues, altlocs, waters/ions).

Optional environment variables: `ROTA_ASSIGN_PRORNA_DIR` (protein-RNA working data for `merge_pipeline.py` and the NtC scripts; default `dataset/`) and `ROTA_ASSIGN_DSSR_CACHE` (directory for cached x3dna-dssr JSON; no caching if unset).

---

## 💻 Installation & Environment Setup

### Option A: Using Conda (Recommended)

1. **Clone the Repository**:
   ```bash
   git clone --recurse-submodules git@github.com:shrikantcombio/rota_assign.git
   cd rota_assign
   ```

2. **Create and Activate Conda Environment**:
   ```bash
   conda env create -f environment.yml
   conda activate rota_assign
   ```

3. **Verify External Binaries**:
   Ensure `naccess`, `dssp` (or `mkdssp`), and `clustalw2` are executable:
   ```bash
   naccess -help
   dssp --version
   clustalw2 -help
   ```

### Option B: Using Pip

```bash
pip install -r requirements.txt
```

---

## 📁 Repository Architecture

```text
rota_assign/
├── dataset/                        # Input PDB structures, interaction files, & results
│   ├── 1ASY.pdb                    # Cleaned bound protein-RNA complex structure
│   ├── 1EOV.pdb                    # Cleaned free/unbound protein structure
│   ├── 1ASY_1EOV_C.pdb             # Sequence-aligned bound complex equal structure
│   ├── 1ASY_1EOV_U.pdb             # Sequence-aligned unbound equal structure
│   └── results/                    # Generated intermediate & final output files
├── macromol_torsion/               # Core torsion calculation engine (Protein & RNA)
│   ├── calculator.py               # ComplexTorsionCalculator master class
│   ├── protein_backbone.py         # Mainchain (Phi, Psi, Omega) dihedrals
│   ├── protein_sidechain.py        # Sidechain (Chi1..Chi5) dihedrals
│   ├── rna_backbone.py             # RNA 7D suite dihedrals (Alpha..Zeta, Delta, Chi)
│   ├── dssp_runner.py              # Secondary Structure assignment (H, S, L)
│   └── naccess_runner.py           # SASA & Interface status assignment (I vs N)
├── rota_lib_dict/                  # Cached JSON Rotamer Library Dictionaries
│   ├── bbd/                        # Pre-built Backbone Dependent JSON dictionaries
│   └── bbi/                        # Pre-built Backbone Independent JSON dictionaries
├── script/                         # High-level pipeline & analysis modules
│   ├── build_rotalib.py            # Rotamer library dict builder (SSD, BBD, BBI) & JSON cache
│   ├── get_rotastate.py            # Rotamer state calculator & SCWRL4 style matching
│   ├── rotaassign.py               # ProteinRotamerAssigner (bound/unbound/both modes)
│   ├── rna_conformer_assigner.py   # StructAnnotate integration (NtC, CANA, confalA)
│   ├── analyze_rotastates.py       # Arginine forks & rotamer shift analyzer
│   ├── prince_subzone.py          # Core / Rim Interface Sub-zoning Assigner (PRince dASA)
│   ├── merge_pipeline.py           # Master end-to-end processing pipeline
│   └── README.md                   # Script package documentation
├── utils/                          # Preprocessing & alignment utilities
│   ├── clean_structure.py          # PDB/CIF StructureCleaner module
│   └── equal_residue_aligner.py    # EqualResidueAligner sequence alignment module
├── library/                        # Rotamer & conformer reference libraries
│   ├── RBPs_BBD_rotamer_lib/       # Backbone Dependent (BBD) rotamer libraries (JSON/CSV)
│   └── RBPs_BBI_rotamer_lib/       # Backbone Independent (BBI) rotamer libraries (JSON/CSV)
├── StructAnnotate/                 # Cloned StructAnnotate package for RNA conformers
├── test/                           # Test suite & verification scripts
│   ├── run_all_tests.py            # Master test runner script
│   ├── test_torsion_calculation.py # Torsion & interface calculation test
│   ├── rota_conf_assign.py         # Rotamer & RNA conformer assignment test
│   └── test_results/               # Directory where test outputs are saved
├── environment.yml                 # Conda environment definition file
├── requirements.txt                # Pip requirements file
└── README.md                       # Main repository documentation
```

---

## 🚀 Quick Start & Full Pipeline Workflow

### Step 1: Structure Cleaning (`utils/clean_structure.py`)
Clean raw PDB/mmCIF structures by removing unwanted heteroatoms, handling alternate locations (altlocs), optionally keeping waters/ions, setting occupancies, and exporting standardized structures:

```python
from pathlib import Path
from utils.clean_structure import StructureCleaner

cleaner = StructureCleaner(
    keep_waters=True,
    keep_ions=True,
    altloc_strategy='A',
    set_occupancy=1.0
)

# Clean bound and unbound structures
out_bound, logs_b = cleaner.clean_structure("dataset/1ASY.cif", "dataset/1ASY.pdb")
out_unbound, logs_u = cleaner.clean_structure("dataset/1EOV.cif", "dataset/1EOV.pdb")
```

---

### Step 2: Equal Residue Sequence Alignment (`utils/equal_residue_aligner.py`)
Align protein chains between bound complex and unbound free protein structures using ClustalW2 to produce pair-aligned, equal-length PDB files (`1ASY_1EOV_C.pdb` and `1ASY_1EOV_U.pdb`):

```python
from utils.equal_residue_aligner import EqualResidueAligner

aligner = EqualResidueAligner()
results = aligner.align_and_extract_equal_residues(
    bound_file="dataset/1ASY.pdb",
    bound_chains="A",
    unbound_file="dataset/1EOV.pdb",
    unbound_chains="A",
    rna_chains="R",
    output_dir="dataset/1ASY_1EOV"
)

# Generated files: dataset/1ASY_1EOV/1ASY_1EOV_C.pdb, dataset/1ASY_1EOV/1ASY_1EOV_U.pdb
```

---

### Step 3: Torsion Angle, DSSP, & NACCESS SASA Calculation (`macromol_torsion/`)
Calculate mainchain ($\phi, \psi, \omega$), sidechain ($\chi_1 \dots \chi_5$), secondary structure (`H`, `S`, `L`), NACCESS SASA interface status (`I` vs `N`), and RNA 7D suite angles:

```python
from pathlib import Path
from macromol_torsion import ComplexTorsionCalculator

calc = ComplexTorsionCalculator()

# Comprehensive complex analysis
results = calc.analyze_complex(
    bound_file=Path("dataset/1ASY_1EOV/1ASY_1EOV_C.pdb"),
    unbound_file=Path("dataset/1ASY_1EOV/1ASY_1EOV_U.pdb"),
    protein_chains="A",
    rna_chains="R",
    orig_bound_file=Path("dataset/1ASY.pdb"),
    orig_unbound_file=Path("dataset/1EOV.pdb"),
    output_file="dataset/results/1ASY_1EOV.csv"
)

df_protein = results['PROTEIN']
df_rna = results['RNA']
```

---

### Step 4: Rotamer & RNA Conformer Assignment (`script/`)
Assign protein rotameric states using **Backbone Dependent (BBD)** rotamer libraries (`library/RBPs_BBD_rotamer_lib/`) and annotate RNA dinucleotide conformers (`NtC`, `CANA`, `confalA`) using **StructAnnotate**:

```python
from script.rotaassign import assign_rotamer_states
from script.rna_conformer_assigner import assign_rna_conformers

# 1. Assign Protein Rotamer States (bound & unbound)
df_protein_assigned = assign_rotamer_states(df_protein, mode='both')

# 2. Assign RNA Dinucleotide Conformers (NtC, CANA, confalA)
df_rna_assigned = assign_rna_conformers(df_rna)
```

**Rotamer names** follow Lovell et al. (2000): every chi up to chi4 is compared with the library, and
rotameric chis are named by letter (`p` ≈ +60°, `t` ≈ 180°, `m` ≈ −60°; capitals `P`/`T`/`M` mark the
off-rotamer bins between them). A terminal sp2 chi (ASN, ASP, GLN, GLU, HIS, PHE, TYR, TRP, and ARG chi4)
is given as the matched library value instead, e.g. ASP `m-20`, GLN `mt-30`; LEU is `tp`, LYS `mttt`.
The assigned name is that of the closest rotamer in the library cell, so a residue whose own
conformation is not populated in that cell receives the nearest populated rotamer.
`get_rotastate.classify_rotamer(aa, chis)` gives the class directly from a residue's own chi angles,
without a library.

**Bound/unbound pairing.** `macromol_torsion.pairing` pairs residues by chain order and sequence
(`pair_equal_residues` for equal-residue files, `pair_residues_across_entries` for different PDB
entries); residue labels are not used because unbound chains are usually relabelled or renumbered.
Angles that cannot be computed (chain termini and breaks, missing atoms, no unbound partner) are NaN.

---

### Step 5: Master Integration Pipeline (`script/merge_pipeline.py`)

Run the master pipeline via CLI for bound-unbound pair structures or single complex mode:

```bash
# Complex + Unbound Protein Mode
python script/merge_pipeline.py --pdb dataset/1ASY.cif --ub_protein dataset/1EOV.cif --interactions dataset/1ASY_rnp_interactions_model1_energy.csv -o dataset/results/

# Complex + Unbound RNA Mode
python script/merge_pipeline.py --pdb dataset/1ASY.cif --ub_RNA dataset/2TRA.cif --interactions dataset/1ASY_rnp_interactions_model1_energy.csv -o dataset/results/
```

---

## 🧪 Testing & Verification

To verify that all dependencies and external binaries are functioning properly, run the test suite:

```bash
# Run master test suite (integration tests + regression unit tests)
python3 test/run_all_tests.py

# Regression unit tests only
python3 test/test_regressions.py
```

Outputs will be saved in `test/test_results/`.

---

## 📜 License & Citation

Developed for PhD Research on Protein-RNA Interactions, Rotameric Transitions, and Conformational State Assignments.
