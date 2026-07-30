# Protein-RNA Complex Rotamer & RNA Conformer Assignment Framework (`rota_assign`)

A comprehensive, modular Python framework for structural analysis of Protein-RNA complexes. The framework automates end-to-end processing—from raw PDB structure cleaning and sequence alignment to protein backbone/sidechain dihedral calculation, Backbone Dependent (BBD) rotamer assignment, StructAnnotate RNA dinucleotide conformer annotation, NAxplorer interaction mapping, and Arginine fork / rotamer shift analysis.

---

## 📁 Repository Architecture

```text
rota_assign/
├── dataset/                        # Input PDB structures, interaction files, & results
│   ├── 1ASY.pdb                    # Cleaned bound protein-RNA complex structure
│   ├── 1EOV.pdb                    # Cleaned free/unbound protein structure
│   ├── 1ASY_1EOV_C.pdb             # Sequence-aligned bound complex equal structure
│   ├── 1ASY_1EOV_U.pdb             # Sequence-aligned unbound equal structure
│   ├── 1ASY_interactions_NAxplorer.csv # NAxplorer interaction dataset
│   └── results/                    # Generated intermediate & final output files
│       ├── 1ASY_protein_torsions.csv
│       ├── 1ASY_protein_rotastates.csv
│       ├── 1ASY_rna_conformers.csv
│       ├── 1ASY_interactions.csv
│       ├── 1ASY_final_output.csv
│       ├── 1ASY_arg_forks_analysis.csv
│       └── 1ASY_rotamer_transitions.csv
├── macromol_torsion/               # Core torsion calculation engine (Protein & RNA)
│   ├── calculator.py               # ComplexTorsionCalculator master class
│   ├── protein_backbone.py         # Mainchain (Phi, Psi, Omega) dihedrals
│   ├── protein_sidechain.py        # Sidechain (Chi1..Chi5) dihedrals
│   ├── rna_backbone.py             # RNA 7D suite dihedrals (Alpha..Zeta, Delta, Chi)
│   ├── dssp_runner.py              # Secondary Structure assignment (H, S, L)
│   └── naccess_runner.py           # SASA & Interface status assignment (I vs N)
├── script/                         # High-level pipeline & analysis modules
│   ├── build_rotalib.py            # Rotamer library dict builder (SSD, BBD, BBI)
│   ├── get_rotastate.py            # Rotamer state calculator & BBD nearest grid lookup
│   ├── rotaassign.py               # ProteinRotamerAssigner (bound/unbound/both modes)
│   ├── rna_conformer_assigner.py   # StructAnnotate integration (NtC, CANA, confalA)
│   ├── analyze_rotastates.py       # Arginine forks & rotamer shift analyzer
│   ├── merge_pipeline.py           # Master end-to-end processing pipeline
│   └── README.md                   # Script package documentation
├── utils/                          # Preprocessing & alignment utilities
│   ├── clean_structure.py          # PDB/CIF StructureCleaner module
│   └── equal_residue_aligner.py    # EqualResidueAligner sequence alignment module
├── library/                        # Rotamer & conformer reference libraries
│   └── RBPs_BBD_rotamer_lib/       # Backbone Dependent (BBD) rotamer libraries
├── rota_lib_dict/                  # Pre-compiled JSON rotamer dictionaries
├── StructAnnotate/                 # Cloned StructAnnotate package for RNA conformers
├── test/                           # Test suite & verification scripts
│   └── rota_conf_assign.py         # End-to-end integration test script
├── environment.yml                 # Conda environment definition
└── README.md                       # Main repository documentation
```

---

## 🚀 Quick Start & Full Pipeline Workflow

### Step 1: Structure Cleaning (`utils/clean_structure.py`)
Clean raw PDB/mmCIF structures by removing unwanted heteroatoms, handling alternate locations (altlocs), optionally keeping waters/ions, setting occupancies, and exporting standardized structures into `dataset/`:

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
out_bound, logs_b = cleaner.clean_structure("path/to/raw_1ASY.pdb", "dataset/1ASY.pdb")
out_unbound, logs_u = cleaner.clean_structure("path/to/raw_1EOV.pdb", "dataset/1EOV.pdb")
```

---

### Step 2: Equal Residue Sequence Alignment (`utils/equal_residue_aligner.py`)
Align protein chains between bound complex and unbound free protein structures to produce pair-aligned, equal-length PDB files (`1ASY_1EOV_C.pdb` and `1ASY_1EOV_U.pdb`):

```python
from utils.equal_residue_aligner import EqualResidueAligner

aligner = EqualResidueAligner()
results = aligner.align_and_extract_equal_residues(
    bound_file="dataset/1ASY.pdb",
    bound_chains="A",
    unbound_file="dataset/1EOV.pdb",
    unbound_chains="A",
    rna_chains="R",
    output_dir="dataset/"
)

# Generated files: dataset/1ASY_1EOV_C.pdb, dataset/1ASY_1EOV_U.pdb
```

---

### Step 3: Torsion Angle, DSSP, & NACCESS SASA Calculation (`macromol_torsion/`)
Calculate mainchain ($\phi, \psi, \omega$), sidechain ($\chi_1 \dots \chi_5$), secondary structure (`H`, `S`, `L`), NACCESS Solvent Accessible Surface Area (SASA) interface status (`I` vs `N`), and RNA 7D suite angles:

```python
from pathlib import Path
from macromol_torsion import ComplexTorsionCalculator

calc = ComplexTorsionCalculator()

# Comprehensive complex analysis
results = calc.analyze_complex(
    bound_file=Path("dataset/1ASY_1EOV_C.pdb"),
    unbound_file=Path("dataset/1ASY_1EOV_U.pdb"),
    protein_chains="A",
    rna_chains="R",
    orig_bound_file=Path("dataset/1ASY.pdb"),
    orig_unbound_file=Path("dataset/1EOV.pdb"),
    output_file="dataset/results/pairs/1ASY_1EOV.csv"
)

df_protein = results['PROTEIN']
df_rna = results['RNA']
```

For detailed documentation on the torsion calculation module, see **[`macromol_torsion/README.md`](file:///home/labuser/Projects/PhD_projects/rota_assign/macromol_torsion/README.md)** (if available) or module docstrings.

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

For comprehensive module API documentation, see **[`script/README.md`](file:///home/labuser/Projects/PhD_projects/rota_assign/script/README.md)**.

---

### Step 5: Interaction Mapping & Master Merging Pipeline (`script/merge_pipeline.py`)
Run the master pipeline to integrate all intermediate outputs and map NAxplorer interaction entries (`dataset/1ASY_interactions_NAxplorer.csv`) into the final merged dataset (`dataset/results/1ASY_final_output.csv`):

```bash
# Execute master pipeline via CLI
python script/merge_pipeline.py --pdb 1ASY
```

**Or execute full integration test script:**
```bash
python test/rota_conf_assign.py
```

---

### Step 6: Arginine Forks & Rotamer Transition Shift Analysis (`script/analyze_rotastates.py`)
Analyze bidentate Arginine forks, RNA target moieties, and rotamer preference shifts upon RNA binding (unbound to bound) categorized by amino acid chemical groups:

```bash
python script/analyze_rotastates.py --input dataset/results/1ASY_final_output.csv -o dataset/results/
```

---

## 📊 Final Dataset Output Schema (`dataset/results/1ASY_final_output.csv`)

The final output dataset contains 17 standardized columns:

| Column Name | Type | Description |
| :--- | :--- | :--- |
| `PDB_ID` | `str` | PDB identifier (e.g. `1ASY`) |
| `AA` | `str` | 3-letter amino acid code (e.g. `ASN`, `LYS`) |
| `AA_atom` | `str` | Participating protein atom (e.g. `ND2`, `NZ`, `N`) |
| `AA_chain` | `str` | Protein chain identifier (e.g. `A`) |
| `AA_pos` | `int` | Protein residue position |
| `Nuc` | `str` | 1-letter nucleotide code (e.g. `U`, `A`, `G`, `C`) |
| `Nuc_atom` | `str` | Participating RNA atom (e.g. `OP1`, `O2'`, `N3`) |
| `Nuc_chain` | `str` | RNA chain identifier (e.g. `R`) |
| `Nuc_pos` | `int` | RNA nucleotide position |
| `hb_dis` | `float` | Hydrogen bond / interaction distance (Å) |
| `NtC` | `str` | StructAnnotate NtC conformer cluster code (e.g. `AA08`, `AA00`) |
| `CANA` | `str` | StructAnnotate CANA category string (e.g. `AAA`, `A-B`) |
| `confalA` | `float` | CONFAL score (0.0 to 100.0) |
| `b_phi` | `float` | Bound protein backbone Phi angle ($\phi$) |
| `b_psi` | `float` | Bound protein backbone Psi angle ($\psi$) |
| `b_ss` | `str` | DSSP secondary structure code (`H`: Helix, `S`: Sheet, `L`: Loop) |
| `b_rotastate` | `str` | Assigned BBD rotamer state (e.g. `t17`, `t-152`, `M145`, `MA`) |

---

## 🧪 Testing & Verification

To run unit and integration test scripts:

```bash
# Run master integration test script
python test/rota_conf_assign.py
```

---

## 📄 References & Further Documentation

- **Script Package Documentation**: [`script/README.md`](file:///home/labuser/Projects/PhD_projects/rota_assign/script/README.md)
- **StructAnnotate Package**: [`StructAnnotate/`](file:///home/labuser/Projects/PhD_projects/rota_assign/StructAnnotate/)
- **Backbone Dependent Rotamer Libraries**: [`library/RBPs_BBD_rotamer_lib/`](file:///home/labuser/Projects/PhD_projects/rota_assign/library/RBPs_BBD_rotamer_lib/)

---

## 📜 License & Citation

Developed for PhD Research on Protein-RNA Interactions and Rotameric/Conformational State Assignments.
