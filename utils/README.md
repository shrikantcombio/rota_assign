# `rota_assign/utils` Utility Scripts & Package Documentation

The `utils/` directory contains tools and modules for structure parsing, cleaning, residue standardization, mmCIF processing, sequence alignment, and equal-residue dataset curation for protein and RNA structures.

---

## Directory Overview

| Script / Directory | Primary Purpose |
| :--- | :--- |
| [`equal_residue_aligner.py`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/equal_residue_aligner.py) | **New**: Aligns bound/unbound structure pairs & extracts equal-residue datasets using ClustalW2. |
| [`clean_structure.py`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/clean_structure.py) | Main unified structure cleaner for PDB and mmCIF files. |
| [`mmcif_clean_reader.py`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/mmcif_clean_reader.py) | Parser for PDB and mmCIF files into standardized atom dictionaries. |
| [`renaming.py`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/renaming.py) | Normalizes mmCIF chain (`auth_asym_id`) and residue numbering (`auth_seq_id`). |
| [`pdb_clean.py`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/pdb_clean.py) | Backward-compatible wrapper for occupancy resolution & residue cleaning. |
| [`pdb_HETATM.py`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/pdb_HETATM.py) | Backward-compatible wrapper for converting modified HETATMs to ATOM records. |
| [`pdb_module.py`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/pdb_module.py) | Legacy backward-compatible wrapper calling `EqualResidueAligner`. |
| [`addtional/`](file:///home/labuser/Projects/PhD_projects/rota_assign/utils/addtional) | Subdirectory containing atomic filtering and occupancy utilities from `pdb-tools`. |

---

## Detailed Script Descriptions

### 1. `equal_residue_aligner.py`
**Bound & Unbound Structure Equal-Residue Aligner**

Replaces `pdb_module.py`. Class `EqualResidueAligner` aligns sequence pairs between bound (complex) and unbound (free) structure files using ClustalW2, filters identical/conserved residue positions and common atoms, and writes curated equal-residue dataset files.

**Capabilities**:
- Extract sequence dictionary per chain for both PDB and mmCIF files.
- Execute ClustalW2 sequence alignment (with fallback alignment).
- Match identical (`*`) and similar (`:`, `.`) residue positions.
- Remove non-common extra atoms between bound and unbound residue pairs.
- Output curated structure pairs into designated dataset directories (`equal_residue_dataset/`).

#### Usage Example:
```python
from utils.equal_residue_aligner import EqualResidueAligner

aligner = EqualResidueAligner(default_output_dir="equal_residue_dataset")

# Align bound and unbound structure pair
results = aligner.align_and_extract_equal_residues(
    bound_file="1ASY.pdb",
    bound_chains="A",
    unbound_file="1EOV.pdb",
    unbound_chains="A",
    output_dir="equal_residue_dataset/1ASY_1EOV"
)

print(results['BOUND_EQUAL'])    # Output path to equal bound structure
print(results['UNBOUND_EQUAL'])  # Output path to equal unbound structure
```

#### CLI Usage:
```bash
python utils/equal_residue_aligner.py --bound 1ASY.pdb --bound-chains A --unbound 1EOV.pdb --unbound-chains A -o equal_residue_dataset/
```

---

### 2. `clean_structure.py`
**Unified Structure Cleaner for PDB & mmCIF Formats**

Exposes `StructureCleaner` and `ModResMap` imported from `/home/labuser/clean_pdb/clean_structure.py`.

---

### 3. `mmcif_clean_reader.py`
**Fast mmCIF & PDB Tokenizer/Parser**

Parses `.pdb` and `.cif`/`.mmcif` files without requiring external software like BioPython. Returns column headers and dictionary rows per atom record.

---

### 4. `renaming.py`
**mmCIF Identifier Standardization Tool**

Standardizes mmCIF files by copying `auth_asym_id` (author chain ID) into `label_asym_id` and `auth_seq_id` (author residue number) into `label_seq_id` while preserving all header metadata.

---

### 5. `pdb_module.py`, `pdb_clean.py` & `pdb_HETATM.py`
**Legacy Backward-Compatibility Wrappers**

- **`pdb_module.py`**: Wraps `EqualResidueAligner` for legacy scripts requiring `get_res_details`, `pdbTofasta`, or `mapping_residue`.
- **`pdb_clean.py`**: Wraps `StructureCleaner` to clean modified residues and resolve alternate location occupancies.
- **`pdb_HETATM.py`**: Wraps `StructureCleaner` to convert modified `HETATM` records to `ATOM` records.

---

- **`equal_residue_aligner.py --bound 1ASY_clean.pdb --bound-chains A --unbound 1EOV_clean.pdb --unbound-chains A -o equal_residue_dataset/`**


## Subdirectory: `utils/addtional/`

Contains standalone atomic manipulation tools derived from the standard `pdb-tools` suite:
- **`pdb_selaltloc.py`**: Selects a single alternate location for atoms with fractional occupancy.
- **`pdb_occ.py`**: Modifies the occupancy column of ATOM/HETATM records.
- **`pdb_selatom.py`**: Filters PDB files to retain only matching atom names.
