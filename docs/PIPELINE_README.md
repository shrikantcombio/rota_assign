# Comprehensive Documentation: `rota_assign` End-to-End Pipeline

## 1. Executive Summary & Capabilities Overview

The **`rota_assign`** framework is an integrated computational suite designed for high-throughput structural analysis, side-chain rotamer assignment, RNA dinucleotide conformer annotation, non-covalent interaction classification, and binding energy integration in protein-RNA complexes and free molecular structures.

### Execution Modes & Dynamic Column Schema

`rota_assign` dynamically scopes column outputs in `1ASY_rotastate_RNA_conf.csv` based on whether an unbound protein partner is passed:

1. **Unbound RNA Partner Alone (`--ub_RNA 2TRA.cif`)**:
   - Command: `python script/merge_pipeline.py --pdb 1ASY.cif --ub_RNA 2TRA.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/`
   - **No unbound protein is loaded**. Unbound protein columns (`u_phi`, `u_psi`, `u_ss`, `u_rotastate`) are **omitted** from `1ASY_rotastate_RNA_conf.csv` (25 columns total).
   - Output: `1ASY_2TRA_RNA_conformers.csv` (contains both `1ASY` complex RNA and `2TRA` unbound RNA dinucleotide conformers, with explicit `PDB_ID` column) and `1ASY_rotastate_RNA_conf.csv` (25 columns).
2. **Unbound Protein Partner Alone (`--ub_protein 1EOV.cif`)**:
   - Command: `python script/merge_pipeline.py --pdb 1ASY.cif --ub_protein 1EOV.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/`
   - Unbound protein features are computed. Unbound protein columns (`u_phi`, `u_psi`, `u_ss`, `u_rotastate`) **are included** in `1ASY_rotastate_RNA_conf.csv` (29 columns total).
   - Output: `1ASY_1EOV_protein_rotastates.csv` and `1ASY_rotastate_RNA_conf.csv` (29 columns).
3. **Dual Unbound Case (`--ub_protein 1EOV.cif --ub_RNA 2TRA.cif`)**:
   - Command: `python script/merge_pipeline.py --pdb 1ASY.cif --ub_protein 1EOV.cif --ub_RNA 2TRA.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/`
   - Outputs:
     - `1ASY_1EOV_protein_rotastates.csv`: Complex protein and unbound protein details including backbone dihedrals ($\phi, \psi$), DSSP 3-state, SASA, and rotastates (`b_rotastate`, `ub_rotastate`).
     - `1ASY_2TRA_RNA_conformers.csv`: Complex RNA (`1ASY`, 147 suites) and unbound RNA (`2TRA`, 70 suites) dinucleotide conformers (`NtC`, `CANA`, `confalA`), 7D suite dihedrals, SASA, with explicit `PDB_ID` column.
     - `1ASY_rotastate_RNA_conf.csv`: Unified 29-column dataset combining bound/unbound rotamer details, SASA, RNA conformers, non-covalent interactions, and binding energy.
4. **Complex Structure Alone (`--pdb 1ASY.cif`)**:
   - Command: `python script/merge_pipeline.py --pdb 1ASY.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/`
   - Processes **ONLY** the complex protein and complex RNA chains present in the input complex file.
   - Output: `1ASY_rotastate_RNA_conf.csv` (25 columns, omitting unbound protein columns) and clean 12-column `1ASY_protein_torsions.csv`.

---

## 2. Structural & Physicochemical Features

1. **Full-Length Macromolecular Dihedral & DSSP Calculation**: Computes backbone dihedrals ($\phi, \psi$), side-chain chi angles ($\chi_1..\chi_5$), and DSSP secondary structure codes using full-length structure files (*Kabsch & Sander 1983*).
2. **Rost & Sander (1993) Secondary Structure Corrections**:
   - Categorized into Helices (`H`), Strands (`S`), and Loops (`L`).
   - Sequence rules applied: `B_` $\rightarrow$ `SS` and `B_B` $\rightarrow$ `LLL`.
3. **Physicochemical Interface & Surface Classification**:
   - **Interface Residues (`I`)**: Classified using NACCESS (*Hubbard & Thornton 1993; Lee & Richards 1971*) if residue SASA decreases by at least **$1.0\ \text{\AA}^2$** ($\Delta\text{SASA} \ge 1.0\ \text{\AA}^2$) upon complex formation (*Barik et al., 2012 PRince; Belapure et al., 2023 IRAA*).
   - **Surface Residues (`S`)**: Defined following *Guharoy et al. (2010)* and *Miller et al. (1987)* if relative SASA in NACCESS is **$> 8.0\%$** from extended Ala-X-Ala tripeptide accessibility.
   - **Non-Interface Surface Residues (`N`)**: Surface residues (relative SASA $> 8.0\%$) excluding interface residues ($\Delta\text{SASA} \ge 1.0\ \text{\AA}^2$).
4. **Protein Rotamer State Assignment**: Matches side-chains against BBD, SSD, or BBI libraries using log-likelihood weighted angular deviation scoring. Assigns ON (`p`, `t`, `m`) and OFF (`P`, `T`, `M`) states, `no_chi` for Ala/Gly, and `NA` for missing states.
5. **RNA Dinucleotide Conformer Assignment**: Annotates RNA 7D suites into 28 NtC classes (`AA00`, `AA07`, `AA10`, `OP15`, etc.), CANA strings, and CONFAL quality scores via StructAnnotate (*Kant et al., 2026, doi:10.64898/2026.05.14.725147*).
6. **Dynamic Schema Format**: Output contains 25 columns when no unbound protein is specified, and 29 columns when an unbound protein structure is provided.

---

## 3. Usage Instructions

### Run Unbound RNA Partner (--ub_RNA 2TRA.cif)
```bash
python script/merge_pipeline.py --pdb 1ASY.cif --ub_RNA 2TRA.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/
```

### Run Unbound Protein Partner (--ub_protein 1EOV.cif)
```bash
python script/merge_pipeline.py --pdb 1ASY.cif --ub_protein 1EOV.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/
```

### Run Dual Unbound Case (1ASY bound + 1EOV free protein + 2TRA free RNA)
```bash
python script/merge_pipeline.py --pdb 1ASY.cif --ub_protein 1EOV.cif --ub_RNA 2TRA.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/
```

### Run Complex Alone Case (1ASY)
```bash
python script/merge_pipeline.py --pdb 1ASY.cif --interactions 1ASY_rnp_interactions_model1_energy.csv -o results/
```

---

## 4. Scientific References

1. **Kabsch, W. and Sander, C. (1983)** Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. *Biopolymers*, 22, 2577–2637.
2. **Rost, B. and Sander, C. (1993)** Prediction of Protein Secondary Structure at Better than 70% Accuracy. *Journal of Molecular Biology*, 232, 584–599.
3. **Hubbard, S.J. and Thornton, J.M. (1993)** NACCESS Computer Program, Department of Biochemistry and Molecular Biology, University College London.
4. **Lee, B. and Richards, F.M. (1971)** The interpretation of protein structures: Estimation of static accessibility. *Journal of Molecular Biology*, 55, 379-IN4.
5. **Barik, A., Mishra, A. and Bahadur, R.P. (2012)** PRince: a web server for structural and physicochemical analysis of protein-RNA interface. *Nucleic Acids Research*, 40, W408–W412.
6. **Belapure, J., Sorokina, M. and Kastritis, P.L. (2023)** IRAA: A statistical tool for investigating a protein–protein interaction interface from multiple structures. *Protein Science*, 32, e4523.
7. **Guharoy, M., Janin, J. and Robert, C.H. (2010)** Side-chain rotamer transitions at protein-protein interfaces. *Proteins: Structure, Function, and Bioinformatics*, 78, 3219–3225.
8. **Miller, S., Janin, J., Lesk, A.M. and Chothia, C. (1987)** Interior and surface of monomeric proteins. *Journal of Molecular Biology*, 196, 641–656.
9. **Kant et al. (2026)** RNA Dinucleotide Conformers and NtC Classification Library. *bioRxiv*, doi:10.64898/2026.05.14.725147.
10. **Kant, S., Nithin, C., Mukherjee, S., Maity, A., & Bahadur, R. P. (2025)** Protein-RNA Docking Benchmark v3.0 Integrated With Binding Affinity. *Proteins: Structure, Function, and Bioinformatics*, doi:10.1002/prot.26825.
