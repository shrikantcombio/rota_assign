"""
merge_pipeline.py
-----------------
Master end-to-end processing pipeline module.
1. Calculates protein & RNA dihedrals, DSSP, and NACCESS SASA -> dataset/results/1ASY_protein_torsions.csv
2. Assigns protein rotamer states using BBD rotamer library (library/RBPs_BBD_rotamer_lib) -> dataset/results/1ASY_protein_rotastates.csv
3. Assigns RNA dinucleotide conformers (NtC, CANA, confalA) using StructAnnotate -> dataset/results/1ASY_rna_conformers.csv
4. Uses interaction details (including NAxplorer dataset/1ASY_interactions_NAxplorer.csv) -> dataset/results/1ASY_interactions.csv
5. Exports final merged dataset into dataset/results/1ASY_final_output.csv matching schema:
   PDB_ID,AA,AA_atom,AA_chain,AA_pos,Nuc,Nuc_atom,Nuc_chain,Nuc_pos,hb_dis,NtC,CANA,confalA,b_phi,b_psi,b_ss,b_rotastate
6. Runs Arginine forks & rotamer transition preference analysis.
"""

from pathlib import Path
import sys
import pandas as pd
import numpy as np

# Path setup
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
RESULTS_DIR = ROOT_DIR / "dataset" / "results"
TEST_DIR = ROOT_DIR / "test"
DATASET_DIR = ROOT_DIR / "dataset"

if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from macromol_torsion import ComplexTorsionCalculator
from script.rotaassign import ProteinRotamerAssigner
from script.rna_conformer_assigner import RNAConformerAssigner
from script.analyze_rotastates import run_full_rotastate_analysis


def map_complex_interactions_and_conformers(df_interactions, df_protein_assigned, df_rna_assigned, target_pdb="1ASY"):
    """
    Maps interactions, energy, bound/unbound rotamers, SASA, and RNA conformers into the exact 29-column schema:
    PDB_ID, AA, AA_atom, AA_chain, AA_pos, Nuc, Nuc_atom, Nuc_chain, Nuc_pos, b_phi, b_psi, b_ss, b_rotastate, u_phi, u_psi, u_ss, u_rotastate, SASA, NtC, CANA, confalA, dist, interaction_type, edge, score, mediator_info, E_int_kcalmol, status, source
    """
    prot_map = {}
    if df_protein_assigned is not None and not df_protein_assigned.empty:
        for idx, row in df_protein_assigned.iterrows():
            if 'LABEL' in row and not pd.isna(row['LABEL']):
                parts = str(row['LABEL']).split()
                chain = parts[1] if len(parts) > 1 else 'A'
                pos = parts[2] if len(parts) > 2 else '0'
            else:
                chain = str(row.get('AA_chain', 'A'))
                pos = str(row.get('AA_pos', '0'))

            u_rota_raw = row.get('ub_rotastate', row.get('u_rotastate', 'NA'))
            u_rota_val = 'NA' if (pd.isna(u_rota_raw) or str(u_rota_raw).strip() in ['', 'nan', 'NaN', 'None']) else str(u_rota_raw).strip()
            
            b_rota_raw = row.get('b_rotastate', row.get('B_ROTASTATE', 'NA'))
            b_rota_val = 'NA' if (pd.isna(b_rota_raw) or str(b_rota_raw).strip() in ['', 'nan', 'NaN', 'None']) else str(b_rota_raw).strip()

            prot_map[(chain, str(pos))] = {
                'b_phi': row.get('B_PHI', row.get('b_phi', 0.0)),
                'b_psi': row.get('B_PSI', row.get('b_psi', 0.0)),
                'b_ss': row.get('B_SS', row.get('b_ss', 'L')),
                'b_rotastate': b_rota_val,
                'u_phi': row.get('U_PHI', row.get('u_phi', 0.0)),
                'u_psi': row.get('U_PSI', row.get('u_psi', 0.0)),
                'u_ss': row.get('U_SS', row.get('u_ss', 'L')),
                'u_rotastate': u_rota_val,
                'SASA': row.get('SASA', row.get('sasa', 'N'))
            }

    rna_map = {}
    if df_rna_assigned is not None and not df_rna_assigned.empty:
        for idx, row in df_rna_assigned.iterrows():
            if 'LABEL' in row and not pd.isna(row['LABEL']):
                parts = str(row['LABEL']).split()
                chain = parts[1] if len(parts) > 1 else 'R'
                pos = parts[2] if len(parts) > 2 else '0'
            else:
                chain = str(row.get('Nuc_chain', 'R'))
                pos = str(row.get('Nuc_pos', '0'))

            rna_map[(chain, str(pos))] = {
                'NtC': row.get('NtC', 'NANT'),
                'CANA': row.get('CANA', 'NAN'),
                'confalA': row.get('confalA', 0.0)
            }

    rows = []
    for idx, r in df_interactions.iterrows():
        pdb_id = str(r.get('pdb4', r.get('PDB_ID', r.get('pdb_id', target_pdb)))).strip().upper()

        aa = str(r.get('prot_resn', r.get('AA', 'GLY'))).strip()
        aa_chain = str(r.get('prot_chain', r.get('AA_chain', 'A'))).strip()
        aa_pos = str(r.get('prot_resid', r.get('AA_pos', '0'))).strip()

        nuc = str(r.get('rna_base', r.get('Nuc', 'U'))).strip()
        nuc_chain = str(r.get('rna_chain', r.get('Nuc_chain', 'R'))).strip()
        nuc_pos = str(r.get('rna_resid', r.get('Nuc_pos', '0'))).strip()

        if 'atoms' in r and not pd.isna(r['atoms']):
            atom_str = str(r['atoms'])
            atom_parts = atom_str.split(':') if ':' in atom_str else [atom_str, '']
            nuc_atom = str(r.get('Nuc_atom', atom_parts[0] if len(atom_parts) > 0 else ''))
            aa_atom = str(r.get('AA_atom', atom_parts[1] if len(atom_parts) > 1 else ''))
        else:
            nuc_atom = str(r.get('Nuc_atom', ''))
            aa_atom = str(r.get('AA_atom', ''))

        dist_val = r.get('dist', r.get('hb_dis', '-'))
        int_type = r.get('interaction', r.get('interaction_type', 'H-BOND'))
        edge_val = r.get('edge', '-')
        score_val = r.get('score', '-')
        mediator_val = r.get('mediator_info', '-')
        e_int_val = r.get('E_int_kcalmol', '-')
        status_val = r.get('status', 'ok')
        source_val = r.get('source', 'NAxplorer_energy')

        p_data = prot_map.get((aa_chain, aa_pos), {
            'b_phi': r.get('b_phi', 0.0),
            'b_psi': r.get('b_psi', 0.0),
            'b_ss': r.get('b_ss', 'L'),
            'b_rotastate': r.get('b_rotastate', 'NA'),
            'u_phi': r.get('u_phi', 0.0),
            'u_psi': r.get('u_psi', 0.0),
            'u_ss': r.get('u_ss', 'L'),
            'u_rotastate': r.get('u_rotastate', 'NA'),
            'SASA': r.get('SASA', 'N')
        })

        r_data = rna_map.get((nuc_chain, nuc_pos), {
            'NtC': r.get('NtC', 'NANT'),
            'CANA': r.get('CANA', 'NAN'),
            'confalA': r.get('confalA', 0.0)
        })

        b_rot_val = p_data['b_rotastate']
        if pd.isna(b_rot_val) or str(b_rot_val).strip() in ['', 'nan', 'NaN', 'None']:
            b_rot_val = 'NA'

        u_rot_val = p_data['u_rotastate']
        if pd.isna(u_rot_val) or str(u_rot_val).strip() in ['', 'nan', 'NaN', 'None']:
            u_rot_val = 'NA'

        rows.append({
            'PDB_ID': pdb_id,
            'AA': aa,
            'AA_atom': aa_atom,
            'AA_chain': aa_chain,
            'AA_pos': aa_pos,
            'Nuc': nuc,
            'Nuc_atom': nuc_atom,
            'Nuc_chain': nuc_chain,
            'Nuc_pos': nuc_pos,
            'b_phi': p_data['b_phi'],
            'b_psi': p_data['b_psi'],
            'b_ss': p_data['b_ss'],
            'b_rotastate': str(b_rot_val).strip(),
            'u_phi': p_data['u_phi'],
            'u_psi': p_data['u_psi'],
            'u_ss': p_data['u_ss'],
            'u_rotastate': str(u_rot_val).strip(),
            'SASA': p_data['SASA'],
            'NtC': r_data['NtC'],
            'CANA': r_data['CANA'],
            'confalA': r_data['confalA'],
            'dist': dist_val,
            'interaction_type': int_type,
            'edge': edge_val,
            'score': score_val,
            'mediator_info': mediator_val,
            'E_int_kcalmol': e_int_val,
            'status': status_val,
            'source': source_val
        })

    df_out = pd.DataFrame(rows)
    has_unbound_protein = (df_protein_assigned is not None and not df_protein_assigned.empty and ('ub_rotastate' in df_protein_assigned.columns or 'U_PHI' in df_protein_assigned.columns))
    
    if has_unbound_protein:
        target_cols = [
            'PDB_ID', 'AA', 'AA_atom', 'AA_chain', 'AA_pos',
            'Nuc', 'Nuc_atom', 'Nuc_chain', 'Nuc_pos',
            'b_phi', 'b_psi', 'b_ss', 'b_rotastate',
            'u_phi', 'u_psi', 'u_ss', 'u_rotastate', 'SASA',
            'NtC', 'CANA', 'confalA',
            'dist', 'interaction_type', 'edge', 'score',
            'mediator_info', 'E_int_kcalmol', 'status', 'source'
        ]
    else:
        target_cols = [
            'PDB_ID', 'AA', 'AA_atom', 'AA_chain', 'AA_pos',
            'Nuc', 'Nuc_atom', 'Nuc_chain', 'Nuc_pos',
            'b_phi', 'b_psi', 'b_ss', 'b_rotastate', 'SASA',
            'NtC', 'CANA', 'confalA',
            'dist', 'interaction_type', 'edge', 'score',
            'mediator_info', 'E_int_kcalmol', 'status', 'source'
        ]
    return df_out[target_cols]


class ComplexMergePipeline:
    """
    Master pipeline for running torsion calculations, rotamer assignment, RNA conformer annotation,
    and merging interaction data into final output datasets under dataset/results/.
    """

    def __init__(self, output_dir=None, lib_type='BBD'):
        self.output_dir = Path(output_dir) if output_dir else RESULTS_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.lib_type = lib_type
        self.rotamer_assigner = ProteinRotamerAssigner(lib_type=lib_type)
        self.rna_assigner = RNAConformerAssigner()
        self.calculator = ComplexTorsionCalculator()

    def process_interactions_file(self, interactions_source, target_pdb="1ASY"):
        """
        Loads interaction details from NAxplorer dataset/1ASY_interactions_NAxplorer.csv,
        1ASY_rnp_interactions_model1_energy.csv, or standard CSV files and writes intermediate dataset/results/1ASY_interactions.csv.
        Populates score from 1ASY_interactions_NAxplorer.csv if missing.
        """
        src_p = Path(interactions_source)
        if not src_p.exists():
            raise FileNotFoundError(f"Interactions source file not found: {src_p}")

        clean_pdb = target_pdb.split('.')[0].split('_')[0].upper()

        try:
            df_raw = pd.read_csv(src_p, sep=',')
            if len(df_raw.columns) == 1 and ('\t' in str(df_raw.columns[0]) or ' ' in str(df_raw.columns[0])):
                df_raw = pd.read_csv(src_p, sep=r'\s+', engine='python')
        except Exception:
            df_raw = pd.read_csv(src_p, sep=r'\s+', engine='python')

        # Load score lookup map from NAxplorer file if available
        nax_map = {}
        nax_file = DATASET_DIR / f"{clean_pdb}_interactions_NAxplorer.csv"
        if not nax_file.exists():
            nax_file = DATASET_DIR / f"{target_pdb}_interactions_NAxplorer.csv"

        if nax_file.exists():
            try:
                try:
                    df_nax = pd.read_csv(nax_file, sep=r'\s+', engine='python')
                except Exception:
                    df_nax = pd.read_csv(nax_file)
                for _, nr in df_nax.iterrows():
                    p_parts = str(nr.get('prot_aa', '')).split(':')
                    r_parts = str(nr.get('rna_nt', '')).split(':')
                    atoms_str = str(nr.get('atoms', ''))
                    if len(p_parts) >= 3 and len(r_parts) >= 3:
                        key = (p_parts[0], str(p_parts[2]), r_parts[0], str(r_parts[2]), atoms_str)
                        nax_map[key] = nr.get('score', '-')
            except Exception:
                pass

        rows = []
        for idx, r in df_raw.iterrows():
            pdb_id = str(r.get('pdb4', r.get('pdb_id', r.get('PDB_ID', clean_pdb)))).split('.')[0].strip().upper()

            # Parse Protein Residue Info
            if 'prot_resn' in r and not pd.isna(r['prot_resn']):
                aa = str(r['prot_resn']).strip()
                aa_chain = str(r.get('prot_chain', 'A')).strip()
                aa_pos = str(r.get('prot_resid', '0')).strip()
            elif 'prot_aa' in r and not pd.isna(r['prot_aa']):
                prot_parts = str(r['prot_aa']).split(':')
                aa_chain = prot_parts[0] if len(prot_parts) > 0 else 'A'
                aa = prot_parts[1] if len(prot_parts) > 1 else 'GLY'
                aa_pos = prot_parts[2] if len(prot_parts) > 2 else '0'
            else:
                aa = str(r.get('AA', 'GLY')).strip()
                aa_chain = str(r.get('AA_chain', 'A')).strip()
                aa_pos = str(r.get('AA_pos', '0')).strip()

            # Parse RNA Residue Info
            if 'rna_base' in r and not pd.isna(r['rna_base']):
                nuc = str(r['rna_base']).strip()
                nuc_chain = str(r.get('rna_chain', 'R')).strip()
                nuc_pos = str(r.get('rna_resid', '0')).strip()
            elif 'rna_nt' in r and not pd.isna(r['rna_nt']):
                rna_parts = str(r['rna_nt']).split(':')
                nuc_chain = rna_parts[0] if len(rna_parts) > 0 else 'R'
                nuc = rna_parts[1] if len(rna_parts) > 1 else 'U'
                nuc_pos = rna_parts[2] if len(rna_parts) > 2 else '0'
            else:
                nuc = str(r.get('Nuc', 'U')).strip()
                nuc_chain = str(r.get('Nuc_chain', 'R')).strip()
                nuc_pos = str(r.get('Nuc_pos', '0')).strip()

            # Parse Atoms
            if 'atoms' in r and not pd.isna(r['atoms']):
                atom_str = str(r['atoms'])
                atom_parts = atom_str.split(':') if ':' in atom_str else [atom_str, '']
                nuc_atom = str(r.get('Nuc_atom', atom_parts[0] if len(atom_parts) > 0 else ''))
                aa_atom = str(r.get('AA_atom', atom_parts[1] if len(atom_parts) > 1 else ''))
            else:
                nuc_atom = str(r.get('Nuc_atom', ''))
                aa_atom = str(r.get('AA_atom', ''))

            # Parse Interaction, Energy, Status details
            dist_val = r.get('dist', r.get('hb_dis', '-'))
            int_type = r.get('interaction', r.get('interaction_type', 'H-BOND'))
            edge_val = r.get('edge', '-')
            score_val = r.get('score', r.get('overlap', '-'))

            if (pd.isna(score_val) or str(score_val).strip() in ['-', '']) and nax_map:
                key = (aa_chain, str(aa_pos), nuc_chain, str(nuc_pos), f"{nuc_atom}:{aa_atom}")
                if key in nax_map:
                    score_val = nax_map[key]

            mediator_val = r.get('mediator_info', '-')
            e_int_val = r.get('E_int_kcalmol', '-')
            status_val = r.get('status', 'ok')
            source_val = r.get('source', 'NAxplorer_energy' if 'E_int_kcalmol' in r else 'NAxplorer')

            rows.append({
                'PDB_ID': pdb_id,
                'AA': aa,
                'AA_atom': aa_atom,
                'AA_chain': aa_chain,
                'AA_pos': aa_pos,
                'Nuc': nuc,
                'Nuc_atom': nuc_atom,
                'Nuc_chain': nuc_chain,
                'Nuc_pos': nuc_pos,
                'dist': dist_val,
                'interaction_type': int_type,
                'edge': edge_val,
                'score': score_val,
                'mediator_info': mediator_val,
                'E_int_kcalmol': e_int_val,
                'status': status_val,
                'source': source_val
            })

        df_int = pd.DataFrame(rows)
        out_int_file = self.output_dir / f"{clean_pdb}_interactions.csv"
        df_int.to_csv(out_int_file, index=False)
        return out_int_file, df_int

    def run_pipeline(self, bound_file=None, unbound_file=None, unbound_rna_file=None, protein_chains="A,B", rna_chains="R,S", interactions_file=None, pdb_id=None, lib_type=None, mode=None, ub_protein=None, ub_RNA=None, ub_rna=None):
        """
        Executes pipeline tailored to user output schema specifications:
        1. Unbound-Unbound Case (bound complex + unbound protein + unbound RNA):
           Output:
           - 1ASY_1EOV_protein_rotastates.csv (complex protein and unbound protein rotastates)
           - 1ASY_2TRA_RNA_conformers.csv (complex RNA and unbound RNA conformers, with PDB_ID column)
           - 1ASY_rotastate_RNA_conf.csv (merged complex interactions, rotamers, and conformers)
        2. Complex alone (e.g. 1ASY):
           Output: 1ASY_rotastate_RNA_conf.csv
        3. Free protein alone (e.g. 1EOV):
           Calculates dihedrals, DSSP, rotastates (NO NACCESS).
           Output: 1EOV_rotastate.csv
        4. Free RNA alone (e.g. 2TRA):
           Calculates RNA dihedrals, NtC conformers.
           Output: 2TRA_RNA_confs.csv
        """
        unbound_file = ub_protein or unbound_file
        unbound_rna_file = ub_RNA or ub_rna or unbound_rna_file
        l_type = lib_type or self.lib_type

        # Determine mode if not explicitly provided
        if mode is None:
            if bound_file and unbound_file and unbound_rna_file:
                mode = 'unbound-unbound'
            elif bound_file:
                mode = 'complex'
            elif unbound_file and not bound_file and not unbound_rna_file:
                mode = 'free-protein'
            elif unbound_rna_file and not bound_file and not unbound_file:
                mode = 'free-rna'
            else:
                mode = 'complex'

        raw_pdb = pdb_id or (Path(bound_file).stem.split('_')[0].split('-')[0] if bound_file else (Path(unbound_file).stem.split('_')[0].split('-')[0] if unbound_file else Path(unbound_rna_file).stem.split('_')[0].split('-')[0]))
        target_pdb = raw_pdb.split('.')[0].split('_')[0].upper()

        # --- MODE 3: FREE PROTEIN ALONE ---
        if mode == 'free-protein':
            unbound_p = Path(unbound_file)
            df_protein = self.calculator.analyze_protein_bound_unbound(
                bound_file=unbound_p,
                unbound_file=unbound_p,
                chains=protein_chains,
                orig_bound_file=unbound_p,
                orig_unbound_file=unbound_p,
                run_dssp=True,
                run_naccess=False
            )
            df_protein_assigned = self.rotamer_assigner.assign_df(df_protein, mode='bound', lib_type=l_type)
            df_protein_assigned['PDB_ID'] = target_pdb
            
            free_prot_file = self.output_dir / f"{target_pdb}_rotastate.csv"
            df_protein_assigned.to_csv(free_prot_file, index=False)
            return {
                'FREE_PROTEIN_ROTASTAT_FILE': free_prot_file,
                'PROTEIN_ROTA_FILE': free_prot_file,
                'FINAL_OUTPUT_FILE': free_prot_file,
                'DF': df_protein_assigned
            }

        # --- MODE 4: FREE RNA ALONE ---
        if mode == 'free-rna':
            unbound_rna_p = Path(unbound_rna_file)
            df_rna = self.calculator.analyze_rna_7d_suite(
                bound_file=unbound_rna_p,
                unbound_file=unbound_rna_p,
                rna_chains=rna_chains,
                complex_file=unbound_rna_p,
                run_naccess=False
            )
            df_rna_assigned = self.rna_assigner.assign_df(df_rna)
            df_rna_assigned['PDB_ID'] = target_pdb

            free_rna_file = self.output_dir / f"{target_pdb}_RNA_confs.csv"
            df_rna_assigned.to_csv(free_rna_file, index=False)
            return {
                'FREE_RNA_CONF_FILE': free_rna_file,
                'RNA_CONF_FILE': free_rna_file,
                'FINAL_OUTPUT_FILE': free_rna_file,
                'DF': df_rna_assigned
            }

        # --- MODES 1 & 2: UNBOUND-UNBOUND & COMPLEX ---
        bound_p = Path(bound_file)
        unbound_p = Path(unbound_file) if unbound_file else None
        unbound_rna_p = Path(unbound_rna_file) if unbound_rna_file else None

        # 1. Analyze protein dihedrals, DSSP & SASA
        df_protein = self.calculator.analyze_protein_bound_unbound(
            bound_file=bound_p,
            unbound_file=unbound_p,
            chains=protein_chains,
            orig_bound_file=bound_p,
            orig_unbound_file=unbound_p,
            run_dssp=True,
            run_naccess=True
        )
        df_protein_assigned = self.rotamer_assigner.assign_df(df_protein, mode='both' if unbound_file else 'bound', lib_type=l_type)

        # 2. Analyze RNA dihedrals & SASA
        df_rna = self.calculator.analyze_rna_7d_suite(
            bound_file=bound_p,
            unbound_file=unbound_rna_p,
            rna_chains=rna_chains,
            complex_file=bound_p,
            run_naccess=True
        )
        df_rna_assigned = self.rna_assigner.assign_df(df_rna)

        # Save intermediate files
        prot_torsion_file = self.output_dir / f"{target_pdb}_protein_torsions.csv"
        df_protein.to_csv(prot_torsion_file, index=False)

        prot_rota_file = self.output_dir / f"{target_pdb}_protein_rotastates.csv"
        df_protein_assigned.to_csv(prot_rota_file, index=False)

        rna_conf_file = self.output_dir / f"{target_pdb}_rna_conformers.csv"
        df_rna_assigned.to_csv(rna_conf_file, index=False)

        # Process Interactions
        if interactions_file:
            int_file_path, df_interactions = self.process_interactions_file(interactions_file, target_pdb=target_pdb)
        else:
            default_energy = DATASET_DIR / f"{target_pdb}_rnp_interactions_model1_energy.csv"
            default_naxplorer = DATASET_DIR / f"{target_pdb}_interactions_NAxplorer.csv"
            if default_energy.exists():
                int_file_path, df_interactions = self.process_interactions_file(default_energy, target_pdb=target_pdb)
            elif default_naxplorer.exists():
                int_file_path, df_interactions = self.process_interactions_file(default_naxplorer, target_pdb=target_pdb)
            else:
                df_interactions = None
                int_file_path = None

        # Build 1ASY_rotastate_RNA_conf.csv (Complex merged output)
        if df_interactions is not None:
            df_merged = map_complex_interactions_and_conformers(
                df_interactions=df_interactions,
                df_protein_assigned=df_protein_assigned,
                df_rna_assigned=df_rna_assigned,
                target_pdb=target_pdb
            )
            complex_merged_file = self.output_dir / f"{target_pdb}_rotastate_RNA_conf.csv"
            df_merged.to_csv(complex_merged_file, index=False)
            
            # Legacy alias 1ASY_final_output.csv
            df_merged.to_csv(self.output_dir / f"{target_pdb}_final_output.csv", index=False)
            analysis_results = run_full_rotastate_analysis(complex_merged_file, output_dir=self.output_dir)
        else:
            complex_merged_file = None
            df_merged = None
            analysis_results = None

        # If Unbound Case: Generate pair files
        prot_pair_file = None
        rna_pair_file = None
        if unbound_file:
            prot_stem = Path(unbound_file).stem.split('.')[0].split('_')[0].split('-')[0]
            prot_pair_file = self.output_dir / f"{target_pdb}_{prot_stem}_protein_rotastates.csv"
            df_protein_assigned.to_csv(prot_pair_file, index=False)

        if unbound_rna_file:
            ub_rna_p = Path(unbound_rna_file)
            rna_stem = ub_rna_p.stem.split('.')[0].split('_')[0].split('-')[0].upper()

            df_rna_bound_copy = df_rna_assigned.copy()
            if 'PDB_ID' not in df_rna_bound_copy.columns:
                df_rna_bound_copy.insert(0, 'PDB_ID', target_pdb)

            # Analyze and assign unbound RNA structure dihedrals & conformers
            try:
                df_rna_ub = self.calculator.analyze_rna_7d_suite(
                    bound_file=ub_rna_p,
                    unbound_file=ub_rna_p,
                    rna_chains=rna_chains,
                    complex_file=ub_rna_p,
                    run_naccess=False
                )
                df_rna_ub_assigned = self.rna_assigner.assign_df(df_rna_ub)
                if 'PDB_ID' not in df_rna_ub_assigned.columns:
                    df_rna_ub_assigned.insert(0, 'PDB_ID', rna_stem)
                df_rna_pair = pd.concat([df_rna_bound_copy, df_rna_ub_assigned], ignore_index=True)
            except Exception:
                df_rna_pair = df_rna_bound_copy

            rna_pair_file = self.output_dir / f"{target_pdb}_{rna_stem}_RNA_conformers.csv"
            df_rna_pair.to_csv(rna_pair_file, index=False)

        if unbound_file or unbound_rna_file:
            prot_stem = Path(unbound_file).stem.split('.')[0].split('_')[0].split('-')[0] if unbound_file else "bound"
            rna_stem = Path(unbound_rna_file).stem.split('.')[0].split('_')[0].split('-')[0] if unbound_rna_file else "bound"
            pair_tag = f"{target_pdb}_{prot_stem}_{rna_stem}"
            pair_dir = self.output_dir / "pairs" / pair_tag
            pair_dir.mkdir(parents=True, exist_ok=True)
            df_protein_assigned.to_csv(pair_dir / f"{pair_tag}_protein_rotastates.csv", index=False)
            if unbound_rna_file and 'df_rna_pair' in locals():
                df_rna_pair.to_csv(pair_dir / f"{pair_tag}_RNA_conformers.csv", index=False)
            else:
                df_rna_bound_copy = df_rna_assigned.copy()
                if 'PDB_ID' not in df_rna_bound_copy.columns:
                    df_rna_bound_copy.insert(0, 'PDB_ID', target_pdb)
                df_rna_bound_copy.to_csv(pair_dir / f"{pair_tag}_RNA_conformers.csv", index=False)

        return {
            'PROTEIN_TORSIONS_FILE': prot_torsion_file,
            'PROTEIN_ROTA_FILE': prot_rota_file,
            'RNA_CONF_FILE': rna_conf_file,
            'INTERACTIONS_FILE': int_file_path,
            'COMPLEX_MERGED_FILE': complex_merged_file,
            'UNBOUND_PROTEIN_ROTA_FILE': prot_pair_file,
            'UNBOUND_RNA_CONF_FILE': rna_pair_file,
            'FINAL_OUTPUT_FILE': complex_merged_file,
            'ANALYSIS': analysis_results
        }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="End-to-End Complex Merge Pipeline")
    parser.add_argument("--pdb", default="1ASY", help="PDB ID to process (default: 1ASY)")
    parser.add_argument("--bound", default=None, help="Bound PDB file path")
    parser.add_argument("--ub_protein", "--unbound", default=None, help="Unbound protein PDB/CIF file path")
    parser.add_argument("--ub_RNA", "--ub_rna", "--unbound-rna", default=None, help="Unbound RNA PDB/CIF file path")
    parser.add_argument("--interactions", default=None, help="Interactions CSV file path")
    parser.add_argument("--mode", default=None, help="Execution mode (unbound-unbound, complex, free-protein, free-rna)")
    parser.add_argument("-o", "--output-dir", default=None, help="Output directory (default: dataset/results)")

    args = parser.parse_args()

    out_dir = Path(args.output_dir) if args.output_dir else RESULTS_DIR
    clean_pdb = args.pdb.split('.')[0].split('_')[0].upper()

    bound_file = Path(args.bound) if args.bound else (DATASET_DIR / f"{clean_pdb}.cif" if (DATASET_DIR / f"{clean_pdb}.cif").exists() else DATASET_DIR / f"{clean_pdb}.pdb")
    ub_protein_file = Path(args.ub_protein) if args.ub_protein else None
    ub_rna_file = Path(args.ub_RNA) if args.ub_RNA else None
    int_file = Path(args.interactions) if args.interactions else DATASET_DIR / f"{clean_pdb}_rnp_interactions_model1_energy.csv"

    pipeline = ComplexMergePipeline(output_dir=out_dir)
    res = pipeline.run_pipeline(
        bound_file=bound_file if bound_file and bound_file.exists() else None,
        ub_protein=ub_protein_file if ub_protein_file and ub_protein_file.exists() else None,
        ub_RNA=ub_rna_file if ub_rna_file and ub_rna_file.exists() else None,
        interactions_file=int_file if int_file and int_file.exists() else None,
        pdb_id=clean_pdb,
        mode=args.mode
    )

    print("Pipeline Execution Completed Successfully!")
    if res.get('PROTEIN_TORSIONS_FILE'):
        print("Intermediate Protein Torsions File:", res['PROTEIN_TORSIONS_FILE'])
    if res.get('PROTEIN_ROTA_FILE'):
        print("Intermediate Protein Rotastates File:", res['PROTEIN_ROTA_FILE'])
    if res.get('RNA_CONF_FILE'):
        print("Intermediate RNA Conformers File:  ", res['RNA_CONF_FILE'])
    if res.get('INTERACTIONS_FILE'):
        print("Intermediate Interactions File:    ", res['INTERACTIONS_FILE'])
    if res.get('COMPLEX_MERGED_FILE'):
        print("Complex Merged File:               ", res['COMPLEX_MERGED_FILE'])
    if res.get('UNBOUND_PROTEIN_ROTA_FILE'):
        print("Unbound Protein Rotastates File:   ", res['UNBOUND_PROTEIN_ROTA_FILE'])
    if res.get('UNBOUND_RNA_CONF_FILE'):
        print("Unbound RNA Conformers File:       ", res['UNBOUND_RNA_CONF_FILE'])
    if res.get('FINAL_OUTPUT_FILE'):
        print("Final Output File:                 ", res['FINAL_OUTPUT_FILE'])
    if res.get('ANALYSIS'):
        print("\nAnalysis Reports Output:")
        print("  Arginine Forks File:             ", res['ANALYSIS']['FORK_FILE'])
        print("  Rotamer Transitions File:        ", res['ANALYSIS']['TRANSITION_FILE'])
        print("  Rotamer-Conformer Cross File:    ", res['ANALYSIS']['PREFERENCE_FILE'])


if __name__ == "__main__":
    main()
