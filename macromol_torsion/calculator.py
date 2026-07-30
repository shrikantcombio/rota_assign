"""
ComplexTorsionCalculator Module
-------------------------------
Master calculator class for analyzing Protein-RNA complex structures and free structures.
Integrates Protein and RNA backbone dihedrals, 7D dinucleotide suites, RNA pseudo-torsions,
sugar pucker pseudorotation, glycosidic Chi angles, DSSP, NACCESS SASA interface/surface filtering,
and supports PDB and mmCIF structure formats.
"""

from pathlib import Path
import tempfile
import pandas as pd
import numpy as np

from macromol_torsion.structure_reader import StructureReader
from macromol_torsion.protein_backbone import ProteinBackboneCalculator
from macromol_torsion.protein_sidechain import ProteinSidechainCalculator
from macromol_torsion.dssp import DsspAnalyzer
from macromol_torsion.sasa import SasaAnalyzer

from macromol_torsion.rna_backbone import RNABackboneCalculator
from macromol_torsion.rna_pseudo import RNAPseudoTorsionCalculator
from macromol_torsion.rna_sugar import RNASugarPuckerCalculator
from macromol_torsion.rna_sidechain import RNASidechainCalculator


class ComplexTorsionCalculator:
    """
    Unified master calculator for Protein-RNA complex structures.
    """

    def __init__(self, dssp_binary=None, naccess_binary=None):
        self.protein_bb_calc = ProteinBackboneCalculator()
        self.protein_sc_calc = ProteinSidechainCalculator()
        self.dssp_analyzer = DsspAnalyzer(dssp_binary=dssp_binary)
        self.sasa_analyzer = SasaAnalyzer(naccess_binary=naccess_binary)

        self.rna_bb_calc = RNABackboneCalculator()
        self.rna_pseudo_calc = RNAPseudoTorsionCalculator()
        self.rna_sugar_calc = RNASugarPuckerCalculator()
        self.rna_sc_calc = RNASidechainCalculator()

    def analyze_protein_bound_unbound(self, bound_file, unbound_file, chains=None, orig_bound_file=None, orig_unbound_file=None, run_dssp=True, run_naccess=True, output_file=None):
        """
        Calculates protein backbone and sidechain dihedrals for bound and unbound structure pairs.
        Uses original cleaned structure files (orig_bound_file / orig_unbound_file) for DSSP and NACCESS SASA.
        """
        bound_p = Path(bound_file)
        unbound_p = Path(unbound_file)

        target_bound_dssp_naccess = Path(orig_bound_file) if orig_bound_file else bound_p
        target_unbound_dssp_naccess = Path(orig_unbound_file) if orig_unbound_file else unbound_p

        pdb_tag = bound_p.stem.split('_')[0]

        # 1. Bound Protein Torsions
        b_phi, b_psi, b_omega, b_tau, b_xi, b_labels = self.protein_bb_calc.calculate_dihedrals(bound_p)
        b_chis = self.protein_sc_calc.calculate_all_chis(bound_p)

        # 2. Unbound Protein Torsions
        u_phi, u_psi, u_omega, u_tau, u_xi, u_labels = self.protein_bb_calc.calculate_dihedrals(unbound_p)
        u_chis = self.protein_sc_calc.calculate_all_chis(unbound_p)

        # Map unbound dihedrals by label & (chain, resSeq)
        u_map = {}
        for i, lbl in enumerate(u_labels):
            parts = lbl.split()
            chain = parts[1] if len(parts) > 1 else 'A'
            seq = parts[2] if len(parts) > 2 else '0'
            res_key = (chain, seq)

            data_item = {
                'PHI': u_phi[i], 'PSI': u_psi[i], 'OMEGA': u_omega[i],
                'TAU': u_tau[i], 'XI': u_xi[i] if i < len(u_xi) else 0.0,
                'CHI1': u_chis['CHI1'].get(lbl, 0.0),
                'CHI2': u_chis['CHI2'].get(lbl, 0.0),
                'CHI3': u_chis['CHI3'].get(lbl, 0.0),
                'CHI4': u_chis['CHI4'].get(lbl, 0.0),
                'CHI5': u_chis['CHI5'].get(lbl, 0.0),
            }
            u_map[lbl] = data_item
            u_map[res_key] = data_item

        # 3. DSSP Secondary Structure on original cleaned structure files
        b_ss_map = {}
        u_ss_map = {}
        if run_dssp:
            try:
                b_ss_map = self.dssp_analyzer.get_secondary_structure(target_bound_dssp_naccess)
                u_ss_map = self.dssp_analyzer.get_secondary_structure(target_unbound_dssp_naccess)
            except Exception:
                pass

        # 4. NACCESS SASA Interface Classification on Protein-RNA Complex vs Protein Alone
        interface_set = set()
        if run_naccess and chains:
            try:
                b_atoms = StructureReader.parse_structure(target_bound_dssp_naccess)
                chains_set = set(chains) if isinstance(chains, (list, tuple, set)) else set(str(chains))
                prot_atoms = [a for a in b_atoms if a.get('chain') in chains_set and a.get('record') in ['ATOM', 'HETATM']]

                from clean_structure import StructureCleaner
                sc = StructureCleaner()
                temp_prot_file = Path(tempfile.gettempdir()) / f"{target_bound_dssp_naccess.stem}_prot_alone.pdb"
                sc.write_pdb(prot_atoms, temp_prot_file)

                interface, non_int, surf, core = self.sasa_analyzer.classify_residues(target_bound_dssp_naccess, temp_prot_file, chains)
                interface_set = set(interface)
            except Exception:
                pass

        records = []
        for i, lbl in enumerate(b_labels):
            parts = lbl.split()
            res_name = parts[0]
            chain = parts[1] if len(parts) > 1 else 'A'
            seq = parts[2] if len(parts) > 2 else '0'
            res_key = (chain, seq)

            b_c1 = b_chis['CHI1'].get(lbl, 0.0)
            b_c2 = b_chis['CHI2'].get(lbl, 0.0)
            b_c3 = b_chis['CHI3'].get(lbl, 0.0)
            b_c4 = b_chis['CHI4'].get(lbl, 0.0)
            b_c5 = b_chis['CHI5'].get(lbl, 0.0)

            u_item = u_map.get(lbl) or u_map.get(res_key, {
                'PHI': 0.0, 'PSI': 0.0, 'OMEGA': 0.0, 'TAU': 0.0, 'XI': 0.0,
                'CHI1': 0.0, 'CHI2': 0.0, 'CHI3': 0.0, 'CHI4': 0.0, 'CHI5': 0.0
            })

            b_ss = b_ss_map.get(lbl) or b_ss_map.get(res_key) or 'L'
            u_ss = u_ss_map.get(lbl) or u_ss_map.get(res_key) or 'L'

            sasa_type = 'I' if (lbl in interface_set or res_key in interface_set) else 'N'

            records.append({
                'PDB': pdb_tag,
                'LABEL': lbl,
                'B_PHI': b_phi[i],
                'B_PSI': b_psi[i],
                'B_OMEGA': b_omega[i],
                'B_CHI1': b_c1,
                'B_CHI2': b_c2,
                'B_CHI3': b_c3,
                'B_CHI4': b_c4,
                'B_CHI5': b_c5,
                'B_SS': b_ss,
                'U_PHI': u_item['PHI'],
                'U_PSI': u_item['PSI'],
                'U_OMEGA': u_item['OMEGA'],
                'U_CHI1': u_item['CHI1'],
                'U_CHI2': u_item['CHI2'],
                'U_CHI3': u_item['CHI3'],
                'U_CHI4': u_item['CHI4'],
                'U_CHI5': u_item['CHI5'],
                'U_SS': u_ss,
                'SASA': sasa_type
            })

        df = pd.DataFrame(records)

        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            sep = '\t' if out_p.suffix.lower() in ['.tsv', '.txt'] else ','
            df.to_csv(out_p, sep=sep, index=False)

        return df

    def analyze_rna_7d_suite(self, bound_file, unbound_file=None, rna_chains=None, complex_file=None, run_naccess=True, output_file=None):
        """
        Calculates RNA 7D dinucleotide suite dihedrals (delm, epsm, zetm, alpha, beta, gamma, delta),
        glycosidic Chi angles, inter-nucleotide distances, and SASA interface classification.
        Supports complex_file (e.g. 1ASY.pdb) for full protein-RNA SASA and RNA calculations.
        Automatically evaluates RNA SASA interface against RNA-alone structure if unbound_file lacks RNA.
        """
        target_bound = complex_file if complex_file else bound_file

        labels, alpha, beta, gamma, delta, epsilon, zeta, delm, epsm, zetm = self.rna_bb_calc.calculate_7d_suite(bound_file)

        if len(labels) == 0 and complex_file:
            target_bound = complex_file
            labels, alpha, beta, gamma, delta, epsilon, zeta, delm, epsm, zetm = self.rna_bb_calc.calculate_7d_suite(target_bound)

        lab_chi, chi1, chi2, nn_dist, cc_dist, mu = self.rna_sc_calc.calculate_sidechain_and_distances(target_bound)

        sc_map = {}
        for i, lbl in enumerate(lab_chi):
            sc_map[lbl] = {
                'CHI1': chi1[i] if i < len(chi1) else 0.0,
                'CHI2': chi2[i] if i < len(chi2) else 0.0,
                'NN': nn_dist[i] if i < len(nn_dist) else 0.0,
                'CC': cc_dist[i] if i < len(cc_dist) else 0.0,
                'MU': mu[i] if i < len(mu) else 0.0
            }

        interface_set = set()
        if run_naccess and rna_chains and isinstance(target_bound, (str, Path)):
            try:
                bound_p = Path(target_bound)
                rna_set = set(rna_chains) if isinstance(rna_chains, (list, tuple, set)) else set(str(rna_chains))

                has_unbound_rna = False
                if unbound_file and Path(unbound_file).exists():
                    u_atoms = StructureReader.parse_structure(unbound_file)
                    has_unbound_rna = any(a.get('chain') in rna_set for a in u_atoms)

                if has_unbound_rna:
                    unbound_p = Path(unbound_file)
                else:
                    b_atoms = StructureReader.parse_structure(bound_p)
                    rna_atoms = [a for a in b_atoms if a.get('chain') in rna_set and a.get('record') in ['ATOM', 'HETATM']]

                    from clean_structure import StructureCleaner
                    sc = StructureCleaner()
                    unbound_p = Path(tempfile.gettempdir()) / f"{bound_p.stem}_rna_alone.pdb"
                    sc.write_pdb(rna_atoms, unbound_p)

                interface, non_int, surf, core = self.sasa_analyzer.classify_residues(bound_p, unbound_p, rna_chains)
                interface_set = set(interface)
            except Exception:
                pass

        records = []
        for i, lbl in enumerate(labels):
            sc = sc_map.get(lbl, {'CHI1': 0.0, 'CHI2': 0.0, 'NN': 0.0, 'CC': 0.0, 'MU': 0.0})
            parts = lbl.split()
            chain = parts[1] if len(parts) > 1 else 'R'
            seq = parts[2] if len(parts) > 2 else '0'
            res_key = (chain, seq)

            sasa_type = 'I' if (lbl in interface_set or res_key in interface_set) else 'N'

            records.append({
                'LABEL': lbl,
                'DELM': delm[i],
                'EPSM': epsm[i],
                'ZETM': zetm[i],
                'ALPHA': alpha[i],
                'BETA': beta[i],
                'GAMMA': gamma[i],
                'DELTA': delta[i],
                'CHI1': sc['CHI1'],
                'CHI2': sc['CHI2'],
                'NN_DIST': sc['NN'],
                'CC_DIST': sc['CC'],
                'MU': sc['MU'],
                'SASA': sasa_type
            })

        df = pd.DataFrame(records)

        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            sep = '\t' if out_p.suffix.lower() in ['.tsv', '.txt'] else ','
            df.to_csv(out_p, sep=sep, index=False)

        return df

    def analyze_rna_full(self, structure_input, run_naccess=False, output_file=None):
        """
        Calculates all RNA parameters for a single structure or atom list:
        6-bond backbone, 7D suite, pseudo-torsions (eta, theta, eta', theta'),
        sugar pucker pseudorotation (phase, amplitude, pucker state), and glycosidic Chi angles.
        """
        labels, alpha, beta, gamma, delta, epsilon, zeta, delm, epsm, zetm = self.rna_bb_calc.calculate_7d_suite(structure_input)

        p_labels, eta, theta, etad, thetad = self.rna_pseudo_calc.calculate_pseudo_torsions(structure_input)
        pseudo_map = {}
        for i, lbl in enumerate(p_labels):
            pseudo_map[lbl] = {
                'ETA': eta[i], 'THETA': theta[i],
                'ETA_PRIME': etad[i], 'THETA_PRIME': thetad[i]
            }

        sugar_results = self.rna_sugar_calc.calculate_pseudorotation(structure_input)
        sugar_map = {r['label']: r for r in sugar_results}

        lab_chi, chi1, chi2, nn_dist, cc_dist, mu = self.rna_sc_calc.calculate_sidechain_and_distances(structure_input)
        sc_map = {}
        for i, lbl in enumerate(lab_chi):
            sc_map[lbl] = {
                'CHI1': chi1[i], 'CHI2': chi2[i],
                'NN': nn_dist[i], 'CC': cc_dist[i], 'MU': mu[i]
            }

        records = []
        for i, lbl in enumerate(labels):
            ps = pseudo_map.get(lbl, {'ETA': 0.0, 'THETA': 0.0, 'ETA_PRIME': 0.0, 'THETA_PRIME': 0.0})
            sg = sugar_map.get(lbl, {'phase': 0.0, 'amplitude': 0.0, 'pucker': 'N/A', 'v0': 0.0, 'v1': 0.0, 'v2': 0.0, 'v3': 0.0, 'v4': 0.0})
            sc = sc_map.get(lbl, {'CHI1': 0.0, 'CHI2': 0.0, 'NN': 0.0, 'CC': 0.0, 'MU': 0.0})

            records.append({
                'LABEL': lbl,
                'ALPHA': alpha[i], 'BETA': beta[i], 'GAMMA': gamma[i], 'DELTA': delta[i], 'EPSILON': epsilon[i], 'ZETA': zeta[i],
                'DELM': delm[i], 'EPSM': epsm[i], 'ZETM': zetm[i],
                'ETA': ps['ETA'], 'THETA': ps['THETA'], 'ETA_PRIME': ps['ETA_PRIME'], 'THETA_PRIME': ps['THETA_PRIME'],
                'SUGAR_PHASE': sg['phase'], 'SUGAR_AMP': sg['amplitude'], 'SUGAR_PUCKER': sg['pucker'],
                'GLYCOSIDIC_CHI': sc['CHI2'], 'MU': sc['MU'], 'CC_DIST': sc['CC'], 'NN_DIST': sc['NN']
            })

        df = pd.DataFrame(records)

        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            sep = '\t' if out_p.suffix.lower() in ['.tsv', '.txt'] else ','
            df.to_csv(out_p, sep=sep, index=False)

        return df

    def analyze_complex(self, bound_file, unbound_file, protein_chains, rna_chains, orig_bound_file=None, orig_unbound_file=None, run_dssp=True, run_naccess=True, output_file=None):
        """
        Calculates both Protein and RNA torsion angles for a bound (protein-RNA complex) and unbound structure pair.
        Uses orig_bound_file (e.g. 1ASY_clean.pdb) and orig_unbound_file (e.g. 1EOV_clean.pdb) for DSSP and NACCESS SASA.
        Saves output as CSV or TSV based on file extension.
        :returns: Dictionary with 'PROTEIN' DataFrame and 'RNA' DataFrame.
        """
        bound_p = Path(bound_file)
        unbound_p = Path(unbound_file)

        orig_b = orig_bound_file if orig_bound_file else bound_p
        orig_u = orig_unbound_file if orig_unbound_file else unbound_p

        # 1. Protein Analysis
        df_protein = self.analyze_protein_bound_unbound(
            bound_file=bound_p,
            unbound_file=unbound_p,
            chains=protein_chains,
            orig_bound_file=orig_b,
            orig_unbound_file=orig_u,
            run_dssp=run_dssp,
            run_naccess=run_naccess
        )

        # 2. RNA Analysis
        df_rna = self.analyze_rna_7d_suite(
            bound_file=bound_p,
            unbound_file=unbound_p,
            rna_chains=rna_chains,
            complex_file=orig_b,
            run_naccess=run_naccess
        )

        if output_file:
            out_p = Path(output_file)
            out_p.parent.mkdir(parents=True, exist_ok=True)
            ext = out_p.suffix.lower()
            sep = '\t' if ext in ['.tsv', '.txt'] else ','

            if ext == '.csv' or '_protein' not in out_p.stem:
                prot_out = out_p.parent / out_p.name
                rna_out = out_p.parent / f"{out_p.stem}_rna{ext}"
            else:
                prot_out = out_p
                rna_out = out_p.parent / f"{out_p.stem}_rna{ext}"

            df_protein.to_csv(prot_out, sep=sep, index=False)
            df_rna.to_csv(rna_out, sep=sep, index=False)

        return {
            'PROTEIN': df_protein,
            'RNA': df_rna
        }
