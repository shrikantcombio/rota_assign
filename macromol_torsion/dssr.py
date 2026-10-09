"""
DSSRParser Module
-----------------
Parser for x3dna-dssr (Dissecting the Spatial Structure of RNA / 3DNA).

Supports:
  1. Direct execution of x3dna-dssr binary on structure files (.pdb, .cif, .mmcif).
  2. Parsing pre-generated DSSR JSON files (.json).
  3. Parsing pre-generated DSSR text output files (.out, .txt).

Literature Reference:
  Lu et al. (2015) "DSSR: an integrated software tool for dissecting the spatial structure of RNA."
  Nucleic Acids Res. 43(21):e142 (doi:10.1093/nar/gkv716).
"""

import os
import re
import json
import shutil
import logging
import tempfile
import subprocess
from pathlib import Path
from functools import cached_property
import pandas as pd

try:
    from naxplorer.parsers.base import CSVParser
except ImportError:
    try:
        from .base import CSVParser
    except ImportError:
        try:
            from parsers.base import CSVParser
        except ImportError:
            class CSVParser:
                @staticmethod
                def extract_pdb_id(filename):
                    m = re.search(r'([0-9][A-Za-z0-9]{3})', Path(filename).stem.upper())
                    return m.group(1) if m else None

logger = logging.getLogger(__name__)


def parse_dssr_nt_id(token: str) -> dict:
    """
    Parses DSSR residue identifier tokens like:
      R.U601
      R.1MG637
      A.A1689^E
      B.5MC649
      A.ARG553

    Returns:
      {'chain': str, 'resname': str, 'resid': int, 'icode': str}
    """
    if not token or not isinstance(token, str):
        return {'chain': 'NA', 'resname': 'NA', 'resid': 0, 'icode': '.'}

    token = token.strip()
    pattern = r'^(?P<chain>[^.]+)\.(?P<resname>[A-Za-z0-9_]+?)(?P<resid>-?\d+)(?:\^(?P<icode>[A-Za-z]))?$'
    m = re.match(pattern, token)
    if m:
        return {
            'chain': m.group('chain'),
            'resname': m.group('resname'),
            'resid': int(m.group('resid')),
            'icode': m.group('icode') or '.'
        }

    return {'chain': 'NA', 'resname': token, 'resid': 0, 'icode': '.'}


def extract_sugar_pucker(summary: str) -> str:
    """Extracts sugar pucker conformer (e.g. ~C3'-endo, ~C2'-endo) from summary text."""
    if not summary:
        return 'NA'
    m = re.search(r'(~?C[1-5]\'-[a-z]+)', summary)
    if m:
        return m.group(1)
    return 'NA'


class DSSRParser:
    """
    Parser class for x3dna-dssr analysis.
    """

    DEFAULT_DSSR_PATHS = [
        Path.home() / "bin" / "x3dna-dssr",
        Path.home() / "x3DNA" / "x3dna-dssr",
        Path("/usr/local/bin/x3dna-dssr"),
        Path("/usr/bin/x3dna-dssr"),
    ]

    def __init__(self, input_file: str | Path, dssr_bin: str | Path = None, quiet: bool = True):
        self.input_file = Path(input_file)
        self.quiet = quiet
        self.pdb_id = CSVParser.extract_pdb_id(self.input_file.name) or "UNKNOWN"
        self.dssr_bin = self._resolve_dssr_binary(dssr_bin)
        self._raw_json = None

    def _resolve_dssr_binary(self, custom_bin: str | Path = None) -> str:
        """Resolves executable path for x3dna-dssr."""
        if custom_bin:
            cb = Path(custom_bin)
            if cb.exists() and os.access(cb, os.X_OK):
                return str(cb)

        for p in self.DEFAULT_DSSR_PATHS:
            if p.exists() and os.access(p, os.X_OK):
                return str(p)

        which_bin = shutil.which("x3dna-dssr")
        if which_bin:
            return which_bin

        return "x3dna-dssr"

    @cached_property
    def raw_json(self) -> dict:
        """Loads or executes DSSR to retrieve complete raw JSON dictionary."""
        if self._raw_json is not None:
            return self._raw_json

        if not self.input_file.exists():
            logger.warning(f"Input file not found: {self.input_file}")
            return {}

        suffix = self.input_file.suffix.lower()

        # Case 1: Input is already a JSON file
        if suffix == '.json':
            try:
                with open(self.input_file, 'r', encoding='utf-8') as f:
                    self._raw_json = json.load(f)
                return self._raw_json
            except Exception as e:
                logger.error(f"Error reading JSON file {self.input_file}: {e}")
                return {}

        # Case 2: Input is a structure file (.pdb, .cif, .mmcif)
        if suffix in ('.pdb', '.cif', '.mmcif'):
            return self._run_dssr_json(self.input_file)

        # Case 3: Input is a text output file (.out, .txt)
        if suffix in ('.out', '.txt'):
            struct_cand = [
                self.input_file.with_suffix('.cif'),
                self.input_file.with_suffix('.pdb'),
                self.input_file.with_suffix('.mmcif')
            ]
            for sc in struct_cand:
                if sc.exists():
                    return self._run_dssr_json(sc)

            return self._parse_text_out(self.input_file)

        return {}

    def _run_dssr_json(self, structure_file: Path) -> dict:
        """
        Executes x3dna-dssr binary with --json --more flags in an isolated temp directory
        to prevent un-prefixed dssr-* scratch file collisions across multiple PDBs.
        """
        stem = structure_file.stem.split('.')[0].split('_')[0].upper()
        dssr_dir = Path("/home/labuser/Projects/results/PRDBv35_clean/proRNA_data/x3dna-dssr")
        cache_file = dssr_dir / f"{stem}_dssr.json"
        if cache_file.exists() and cache_file.stat().st_size > 0:
            try:
                with open(cache_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    self._raw_json = data
                    return data
            except Exception:
                pass

        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = [self.dssr_bin, f"-i={structure_file.resolve()}", "--json", "--more"]
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=tmpdir)
                data = json.loads(res.stdout)
                self._raw_json = data
                if dssr_dir.exists():
                    try:
                        with open(cache_file, 'w', encoding='utf-8') as f_out:
                            json.dump(data, f_out)
                    except Exception:
                        pass
                return data
            except (subprocess.CalledProcessError, json.JSONDecodeError, FileNotFoundError) as e:
                logger.error(f"Failed to execute x3dna-dssr on {structure_file}: {e}")
                return {}

    def _parse_text_out(self, out_file: Path) -> dict:
        """Fallback parser for standard DSSR text output files (.out)."""
        pairs = []
        try:
            with open(out_file, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()

            in_pairs = False
            for line in lines:
                if 'List of' in line and 'base pair' in line:
                    in_pairs = True
                    continue
                if in_pairs:
                    if not line.strip() or line.startswith('*****'):
                        in_pairs = False
                        continue
                    parts = line.split()
                    if len(parts) >= 6 and parts[0].isdigit():
                        nt1 = parts[1]
                        nt2 = parts[2]
                        bp = parts[3]
                        name = parts[4]
                        lw = parts[6] if len(parts) >= 7 else 'cWW'
                        pairs.append({
                            'index': int(parts[0]),
                            'nt1': nt1,
                            'nt2': nt2,
                            'bp': bp,
                            'name': name,
                            'LW': lw
                        })
        except Exception as e:
            logger.error(f"Error parsing text file {out_file}: {e}")

        return {'pairs': pairs, 'num_pairs': len(pairs)}

    # ---------------------------------------------------------------------------
    # DataFrames & Extraction Methods
    # ---------------------------------------------------------------------------

    @cached_property
    def nts_df(self) -> pd.DataFrame:
        """DataFrame of nucleotides and torsion angles."""
        nts = self.raw_json.get('nts', [])
        if not nts:
            return pd.DataFrame()

        rows = []
        for nt in nts:
            nt_id = nt.get('nt_id', '')
            parsed_id = parse_dssr_nt_id(nt_id)
            summary_txt = nt.get('summary', '')

            row = {
                'pdb_id': self.pdb_id,
                'nt_id': nt_id,
                'chain': parsed_id['chain'],
                'resname': parsed_id['resname'],
                'resid': parsed_id['resid'],
                'icode': parsed_id['icode'],
                'dbn': nt.get('dbn', '.'),
                'summary': summary_txt,
                'alpha': nt.get('alpha'),
                'beta': nt.get('beta'),
                'gamma': nt.get('gamma'),
                'delta': nt.get('delta'),
                'epsilon': nt.get('epsilon'),
                'zeta': nt.get('zeta'),
                'epsilon_zeta': nt.get('epsilon_zeta'),
                'chi': nt.get('chi'),
                'glyco_bond': nt.get('glyco_bond', 'NA'),
                'bb_type': nt.get('bb_type', 'NA'),
                'sugar_pucker': nt.get('sugar_pucker') or extract_sugar_pucker(summary_txt)
            }
            rows.append(row)

        return pd.DataFrame(rows)

    @cached_property
    def pairs_df(self) -> pd.DataFrame:
        """DataFrame of base pairs and helical parameters."""
        pairs = self.raw_json.get('pairs', [])
        if not pairs:
            return pd.DataFrame()

        rows = []
        for bp in pairs:
            nt1 = bp.get('nt1', '')
            nt2 = bp.get('nt2', '')
            p1 = parse_dssr_nt_id(nt1)
            p2 = parse_dssr_nt_id(nt2)

            row = {
                'pdb_id': self.pdb_id,
                'bp_index': bp.get('index'),
                'nt1_id': nt1,
                'nt1_chain': p1['chain'],
                'nt1_resname': p1['resname'],
                'nt1_resid': p1['resid'],
                'nt1_icode': p1['icode'],
                'nt2_id': nt2,
                'nt2_chain': p2['chain'],
                'nt2_resname': p2['resname'],
                'nt2_resid': p2['resid'],
                'nt2_icode': p2['icode'],
                'bp': bp.get('bp', ''),
                'name': bp.get('name', ''),
                'saenger': bp.get('Saenger', 'NA'),
                'lw': bp.get('LW', 'NA'),
                'dssr': bp.get('DSSR', 'NA'),
                'shear': bp.get('shear'),
                'stretch': bp.get('stretch'),
                'stagger': bp.get('stagger'),
                'buckle': bp.get('buckle'),
                'propeller': bp.get('propeller'),
                'opening': bp.get('opening')
            }
            rows.append(row)

        return pd.DataFrame(rows)

    @cached_property
    def stems_df(self) -> pd.DataFrame:
        """DataFrame of RNA stems."""
        stems = self.raw_json.get('stems', [])
        if not stems:
            return pd.DataFrame()

        rows = []
        for s in stems:
            rows.append({
                'pdb_id': self.pdb_id,
                'stem_index': s.get('index'),
                'helix_index': s.get('helix_index'),
                'num_pairs': s.get('num_pairs'),
                'strand1': s.get('strand1', ''),
                'strand2': s.get('strand2', ''),
                'bp_type': s.get('bp_type', ''),
                'helix_form': s.get('helix_form', ''),
                'helical_rise': s.get('helical_rise'),
                'helical_radius': s.get('helical_radius')
            })
        return pd.DataFrame(rows)

    @cached_property
    def helices_df(self) -> pd.DataFrame:
        """DataFrame of RNA helices."""
        helices = self.raw_json.get('helices', [])
        if not helices:
            return pd.DataFrame()

        rows = []
        for h in helices:
            rows.append({
                'pdb_id': self.pdb_id,
                'helix_index': h.get('index'),
                'num_stems': h.get('num_stems'),
                'num_pairs': h.get('num_pairs'),
                'strand1': h.get('strand1', ''),
                'strand2': h.get('strand2', ''),
                'bp_type': h.get('bp_type', ''),
                'helix_form': h.get('helix_form', ''),
                'helical_rise': h.get('helical_rise'),
                'helical_radius': h.get('helical_radius')
            })
        return pd.DataFrame(rows)

    @cached_property
    def hairpins_df(self) -> pd.DataFrame:
        """DataFrame of hairpin loops."""
        hairpins = self.raw_json.get('hairpins', [])
        if not hairpins:
            return pd.DataFrame()

        rows = []
        for hp in hairpins:
            rows.append({
                'pdb_id': self.pdb_id,
                'index': hp.get('index'),
                'type': hp.get('type', 'hairpin loop'),
                'num_nts': hp.get('num_nts'),
                'nts_short': hp.get('nts_short', ''),
                'summary': hp.get('summary', '')
            })
        return pd.DataFrame(rows)

    @cached_property
    def junctions_df(self) -> pd.DataFrame:
        """DataFrame of multi-way junctions."""
        junctions = self.raw_json.get('junctions', [])
        if not junctions:
            return pd.DataFrame()

        rows = []
        for jc in junctions:
            rows.append({
                'pdb_id': self.pdb_id,
                'index': jc.get('index'),
                'type': jc.get('type', ''),
                'num_nts': jc.get('num_nts'),
                'num_stems': jc.get('num_stems'),
                'nts_short': jc.get('nts_short', '')
            })
        return pd.DataFrame(rows)

    @cached_property
    def hbonds_df(self) -> pd.DataFrame:
        """DataFrame of hydrogen bonds."""
        hbonds = self.raw_json.get('hbonds', [])
        if not hbonds:
            return pd.DataFrame()

        rows = []
        for hb in hbonds:
            atom1 = hb.get('atom1_id', '')
            atom2 = hb.get('atom2_id', '')
            rows.append({
                'pdb_id': self.pdb_id,
                'hbond_index': hb.get('index'),
                'atom1_id': atom1,
                'atom2_id': atom2,
                'distance': hb.get('distance'),
                'atom_pair': hb.get('atom_pair', ''),
                'residue_pair': hb.get('residue_pair', ''),
                'donAcc_type': hb.get('donAcc_type', '')
            })
        return pd.DataFrame(rows)

    @cached_property
    def multiplets_df(self) -> pd.DataFrame:
        """DataFrame of base multiplets (triplets/quadruplexes)."""
        multiplets = self.raw_json.get('multiplets', [])
        if not multiplets:
            return pd.DataFrame()

        rows = []
        for mp in multiplets:
            rows.append({
                'pdb_id': self.pdb_id,
                'index': mp.get('index'),
                'num_nts': mp.get('num_nts'),
                'nts_short': mp.get('nts_short', ''),
                'nts_long': mp.get('nts_long', '')
            })
        return pd.DataFrame(rows)

    @cached_property
    def dbn_string(self) -> str:
        """Returns secondary structure dot-bracket notation string."""
        dbn = self.raw_json.get('dbn', {})
        if isinstance(dbn, dict):
            if 'sstr' in dbn:
                return dbn['sstr']
            if 'all_chains' in dbn and isinstance(dbn['all_chains'], dict):
                return dbn['all_chains'].get('sstr', '')
            return dbn.get('sstruc', '')
        return str(dbn)

    # ---------------------------------------------------------------------------
    # Getter Aliases
    # ---------------------------------------------------------------------------

    def get_base_pairs(self) -> pd.DataFrame:
        return self.pairs_df

    def get_nucleotides(self) -> pd.DataFrame:
        return self.nts_df

    def get_stems(self) -> pd.DataFrame:
        return self.stems_df

    def get_helices(self) -> pd.DataFrame:
        return self.helices_df

    def get_hairpins(self) -> pd.DataFrame:
        return self.hairpins_df

    def get_junctions(self) -> pd.DataFrame:
        return self.junctions_df

    def get_hbonds(self) -> pd.DataFrame:
        return self.hbonds_df

    def get_multiplets(self) -> pd.DataFrame:
        return self.multiplets_df

    def get_secondary_structure(self) -> dict:
        """
        Extracts nucleotide-level secondary structure annotations:
        - Dot-bracket notation (dbn)
        - Secondary structure element type (ss_type): 'stem', 'hairpin', 'junction',
          'internal-loop', 'bulge', 'single-stranded', 'isolated-pair', 'pseudoknot', 'unpaired'
        - Single-letter SS code (ss_code): 'S', 'H', 'J', 'I', 'B', 'U', 'P', 'K', 'L'
        - Base pairing partner (paired_with) and base pair classification (bp_type)
        - DSSR summary descriptors

        Indexed by multiple convenient keys:
        - DSSR nt_id (e.g. 'R.U601')
        - Tuple (chain, str(resid)) (e.g. ('R', '601'))
        - Tuple (chain, int(resid)) (e.g. ('R', 601))
        - Residue label (e.g. 'U   R  601')
        """
        pairs = self.raw_json.get('pairs', [])
        pair_map = {}
        for bp in pairs:
            nt1 = bp.get('nt1', '')
            nt2 = bp.get('nt2', '')
            bp_name = bp.get('name', '.')
            bp_seq = bp.get('bp', '.')
            lw = bp.get('LW', '.')
            dssr_lw = bp.get('DSSR', '.')
            if nt1:
                pair_map[nt1] = {
                    'partner': nt2,
                    'bp': bp_seq,
                    'name': bp_name,
                    'lw': lw,
                    'dssr': dssr_lw
                }
            if nt2:
                pair_map[nt2] = {
                    'partner': nt1,
                    'bp': bp_seq,
                    'name': bp_name,
                    'lw': lw,
                    'dssr': dssr_lw
                }

        nts = self.raw_json.get('nts', [])
        ss_dict = {}
        for nt in nts:
            nt_id = nt.get('nt_id', '')
            parsed_id = parse_dssr_nt_id(nt_id)
            chain = parsed_id['chain']
            resname = parsed_id['resname']
            resid = parsed_id['resid']
            icode = parsed_id['icode']

            dbn = nt.get('dbn', '.')
            summary = nt.get('summary', '')

            # Classify secondary structure element motif:
            if dbn in ['(', ')', '[', ']', '{', '}', '<', '>']:
                if 'isolated-canonical' in summary or 'isolated' in summary:
                    ss_type = 'isolated-pair'
                    ss_code = 'P'
                elif 'pseudoknotted' in summary:
                    ss_type = 'pseudoknot'
                    ss_code = 'K'
                else:
                    ss_type = 'stem'
                    ss_code = 'S'
            else:
                if 'hairpin-loop' in summary:
                    ss_type = 'hairpin'
                    ss_code = 'H'
                elif 'junction-loop' in summary:
                    ss_type = 'junction'
                    ss_code = 'J'
                elif 'internal-loop' in summary or 'iloop' in summary:
                    ss_type = 'internal-loop'
                    ss_code = 'I'
                elif 'bulge' in summary:
                    ss_type = 'bulge'
                    ss_code = 'B'
                elif 'ss-non-loop' in summary or 'single-stranded' in summary:
                    ss_type = 'single-stranded'
                    ss_code = 'U'
                else:
                    ss_type = 'unpaired'
                    ss_code = 'L'

            p_info = pair_map.get(nt_id, {})
            paired_with = p_info.get('partner', '.')
            bp_type = p_info.get('name', '.')

            rec = {
                'nt_id': nt_id,
                'chain': chain,
                'resname': resname,
                'resid': resid,
                'icode': icode,
                'dbn': dbn,
                'ss_type': ss_type,
                'ss_code': ss_code,
                'paired_with': paired_with,
                'bp_type': bp_type,
                'summary': summary
            }

            lbl = f"{resname:<3s} {chain:1s} {str(resid):>4s}".strip()

            ss_dict[nt_id] = rec
            ss_dict[(chain, str(resid))] = rec
            ss_dict[(chain, int(resid))] = rec
            ss_dict[lbl] = rec

        return ss_dict

    def summary(self) -> dict:
        """Summary metrics dictionary."""
        return {
            'pdb_id': self.pdb_id,
            'num_nts': len(self.nts_df),
            'num_pairs': len(self.pairs_df),
            'num_stems': len(self.stems_df),
            'num_helices': len(self.helices_df),
            'num_hairpins': len(self.hairpins_df),
            'num_junctions': len(self.junctions_df),
            'num_hbonds': len(self.hbonds_df),
            'dbn': self.dbn_string
        }


class DSSRAnalyzer:
    """
    Interface to DSSR (x3dna-dssr) for RNA secondary structure assignment.
    Follows the DsspAnalyzer design pattern in macromol_torsion.
    """

    def __init__(self, dssr_binary=None):
        self.dssr_bin = dssr_binary

    def analyze(self, structure_file: str | Path) -> DSSRParser:
        return DSSRParser(structure_file, dssr_bin=self.dssr_bin)

    def get_secondary_structure(self, structure_file: str | Path) -> dict:
        parser = self.analyze(structure_file)
        return parser.get_secondary_structure()


def parse_dssr(input_file: str | Path, dssr_bin: str | Path = None) -> DSSRParser:
    """Convenience function to parse structure or DSSR output with DSSRParser."""
    parser = DSSRParser(input_file, dssr_bin=dssr_bin)
    return parser
