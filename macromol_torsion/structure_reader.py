"""
StructureReader Module
-----------------------
Unified parser for PDB and mmCIF structure files using mmcif_clean_reader.py.
"""

from pathlib import Path
import sys
import os
import tempfile

ROTAMER_LIB_SCRIPT = Path('/home/labuser/Projects/PhD_projects/RBPs_rotamer_lib/script')
if ROTAMER_LIB_SCRIPT.exists() and str(ROTAMER_LIB_SCRIPT) not in sys.path:
    sys.path.insert(0, str(ROTAMER_LIB_SCRIPT))

try:
    from utils.mmcif_clean_reader import parse_file as parse_cif_or_pdb
except ImportError:
    parse_cif_or_pdb = None


class StructureReader:
    """
    Parser for PDB and mmCIF files. Standardizes atom records into dictionaries.
    """

    @classmethod
    def parse_structure(cls, file_path):
        """
        Parses a PDB or CIF file and returns a list of standardized atom dictionaries.
        """
        file_p = Path(file_path)

        if parse_cif_or_pdb is not None:
            cols, rows = parse_cif_or_pdb(str(file_p), all_models=False)
            atoms = []
            for r in rows:
                rec = r.get('_atom_site.group_PDB', 'ATOM')
                aname = r.get('_atom_site.auth_atom_id', r.get('_atom_site.label_atom_id', '')).strip()
                resname = r.get('_atom_site.auth_comp_id', r.get('_atom_site.label_comp_id', '')).strip()
                chain = r.get('_atom_site.auth_asym_id', r.get('_atom_site.label_asym_id', 'A')).strip()
                res_seq = r.get('_atom_site.auth_seq_id', r.get('_atom_site.label_seq_id', '0')).strip()
                icode = r.get('_atom_site.pdbx_PDB_ins_code', '.').strip()
                if icode in ['.', '?', '']:
                    icode = ' '

                try:
                    x = float(r.get('_atom_site.Cartn_x', 0.0))
                    y = float(r.get('_atom_site.Cartn_y', 0.0))
                    z = float(r.get('_atom_site.Cartn_z', 0.0))
                except (ValueError, TypeError):
                    continue

                try:
                    occ = float(r.get('_atom_site.occupancy', 1.0))
                except (ValueError, TypeError):
                    occ = 1.0

                elem = r.get('_atom_site.type_symbol', '').strip()

                # Format standard residue label matching legacy code (e.g. 'ARG A  18' or '  A B  12')
                res_label = f"{resname:<3s} {chain:1s} {res_seq:>4s}".strip()

                atoms.append({
                    'record': rec,
                    'atom_name': aname,
                    'resname': resname,
                    'chain': chain,
                    'resSeq': res_seq,
                    'icode': icode,
                    'res_label': res_label,
                    'xyz': [x, y, z],
                    'elem': elem,
                    'occupancy': occ,
                    'raw_line': r
                })
            return atoms

        else:
            return cls._parse_pdb_lines_fallback(file_p)

    @classmethod
    def _parse_pdb_lines_fallback(cls, file_path):
        atoms = []
        with open(file_path, 'r') as f:
            for line in f:
                rec = line[0:6].strip()
                if rec in ['ATOM', 'HETATM']:
                    aname = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21:22].strip() or 'A'
                    res_seq = line[22:26].strip()
                    icode = line[26:27].strip() or ' '
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        continue

                    res_label = f"{resname:<3s} {chain:1s} {res_seq:>4s}".strip()
                    atoms.append({
                        'record': rec,
                        'atom_name': aname,
                        'resname': resname,
                        'chain': chain,
                        'resSeq': res_seq,
                        'icode': icode,
                        'res_label': res_label,
                        'xyz': [x, y, z],
                        'elem': line[76:78].strip() or (aname[0] if aname else 'C'),
                        'occupancy': 1.0
                    })
        return atoms

    @classmethod
    def read_pdb_lines(cls, file_path):
        file_p = Path(file_path)
        if file_p.suffix.lower() == '.pdb':
            with open(file_p, 'r') as f:
                return f.readlines()

        atoms = cls.parse_structure(file_p)
        pdb_lines = []
        for i, a in enumerate(atoms, start=1):
            line = (
                f"{a['record']:<6s}{i:>5d} {a['atom_name']:<4s} {a['resname']:<3s} "
                f"{a['chain']:1s}{int(a['resSeq']):>4d}{a['icode']:1s}   "
                f"{a['xyz'][0]:8.3f}{a['xyz'][1]:8.3f}{a['xyz'][2]:8.3f}"
                f"{a['occupancy']:6.2f} 0.00          {a['elem']:>2s}\n"
            )
            pdb_lines.append(line)
        return pdb_lines

    @classmethod
    def write_temp_pdb(cls, file_path, output_dir=None):
        file_p = Path(file_path)
        if file_p.suffix.lower() == '.pdb':
            return file_p

        lines = cls.read_pdb_lines(file_p)
        out_dir = Path(output_dir) if output_dir else Path(tempfile.gettempdir())
        temp_pdb = out_dir / f"{file_p.stem}_temp.pdb"
        with open(temp_pdb, 'w') as f:
            f.writelines(lines)
        return temp_pdb
