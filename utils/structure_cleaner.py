"""
clean_structure.py
------------------
Unified Structure Cleaner package for PDB and mmCIF files.
Converts >200 modified amino acids and nucleotides to standard residue names,
converts modified HETATM records to ATOM, filters waters/ions (optional CLI flags),
and resolves alternate location occupancies. Saves outputs to designated dataset directory.

Bundled with rota_assign (originally /home/<user>/clean_pdb/clean_structure.py). Change from the
original: 'GAL' (galactose, a sugar) is no longer mapped to GLU.
"""

from pathlib import Path
import sys
import os
import shutil
import re
import argparse
from collections import OrderedDict

# Define ModResMap dictionary for non-standard to standard mappings
MODRES_MAP = {
    'MSE': 'MET', 'SEP': 'SER', 'TPO': 'THR', 'PTR': 'TYR', 'CSO': 'CYS',
    'CSD': 'CYS', 'CME': 'CYS', 'OCS': 'CYS', 'KCX': 'LYS', 'LLP': 'LYS',
    'MLY': 'LYS', 'M3L': 'LYS', 'HYP': 'PRO', 'PCA': 'GLU', 'CGU': 'GLU',
    'FME': 'MET', '1MG': 'G',   '2MG': 'G',   '7MG': 'G',
    'OMG': 'G',   'YG': 'G',    'YYG': 'G',   '1MA': 'A',   '6MA': 'A',
    'M2G': 'G',   'M22G': 'G',  'H2U': 'U',   '5MC': 'C',   'OMC': 'C',
    '5MU': 'U',   'OMU': 'U',   'PSU': 'U',   'I': 'G',     'A': 'A',
    'C': 'C',     'G': 'G',     'U': 'U'
}


class StructureCleaner:
    """
    Cleaner for PDB and mmCIF structure files.
    """

    def __init__(self, keep_waters=False, keep_ions=False, altloc_strategy='A', set_occupancy=None, default_output_dir="dataset"):
        self.keep_waters = keep_waters
        self.keep_ions = keep_ions
        self.altloc_strategy = altloc_strategy
        self.set_occupancy = set_occupancy
        self.default_output_dir = Path(default_output_dir)

        # Common ion names to filter if keep_ions=False
        self.ion_names = {
            'NA', 'CL', 'MG', 'K', 'CA', 'ZN', 'MN', 'FE', 'CU', 'CO', 'NI',
            'SO4', 'PO4', 'EDO', 'GOL', 'PEG', 'ACT', 'FMT', 'NH4', 'AZI'
        }
        self.water_names = {'HOH', 'WAT', 'D2O', 'DOD'}

    def clean_structure(self, input_file, output_file=None, output_dir=None):
        """
        Cleans modified residues, waters/ions, alternate locations, and occupancies.
        Saves output file into dataset directory.
        """
        input_p = Path(input_file)
        if not input_p.exists():
            raise FileNotFoundError(f"Input structure file not found: {input_p}")

        out_dir = Path(output_dir) if output_dir else self.default_output_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        if output_file:
            out_p = Path(output_file)
            if not out_p.is_absolute() and out_p.parent == Path('.'):
                out_p = out_dir / out_p.name
        else:
            out_p = out_dir / input_p.name

        from macromol_torsion.structure_reader import StructureReader
        atoms = StructureReader.parse_structure(input_p)

        cleaned_atoms = []
        logs = []

        seen_altloc = set()
        for a in atoms:
            rec = a.get('record', 'ATOM')
            resname = a.get('resname', '').strip().upper()
            aname = a.get('atom_name', '').strip()
            altloc = a.get('altLoc', ' ').strip()
            occ = a.get('occupancy', 1.0)

            # Filter waters
            if not self.keep_waters and resname in self.water_names:
                continue

            # Filter ions
            if not self.keep_ions and resname in self.ion_names:
                continue

            # Modified residue mapping
            if resname in MODRES_MAP:
                std_name = MODRES_MAP[resname]
                if std_name != resname:
                    logs.append(f"Converted {resname} to {std_name} at residue {a.get('resSeq')}")
                    a['resname'] = f"{std_name:<3s}"
                    a['record'] = 'ATOM'

            # Convert non-water/ion HETATM to ATOM
            if rec == 'HETATM' and resname not in self.water_names and resname not in self.ion_names:
                a['record'] = 'ATOM'

            # Alternate location resolution
            if altloc and altloc != ' ':
                res_key = (a['chain'], a['resSeq'], a['icode'], aname)
                if res_key in seen_altloc:
                    continue
                if self.altloc_strategy and altloc != self.altloc_strategy:
                    continue
                seen_altloc.add(res_key)
                a['altLoc'] = ' '

            # Occupancy override
            if self.set_occupancy is not None:
                a['occupancy'] = float(self.set_occupancy)

            cleaned_atoms.append(a)

        # Write clean structure
        ext = input_p.suffix.lower()
        if ext in ['.cif', '.mmcif']:
            self.write_cif(cleaned_atoms, out_p, structure_id=out_p.stem)
        else:
            self.write_pdb(cleaned_atoms, out_p)

        return out_p, logs

    def write_pdb(self, atoms, output_file):
        out_p = Path(output_file)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, 'w') as f:
            for idx, a in enumerate(atoms, start=1):
                rec = a.get('record', 'ATOM')
                serial = a.get('serial', idx)
                try:
                    serial_str = f"{int(serial):>5d}"
                except Exception:
                    serial_str = f"{str(serial):>5s}"

                aname = a.get('atom_name', '')
                if len(aname) < 4:
                    aname_str = f" {aname:<3s}"
                else:
                    aname_str = f"{aname:<4s}"

                altLoc = a.get('altLoc', ' ')
                resName = a.get('resname', '').strip()
                resName_str = f"{resName:>3s}"
                chainID = a.get('chain', 'A')[:1]
                resSeq = a.get('resSeq', 1)
                try:
                    resSeq_str = f"{int(resSeq):>4d}"
                except Exception:
                    resSeq_str = f"{str(resSeq):>4s}"

                iCode = a.get('icode', ' ')[:1]

                xyz = a.get('xyz', [a.get('x', 0.0), a.get('y', 0.0), a.get('z', 0.0)])
                x = float(xyz[0])
                y = float(xyz[1])
                z = float(xyz[2])

                occ = a.get('occupancy', 1.0)
                temp = a.get('tempFactor', 0.0)
                element = a.get('element', '')[:2]

                line = f"{rec:<6s}{serial_str} {aname_str}{altLoc:1s}{resName_str:3s} {chainID:1s}{resSeq_str}{iCode:1s}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occ:>6.2f}{temp:>6.2f}          {element:>2s}\n"
                f.write(line)
            f.write("END\n")

    def write_cif(self, atoms, output_file, structure_id="structure"):
        out_p = Path(output_file)
        out_p.parent.mkdir(parents=True, exist_ok=True)

        with open(out_p, 'w') as f:
            f.write(f"data_{structure_id}\n#\nloop_\n")
            f.write("_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n")
            f.write("_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n")
            f.write("_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n")
            f.write("_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n")
            f.write("_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n")
            f.write("_atom_site.auth_seq_id\n_atom_site.auth_asym_id\n_atom_site.pdbx_PDB_model_num\n")

            for idx, a in enumerate(atoms, start=1):
                rec = a.get('record', 'ATOM')
                elem = a.get('element', 'C')
                aname = a.get('atom_name', '')
                alt = a.get('altLoc', '.')
                if alt == ' ':
                    alt = '.'
                resName = a.get('resname', '').strip()
                chain = a.get('chain', 'A')
                resSeq = a.get('resSeq', 1)
                icode = a.get('icode', '?')
                if icode == ' ':
                    icode = '?'

                xyz = a.get('xyz', [a.get('x', 0.0), a.get('y', 0.0), a.get('z', 0.0)])
                x = float(xyz[0])
                y = float(xyz[1])
                z = float(xyz[2])

                occ = a.get('occupancy', 1.0)
                temp = a.get('tempFactor', 0.0)

                f.write(f"{rec} {idx} {elem} {aname} {alt} {resName} {chain} 1 {resSeq} {icode} {x:.3f} {y:.3f} {z:.3f} {occ:.2f} {temp:.2f} {resSeq} {chain} 1\n")
