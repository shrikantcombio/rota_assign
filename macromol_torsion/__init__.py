"""
Macromolecular Torsion Angle Calculation Package (macromol_torsion)
=====================================================================
Unified package for calculating Protein and RNA torsion angles,
RNA 7D dinucleotide suites, RNA pseudo-torsions, sugar pucker pseudorotations,
glycosidic side-chains, DSSP secondary structure, NACCESS SASA interface classification,
and comparative analysis of Protein-RNA complex structures (PDB & CIF).
"""

from macromol_torsion.geometry import PDBGeometry
from macromol_torsion.structure_reader import StructureReader
from macromol_torsion.protein_backbone import ProteinBackboneCalculator
from macromol_torsion.protein_sidechain import ProteinSidechainCalculator
from macromol_torsion.dssp import DsspAnalyzer
from macromol_torsion.sasa import SasaAnalyzer

from macromol_torsion.rna_backbone import RNABackboneCalculator
from macromol_torsion.rna_pseudo import RNAPseudoTorsionCalculator
from macromol_torsion.rna_sugar import RNASugarPuckerCalculator
from macromol_torsion.rna_sidechain import RNASidechainCalculator

from macromol_torsion.calculator import ComplexTorsionCalculator

__all__ = [
    "PDBGeometry",
    "StructureReader",
    "ProteinBackboneCalculator",
    "ProteinSidechainCalculator",
    "DsspAnalyzer",
    "SasaAnalyzer",
    "RNABackboneCalculator",
    "RNAPseudoTorsionCalculator",
    "RNASugarPuckerCalculator",
    "RNASidechainCalculator",
    "ComplexTorsionCalculator",
]
