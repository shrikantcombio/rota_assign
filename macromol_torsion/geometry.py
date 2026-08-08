"""
PDBGeometry Module
------------------
Provides 3D vector geometry calculations: distance, dot product, cross product,
bond angles, and signed dihedral torsion angles.
"""

import math
from math import sqrt, pi, acos, degrees, atan2


class PDBGeometry:
    @staticmethod
    def distance(a, b):
        """Calculates Euclidean distance between 3D points a and b."""
        return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))

    @staticmethod
    def sq_distance(a, b):
        """Calculates squared distance between 3D points a and b."""
        return sum((a[i] - b[i]) ** 2 for i in range(3))

    @staticmethod
    def dot_product(a, b):
        """Calculates dot product of 3D vectors a and b."""
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    @staticmethod
    def cross_product(a, b):
        """Calculates cross product of 3D vectors a and b."""
        return [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ]

    @staticmethod
    def vector_magnitude(a):
        """Calculates magnitude of a 3D vector."""
        return math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)

    @classmethod
    def get_angle(cls, a, b):
        """Calculates angle between two 3D vectors in degrees."""
        mag_a = cls.vector_magnitude(a)
        mag_b = cls.vector_magnitude(b)
        if mag_a == 0.0 or mag_b == 0.0:
            return 0.0

        dot = cls.dot_product(a, b)
        cos_theta = max(-1.0, min(1.0, dot / (mag_a * mag_b)))
        return math.degrees(math.acos(cos_theta))

    @classmethod
    def get_dihedral(cls, atm1, atm2, atm3, atm4):
        """
        Calculates signed dihedral angle between four atom coordinates (atm1, atm2, atm3, atm4).
        Follows standard IUPAC sign convention (-180.0 to 180.0 degrees).
        Matches Bio.PDB.calc_dihedral exactly.
        """
        b1 = [atm2[i] - atm1[i] for i in range(3)]
        b2 = [atm3[i] - atm2[i] for i in range(3)]
        b3 = [atm4[i] - atm3[i] for i in range(3)]

        n1 = cls.cross_product(b1, b2)
        n2 = cls.cross_product(b2, b3)

        angle = cls.get_angle(n1, n2)

        if cls.dot_product(cls.cross_product(n1, n2), b2) < 0:
            angle = -angle

        return round(float(angle), 2)

    @classmethod
    def tau_angle(cls, n_coord, ca_coord, c_coord):
        """
        Calculates tau bond angle N-CA-C in degrees.
        """
        a = [n_coord[i] - ca_coord[i] for i in range(3)]
        b = [c_coord[i] - ca_coord[i] for i in range(3)]
        return round(float(cls.get_angle(a, b)), 2)
