from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from roof_shear_v316 import (band_intervals, digital_trapezoid, folded_example,
                              gate_certificate, pt, terrain_certificate)


class RoofShearTests(unittest.TestCase):
    def test_contained_optimum_against_interval_enumeration(self):
        for heights in product(range(1, 5), repeat=5):
            cert = terrain_certificate(heights)
            exact = max(min(j-i, min(heights[i:j]))
                        for i in range(5) for j in range(i+1, 6))
            self.assertEqual(cert["contained_side"], exact)
            self.assertGreaterEqual(cert["boundary_side"], exact)

    def test_sharp_and_disconnected_height_sections(self):
        self.assertEqual(terrain_certificate([3]*8)["boundary_side"], 3)
        self.assertEqual(terrain_certificate([1,1,1,2,1,2]),
                         {"contained_side":1, "boundary_side":2, "base_left":3})
        for invalid in ([], [1,0], [-1], [1.5]):
            with self.assertRaises(ValueError):
                terrain_certificate(invalid)

    def test_clips_low_band_interior_of_long_edges(self):
        got = band_intervals([pt(0,1), pt(10,6), pt(8,2)], Q(3))
        self.assertEqual(got, [(0,4),(8,Q(17,2))])

    def test_folded_base_independent_census(self):
        cert = folded_example()
        self.assertEqual(cert["side_floor"], "7/4")
        self.assertEqual(cert["regular_marked_squares"], 1)
        self.assertGreater(Q(cert["squared_sides"][0]), 25)

    def test_mesh_bound(self):
        for n in (4,5,11,32):
            upper, base = digital_trapezoid(n)
            cert = gate_certificate(upper, base, Q(1,2*n))
            self.assertEqual(Q(cert["side_floor"]), 5-Q(5,2*n))

    def test_gate_failures(self):
        u = [pt(1,2),pt(2,5),pt(8,5),pt(9,2)]
        b = [pt(-1,0),pt(12,0)]
        self.assertEqual(gate_certificate(u,b,Q(0))["side_floor"], "2")
        with self.assertRaisesRegex(ValueError,"endpoints"):
            gate_certificate(u,[pt(1,0),pt(12,0)],Q(1,4))
        with self.assertRaisesRegex(ValueError,"strip"):
            gate_certificate(u,[pt(-1,1),pt(12,0)],Q(1,4))
        with self.assertRaisesRegex(ValueError,"two-gate"):
            gate_certificate([pt(1,2),pt(5,Q(9,4)),pt(9,2)],b,Q(1,4))
