from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import census_polygon, pt
from scale_transfer_v310 import (
    audit_switch, flux_cover, integration_example, multi_patch_cover,
    parity_intervals, seed, stretch, stretch_census, triple_bound,
)


class ScaleTransferTests(unittest.TestCase):
    def test_continuous_parameter_audits(self):
        self.assertEqual(stretch_census(), 13)
        audit_switch()

    def test_independent_orbit_census(self):
        for d in (Q(1, 100), Q(1, 1000)):
            shear, = census_polygon(seed(d)).squares
            vertical, = census_polygon(stretch(d)).squares
            self.assertGreater(shear.side_squared, Q(39, 100))
            self.assertLess(sum(p.x for p in shear.vertices)/4, Q(47, 100))
            self.assertEqual(vertical.side_squared, (9*(1+d)/(17+9*d))**2)
            self.assertEqual(sum(p.x for p in vertical.vertices)/4, Q(1, 2))

    def test_odd_profile_counts_multiplicity(self):
        sizes = [Q(2), Q(2), Q(3), Q(5)]
        self.assertEqual(parity_intervals(sizes, (Q(1), Q(4))), [(Q(3), Q(4))])

    def test_dependent_triple_and_small_sizes(self):
        aligned = (((1, 0), 1), ((0, 1), 1), ((1, 1), 1))
        self.assertIsNone(triple_bound(aligned, pt(0, 0), Q(1, 10)))
        self.assertEqual(flux_cover([], pt(0, 0), Q(1, 4), (0, 1)), [(0, Q(1, 2))])

    def test_two_box_cover_includes_diagonal_and_side(self):
        boxes = [(pt(0, 0), Q(0)), (pt(3, 4), Q(0))]
        self.assertEqual(multi_patch_cover([], boxes, (1, 30)),
                         [(Q(25, 2), Q(25, 2)), (25, 25)])

    def test_triangle_replacements(self):
        data = integration_example()
        self.assertLess(Q(data["lost_measure_upper"]), Q(1, 2000))
        self.assertGreater(Q(data["squared_side_floor"]), Q(29, 10))
        lo, hi = map(Q, data["four_patch_squared_side_enclosure"])
        self.assertLess(Q(2937, 1000), lo)
        self.assertTrue(lo < Q(144, 49) < hi)
        self.assertLess(hi, Q(2941, 1000))
        moved = Q(data["four_contact_replacement_squared_side"])
        self.assertNotEqual(moved, Q(144, 49))
        self.assertTrue(lo <= moved <= hi)
