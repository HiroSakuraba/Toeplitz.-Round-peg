import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from liouville_primitive_experiments_v276 import (
    action_weight, block_action, block_sizes, block_vertices, closed_truncation,
    signed_action,
)


class PrimitiveExperimentsTests(unittest.TestCase):
    def test_block_ranges_agree_with_independent_edge_summation(self):
        for scale in range(5):
            for counterwind in (False, True):
                points = list(block_vertices(scale, counterwind))
                values = [Q(0)]
                for a, b in zip(points, points[1:]):
                    values.append(values[-1] + (a[0]*b[1]-a[1]*b[0])/2)
                self.assertEqual(block_action(scale, counterwind),
                                 (values[-1], min(values), max(values)))

    def test_subinterval_weight(self):
        _, burst = block_sizes(3)
        points = list(block_vertices(3, False))
        self.assertEqual(action_weight(3, burst+1, 2*burst),
                         signed_action(points[burst+1:2*burst+1]))

    def test_shrinking_blocks_have_very_different_action_scales(self):
        for scale in range(11):
            self.assertGreaterEqual(block_action(scale, False)[2], Q(1, 2))
            _, low, high = block_action(scale, True)
            self.assertGreaterEqual(low, 0)
            self.assertLessEqual(high, Q(2, 2**scale))

    def test_both_closed_truncations_have_fixed_area(self):
        for blocks in range(1, 6):
            for counterwind in (False, True):
                self.assertEqual(signed_action(closed_truncation(blocks, counterwind)), 1)

    def test_strict_radial_coordinate_and_exact_splices(self):
        for scale in range(4):
            for counterwind in (False, True):
                points = list(block_vertices(scale, counterwind))
                radii = [abs(x)+abs(y) for x, y in points]
                self.assertTrue(all(a > b for a, b in zip(radii, radii[1:])))
                self.assertEqual(points[-1], (Q(1, 2**(scale+1)), 0))
