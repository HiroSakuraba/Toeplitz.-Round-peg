from fractions import Fraction as Q
from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from asymmetric_spiral_v374 import asymmetric_drawing  # noqa: E402
from general_comb_v380 import (  # noqa: E402
    check_theorem, gap_criticality, polygon_far_sides, radius, search)
from wild_point_parity_v370 import sector_counts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Theorem(unittest.TestCase):
    def test_two_quarter_turn_sides_on_smooth_families(self):
        for kappa, omega, expect_aligned in ((0.0, 1.0, False), (0.3, 1.0, True)):
            sq = search(kappa, omega, n_random=1500, seed=4)
            self.assertTrue(sq)
            far = [check_theorem(v, Ths, kappa, omega) for v, Ths in sq]
            self.assertGreaterEqual(min(far), 2)
            self.assertEqual(any(f == 2 for f in far), expect_aligned)
            for (v, Ths), f in zip(sq, far):
                if f == 2:
                    side_over_gap, crit = gap_criticality(v, Ths, kappa, omega)
                    self.assertLess(crit, 0.01)
                    self.assertLess(abs(side_over_gap - 1), 0.05)

    def test_on_a_polygon_drawing(self):
        V = asymmetric_drawing(4, Q(1))
        _, _, sq = sector_counts(V, Q(1, 20))
        far = [polygon_far_sides(V, e) for e in sq.values()]
        self.assertGreaterEqual(min(far), 2)

    def test_family_parameters(self):
        Th = np.linspace(3, 400, 100000)
        for kappa, omega, backs_up in ((0.3, 1.0, False), (0.5, 3.0, True)):
            self.assertGreater((radius(Th, kappa, omega) - radius(Th + pi, kappa, omega)).min(), 0)  # Jordan
            self.assertEqual(bool((np.diff(radius(Th, kappa, omega)) > 0).any()), backs_up)

    def test_recorded(self):
        data = json.loads((ROOT / 'results/v380/general_comb_check.json').read_text())
        for key, s in data['summary'].items():
            self.assertNotIn('0', s['far_side_counts'])
            self.assertNotIn('1', s['far_side_counts'])
            self.assertGreater(s['min_side_over_quarter_gap'], 0.999)      # sampling tolerance for G
            has_aligned = '2' in s['far_side_counts']
            self.assertEqual(has_aligned, s['gap_critical_points'] > 0)
        for key, c in data['aligned_vs_gap_criticality'].items():
            self.assertLess(c['crit_max'], 0.01)
        for key, counts in data['polygon_drawings_far_sides'].items():
            self.assertTrue(all(int(k) >= 2 for k in counts))


if __name__ == '__main__':
    unittest.main()
