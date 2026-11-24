from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from invisibility_v377 import arc_function, near_squares, turn_spacing  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class NearSquares(unittest.TestCase):
    def test_free_corner_within_one_spacing(self):
        s = 0.01
        for lam in (0, 1):
            sq = near_squares(s, lam, n_theta=8000)
            self.assertTrue(sq)
            d = turn_spacing(s, lam)
            A, th, res = min(sq, key=lambda q: abs(q[0]))
            self.assertLess(abs(A), 0.6 * d)
            self.assertLess(res, 1e-9)
            # it is a square of side s whose corners B, C, D lie on the arc (arc function ~ 0)
            e = np.exp(1j * th)
            for z in (A + s * e, A + (1 + 1j) * s * e, A + 1j * s * e):
                self.assertLess(abs(arc_function(z, lam)), 1e-9)
                self.assertGreater(abs(z) ** -2, 3.0)   # inside the arc (Theta >= theta0)

    def test_arc_function_vanishes_on_both_arms(self):
        for lam in (0, 1):
            for Th in (50.0, 400.0):
                r = Th ** -0.5
                ph = lam * 0.4 * np.sin(5 * np.log(1 / r))
                for z in (r * np.exp(1j * Th), r * np.exp(1j * (Th + pi + ph))):
                    self.assertLess(abs(arc_function(z, lam)), 1e-9)

    def test_recorded(self):
        data = json.loads((ROOT / 'results/v377/invisibility_check.json').read_text())
        ratios = [x for row in data['near_squares'] for x in row['free_corner_over_spacing']]
        ratios += data['near_squares_neighbour_cell']['free_corner_over_spacing']
        self.assertTrue(ratios)
        self.assertLess(max(ratios), 1.42)   # below the sketch's bound sqrt(2) for g = pi
        scan = [r for r in data['depth_scan'] if isinstance(r['d5'], list)]
        self.assertEqual(len(scan), 24)
        self.assertEqual(sum(r['d5'][0] != r['d7'][0] for r in scan), 9)


if __name__ == '__main__':
    unittest.main()
