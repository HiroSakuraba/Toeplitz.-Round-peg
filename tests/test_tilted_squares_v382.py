from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from smooth_squares_v379 import classify  # noqa: E402
from tilted_squares_v382 import arc_residual, h_value, predicted_radius, scan  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Branches(unittest.TestCase):
    def test_every_integer_in_the_range_is_a_square(self):
        for n1, n3 in ((2, 3), (-2, 3), (4, 5)):
            Ths, vals, sq = scan(n1, n3, 6.0, 600.0, steps=1500, rmax=3 ** -0.5)
            v = [x for x in vals if x is not None]
            d = np.diff(v)
            self.assertTrue((d >= -1e-12).all() or (d <= 1e-12).all())        # h is monotone
            lo, hi = sorted((min(v), max(v)))
            want = {j for j in range(int(np.ceil(lo)), int(np.floor(hi)) + 1) if j != 0}
            self.assertEqual({s['j'] for s in sq}, want)
            for s in sq:
                self.assertLess(arc_residual(s['corners']), 1e-10)
                self.assertGreater(s['j'] * n1 * n3, 0)
                Z = s['corners']
                c = classify(np.array([Z[0].real, Z[0].imag, np.angle((Z[1] - Z[0]) / s['side']), s['side']]))
                self.assertEqual(c[:3], (n1, n3, s['j']))                       # v379's labelling agrees

    def test_h_tends_to_zero_deep_in(self):
        v, _ = h_value(2000.0, 3, 3)
        self.assertGreater(v, 0)
        self.assertLess(v, 0.2)

    def test_recorded(self):
        data = json.loads((ROOT / 'results/v382/tilted_check.json').read_text())
        self.assertEqual(len(data['branches']), 252)
        self.assertTrue(all(b['monotone'] for b in data['branches']))
        sq = data['squares']
        self.assertEqual(len(sq), 807)
        self.assertLess(max(s['residual'] for s in sq), 1e-10)
        self.assertTrue(all(s['j'] * s['n1'] * s['n3'] > 0 for s in sq))
        acc = data['prediction_accuracy']
        self.assertLess(acc['opposite_sign t/r in [0,0.1)']['median'], 0.02)
        self.assertGreater(acc['same_sign t/r in [0,0.1)']['median'], 0.05)
        self.assertEqual(data['v379_found_on_branches']['at_max_ratio_0_32'], data['v379_squares_in_range'])


if __name__ == '__main__':
    unittest.main()
