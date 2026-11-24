from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from comb_law_v375 import confinement_radius  # noqa: E402
from smooth_squares_v379 import (  # noqa: E402
    aligned_near_squares, classify, corners, min_side, outer_radius, phase_rel, predicted_radius, snap,
    solve_square)

ROOT = Path(__file__).resolve().parents[1]


class Rule(unittest.TestCase):
    def test_predicted_square_is_found(self):
        # (n1, n3, j) = (2, 3, 1): predicted centre radius 0.1629, side about 0.0245
        r = predicted_radius(2, 3, 1)
        self.assertAlmostEqual(r, 0.1629, places=3)
        hits = []
        D = pi * 0.5 * r ** 3
        for sgn in (1, -1):
            beta = sgn * np.arctan2(2, 3)
            for psi in np.linspace(0, 2 * pi, 48, endpoint=False):
                A = snap(r * np.exp(1j * psi))
                er = A / abs(A)
                d = np.cos(beta) * 1j * er + np.sin(beta) * er
                v = solve_square(np.array([A.real, A.imag, np.angle(d), D * 13 ** 0.5]))
                if v is not None:
                    hits.append(classify(v))
        self.assertTrue(any(abs(c[4] / r - 1) < 0.02 and c[2] * c[0] * c[1] > 0 for c in hits))

    def test_recorded_smooth_squares(self):
        # no aligned squares, sign rule, outer radius, and accuracy for small squares
        data = json.loads((ROOT / 'results/v379/smooth_comb_check.json').read_text())
        rows = data['smooth_squares']
        self.assertGreaterEqual(len(rows), 150)
        self.assertTrue(all(r['j'] != 0 and r['n1'] * r['n3'] != 0 for r in rows))
        self.assertTrue(all(r['j'] * r['n1'] * r['n3'] > 0 for r in rows))
        self.assertLess(max(r['centre'] / r['outer'] for r in rows), 1.01)
        self.assertLess(max(r['min_corner'] / confinement_radius(r['side']) for r in rows), 0.65)
        small = [r for r in rows if r['side'] / r['centre'] < 0.15]
        self.assertLess(np.median([abs(r['centre'] / r['pred'] - 1) for r in small]), 0.05)

    def test_outer_radius_inverts_min_side(self):
        for r in (0.1, 0.2, 0.3):
            self.assertAlmostEqual(outer_radius(min_side(r)), r, places=12)


class Aligned(unittest.TestCase):
    def test_mismatch_keeps_its_sign(self):
        for sb, sd in ((1, 1), (1, -1)):
            ms = [m for TA in np.linspace(10, 60, 60) for t, m in aligned_near_squares(TA, 1, sb, sd) if t < 0.08]
            self.assertTrue(ms)
            self.assertTrue(all(m > 0 for m in ms) or all(m < 0 for m in ms))
            self.assertLess(max(abs(m) for m in ms), 0.5)

    def test_recorded(self):
        data = json.loads((ROOT / 'results/v379/smooth_comb_check.json').read_text())
        for row in data['aligned_near_squares']:
            self.assertTrue(row['m_min'] > 0 or row['m_max'] < 0)
            self.assertLess(max(abs(row['m_min']), abs(row['m_max'])), 0.5)
        for s, row in data['drawing_kinks'].items():
            self.assertEqual(row['0.12'], row['total'])


if __name__ == '__main__':
    unittest.main()
