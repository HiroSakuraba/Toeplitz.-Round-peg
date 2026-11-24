from fractions import Fraction as Q
from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from general_comb_v380 import corners, radius, solve_square  # noqa: E402
from two_curve_v381 import (  # noqa: E402
    Graph, H, hypothesis_margin, pass_frame, square, zeros)

ROOT = Path(__file__).resolve().parents[1]


def graph(f, lo=-3.0, hi=3.0, n=60001):
    X = np.linspace(lo, hi, n)
    return Graph(X, f(X))


class Reduction(unittest.TestCase):
    def test_flat_floor_square_is_centred_on_the_gap_maximum(self):
        # f1 = 0, f2 = 1/10 (1 + cos x)/2 + 1/20: with s = 0, R = 0 and squares solve g(x1 + u) = g(x1)
        f1 = graph(lambda x: 0 * x)
        f2 = graph(lambda x: 0.05 + 0.05 * (1 + np.cos(x)))
        zs = zeros(f1, f2, np.linspace(-1, 1, 401), 0.5)
        self.assertEqual(len(zs), 1)
        x1, u = zs[0]
        self.assertAlmostEqual(x1, -u / 2, places=6)
        self.assertAlmostEqual(u, 0.05 + 0.05 * (1 + np.cos(u / 2)), places=6)
        P = square(f1, f2, x1, u)
        self.assertAlmostEqual(P[2][1], f2(P[2][0]), places=8)
        self.assertAlmostEqual(P[3][1], f2(P[3][0]), places=8)

    def test_localization_bound_on_curved_curves(self):
        f1 = graph(lambda x: 0.08 * x * x)
        f2 = graph(lambda x: 0.3 + 0.1 * x * x + 0.02 * np.sin(3 * x))
        for x1, u in zeros(f1, f2, np.linspace(-1.5, 1.2, 801), 0.6):
            dg, B = hypothesis_margin(f1, f2, x1, 0.6)
            h, _ = H(f1, f2, x1, 0.6)
            self.assertLessEqual(abs(dg), B + abs(h) + 1e-12)


class Spirals(unittest.TestCase):
    def test_square_at_a_stationary_gap_and_none_on_the_plain_spiral(self):
        kappa, omega = 0.3, 1.0
        rf = lambda T: radius(T, kappa, omega)
        T = np.linspace(10.5, 12.0, 30001)          # contains one interior gap minimum (near 11.18)
        g = rf(T) - rf(T + pi)
        Tc = T[np.argmin(g)]
        f1, f2, P = pass_frame(rf, Tc, 0.0)
        g0 = float(f2(0) - f1(0))
        zs = zeros(f1, f2, np.linspace(-3 * g0, 3 * g0, 241), 2.5 * g0)
        self.assertEqual(len(zs), 1)
        Z = [P(q) for q in square(f1, f2, *zs[0])]
        e = (Z[0] - Z[1]) / abs(Z[0] - Z[1])
        out = solve_square(np.array([Z[1].real, Z[1].imag, np.angle(e), abs(Z[0] - Z[1])]), kappa, omega)
        self.assertIsNotNone(out)
        self.assertLess(abs(corners(out[0])[0] - Z[1]), 1e-6)
        rp = lambda T: radius(T, 0.0, 1.0)
        f1, f2, P = pass_frame(rp, 13.0, 0.0)
        g0 = float(f2(0) - f1(0))
        self.assertEqual(zeros(f1, f2, np.linspace(-3 * g0, 3 * g0, 241), 2.5 * g0), [])

    def test_recorded(self):
        data = json.loads((ROOT / 'results/v381/two_curve_check.json').read_text())
        for key, r in data['existence_at_stationary_gaps'].items():
            self.assertEqual(r['zeros'], r['polished'])
            self.assertGreaterEqual(r['windows_with_zero'], r['certified_minima'] + r['certified_maxima'])
            if r['stationary'] == 0:
                self.assertEqual(r['zeros'], 0)
        for key, r in data['localization_on_v380_aligned_squares'].items():
            self.assertEqual(r['localization_holds'], r['aligned'])
            self.assertLess(r['max_dg_minus_residual_over_bound'], 1 + 1e-8)
        for key, r in data['polygon_drawings_H_zeros_vs_exact'].items():
            self.assertEqual(r['aligned_exact'], r['exact_matched'])
            self.assertEqual(r['predicted_unmatched'], 0)


if __name__ == '__main__':
    unittest.main()
