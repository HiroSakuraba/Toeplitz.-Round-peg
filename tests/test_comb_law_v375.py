from fractions import Fraction as Q
from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from asymmetric_spiral_v374 import asymmetric_drawing  # noqa: E402
from comb_law_v375 import check_drawing, confinement_radius, hypotheses, spacing_bound  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Bound(unittest.TestCase):
    def test_spacing_bound_is_increasing_and_inverted(self):
        rs = [0.05 * k for k in range(1, 12)]
        ds = [spacing_bound(r) for r in rs]
        self.assertEqual(ds, sorted(ds))
        for t in (1e-4, 0.01, 0.05):
            self.assertAlmostEqual(spacing_bound(confinement_radius(t)), t, places=12)

    def test_asymptotic_comb_radius(self):
        for g in (pi, pi - 0.4):
            t = 1e-8
            self.assertAlmostEqual(confinement_radius(t, g=g) / (t / (0.5 * g)) ** (1 / 3), 1, places=3)

    def test_hypotheses(self):
        self.assertTrue(hypotheses(0.05, confinement_radius(0.05), L=0)['ok'])
        self.assertTrue(hypotheses(0.05, confinement_radius(0.05, g=pi - 0.4), L=2)['ok'])
        self.assertFalse(hypotheses(0.05, 0.1, L=2)['ok'])


class LemmaA(unittest.TestCase):
    def test_non_neighbours_on_the_continuous_arc_are_spaced(self):
        for lam in (0, 1):
            g = pi - 0.4 * lam
            Th = np.linspace(3.0, 3.0 + 8 * pi, 2500)
            r = Th ** -0.5
            Z = np.r_[r * np.exp(1j * Th), r * np.exp(1j * (Th + pi + lam * 0.4 * np.sin(5 * np.log(1 / r))))]
            T, arm, R = np.r_[Th, Th], np.r_[np.zeros(len(Th)), np.ones(len(Th))], np.abs(Z)
            for i in range(0, len(Z), 7):
                d = np.abs(Z - Z[i])
                m = ~((arm == arm[i]) & (np.abs(T - T[i]) <= pi / 2)) & (d < 0.08)
                if m.any():
                    self.assertGreater((d[m] / spacing_bound(np.minimum(R, R[i])[m], 0.5, g)).min(), 1)


class Drawings(unittest.TestCase):
    def test_every_small_square_has_a_corner_inside_the_bound(self):
        for lam in (0, 1):
            rows = check_drawing(asymmetric_drawing(4, Q(lam)), Q(1, 20), g=pi - 0.4 * lam)
            self.assertTrue(all(m <= rho for _, m, rho in rows))
            top = max(rows, key=lambda row: row[1])
            self.assertTrue(0.85 < top[1] / top[2] < 0.95)

    def test_recorded_check(self):
        data = json.loads((ROOT / 'results/v375/comb_bound_check.json').read_text())
        self.assertEqual(sum(r['squares'] for r in data['drawings']), 464)
        for r in data['drawings']:
            self.assertEqual(r['violations'], 0)
            self.assertTrue(r['hypotheses_ok'])
            self.assertLess(r['worst_ratio'], 1)


if __name__ == '__main__':
    unittest.main()
