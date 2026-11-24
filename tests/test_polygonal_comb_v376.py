from fractions import Fraction as Q
from math import pi, sqrt
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from asymmetric_spiral_v374 import assert_jordan, asymmetric_drawing  # noqa: E402
from comb_law_v375 import spacing_bound  # noqa: E402
from polygonal_comb_v376 import (  # noqa: E402
    analyse_square, curvature_factor, drawing_from_thetas, edge_data, gamma, refined_thetas)
from wild_point_parity_v370 import sector_counts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Theorem(unittest.TestCase):
    def test_stars_and_bound_on_depth_four(self):
        for lam in (0, 1):
            V = asymmetric_drawing(4, Q(lam))
            edges, nu = edge_data(V, 4, lam)
            self.assertLess(nu, 1e-7)
            _, _, sq = sector_counts(V, Q(1, 20))
            applied = 0
            for e in sq.values():
                r = analyse_square(e, edges, nu, V, lam)
                if r is None:
                    continue
                # every corner is within h of its star on the smooth arc
                for p, q in zip(e[1], e[2]):
                    k = int(q)
                    arm, T0, T1 = edges[k]
                    z = gamma(arm, T0 + float(q - k) * (T1 - T0), lam)
                    self.assertLessEqual(abs(complex(float(p[0]), float(p[1])) - z), r['h'])
                self.assertTrue(r['holds'])
                applied += r['ok']
            self.assertEqual(applied, {0: 56, 1: 34}[lam])

    def test_uniform_drawing_is_reproduced(self):
        T = [3.0 + 4 * pi * k / 48 for k in range(49)]
        self.assertEqual(drawing_from_thetas(T, 1), asymmetric_drawing(4))


class Mesh(unittest.TestCase):
    def ratio(self, lam, Ta, Tb):
        r = Ta ** -0.5
        h = (Tb - Ta) ** 2 / 8 * r * curvature_factor(0.5 * r * r, 2 * lam, lam * sqrt(104)) + 1e-7
        return 2 * h / spacing_bound(Tb ** -0.5, 0.5, pi - 0.4 * lam)

    def test_uniform_mesh_loses_the_criterion(self):
        for lam, depth in ((0, 28), (1, 24)):
            first = next(d for d in range(4, 80)
                         if self.ratio(lam, 3 + d * pi, 3 + d * pi + pi / 12) >= 1)
            self.assertEqual(first, depth)

    def test_refined_mesh_keeps_it_and_stays_simple(self):
        T = refined_thetas(100)
        for lam, cap in ((0, 0.07), (1, 0.15)):
            self.assertLess(max(self.ratio(lam, a, b) for a, b in zip(T, T[1:])), cap)
        assert_jordan(drawing_from_thetas(refined_thetas(75), 1))

    def test_recorded_check(self):
        data = json.loads((ROOT / 'results/v376/polygonal_comb_check.json').read_text())
        for r in data['drawings']:
            self.assertTrue(r['holds_all'])
            self.assertEqual(r['off_arms'], 2)
        self.assertEqual(sum(r['hypotheses_ok'] for r in data['drawings']), 414)


if __name__ == '__main__':
    unittest.main()
