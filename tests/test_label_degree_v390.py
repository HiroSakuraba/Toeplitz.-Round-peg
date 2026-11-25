from pathlib import Path
import json
import sys
import unittest
from math import pi

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from general_comb_v380 import corners, search  # noqa: E402
from label_degree_v390 import (anchored_corners, anchored_from_square, anchored_jacobian_sign,  # noqa: E402
                               anchored_map, chart_check, partner_anchored_label)
from two_turn_squares_v383 import window  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ZERO = lambda r: 0.0  # noqa: E731


class AnchoredLabels(unittest.TestCase):
    def setUp(self):
        self.sq = search(0.0, 1.0, n_random=500, seed=7)
        self.assertTrue(self.sq)

    def test_coordinates_reproduce_the_square(self):
        for v, Ths in self.sq:
            Z = corners(v)
            _, k, _ = window(Z, Ths)
            for i in range(4):
                x, lab = anchored_from_square(Z, list(map(int, k)), i)
                W = anchored_corners(x, lab[0], 0.0, ZERO)
                self.assertTrue(np.allclose(W, [Z[(i + q) % 4] for q in range(4)], atol=1e-12))

    def test_zeros_and_local_degree(self):
        for v, Ths in self.sq:
            Z = corners(v)
            _, k, _ = window(Z, Ths)
            for i in range(4):
                x, lab = anchored_from_square(Z, list(map(int, k)), i)
                self.assertLess(np.max(np.abs(anchored_map(x, lab, 0.0, ZERO))), 1e-9)
                self.assertNotEqual(anchored_jacobian_sign(x, lab, 0.0, ZERO)[0], 0)

    def test_partners_have_other_labels_and_are_zeros(self):
        for v, Ths in self.sq:
            Z = corners(v)
            _, k, _ = window(Z, Ths)
            x, lab = anchored_from_square(Z, list(map(int, k)), 0)
            plab = partner_anchored_label(lab)
            self.assertNotEqual(plab, lab)
            self.assertLess(np.max(np.abs(anchored_map(x, plab, 0.0, ZERO))), 1e-9)
            P = anchored_corners(x, plab[0], 0.0, ZERO)
            self.assertTrue(np.allclose(P, [-z for z in Z], atol=1e-12))     # it is the half-turn image

    def test_no_branch_cut(self):
        c = chart_check()
        self.assertLess(c["anchored_jump"], 1e-6)
        self.assertAlmostEqual(c["absolute_jump"], 2 * pi, places=6)

    def test_continuous_in_lift_through_odd_multiples_of_pi(self):
        lab = (0, 2, 3, 1)
        w = 0.05 * (1 - 1j) / 2
        T = np.linspace(41 * pi - 1e-3, 41 * pi + 1e-3, 201)
        vals = np.array([anchored_map(np.array([t, w.real, w.imag]), lab, 0.0, ZERO) for t in T])
        self.assertLess(np.max(np.abs(np.diff(vals, axis=0))), 1e-3)


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v390" / "label_check.json").read_text())
        p = d["plain"]
        self.assertEqual(p["configurations"], 4 * p["squares"])
        self.assertEqual(p["partner_same_label"], 0)
        self.assertEqual(p["partner_is_zero_of_partner_label"], p["configurations"])
        self.assertEqual(p["nondegenerate"], p["configurations"])
        self.assertLess(p["max_residual"], 1e-9)
        u = d["uniqueness"]
        self.assertEqual(u["labels_shared_by_two_squares"], 0)
        self.assertEqual(u["anchored_labels"], 4 * u["squares"])
        self.assertEqual(u["on_branch"], u["squares"])


if __name__ == "__main__":
    unittest.main()
