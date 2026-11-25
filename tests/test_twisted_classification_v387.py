from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from general_comb_v380 import corners, search  # noqa: E402
from twisted_classification_v387 import (  # noqa: E402
    W1_of, aligned_bounds, aligned_margin_leading, arm_two_curvature_ok, default_phi_derivs, default_profile,
    diagonal_ok, on_branch, ratio_limit_211, three_one_inside_threshold, three_one_outside_ok)

ROOT = Path(__file__).resolve().parents[1]
P = default_profile()


class Criteria(unittest.TestCase):
    def test_window_criteria(self):
        W = 0.3
        self.assertTrue(diagonal_ok(W, P["Phi"], P["L"]))
        self.assertTrue(three_one_outside_ok(W, P["Phi"], P["L"]))
        self.assertAlmostEqual(three_one_inside_threshold(W, P["Phi"], P["L"]), 6.41, places=1)
        self.assertGreater(ratio_limit_211(W, P["Phi"], P["L"]), 1.19)
        self.assertGreater(W1_of(W, P["L"]), W)

    def test_aligned_beyond_threshold(self):
        for n in (1, 2, 3, 7):
            cases = ("outer", "inner") if n % 2 else ("both", "none")
            for c in cases:
                self.assertGreater(aligned_margin_leading(n, 0.3, P["Phi"], P["L"], 0.5, c), 0)
                for Lf in np.geomspace(6.41, 1e7, 200):
                    lo, hi = aligned_bounds(Lf, n, 0.3, P["Phi"], P["L"], 0.5, c)
                    self.assertGreater(lo, hi)

    def test_steep_twist_breaks_the_margin(self):
        self.assertLess(aligned_margin_leading(1, 0.3, P["Phi"], 6.0), 0)

    def test_arm_two_convex(self):
        phi, dphi, d2phi = default_phi_derivs()
        for T in np.geomspace(3, 1e7, 3000):
            self.assertTrue(arm_two_curvature_ok(T, phi, dphi, d2phi))


class Completeness(unittest.TestCase):
    def test_search_squares_lie_on_branches(self):
        sq = search(0.0, 1.0, n_random=600, seed=7)
        self.assertTrue(sq)
        for v, _ in sq:
            self.assertTrue(on_branch(corners(v))[1])


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v387" / "classification_check.json").read_text())
        self.assertTrue(all(v == 0 for v in d["labelled_lam1"].values()))
        self.assertTrue(set(d["random_lam1"]["counts"]) <= {"1-1-1-1", "2-1-1 diagonal"})
        c = d["v382_completeness"]
        self.assertEqual(c["squares"], c["on_branch"])
        self.assertLess(d["criteria"]["0.3"]["aligned_from"], 6.41)


if __name__ == "__main__":
    unittest.main()
