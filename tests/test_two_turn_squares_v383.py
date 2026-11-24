from math import pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from general_comb_v380 import radius, search  # noqa: E402
from two_turn_squares_v383 import (  # noqa: E402
    classify_found, cond_large_kappa, cond_small_kappa, consts, exact_q_gap, flag_residual, labelled_search,
    p_bound, ratio_limit_211, side_limit_211, split_type, threshold, x_lower)

ROOT = Path(__file__).resolve().parents[1]


class Identities(unittest.TestCase):
    def test_british_flag_on_random_rectangles(self):
        rng = np.random.default_rng(0)
        for _ in range(200):
            A = complex(*rng.normal(size=2))
            e = np.exp(1j * rng.uniform(0, 2 * pi))
            w, h = rng.uniform(0.1, 2, size=2)
            Z = [A, A + w * e, A + w * e + 1j * h * e, A + 1j * h * e]
            self.assertLess(abs(flag_residual(Z)), 1e-12)

    def test_angle_ratio_lower_bound(self):
        rng = np.random.default_rng(1)
        W = 1.0
        g, T = consts(W)
        for _ in range(5000):
            p, q = rng.uniform(1e-4, W, size=2)        # |phi| and |phi'|
            s_true = np.sin(q / 2) / np.sin(p / 2)
            s = s_true * rng.uniform(0.5, 1.0)          # any valid lower bound for the sine ratio
            self.assertGreaterEqual(q / p, x_lower(s, T, g) - 1e-12)

    def test_split_types(self):
        self.assertEqual(split_type((3, 3, 3, 3)), "4")
        self.assertEqual(split_type((3, 3, 4, 3)), "3-1")
        self.assertEqual(split_type((3, 4, 3, 4)), "2-2 diagonal")
        self.assertEqual(split_type((4, 3, 3, 4)), "2-2 adjacent")
        self.assertEqual(split_type((3, 5, 3, 4)), "2-1-1 diagonal")
        self.assertEqual(split_type((3, 3, 5, 4)), "2-1-1 adjacent")
        self.assertEqual(split_type((1, 2, 3, 4)), "1-1-1-1")


class Thresholds(unittest.TestCase):
    def test_half_spiral_threshold(self):
        th = threshold(0.5)
        self.assertLess(th["adjacent"], 2.5)
        self.assertAlmostEqual(th["three_one"], 6.779, places=2)
        self.assertAlmostEqual(th["overall"], th["three_one"])

    def test_conditions_hold_beyond_threshold(self):
        for a in (0.25, 0.5, 1.0, 2.0):
            th = threshold(a)
            k1 = th["kappa1"]
            W = 1.0
            for L in np.geomspace(th["adjacent"] * 1.001, 1e7, 300):
                self.assertTrue(cond_small_kappa(L, a, W, k1))
                self.assertTrue(cond_large_kappa(L, a, W, k1))
                self.assertLess(p_bound(L, a, W), 1)
            for L in np.geomspace(th["three_one"] * 1.001, 1e7, 300):
                self.assertGreater((pi - W) * (L / (L + pi)) ** (2 * a + 1), W)

    def test_adjacent_contradiction_directly(self):
        """The best geometric lower bound for Q beats the spiral's upper bound at every floor and pass gap."""
        for a in (0.5, 1.0):
            L0 = threshold(a)["adjacent"]
            for L in np.geomspace(L0 * 1.001, 1e6, 120):
                for n in (1, 2, 3, 5, 10, 50, 1000):
                    lo, hi = exact_q_gap(L, n, a)
                    self.assertGreater(lo, hi)

    def test_211_limits(self):
        self.assertAlmostEqual(ratio_limit_211(0.5), ((pi - 1) / 1) ** 0.25)
        self.assertAlmostEqual(side_limit_211(0.5), 0.1342, places=3)


class Searches(unittest.TestCase):
    def test_plain_spiral_has_no_two_pass_squares(self):
        out = labelled_search(lambda T: radius(T, 0.0, 1.0), 0.5, n_seeds=40, seed=3)
        self.assertEqual({k: len(v["solutions"]) for k, v in out.items()}, {k: 0 for k in out})

    def test_wobbling_spiral_has_aligned_squares(self):
        out = labelled_search(lambda T: radius(T, 0.5, 3.0), 0.5, n_seeds=60, seed=3,
                              patterns={"2-2 adjacent": (0, 0, 1, 1)})
        self.assertGreater(len(out["2-2 adjacent"]["solutions"]), 0)

    def test_random_search_patterns_and_bounds(self):
        plain = classify_found(search(0.0, 1.0, n_random=800, seed=5), 0.5, 0.0, 1.0)
        self.assertTrue(set(plain["counts"]) <= {"1-1-1-1", "2-1-1 diagonal"})
        self.assertLessEqual(plain["W_max"], 1.0)
        self.assertGreater(plain["floor_min"], threshold(0.5)["overall"])
        wob = classify_found(search(0.5, 3.0, n_random=800, seed=5), 0.5, 0.5, 3.0)
        self.assertGreater(wob["aligned"], 0)
        self.assertEqual(wob["aligned_inside_bounds"], wob["aligned"])


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v383" / "two_turn_check.json").read_text())
        plain = d["random_search"]["0.0,1.0"]
        self.assertTrue(set(plain["counts"]) <= {"1-1-1-1", "2-1-1 diagonal"})
        for fam in d["random_search"].values():
            self.assertEqual(fam["aligned_inside_bounds"], fam["aligned"])
            self.assertLess(fam["flag_residual_max"], 1e-13)
        self.assertTrue(all(v["solutions"] == 0 for v in d["labelled_search"]["0.0,1.0"].values()))
        self.assertGreater(d["labelled_search"]["0.5,3.0"]["2-2 adjacent"]["solutions"], 0)
        self.assertEqual(set(d["v382_patterns"]), {"1-1-1-1", "2-1-1 diagonal"})
        for x in d["adjacent_margin"]:
            self.assertGreater(x["q_lower"], x["q_upper"])


if __name__ == "__main__":
    unittest.main()
