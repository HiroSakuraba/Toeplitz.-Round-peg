import json
import sys
import unittest
from math import log, pi, sin
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import pinned_squares_v391 as P  # noqa: E402
from two_turn_squares_v383 import window  # noqa: E402

DATA = ROOT / "results" / "v391" / "pinned_check.json"


class TestBand(unittest.TestCase):
    def test_upper_wall_for_small_and_large_labels(self):
        for lab in [(1, 3, 1), (2, 3, 1), (5, 6, 1), (-6, 7, -1), (-8, -6, 1), (3, 7, 3)]:
            b = P.pinned_band(*lab, 0.8, 4.0, grid=800)
            self.assertTrue(b["hi_wall"], lab)
            self.assertLess(b["hi"], 20 * b["centre"])

    def test_big_label_band_is_closed_and_contains_resonance(self):
        for Phi in (0.4, 1.2):
            b = P.pinned_band(12, -13, -1, Phi, 5 * Phi, grid=800)
            self.assertTrue(b["closed"])
            self.assertEqual(b["unknown_in_band"], 0)
            self.assertLess(b["lo"], b["centre"])
            self.assertGreater(b["hi"], b["centre"])

    def test_band_grows_with_twist(self):
        b1 = P.pinned_band(12, -13, -1, 0.4, 2.0, grid=800)
        b2 = P.pinned_band(12, -13, -1, 1.2, 6.0, grid=800)
        self.assertLessEqual(b2["lo"], b1["lo"])
        self.assertGreaterEqual(b2["hi"], b1["hi"])

    def test_even_level_is_never_pinned_from_above_by_parity(self):
        # j = 0 (even arm count) leaves no forced gap pi |j|: the numerator bound is not positive
        self.assertEqual(P.band_status(400.0, 8, 8, 0, 0.4, 2.0), "unknown")


class TestSideBound(unittest.TestCase):
    def test_supersolution(self):
        q, u = P.side_bound(300.0, 6, 7, 0.4, 2.0)
        self.assertLessEqual(P.side_map(q, 300.0, 6, 7, 0.4, 2.0), q)
        self.assertGreater(u, q / np.sqrt(2))

    def test_out_of_range_is_unknown(self):
        self.assertIsNone(P.side_bound(10.0, 6, 7, 0.4, 2.0))


class TestBranchSafeResidual(unittest.TestCase):
    def test_no_jump_across_negative_axis(self):
        # a small square whose centre sits just above / below the negative real axis
        ks = [100, 100, 100, 100]
        phi = lambda r: 0.0
        vals_w, vals_r = [], []
        t, ang = 1e-4, 0.3
        off = (1 + 1j) * t * np.exp(1j * ang) / 2   # centre minus corner A
        for y in (1e-7, -1e-7):
            A = complex(-0.05, y) - off
            v = np.array([A.real, A.imag, ang, t])
            Z = P.corners_of(v)
            th, _, _ = window(Z, [0, 0, 0, 0])
            vals_w.append(np.array([abs(z) ** -2 - t - pi * k for z, t, k in zip(Z, th, ks)]))
            vals_r.append(P.label_residual(v, ks, 0.0, phi, pi))
        self.assertGreater(np.max(np.abs(vals_w[0] - vals_w[1])), 6.0)   # v383 window wraps
        self.assertLess(np.max(np.abs(vals_r[0] - vals_r[1])), 1e-2)     # unwrapped residual does not


class TestResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads(DATA.read_text())

    def test_survival(self):
        for m, r in self.d["survival"].items():
            self.assertEqual(r["squares"], 156)
            self.assertGreaterEqual(r["reached"], 150, m)
            self.assertEqual(r["degree_kept"], r["reached"], m)
            self.assertEqual(r["below_hi"], r["squares"], m)       # no square approaches the point
            self.assertLessEqual(r["max_side_ratio_over_estimate"], 1.0, m)

    def test_closure(self):
        for Phi, r in self.d["closure"].items():
            self.assertEqual(r["hi_walls"], r["labels"], Phi)
            self.assertIsNotNone(r["N_all_close"], Phi)
        for Phi, rows in self.d["large_labels"].items():
            self.assertTrue(all(rows.values()), Phi)

    def test_one_continuation(self):
        row = next(r for r in self.d["survival"]["1"]["rows"] if r["reached"])
        self.assertIsNotNone(row["degree_end"])
        self.assertNotEqual(row["degree_end"], 0)


if __name__ == "__main__":
    unittest.main()
