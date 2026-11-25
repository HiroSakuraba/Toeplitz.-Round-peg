import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import pinned_squares_v391 as P  # noqa: E402
from mpmath import iv  # noqa: E402

DATA = ROOT / "results" / "v391" / "pinned_check.json"


class TestSampledBand(unittest.TestCase):
    def test_walls_are_excluded_samples_outside_the_allowed_run(self):
        b = P.pinned_band(12, -13, -1, 0.4, 2.0, grid=800)
        self.assertLess(b["lo_wall_lift"], b["allowed_lo"])
        self.assertGreater(b["hi_wall_lift"], b["allowed_hi"])
        for T in (b["lo_wall_lift"], b["hi_wall_lift"]):
            self.assertEqual(P.band_status(T, 12, -13, -1, 0.4, 2.0), "excl")
        for T in (b["allowed_lo"], b["allowed_hi"], b["centre"]):
            self.assertNotEqual(P.band_status(T, 12, -13, -1, 0.4, 2.0), "excl")

    def test_even_level_forces_nothing(self):
        self.assertEqual(P.band_status(400.0, 8, 8, 0, 0.4, 2.0), "unknown")

    def test_sign_compatibility(self):
        self.assertTrue(P.compatible(5, 6, 1))
        self.assertTrue(P.compatible(-6, 7, -1))
        self.assertFalse(P.compatible(5, 6, -1))


class TestCertificates(unittest.TestCase):
    def test_big_label_region_and_tail(self):
        for Phi in (0.4, 1.2):
            b = P.certified_band(12, -13, -1, Phi, 5 * Phi, grid=800)
            self.assertTrue(b["certified_closed"], Phi)
            reg = b["region"]
            for k in ("lo_face", "hi_face", "outer_side_face", "inner_side_face"):
                self.assertTrue(reg[k], (Phi, k))
            self.assertLess(reg["q_in"], reg["q_hat"])
            self.assertLess(reg["q_hat"], P.Q_MAX)
            self.assertIsNotNone(b["tail_from"])
            self.assertTrue(P.certified_tail(b["tail_from"], 12, -13, -1, Phi, 5 * Phi))

    def test_resonance_is_not_excluded(self):
        b = P.pinned_band(12, -13, -1, 0.4, 2.0, grid=800)
        self.assertFalse(P.certified_excluded(b["centre"], 12, -13, -1, 0.4, 2.0))

    def test_small_label_has_a_certified_tail(self):
        for lab in [(1, 3, 1), (2, 3, 1), (5, 6, 1), (-8, -6, 1), (3, 7, 3)]:
            b = P.certified_band(*lab, 0.8, 4.0, grid=800)
            self.assertIsNotNone(b["tail_from"], lab)

    def test_interval_side_bound_dominates_float(self):
        q_f = P.side_bound(300.0, 6, 7, 0.4, 2.0)
        q_c = P.certified_side_bound(iv.mpf(300.0), 6, 7, 0.4, 2.0)
        self.assertIsNotNone(q_c)
        self.assertGreaterEqual(q_c, q_f)
        g = P.iv_side_map(iv.mpf(q_c), iv.mpf(300.0), 6, 7, 0.4, 2.0)
        self.assertLessEqual(float(g.b), q_c)

    def test_failed_certificate_is_not_exclusion(self):
        self.assertIsNone(P.certified_side_bound(iv.mpf(10.0), 6, 7, 0.4, 2.0))
        self.assertFalse(P.certified_excluded(10.0, 6, 7, 1, 0.4, 2.0))


class TestContinuation(unittest.TestCase):
    """Regenerates real continuations instead of reading stored results."""

    def test_two_squares_survive(self):
        starts = P.survival_starts()
        self.assertEqual(len(starts), 156)
        for (x0, label, _), m in ((starts[0], 1), (starts[len(starts) // 2], 3)):
            self.assertTrue(P.odd_arm_count(label))
            row = P.survival_row(x0, label, m)
            self.assertTrue(row["reached"])
            self.assertNotEqual(row["degree_start"], 0)
            self.assertEqual(row["degree_end"], row["degree_start"])
            self.assertLess(row["residual_end"], 1e-10)


class TestResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads(DATA.read_text())

    def test_survival(self):
        for m, r in self.d["survival"].items():
            self.assertEqual(r["squares"], 156, m)
            self.assertEqual(r["reached"], r["squares"], m)
            self.assertEqual(r["degree_kept"], r["squares"], m)
            self.assertEqual(r["certified_tails"], r["squares"], m)
            self.assertEqual(r["below_tail"], r["squares"], m)
            self.assertEqual(r["inside_certified_region"], r["certified_closed"], m)
            self.assertLessEqual(r["max_side_over_bound"], 1.0, m)
            self.assertLess(r["max_residual_end"], 1e-10, m)
            for row in r["rows"]:
                self.assertNotEqual(row["degree_start"], 0)
                self.assertEqual(row["degree_end"], row["degree_start"])
                self.assertTrue(np.isfinite([row["lift_min"], row["lift_max"], row["q_min"], row["q_max"]]).all())

    def test_closure(self):
        for Phi, r in self.d["closure"].items():
            self.assertEqual(r["certified_tails"], r["labels"], Phi)
            self.assertIsNotNone(r["N_all_certified"], Phi)
        for Phi, rows in self.d["large_labels"].items():
            self.assertTrue(all(rows.values()), Phi)


if __name__ == "__main__":
    unittest.main()
