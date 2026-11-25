import json
import sys
import unittest
from pathlib import Path

import numpy as np
from mpmath import iv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import uniform_closure_v392 as U  # noqa: E402

DATA = ROOT / "results" / "v392" / "uniform_check.json"


class TestCertificates(unittest.TestCase):
    def test_all_levels_constant_is_certified_and_minimal_for_the_method(self):
        self.assertTrue(U.certify_levels(0.4, 2.0, 11, 1)["certified"])
        self.assertFalse(U.certify_levels(0.4, 2.0, 10, 1)["certified"])

    def test_uniform_family_and_single(self):
        self.assertTrue(U.certify_family(0.4, 2.0, 1, 1 / 11, 1 / 11, uniform=True)["certified"])
        self.assertTrue(U.certify_family(0.4, 2.0, 1, 1 / 10, 1 / 11)["certified"])          # m = 10, every M >= 11
        self.assertTrue(U.certify_family(0.4, 2.0, 1, 1 / 10, 1 / 10, single=True)["certified"])

    def test_small_label_is_not_certified(self):
        self.assertIsNone(U.label_ok(0.4, 2.0, 1, 9, 9))

    def test_even_level_is_refused(self):
        self.assertFalse(U.certify_family(0.4, 2.0, 2, 1 / 50, 1 / 50, uniform=True)["certified"])

    def test_region_contents(self):
        c = U.certify_levels(1.2, 6.0, 16, 1)
        self.assertTrue(c["certified"])
        self.assertLess(c["s_lo"], 1.0)
        self.assertGreater(c["s_hi"], 1.0)
        self.assertLess(c["Q_hat"] / 16, U.Q_MAX)
        self.assertGreater(c["side_floor"], 0)


class TestMonotonicity(unittest.TestCase):
    """The uniform arguments rely on these monotonicities; spot-check them on interval enclosures."""

    def test_side_map_increases_with_eps_and_mu_and_decreases_with_s(self):
        q = iv.mpf(0.05)
        base = U.side_G(q, iv.mpf(0.05), iv.mpf(0.03), iv.mpf(1.0), 0.4, 2.0, 1).b
        self.assertGreater(U.side_G(q, iv.mpf(0.06), iv.mpf(0.03), iv.mpf(1.0), 0.4, 2.0, 1).b, base)
        self.assertGreater(U.side_G(q, iv.mpf(0.05), iv.mpf(0.04), iv.mpf(1.0), 0.4, 2.0, 1).b, base)
        self.assertLess(U.side_G(q, iv.mpf(0.05), iv.mpf(0.03), iv.mpf(1.2), 0.4, 2.0, 1).b, base)

    def test_scaled_map_increases_with_eps(self):
        Q = 0.6
        vals = [float(U.side_G(iv.mpf(e * Q), iv.mpf(e), iv.mpf(e), iv.mpf(0.8), 0.4, 2.0, 1).b) / e
                for e in (0.01, 0.03, 0.06, 0.09)]
        self.assertTrue(all(x < y for x, y in zip(vals, vals[1:])))

    def test_floor_decreases_with_eps_and_s(self):
        f = [U.side_floor_scaled(e, e, 1.5, 0.6, 0.4, 2.0, 1) for e in (0.02, 0.05, 0.08)]
        self.assertTrue(all(x > y for x, y in zip(f, f[1:])))
        self.assertGreater(U.side_floor_scaled(0.05, 0.05, 1.2, 0.6, 0.4, 2.0, 1),
                           U.side_floor_scaled(0.05, 0.05, 1.6, 0.6, 0.4, 2.0, 1))

    def test_levels_form_matches_per_level_form(self):
        j, m, M, s = 3, 40, 70, 0.9
        a = U.identity_check(iv.mpf(1 / m), iv.mpf(1 / M), s, 0.02, 0.4, 2.0, j, (1, -1))
        b = U.identity_levels(iv.mpf(j / m), iv.mpf(j / M), iv.mpf(1 / j), s, 0.02, 0.4, 2.0, (1, -1))
        pj = np.pi * j
        self.assertAlmostEqual(float(a["R_lo"]) / pj, float(b["R_lo"]), places=9)
        self.assertAlmostEqual(float(a["N_hi"]) / pj, float(b["N_hi"]), places=9)


class TestConcreteLabels(unittest.TestCase):
    def test_v391_agrees_at_the_faces(self):
        cert = U.certify_levels(0.8, 4.0, 13, 1)
        rows = U.cross_check(cert, [(13, 13, 1), (20, 57, 1), (40, 41, 3)])
        self.assertTrue(all(r["v391_agrees"] for r in rows))

    def test_big_label_square_stays_inside(self):
        cert = U.certify_levels(1.2, 6.0, 16, 1)
        rows = U.big_label_survival({"1.2": cert}, [(-20, -35, 1)], ms=(3,))
        r = rows[0]
        self.assertTrue(r["found"] and r["covered"] and r["reached"] and r["inside"])
        self.assertNotEqual(r["degree_start"], 0)
        self.assertEqual(r["degree_end"], r["degree_start"])


class TestResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads(DATA.read_text())

    def test_levels(self):
        want = {"0.4": 11, "0.8": 13, "1.2": 16}
        for P, c in want.items():
            self.assertEqual(self.d["levels"][P]["1"]["c"], c)
            for j0, cert in self.d["levels"][P].items():
                self.assertTrue(cert["certified"])

    def test_per_level_thresholds(self):
        want = {"0.4": 10, "0.8": 12, "1.2": 14}
        for P, t in want.items():
            self.assertEqual(self.d["per_level"][P]["1"]["threshold"], t)

    def test_cross_check_and_survival(self):
        for P, rows in self.d["cross_check"].items():
            self.assertTrue(all(r["v391_agrees"] for r in rows), P)
        for r in self.d["survival"]:
            self.assertTrue(r["found"])
            self.assertTrue(r["reached"] and r["inside"], r)
            self.assertNotEqual(r["degree_start"], 0)
            self.assertEqual(r["degree_end"], r["degree_start"])


if __name__ == "__main__":
    unittest.main()
