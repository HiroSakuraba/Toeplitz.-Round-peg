import json
import sys
import unittest
from math import pi
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import contact_twist_v394 as C  # noqa: E402

DATA = ROOT / "results" / "v394" / "contact_check.json"


class TestModel(unittest.TestCase):
    def test_arms_meet_exactly_at_twist_pi(self):
        for r in (0.3, 0.05, 0.002):
            p1, p2 = C.arm_points(r, pi)
            self.assertLess(abs(p1 - p2), 1e-9 * r)
            p1, p2 = C.arm_points(r, pi - 0.1)
            self.assertGreater(abs(p1 - p2), 0.09 * r)
            p1, p2 = C.arm_points(r, -pi)
            self.assertLess(abs(p1 - p2), 1e-9 * r)


class TestCertificates(unittest.TestCase):
    def test_plain_tail(self):
        s_tail, ok = C.plain_tail()
        self.assertTrue(ok)
        self.assertAlmostEqual(s_tail, 1.5, places=6)

    def test_entry_near_contact(self):
        c = C.best_c(pi - 0.3, 5)
        self.assertIsNotNone(c)
        cert = C.certify_entry(pi - 0.3, 5, c)
        cert["c"] = c
        self.assertTrue(C.entry_ok(cert))
        self.assertGreaterEqual(cert["s_lo"], C.S_LO_MIN)
        self.assertFalse(C.certify_entry(pi - 0.3, 5, c - 1)["certified"])

    def test_annulus(self):
        cert = C.certify_entry(pi - 0.5, 5, 60)
        self.assertTrue(cert["certified"])
        r_in, r_out = C.annulus(cert, 60, 61, 1)
        r2_in, r2_out = C.annulus(cert, 120, 122, 1)
        self.assertLess(r_in, r_out)
        self.assertAlmostEqual(r_in / r2_in, 2.0, delta=0.05)     # radii scale like 1 / m


class TestResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads(DATA.read_text())

    def test_degree_for_all_regions(self):
        self.assertTrue(self.d["degree_for_all_regions"])
        for j, c in self.d["coarse"].items():
            self.assertTrue(c["tail_certified"])
            self.assertEqual(len(c["patterns"]), 4)

    def test_table(self):
        rows = self.d["table"]
        self.assertTrue(all(r["c"] and r["degree_applies"] for r in rows))
        by = {(r["delta"], r["L"]): r["c"] for r in rows}
        for (d, L), c in by.items():
            for (d2, L2), c2 in by.items():
                if d2 <= d and L2 >= L:
                    self.assertGreaterEqual(c2, c - 1, ((d, L), (d2, L2)))   # harder bounds never need fewer

    def test_examples(self):
        for ex in self.d["examples"]:
            self.assertTrue(ex["reached"] and ex["inside"], ex)
            self.assertNotEqual(ex["degree_start"], 0)
            self.assertEqual(ex["degree_end"], ex["degree_start"])


if __name__ == "__main__":
    unittest.main()
