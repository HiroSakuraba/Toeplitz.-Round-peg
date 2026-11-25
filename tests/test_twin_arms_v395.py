import json
import sys
import unittest
from math import pi
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import twin_arms_v395 as T  # noqa: E402
from uniform_closure_v392 import certify_levels, plain_square  # noqa: E402

DATA = ROOT / "results" / "v395" / "twin_check.json"


class TestTwinLemma(unittest.TestCase):
    def test_numerical_illustrations(self):
        t = T.twin_checks()
        self.assertIsNotNone(t["opposite_winding_meets_at"])
        self.assertIsNotNone(t["winding_vs_straight_meets_at"])
        self.assertIsNotNone(t["unequal_rates_meet_at"])
        self.assertIsNone(t["twins_meet_at"])
        lo, hi = t["twins_relative_twist_range"]
        self.assertGreater(lo, -pi)
        self.assertLess(hi, pi)


class TestDeformedMap(unittest.TestCase):
    def test_reduces_to_v390_without_deformation(self):
        from label_degree_v390 import anchored_map
        x, lab = plain_square(20, 21, 1)
        zero = lambda r: 0.0
        a = T.label_map2(x, lab, 0.0, T.example_e1, T.example_phi)
        b = anchored_map(x, lab, 0.0, zero)
        self.assertLess(max(abs(a - b)), 1e-9)

    def test_constant_common_rotation_changes_nothing(self):
        x, lab = plain_square(20, 21, 1)
        zero = lambda r: 0.0
        a = T.label_map2(x, lab, 1.0, lambda r: 0.0, zero)
        b = T.label_map2(x, lab, 1.0, lambda r: 123.4, zero)
        self.assertLess(max(abs(a - b)), 1e-9)

    def test_zeros_are_inscribed_squares(self):
        x0, lab = plain_square(60, 64, 1)
        path = T.continue2(x0, lab, T.example_e1, T.example_phi, steps=200)
        self.assertEqual(path[-1][0], 1.0)
        self.assertLess(T.on_curve_error(path[-1][1], lab, 1.0, T.example_e1, T.example_phi), 1e-8)


class TestCertificates(unittest.TestCase):
    def test_common_slope_costs_labels(self):
        self.assertTrue(certify_levels(pi - 0.3, 5, 32, 1)["certified"])
        self.assertTrue(certify_levels(pi - 0.3, (5, 0), 32, 1)["certified"])
        self.assertFalse(certify_levels(pi - 0.3, (5, 1), 32, 1)["certified"])

    def test_example_bounds(self):
        Phi, L, L1 = T.example_bounds(0.005)
        self.assertEqual((Phi, L), (2.0, 8.0))
        self.assertLess(L1, 4.0)


class TestResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads(DATA.read_text())

    def test_table(self):
        rows = self.d["table"]
        self.assertTrue(all(r["c"] and r["degree_applies"] for r in rows))
        by = {(r["delta"], r["L"], r["L1"]): r["c"] for r in rows}
        for (d, L, L1), c in by.items():
            for (d2, L2, L12), c2 in by.items():
                if d2 <= d and L2 >= L and L12 >= L1:
                    self.assertGreaterEqual(c2, c - 1)

    def test_examples(self):
        for e in self.d["examples"]:
            self.assertTrue(e["covered"], e)
            self.assertTrue(e["reached"] and e["inside"], e)
            self.assertLess(e["on_curve_error"], 1e-8)
            self.assertNotEqual(e["degree_start"], 0)
            self.assertEqual(e["degree_end"], e["degree_start"])


if __name__ == "__main__":
    unittest.main()
