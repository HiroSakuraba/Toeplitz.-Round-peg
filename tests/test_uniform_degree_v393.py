import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import uniform_degree_v393 as G  # noqa: E402
from uniform_degree_v393 import I  # noqa: E402

DATA = ROOT / "results" / "v393" / "degree_check.json"
CERTS = json.loads((ROOT / "results" / "v392" / "uniform_check.json").read_text())["levels"]


class TestScaledMap(unittest.TestCase):
    def test_vanishes_at_a_plain_square(self):
        from uniform_closure_v392 import plain_square
        for lab in [(12, 13, 1), (-16, 23, -1)]:
            x, _ = plain_square(*lab)
            (s, U, V), p = G.scaled_coords(*lab, x)
            P = G.psi_box([I(s), I(U), I(V)], (I(p["eps"]), I(p["r"]), p["sig1"], p["sig3"], p["b_is_min"], lab[2]))
            self.assertLess(max(abs(float(v.mid)) for v in P), 1e-10)
            self.assertLess(max(abs(s - 1), abs(U - 1), abs(V - 1)), 0.1)

    def test_limit_system(self):
        rng = np.random.default_rng(3)
        for _ in range(20):
            s, U, V = rng.uniform(0.5, 1.5, 3)
            for sig1, sig3, bm in G.patterns(1):
                P = G.psi_box([I(s), I(U), I(V)], (I(0), I(0.7), sig1, sig3, bm, 1))
                want = (s * U - 1, s * V - 1, sig1 * sig3 * s * U * V - 1)
                for p, w in zip(P, want):
                    self.assertAlmostEqual(float(p.mid), w, places=12)

    def test_jacobian_matches_finite_differences(self):
        par = (I(0.05), I(0.8), 1, 1, True, 1)
        z = np.array([1.02, 0.97, 0.99])
        J = np.array([[float(e.mid) for e in row] for row in G.jac_box([I(t) for t in z], par)])
        h = 1e-6
        for k in range(3):
            dz = np.zeros(3)
            dz[k] = h
            f1 = np.array([float(v.mid) for v in G.psi_box([I(t) for t in z + dz], par)])
            f0 = np.array([float(v.mid) for v in G.psi_box([I(t) for t in z - dz], par)])
            np.testing.assert_allclose(J[:, k], (f1 - f0) / (2 * h), atol=1e-6)


class TestCertificates(unittest.TestCase):
    def test_krawczyk_and_local_localization_on_a_piece(self):
        star = [I(1 - G.B_STAR, 1 + G.B_STAR)] * 3
        leaves = G.certify_piece(I(0, 1 / 88), I(0.5, 0.5625), 1, 1, True, 1, star)
        self.assertIsNotNone(leaves)
        self.assertEqual({lf["det_sign"] for lf in leaves}, {-1})

    def test_krawczyk_fails_on_a_box_without_the_zero(self):
        par = (I(0, 1 / 88), I(0.5, 0.5625), 1, 1, True, 1)
        self.assertFalse(G.krawczyk([I(1.3, 1.4), I(0.9, 1.1), I(0.9, 1.1)], par)[0])

    def test_star_inside_v392_region(self):
        for P, c in CERTS.items():
            self.assertTrue(G.star_inside_region(c["1"], 2.0, 1, c["1"]["c"])["ok"], P)


class TestResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads(DATA.read_text())

    def test_certified(self):
        self.assertTrue(all(self.d["certified"].values()))
        self.assertEqual(sorted(self.d["levels"]), ["-1", "1"])
        for j, lv in self.d["levels"].items():
            self.assertEqual(len(lv["fine"]), 4)
            for k, v in lv["fine"].items():
                self.assertTrue(v["certified"], (j, k))
                self.assertIn(v["det_sign"], (-1, 1))
                self.assertTrue(all(lf["det_sign"] == v["det_sign"] for lf in v["leaves"]))
            for P, co in lv["coarse"].items():
                self.assertTrue(all(c["certified"] for c in co.values()), (j, P))
        for P, st in self.d["star_inside"].items():
            self.assertTrue(all(x["ok"] for x in st.values()), P)

    def test_validation(self):
        by_pattern = {}
        for row in self.d["validation"]:
            self.assertTrue(row["in_certified_box"], row)
            self.assertNotEqual(row["local_degree"], 0)
            by_pattern.setdefault(row["pattern"], set()).add(row["local_degree"])
        self.assertTrue(all(len(v) == 1 for v in by_pattern.values()))


if __name__ == "__main__":
    unittest.main()
