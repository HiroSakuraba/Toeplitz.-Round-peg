import json
import sys
import unittest
from pathlib import Path

import numpy as np
from mpmath import mp, mpf, quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

import any_exponent_v396 as X  # noqa: E402
import uniform_degree_v393 as G  # noqa: E402
from uniform_degree_v393 import D, I  # noqa: E402

DATA = ROOT / "results" / "v396" / "exponent_check.json"


class TestRemainders(unittest.TestCase):
    def test_K_G_H_enclose_exact_values(self):
        mp.dps = 30
        for p in (0.6, 1.0, 2.0):
            for rho in (-0.05, 0.01, 0.08):
                K = quad(lambda t: 2 * (1 - t) * (1 + t * rho) ** (-p - 2), [0, 1])
                Gx = (K - 1) / rho
                H = (Gx + (p + 2) / mpf(3)) / rho
                r = D.var(I(rho), 0)
                for fn, exact in ((G.d_K, K), (G.d_G, Gx), (G.d_H, H)):
                    v = fn(r, I(p)).v
                    self.assertTrue(v.a <= exact <= v.b, (fn.__name__, p, rho))
                # (1 + rho)^(-p) = 1 - p rho + p(p+1)/2 rho^2 K exactly
                self.assertAlmostEqual(float((1 + rho) ** (-p)), float(1 - p * rho + p * (p + 1) / 2 * rho ** 2 * K),
                                       places=12)

    def test_derivative_enclosures(self):
        mp.dps = 30
        p, rho, h = 1.3, 0.04, mpf("1e-12")
        Kf = lambda x: quad(lambda t: 2 * (1 - t) * (1 + t * x) ** (-p - 2), [0, 1])
        dK = (Kf(rho + h) - Kf(rho - h)) / (2 * h)
        g = G.d_K(D.var(I(rho), 0), I(p)).g[0]
        self.assertTrue(g.a <= dK <= g.b)


class TestGeneralMap(unittest.TestCase):
    def test_encloses_the_closed_form_at_one_half(self):
        rng = np.random.default_rng(1)
        for _ in range(15):
            s, U, V = rng.uniform(0.6, 1.4, 3)
            e, r = rng.uniform(0, 0.09), rng.uniform(0.5, 1)
            for pat in G.patterns(1) + G.patterns(-1):
                j = 1 if pat[0] * pat[1] > 0 else -1
                A = G.psi_box([I(s), I(U), I(V)], (I(e), I(r), *pat, j))
                B = G.psi_box([I(s), I(U), I(V)], (I(e), I(r), *pat, j, I(0.5)))
                for x, y in zip(A, B):
                    self.assertTrue(y.a <= x.mid <= y.b)

    def test_limit_system_for_any_exponent(self):
        for a in (0.3, 1.0, 3.0):
            for s, U, V in ((1.1, 0.9, 1.05), (0.8, 1.3, 0.95)):
                for pat in G.patterns(1):
                    P = G.psi_box([I(s), I(U), I(V)], (I(0), I(0.7), *pat, 1, I(a)))
                    want = (s * U - 1, s * V - 1, pat[0] * pat[1] * s * U * V - 1)
                    for x, w in zip(P, want):
                        self.assertAlmostEqual(float(x.mid), w, places=10)

    def test_jacobian_matches_finite_differences(self):
        par = (I(0.02), I(0.8), 1, 1, True, 1, I(1.6))
        z = np.array([1.01, 0.98, 0.99])
        J = np.array([[float(e.mid) for e in row] for row in G.jac_box([I(t) for t in z], par)])
        h = 1e-6
        for k in range(3):
            dz = np.zeros(3)
            dz[k] = h
            f1 = np.array([float(v.mid) for v in G.psi_box([I(t) for t in z + dz], par)])
            f0 = np.array([float(v.mid) for v in G.psi_box([I(t) for t in z - dz], par)])
            np.testing.assert_allclose(J[:, k], (f1 - f0) / (2 * h), atol=1e-4)

    def test_vanishes_at_plain_squares_for_other_exponents(self):
        from uniform_closure_v392 import plain_square, resonance
        for a, (n1, n3, j) in ((1.5, (30, 31, 1)), (0.3, (40, 70, 1)), (3.0, (45, -60, -1))):
            x, _ = plain_square(n1, n3, j, a)
            self.assertIsNotNone(x)
            m, M = min(abs(n1), abs(n3)), max(abs(n1), abs(n3))
            eps, r = 1 / m, m / M
            s1, s3 = (1 if n1 > 0 else -1), (1 if n3 > 0 else -1)
            b_is_min = abs(n1) == m
            t1, t3 = (r, 1.0) if b_is_min else (1.0, r)
            c = a * abs(j) / (2 * a + 1)
            sc = (x[0] / resonance(n1, n3, j, a), x[1] / (-eps * s1 * c * t1), x[2] / (eps * s3 * c * t3))
            P = G.psi_box([I(v) for v in sc], (I(eps), I(r), s1, s3, b_is_min, j, I(a)))
            for v in P:
                self.assertTrue(v.a - 1e-9 <= 0 <= v.b + 1e-9)


class TestCertificates(unittest.TestCase):
    def test_small_fine_piece(self):
        star = [I(1 - G.B_STAR, 1 + G.B_STAR)] * 3
        lv = X.certify_piece(I(0, 1 / 120), I(0.9, 1.0), I(1.2, 1.3), 1, 1, True, 1, star)
        self.assertIsNotNone(lv)
        self.assertEqual({x["det_sign"] for x in lv}, {-1})

    def test_plain_tail_for_an_exponent_range(self):
        self.assertIsNotNone(X.plain_tail_a(I(0.7, 1.2)))

    def test_tight_fallback_for_one_plus_rho(self):
        # here the form ec (2u + ec |w|^2) gives an enclosure of 1 + rho reaching below 0 (true minimum about 0.18);
        # the fallback |1 + w|^2 must be used, and must still enclose point values
        u, v, ec, p = I(-39.6, -39.5), I(17.5, 17.6), I(0.006, 0.015), I(1 / 6, 1 / 4)
        self.assertLessEqual(float((1 + ec * (2 * u + ec * (u * u + v * v))).a), 0)
        n2, T0, T1 = G.t_terms(D(u), D(v), D(ec), p)
        for e in (0.006, 0.01, 0.015):
            for pp in (1 / 6, 0.2, 1 / 4):
                _, T0p, T1p = G.t_terms(D(I(-39.55)), D(I(17.55)), D(I(e)), I(pp))
                self.assertTrue(T1p.v.a >= T1.v.a and T1p.v.b <= T1.v.b)

    def test_box_too_large_for_the_remainder_is_split(self):
        # on this box even |1 + w|^2 has an enclosure reaching 0 (1 + ec u straddles 0); localization must split it,
        # not fail
        aI = I(0.4, 0.7)
        par = (I(0, 1 / 30), I(0.5, 1), 1, 1, True, 1, aI)
        region = [I(0.25, 1.3), I(-400, 400), I(-400, 400)]
        with self.assertRaises(ValueError):
            G.psi_box(region, par)
        ok, n = G.localize(region, [I(0.85, 1.15)] * 3, par, max_boxes=1)
        self.assertFalse(ok)
        self.assertEqual(n, 2)


class TestResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = json.loads(DATA.read_text())

    def test_degree_certified_for_every_piece(self):
        self.assertEqual(len(self.d["degree_certified"]), 5)
        self.assertTrue(all(self.d["degree_certified"].values()))

    def test_table(self):
        self.assertTrue(all(r["c"] and r["degree_applies"] for r in self.d["table"]))

    def test_validation(self):
        for v in self.d["validation"]:
            self.assertTrue(v["found"], v)
            self.assertTrue(v["in_certified_box"], v)
            self.assertEqual(v["local_degree"], v["certified_sign"], v)


if __name__ == "__main__":
    unittest.main()
