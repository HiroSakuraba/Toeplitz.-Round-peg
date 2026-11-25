from math import log, pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from symmetric_window_v385 import (  # noqa: E402
    best_cutoff, d_g, eta_star, eta_star_leading, lattice_check, lemma_lhs, phase_defect, window_inner_radius)

ROOT = Path(__file__).resolve().parents[1]


class Lemma(unittest.TestCase):
    def test_eta_star_solves_the_inequality(self):
        for s in (0.05, 0.01, 0.002):
            e = eta_star(s)
            self.assertLess(lemma_lhs(e * (1 - 1e-9), s), pi / 2)
            self.assertGreaterEqual(lemma_lhs(e * (1 + 1e-6), s), pi / 2)

    def test_leading_order(self):
        for s in (0.02, 0.01, 0.005):
            self.assertLess(abs(eta_star(s) / eta_star_leading(s) - 1), 2e-3)
        self.assertLess(eta_star(0.01, Phi=0.4), eta_star(0.01))
        self.assertLess(eta_star(0.01, nu=0.1), eta_star(0.01))

    def test_bound_against_lattice(self):
        for lam in (0, 1):
            rows = lattice_check(np.geomspace(0.004, 0.025, 8), lam)
            for r in rows:
                self.assertGreaterEqual(r["ratio"], 1.0)

    def test_phase_defect(self):
        r, h = 0.01, 1e-9
        self.assertAlmostEqual(phase_defect(h, r) / (2 * h * r ** -3), 1.0, places=3)


class Window(unittest.TestCase):
    def test_window_radius(self):
        prev = None
        for eps in (0.3, 0.1, 0.03, 0.01, 1e-3, 1e-5):
            rho, s = window_inner_radius(eps)
            self.assertAlmostEqual(s, d_g(eps / 8))
            self.assertLess(rho, s / 2)
            ratio = log(rho) / log(eps)
            if prev is not None:
                self.assertLess(ratio, prev)        # tends down to the exponent ((1+a)/a)^2 = 9
            self.assertGreater(ratio, 9)
            prev = ratio
        self.assertAlmostEqual(best_cutoff(0.1), d_g(0.1 / 8))


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v385" / "window_check.json").read_text())
        for k in ("lattice_lam0", "lattice_lam1", "solver_lam0", "solver_lam1"):
            self.assertGreaterEqual(d[k]["min_ratio"], 1.0)
        self.assertLess(d["lattice_lam0"]["min_ratio"], 1.01)     # the bound is attained


if __name__ == "__main__":
    unittest.main()
