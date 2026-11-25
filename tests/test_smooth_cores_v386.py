from math import exp, pi
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from smooth_cores_v386 import (  # noqa: E402
    c2, confinement_terms, default_twist, effective_bound_holds, r_smooth, sector_two_arm, selection_data,
    twisted_search, twisted_selection)
from tilted_squares_v382 import arm_point, three_corners  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Constants(unittest.TestCase):
    def test_c2(self):
        self.assertAlmostEqual(c2(0.5), 0.5 * (2 * pi / 2) ** 0.5)
        t = 0.01
        self.assertAlmostEqual(c2() * r_smooth(t) ** 2, t)


class PlainSpiral(unittest.TestCase):
    def test_identity_and_confinement_on_v382_squares(self):
        d = json.loads((ROOT / "results" / "v382" / "tilted_check.json").read_text())
        n = 0
        for s in d["squares"][::25]:
            sol = three_corners(s["Theta_A"], s["n1"], s["n3"])
            if sol is None:
                continue
            th, t = sol
            A = arm_point(s["Theta_A"])
            e = np.exp(1j * th)
            Z = [A, A + t * e, A + (1 + 1j) * t * e, A + 1j * t * e]
            r = selection_data(Z, [abs(z) ** -2 for z in Z])
            self.assertLess(abs(r["flag"]), 1e-14)
            self.assertGreaterEqual(r["t"] / (c2() * r["rmin"] ** 2), 1.0)
            if r["t"] / r["absZ"] <= 0.2:
                self.assertTrue(r["inside"])
                self.assertTrue(confinement_terms(r)[3])
                n += 1
        self.assertGreater(n, 10)


class Twisted(unittest.TestCase):
    def test_crossed_square_at_the_comb_scale(self):
        r0 = exp(-4 * pi / 5)                       # a zero of the twist 0.4 sin(5 log 1/r)
        found = twisted_search([r0 * exp(-0.0405)], [(-1, -1)])
        self.assertEqual(len(found), 1)
        Z, ks = found[0]
        d = twisted_selection(Z, ks, default_twist)
        self.assertTrue(d["inside"])
        self.assertEqual(d["j"], 0)
        self.assertEqual(d["arm_two"], "diagonal")
        self.assertLess(d["ratio"], 0.25)           # far below the plain spiral's smooth bound (ratio >= 1)
        self.assertLess(d["t"] / (pi * 0.5 * d["rmin"] ** 3), 1.6)      # side about sqrt2 turn spacings
        self.assertEqual(sector_two_arm(ks, Z), "crossed")
        self.assertTrue(effective_bound_holds(d, Z))


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v386" / "smooth_check.json").read_text())
        v = d["v382"]
        self.assertEqual(v["identity_inside"], v["small_ratio"])
        self.assertEqual(v["alternative_holds"], v["small_ratio"])
        self.assertGreaterEqual(v["min_t_over_c2r2"], 1.0)
        self.assertLess(v["min_t_over_c2r2"], 1.1)
        p = d["v380_plain"]
        self.assertEqual(p["alternative_holds"], p["small_ratio"])
        tw = d["twisted"]
        self.assertEqual(tw["identity_inside"], tw["squares"])
        self.assertEqual(tw["bound_holds"], tw["squares"])
        self.assertEqual(tw["sectors"], {"crossed": tw["squares"]})


if __name__ == "__main__":
    unittest.main()
