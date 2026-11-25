from math import pi
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from slow_twist_v388 import best_cutoff, c_L, r_graceful, rho_max  # noqa: E402
from smooth_cores_v386 import r_smooth  # noqa: E402
from symmetric_window_v385 import window_inner_radius  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Radii(unittest.TestCase):
    def test_graceful_radius(self):
        t = 1e-4
        self.assertAlmostEqual(r_graceful(t), r_smooth(t))
        self.assertGreater(r_graceful(t, Phi=0.4), r_smooth(t))
        self.assertGreater(r_graceful(t, L=100.0), r_graceful(t))
        self.assertAlmostEqual(c_L(0.5), 2 * pi ** 2 / (2 * 2 ** 0.5))

    def test_cutoff_and_window(self):
        for eps in (0.1, 0.01, 1e-4):
            s = best_cutoff(eps, Phi=0.4, L=0.5)
            self.assertLessEqual(8 * r_graceful(s, Phi=0.4, L=0.5), eps * (1 + 1e-9))
            self.assertGreater(8 * r_graceful(s * 1.001, Phi=0.4, L=0.5), eps)
            r_smooth_case, _ = rho_max(eps, Phi=0.4, L=0.5)
            r_poly, _ = window_inner_radius(eps)
            self.assertGreater(r_smooth_case, r_poly)
            self.assertGreaterEqual(rho_max(eps, Phi=0.4, L=0.01)[0], r_smooth_case)


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v388" / "slow_twist_check.json").read_text())
        for lam, v in d["graceful"].items():
            self.assertGreater(v["graceful"], 0)
            self.assertLessEqual(v["max_ratio_graceful"], 1.0)
            self.assertLessEqual(v["max_plain_ratio_graceful"], 1.0)
            self.assertGreater(v["max_plain_ratio_crossed"], 1.0)
        rows = [r for r in d["exponents"] if r["eps"] == 1e-05 and r["L"] == 0.0]
        self.assertLess(rows[0]["log_ratio"], 7.5)


if __name__ == "__main__":
    unittest.main()
