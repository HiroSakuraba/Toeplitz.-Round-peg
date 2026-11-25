from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from general_comb_v380 import corners, search  # noqa: E402
from label_degree_v390 import jacobian_sign, partner_label, phase_map  # noqa: E402
from two_turn_squares_v383 import window  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Labels(unittest.TestCase):
    def setUp(self):
        self.sq = search(0.0, 1.0, n_random=500, seed=7)
        self.assertTrue(self.sq)

    def test_zeros_and_local_degree(self):
        for v, Ths in self.sq:
            Z = corners(v)
            th, k, _ = window(Z, Ths)
            ks = list(map(int, k))
            w = np.array([Z[0].real, Z[0].imag, np.angle(Z[1] - Z[0]), abs(Z[1] - Z[0])])
            self.assertLess(np.max(np.abs(phase_map(w, ks, 0.0, lambda r: 0.0))), 1e-8)
            self.assertNotEqual(jacobian_sign(w, ks, 0.0, lambda r: 0.0)[0], 0)

    def test_partners_have_other_labels(self):
        for v, Ths in self.sq:
            Z = corners(v)
            th, k, _ = window(Z, Ths)
            ks = list(map(int, k))
            pl = partner_label(Z, ks)
            self.assertNotEqual(pl, ks)
            self.assertEqual([x % 2 for x in pl], [(x + 1) % 2 for x in ks])     # the half-turn swaps the arms


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v390" / "label_check.json").read_text())
        p = d["plain"]
        self.assertEqual(p["partner_same_label"], 0)
        self.assertEqual(p["nondegenerate"], p["squares"])
        u = d["uniqueness"]
        self.assertEqual(u["distinct_labels"], u["squares"])
        self.assertEqual(u["on_branch"], u["squares"])


if __name__ == "__main__":
    unittest.main()
