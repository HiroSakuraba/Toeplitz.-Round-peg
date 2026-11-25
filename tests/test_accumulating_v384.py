from fractions import Fraction as Q
from math import pi, sqrt
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from accumulating_v384 import (  # noqa: E402
    STRIP, cluster_split, default_twist, drawing, large_parity, separated, sparse_range_check, strip_flattenable)

ROOT = Path(__file__).resolve().parents[1]


class Geometry(unittest.TestCase):
    def test_three_corners_never_fit_in_the_strip(self):
        rng = np.random.default_rng(0)
        s = 1.0
        eta = 0.999 * STRIP * s
        for _ in range(20000):
            A = complex(rng.uniform(-2, 2), rng.uniform(-eta, eta))
            e = np.exp(1j * rng.uniform(0, 2 * pi))
            Z = [A, A + s * e, A + (1 + 1j) * s * e, A + 1j * s * e]
            self.assertLessEqual(sum(abs(z.imag) <= eta for z in Z), 2)

    def test_sparse_hypotheses_below_threshold(self):
        self.assertEqual(sparse_range_check(np.geomspace(1e-7, 0.0303, 1500)), [])
        self.assertTrue(sparse_range_check([0.031, 0.04, 0.05]))

    def test_splits(self):
        self.assertEqual(strip_flattenable("sparse", list(range(2, 8)), Q(1, 60)), ([4, 5, 6, 7], [2, 3]))
        self.assertEqual(strip_flattenable("sparse", list(range(2, 8)), Q(1, 480)), ([6, 7], [2, 3, 4, 5]))
        big, cluster = cluster_split("self-similar", range(2, 9), Q(1, 64))
        self.assertEqual((big, cluster), ([2, 3], [4, 5, 6, 7, 8]))
        self.assertTrue(separated("self-similar", big, Q(1, 64)))
        self.assertEqual(strip_flattenable("self-similar", cluster, Q(1, 64)), ([6, 7, 8], [4, 5]))


class Counter(unittest.TestCase):
    def test_tiny_spirals_are_counted(self):
        """Regression for the absolute determinant cutoff: a spiral 16 times smaller has the same squares."""
        out = []
        for k in (5, 7):
            p, c, _ = large_parity(drawing("sparse", k, k, 2, 12, default_twist), Q(1, 60))
            out.append((p, c))
        self.assertEqual(out[0], out[1])
        self.assertEqual(out[0][0], 1)          # Lemma 384.2: a spiral inside the strip cannot change L

    def test_scale_invariance(self):
        a = large_parity(drawing("sparse", 3, 3, 2, 12, default_twist), Q(1, 60))
        b = large_parity(drawing("sparse", 4, 4, 2, 12, default_twist), Q(1, 240))
        self.assertEqual(a[:2], b[:2])


class Results(unittest.TestCase):
    def setUp(self):
        self.d = json.loads((ROOT / "results" / "v384" / "trials.json").read_text())

    def test_sparse_additivity(self):
        sparse = [r for r in self.d if r["model"] == "sparse"]
        self.assertEqual(len(sparse), 8)
        for r in sparse:
            self.assertTrue(r["separated"])
            self.assertEqual(r["full"]["parity"], r["additive_prediction"])
            self.assertEqual(r["full"]["parity"], r["kept_only"]["parity"])
            for k, v in r["singles"].items():
                self.assertEqual(r["full"]["domains"].get(k), v["domains"].get(k))
            self.assertNotIn("straddling", r["full"]["domains"])
            if not r["twisted"]:
                self.assertEqual(r["full"]["parity"], 1)
                for v in r["full"]["domains"].values():
                    self.assertEqual(v["aligned"], v["reversed"])

    def test_sparse_recursion(self):
        """L(s/4) = L(s) + defect of the largest spiral at s/4 (Corollary 384.6), on the twisted drawings."""
        tw = {r["s"]: r for r in self.d if r["model"] == "sparse" and r["twisted"]}
        for a, b in (("1/60", "1/240"), ("1/120", "1/480")):
            defect = 1 + tw[b]["singles"]["2"]["parity"]
            self.assertEqual(tw[b]["full"]["parity"], (tw[a]["full"]["parity"] + defect) % 2)

    def test_self_similar_records(self):
        ss = [r for r in self.d if r["model"] == "self-similar"]
        self.assertTrue(ss)
        for r in ss:
            if r.get("kind") == "interaction":
                self.assertEqual(r["joint"]["parity"], r["additive_prediction"])
            else:
                self.assertEqual(r["full"]["parity"], r["additive_prediction"])


if __name__ == "__main__":
    unittest.main()
