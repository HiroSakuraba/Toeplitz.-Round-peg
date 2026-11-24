from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from asymmetric_spiral_v374 import asymmetric_drawing  # noqa: E402
from surgery_v378 import (  # noqa: E402
    _edges, graceful, lattice_count, lattice_sector, near_squares, predicted_count, surgery_trial)

ROOT = Path(__file__).resolve().parents[1]
DEN = 2 ** 30


def rat(x):
    return Q(round(x * DEN), DEN)


class Surgery(unittest.TestCase):
    def bump_over(self, V, s, k, picks, extra=0.6):
        P0, P1 = _edges(V)
        a, w = P0[k], P1[k] - P0[k]
        L = np.hypot(*w)
        ts = [((q['A'] - a) @ w) / (w @ w) for q in picks]
        feet = [a + t * w for t in ts]
        nrm = (picks[0]['A'] - feet[0]) / np.hypot(*(picks[0]['A'] - feet[0]))
        hgt = max(np.hypot(*(q['A'] - f)) for q, f in zip(picks, feet)) * (1 + extra) + 1e-4
        t1, t2 = min(ts) - 0.002 / L, max(ts) + 0.002 / L
        foot = a + (t1 + t2) / 2 * w
        apex = (rat(foot[0] + hgt * nrm[0]), rat(foot[1] + hgt * nrm[1]))
        return surgery_trial(V, s, k, rat(t1), rat(t2), apex)

    def candidates(self, V, s, k):
        P0, P1 = _edges(V)
        a, w = P0[k], P1[k] - P0[k]
        out = []
        for q in near_squares(V, float(s), tuple((P0[k] + P1[k]) / 2), 0.024):
            t = ((q['A'] - a) @ w) / (w @ w)
            if 0.1 < t < 0.9 and 1e-4 < np.hypot(*(q['A'] - (a + t * w))) < 0.015:
                out.append((q, graceful(q, k + t, V)))
        return out

    def test_single_graceful_and_single_crossed(self):
        s = Q(1, 20)
        V = asymmetric_drawing(4, Q(1))
        # edge 27: a graceful near-square close to the curve; edge 37: a crossed one
        (q27, sec27), = [c for c in self.candidates(V, s, 27) if c[1] != 'crossed'][:1]
        r = self.bump_over(V, s, 27, [q27])
        self.assertEqual((r['swept_graceful'], r['actual_flip']), (1, 1))
        crossed = [c for c in self.candidates(V, s, 37) if c[1] == 'crossed']
        r = self.bump_over(V, s, 37, [crossed[0][0]])
        self.assertEqual((r['swept'], r['swept_graceful'], r['actual_flip']), (1, 0, 0))

    def test_recorded_trials(self):
        data = json.loads((ROOT / 'results/v378/surgery_and_lattice.json').read_text())
        trials = [t for t in data['surgery_trials'] if 'actual_flip' in t]
        self.assertEqual(len(trials), 42)
        self.assertTrue(all(t['actual_flip'] == t['predicted_flip'] for t in trials))
        kinds = {t['kind'] for t in trials}
        self.assertEqual(kinds, {'targeted', 'random', 'crossed single', 'pair'})
        self.assertTrue(any(t['swept'] == 3 for t in trials))


class Lattice(unittest.TestCase):
    def test_count_and_sector_rule(self):
        for lam in (0, 1):
            f = lattice_count(0.01, 16, lam)
            self.assertLess(abs(len(f) - predicted_count(16)), 4)
            for arms, z in f:
                self.assertEqual(lattice_sector(arms, z, 0.01, lam) == 'graceful', arms[0] != arms[1])

    def test_agrees_with_v377_solver(self):
        f = lattice_count(0.01, 2.5, 0)
        self.assertEqual(sorted(round(np.hypot(*z) / np.pi, 3) for _, z in f), [0.38, 0.38])

    def test_recorded_lattice(self):
        data = json.loads((ROOT / 'results/v378/surgery_and_lattice.json').read_text())
        for row in data['lattice']:
            self.assertTrue(row['graceful_iff_BD_on_different_arms'])
            if row['M'] >= 16:
                self.assertLess(abs(row['count'] / row['predicted'] - 1), 0.03)


if __name__ == '__main__':
    unittest.main()
