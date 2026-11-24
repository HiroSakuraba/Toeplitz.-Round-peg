from fractions import Fraction as Q
from pathlib import Path
import json
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import pt, validate_jordan_polygon  # noqa: E402
from twist_stability_v371 import large_graceful_parity, twisted  # noqa: E402
from wild_point_parity_v370 import drawing, sector_counts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Twisted(unittest.TestCase):
    def test_twisted_drawings_at_depth_four(self):
        expected = {0.0: (29, 29, 1), 0.01: (29, 23, 1), 0.1: (28, 33, 0)}
        for tw, (al, rev, parity) in expected.items():
            V, _ = twisted(4, tw)
            validate_jordan_polygon([pt(*p) for p in V])
            par, counts = large_graceful_parity(V, Q(1, 20))
            self.assertEqual((counts['aligned'], counts['reversed'], par), (al, rev, parity))

    def test_recorded_scan(self):
        rows = json.loads((ROOT / 'results/v371/twist_scan.json').read_text())['rows']
        by = {(r['twist'], r['depth']): r for r in rows}
        for r in rows:
            self.assertEqual(r['degenerate'], 0)
        for d in range(2, 8):
            sym = by[(0.001, d)]
            self.assertEqual(sym['1/20']['aligned'], sym['1/20']['reversed'])
            for s in ('1/20', '1/40'):
                c = by[(0.01, d)][s]
                self.assertEqual((c['aligned'] + c['reversed']) % 2, 0)
        # twist 0.1 at cutoff 1/20: the same aligned square is missing at depths 3 to 6
        for d in range(3, 7):
            self.assertEqual(by[(0.1, d)]['1/20']['aligned'], by[(0.001, d)]['1/20']['aligned'] - 1)


class Stability(unittest.TestCase):
    def test_small_vertex_perturbations_keep_the_large_graceful_count_odd(self):
        V0, _ = drawing(0.5, 3.0, 2, 12)
        rng = random.Random(3)
        for _ in range(4):
            V = [(x + Q(rng.randint(-50, 50), 10**6), y + Q(rng.randint(-50, 50), 10**6)) for x, y in V0]
            validate_jordan_polygon([pt(*p) for p in V])
            total, degenerate, _ = sector_counts(V, Q(100))
            self.assertEqual(degenerate, [])
            self.assertEqual((total['aligned'] % 2, total['reversed'] % 2), (1, 0))
            par, _ = large_graceful_parity(V, Q(1, 10))
            self.assertEqual(par, 1)


if __name__ == '__main__':
    unittest.main()
