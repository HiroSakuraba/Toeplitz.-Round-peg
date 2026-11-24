from fractions import Fraction as Q
from pathlib import Path
import json
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from certificate_margin_v368 import (  # noqa: E402
    margin_bracket, margin_upper, registered_distance_sq, sqrt_hi, sqrt_lo,
)
from refinement_stress_v367 import bump_across  # noqa: E402
from sweep_search_v366 import RegisteredPolygon, certificate, locus, search, symmetric_circle  # noqa: E402
from three_patch_sweep_v364 import target_roof  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DECAGON = [tuple(map(Q, v)) for v in json.loads((ROOT / 'research/v247-refined-ten-gon.json').read_text())['vertices']]
bins = lambda *b, span=1: [(Q(x, 16), Q(span, 16)) for x in b]
UNION = [(Q(6, 16), Q(2, 16)), (Q(1, 16), Q(3, 16)), (Q(10, 16), Q(3, 16))]


class Brackets(unittest.TestCase):
    def test_tight_brackets(self):
        for poly, arcs, orient, low, high in (
                (RegisteredPolygon(target_roof()), bins(1, 7, 11), 1, Q(963, 100000), Q(964, 100000)),
                (RegisteredPolygon(DECAGON), bins(7, 3, 10), -1, Q(1077, 100000), Q(1078, 100000)),
                (RegisteredPolygon(DECAGON), UNION, -1, Q(2776, 100000), Q(2779, 100000))):
            lo, hi, _ = margin_bracket(poly, arcs, orient, 2000)
            self.assertLessEqual(lo, hi)
            self.assertLessEqual(hi, margin_upper(poly, arcs, orient))
            self.assertGreater(lo, low)
            self.assertLess(hi, high)

    def test_square_roots_and_registered_distance(self):
        self.assertLessEqual(sqrt_lo(2) ** 2, 2)
        self.assertGreaterEqual(sqrt_hi(2) ** 2, 2)
        self.assertEqual(sqrt_lo(Q(9, 4)), Q(3, 2))
        a = RegisteredPolygon([(0, 0), (4, 0), (4, 4), (0, 4)])
        b = RegisteredPolygon([(0, 0), (4, 0), (4, 4), (Q(-1, 2), 4)])
        self.assertEqual(registered_distance_sq(a, b), Q(1, 4))


class Tolerance(unittest.TestCase):
    def test_random_redraws_below_the_margin_keep_the_degree(self):
        poly = RegisteredPolygon(DECAGON)
        arcs = bins(7, 3, 10)
        lo, _, _ = margin_bracket(poly, arcs, -1, 1500)
        rng = random.Random(7)
        kept = 0
        for _ in range(12):
            m = lo * Q(99, 100) / 2
            V = [(x + Q(rng.randint(-100, 100), 100) * m, y + Q(rng.randint(-100, 100), 100) * m) for x, y in DECAGON]
            moved = RegisteredPolygon(V)
            self.assertLess(registered_distance_sq(poly, moved), lo * lo)
            c = certificate(moved, *arcs, -1)
            kept += c is not None and c['degree'] == 1
        self.assertEqual(kept, 12)

    def test_v367_killing_bump_exceeds_the_union_margin(self):
        data = json.loads((ROOT / 'results/v367/decagon_greedy_push.json').read_text())
        final = RegisteredPolygon([tuple(map(Q, v)) for v in data['final_vertices']], data['final_times'])
        lo, hi, _ = margin_bracket(final, UNION, -1, 1500)
        D = [e[3] for seg in locus(final, *bins(6, 1, 12), -1) for e in seg
             if 'A0' in e[4] and seg[0][:3] != seg[1][:3]][0]
        killed = bump_across(final, D)
        self.assertEqual(certificate(final, *UNION, -1)['degree'], 1)
        self.assertEqual(certificate(killed, *UNION, -1)['degree'], 0)
        self.assertGreater(registered_distance_sq(final, killed), (16 * hi) ** 2)


class CircleRatio(unittest.TestCase):
    def test_margin_is_a_small_fraction_of_the_bump(self):
        for b in (Q(1, 64), Q(1, 256), Q(1, 1024)):
            poly = RegisteredPolygon(symmetric_circle(2, bump=b))
            res = search(poly, bins=16, spans=(1, 2))
            best = max(margin_bracket(poly, bins(*c['arcs'], span=c['span']), c['orient'], 800)[1] for c in res)
            # the bumped vertex is at distance b from the unit circle, so every registration has delta >= b
            self.assertLess(best, b / 10)


if __name__ == '__main__':
    unittest.main()
