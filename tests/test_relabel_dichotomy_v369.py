from fractions import Fraction as Q
from itertools import product
from math import atan2
from pathlib import Path
import json
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import PositiveDimensionalSquares, census_polygon, pt  # noqa: E402
from relabel_dichotomy_v369 import (  # noqa: E402
    L, T, big_gap_label, gaps, is_row, residual, same_row, separated_pair, shift, small_square_scale_sq,
)
from sweep_search_v366 import RegisteredPolygon, symmetric_circle  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DECAGON = [tuple(map(Q, v)) for v in json.loads((ROOT / 'research/v247-refined-ten-gon.json').read_text())['vertices']]


class Equivariance(unittest.TestCase):
    def test_residual_of_the_relabeled_row_is_linear(self):
        rng = random.Random(1)
        for sigma in (1, -1):
            for _ in range(300):
                P = tuple((Q(rng.randint(-60, 60), 7), Q(rng.randint(-60, 60), 7)) for _ in range(4))
                r = residual(P, sigma)
                self.assertEqual(residual(shift(P), sigma), L(r, sigma))
                x = r
                for _ in range(4):
                    x = L(x, sigma)
                self.assertEqual(x, r)

    def test_relabeling_is_free_and_shifts_the_big_gap(self):
        rng = random.Random(2)
        for _ in range(500):
            t = sorted(Q(rng.randint(0, 96), 96) for _ in range(4))
            row = tuple(t)
            if not is_row(row):
                continue
            y = row
            for k in range(1, 4):
                y = T(y)
                self.assertFalse(same_row(row, y))
            self.assertTrue(same_row(row, T(T(T(T(row))))))
            lab = big_gap_label(row)
            if lab is not None:
                self.assertEqual(big_gap_label(T(row)), (lab - 1) % 4)

    def test_balanced_rows_have_a_pair_one_third_apart(self):
        n = 24
        for g in product(range(n + 1), repeat=3):
            g4 = n - sum(g)
            if g4 < 0 or max(g + (g4,)) * 2 > n:
                continue
            row = (Q(0), Q(g[0], n), Q(g[0] + g[1], n), Q(g[0] + g[1] + g[2], n))
            _, d = separated_pair(row)
            self.assertGreaterEqual(d, Q(1, 3))
        # sharpness: three equally spaced positions with one repeated
        self.assertEqual(separated_pair((Q(0), Q(1, 3), Q(1, 3), Q(2, 3)))[1], Q(1, 3))


class Examples(unittest.TestCase):
    def test_symmetric_octagon_has_an_invariant_family(self):
        poly = RegisteredPolygon(symmetric_circle(2))
        for k in range(16):
            t = Q(k, 61)
            row = (t, t + Q(1, 4), t + Q(1, 2), t + Q(3, 4))
            pts = tuple(poly.point(x) for x in row)
            self.assertEqual(residual(pts), ((0, 0), (0, 0)))
            self.assertTrue(same_row(T(row), tuple(x + Q(1, 4) for x in row)))
            self.assertEqual(gaps(row), (Q(1, 4),) * 4)

    def test_decagon_has_one_isolated_square_and_no_small_squares(self):
        census = census_polygon([pt(*v) for v in DECAGON])
        self.assertEqual((census.count, census.graceful_count), (1, 1))
        tau2 = small_square_scale_sq(RegisteredPolygon(DECAGON))
        self.assertTrue(all(2 * s.side_squared >= tau2 for s in census.squares))

    def test_small_square_scale_on_random_star_polygons(self):
        rng = random.Random(5)
        checked = 0
        while checked < 6:
            n = rng.randint(5, 6)
            pts = set()
            while len(pts) < n:
                pts.add((rng.randint(-9, 9), rng.randint(-9, 9)))
            pts = sorted(pts, key=lambda p: atan2(p[1] - Q(1, 3), p[0] - Q(1, 7)))
            try:
                poly = RegisteredPolygon(pts)
                census = census_polygon([pt(*p) for p in pts])
            except (ValueError, PositiveDimensionalSquares):
                continue
            tau2 = small_square_scale_sq(poly)
            self.assertTrue(all(2 * s.side_squared >= tau2 for s in census.squares))
            checked += 1


if __name__ == '__main__':
    unittest.main()
