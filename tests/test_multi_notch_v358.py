from fractions import Fraction as Q
from pathlib import Path
import json
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from multi_notch_v358 import (  # noqa: E402
    NotchedRectangle, certify, closed_form_components, contact_squares, covering_theorem,
    inradius_candidates, largest_axis_square, radius_bracket, reaches_disk_target, redraw_tolerance,
    slide_parity, traced_component,
)
from convex_l_disk_v344_v345 import OrthogonalL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def random_shape(rng, k, size=16, den=1):
    while True:
        H = Q(rng.randint(3 * den, size * den), den)
        W = Q(rng.randint(int(H * den) + 1, 3 * size * den), den)
        pts = sorted({Q(rng.randint(1, int(W * den) - 1), den) for _ in range(2 * k)})
        if len(pts) == 2 * k:
            cuts = [(pts[2 * i], pts[2 * i + 1], Q(rng.randint(1, int(H * den) - 1), den)) for i in range(k)]
            return NotchedRectangle.make(W, H, cuts)


class ExactInradius(unittest.TestCase):
    def test_candidates_contain_the_inradius(self):
        rng = random.Random(3)
        for _ in range(150):
            shape = random_shape(rng, rng.randint(1, 5), den=rng.choice((1, 2, 3)))
            lo, hi = radius_bracket(shape, 70)
            if lo == hi:
                continue
            inside = [c for c in inradius_candidates(shape) if c.compare(lo) >= 0 and c.compare(hi) <= 0]
            self.assertTrue(inside, shape)

    def test_exact_ties(self):
        # integer shapes whose best odd component has squared side exactly 2 r^2
        for args, side_sq in (((10, 9, [(1, 2, 4), (5, 9, 8)]), Q(32)),
                              ((10, 9, [(1, 3, 3), (4, 7, 7)]), Q(49, 2)),
                              ((10, 9, [(1, 3, 3), (4, 8, 7)]), Q(49, 2))):
            shape = NotchedRectangle.make(*args)
            self.assertIs(reaches_disk_target(shape, side_sq), True)
            self.assertIs(reaches_disk_target(shape, side_sq - Q(1, 10**12)), False)

    def test_one_notch_agrees_with_the_comparison_L(self):
        # v351: one notch has the inradius of the comparison L built on its wider tower
        rng = random.Random(5)
        for _ in range(120):
            shape = random_shape(rng, 1)
            (a, b, k), = shape.cuts
            wide = max(a, shape.W - b)
            L = OrthogonalL(shape.W, shape.H, wide, k)
            for side_sq in (Q(rng.randint(1, 4 * int(shape.H) ** 2), 4), shape.H ** 2 / 2, 2 * (wide / 2) ** 2):
                self.assertEqual(reaches_disk_target(shape, side_sq), L.reaches_disk_target(side_sq))


class ClosedForms(unittest.TestCase):
    def test_every_closed_form_matches_the_independent_trace(self):
        rng = random.Random(7)
        seen = set()
        for n in range(60):
            shape = random_shape(rng, 1 + n % 5, size=14)
            if shape.H >= shape.W:
                continue
            polygon = shape.polygon()
            for comp in closed_form_components(shape):
                floor, deg, _ = traced_component(polygon, comp.seed)
                self.assertEqual(deg % 2 == 1, comp.odd, (shape, comp))
                if comp.odd:
                    self.assertEqual(floor, comp.side_sq, (shape, comp))
                seen.add(comp.kind.replace('-mirror', ''))
        self.assertGreaterEqual(seen, {'wall-pair', 'wall-floor', 'floor-wall', 'roof-slide', 'floor-slide',
                                       'anchored-floor-wall', 'anchored-floor-roof', 'anchored-wall-wall'})

    def test_slide_parity_table(self):
        # uncut slide: A leaves on the left, B on the right
        self.assertTrue(slide_parity(False, 'down', True, False, False, True))
        self.assertFalse(slide_parity(False, 'down', True, True, False, False))
        self.assertFalse(slide_parity(False, 'down', False, False, True, True))
        self.assertTrue(slide_parity(False, 'down', True, True, True, True))
        # a single square where A arrives from the west and B leaves to the east
        self.assertFalse(slide_parity(True, 'down', False, True, True, False))
        self.assertTrue(slide_parity(False, 'up', True, False, False, True))

    def test_slides_with_ends_at_notch_corners(self):
        # 28 x 17: the slide x in [0, 1] runs from the box corner to the corner (18, 17) of a notch
        shape = NotchedRectangle.make(28, 17, [(11, 13, 16), (18, 24, 5)])
        slides = [c for c in closed_form_components(shape) if c.kind == 'roof-slide' and c.detail[1] == 0]
        self.assertEqual([(c.detail, c.odd) for c in slides], [((17, 0, 1), True)])
        self.assertEqual(traced_component(shape.polygon(), slides[0].seed)[:2], (289, 1))


class Certificates(unittest.TestCase):
    def check(self, shape):
        cert = certify(shape)
        self.assertIsNotNone(cert, shape)
        if cert.level == 'theorem':
            self.assertEqual(cert.kind, covering_theorem(shape))
            return cert
        floor, deg, _ = traced_component(shape.polygon(), cert.seed)
        self.assertEqual(deg % 2, 1)
        self.assertGreaterEqual(floor, cert.side_sq)
        self.assertTrue(reaches_disk_target(shape, cert.side_sq))
        return cert

    def test_random_rational_shapes_with_up_to_five_notches(self):
        rng = random.Random(11)
        levels = set()
        for n in range(80):
            levels.add(self.check(random_shape(rng, 1 + n % 5, den=rng.choice((1, 2, 3)))).level)
        self.assertGreaterEqual(levels, {'theorem', 'closed-form'})

    def test_shapes_that_need_a_traced_component(self):
        # exact ties at notch corners, and slides that turn at a box corner along a wall at distance H
        for args in ((10, 9, [(1, 2, 4), (5, 9, 8)]), (8, 4, [(1, 2, 3), (5, 6, 1)]),
                     (42, 9, [(4, 5, 7), (9, 36, 3)]), (31, 18, [(14, 18, 4)]),
                     (15, 5, [(1, 2, 4), (3, 4, 3), (5, 11, 1)])):
            self.assertEqual(self.check(NotchedRectangle.make(*args)).level, 'traced')

    def test_transverse_search_agrees_with_closed_forms(self):
        rng = random.Random(13)
        for n in range(10):
            shape = random_shape(rng, 1 + n % 3, size=12)
            if shape.H >= shape.W:
                continue
            strict = {sq for _, sq in contact_squares(shape.polygon())}
            for comp in closed_form_components(shape):
                if 'slide' not in comp.kind:
                    self.assertIn(min(comp.seed[j:] + comp.seed[:j] for j in range(4)), strict)


class FiveNotchExample(unittest.TestCase):
    shape = NotchedRectangle.make(46, 10, [(7, 12, 8), (13, 29, 5), (30, 31, 3), (33, 38, 2), (40, 43, 9)])

    def test_inradius_and_largest_square(self):
        lo, hi = radius_bracket(self.shape, 80)
        r = [c for c in inradius_candidates(self.shape) if c.compare(lo) >= 0 and c.compare(hi) <= 0]
        self.assertEqual({(c.p / c.m, c.s / c.m, c.q) for c in r}, {(15, Q(-1, 2), 448)})  # 15 - 4 sqrt 7
        self.assertEqual(largest_axis_square(self.shape), 8)
        self.assertIsNone(covering_theorem(self.shape))

    def test_redraw_tolerance(self):
        wall_pair = ((Q(29), Q(0)), (Q(38), Q(0)), (Q(38), Q(9)), (Q(29), Q(9)))
        tilted = ((Q(8), Q(0)), (Q(10), Q(8)), (Q(2), Q(10)), (Q(0), Q(2)))
        self.assertEqual(redraw_tolerance(self.shape.polygon(), wall_pair, Q(1, 32)), (81, Q(13, 16), 1, 3))
        self.assertEqual(redraw_tolerance(self.shape.polygon(), tilted, Q(1, 32)), (68, Q(13, 16), 2, 3))
        for square in (wall_pair, tilted):
            self.assertEqual(traced_component(self.shape.polygon(), square)[:2],
                             ((square[1][0] - square[0][0]) ** 2 + (square[1][1] - square[0][1]) ** 2, 1))
        # 9 - 13/16 = 131/16 and sqrt(68) - 13/16 > 37/5 both clear sqrt(2) r
        self.assertTrue(reaches_disk_target(self.shape, Q(131, 16) ** 2))
        self.assertGreater(Q(68), (Q(37, 5) + Q(13, 16)) ** 2)
        self.assertTrue(reaches_disk_target(self.shape, Q(37, 5) ** 2))


class SweepRecords(unittest.TestCase):
    def test_recorded_sweeps_have_no_misses(self):
        files = sorted((ROOT / 'results' / 'v358').glob('*.json'))
        self.assertTrue(files)
        for f in files:
            data = json.loads(f.read_text())
            self.assertEqual(data['misses'], [], f.name)
            self.assertEqual(sum(data['levels'].values()), data['shapes'])


if __name__ == '__main__':
    unittest.main()
