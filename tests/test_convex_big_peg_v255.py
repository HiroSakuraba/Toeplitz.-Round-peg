"""v255: convex sharp Big Peg, hidden-corner lemma, Pochuev's hidden corner.

Exact tests use only the standard library. Numeric tests need numpy, scipy
(and cma for none of them); they are skipped when those are missing unless
REQUIRE_NUMERIC_DEPS=1, in which case a missing dependency is a failure.
"""
import os
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from polygon_orbit_census_v223 import pt, census_polygon  # noqa: E402
import convex_big_peg_v255 as v  # noqa: E402

try:
    import numpy  # noqa: F401
    import scipy  # noqa: F401
    HAVE_NUMERIC = True
except ImportError:
    HAVE_NUMERIC = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"


class HiddenCornerDecision(unittest.TestCase):
    def test_square_inscribed_in_disk_has_no_hidden_corner(self):
        for u in (pt(1, 0), pt(Q(3, 5), Q(4, 5)), pt(Q(5, 13), Q(12, 13))):
            sq = v.square_from(pt(0, 0), u)
            self.assertEqual(v.hidden_corners(sq, Q(1)), [])

    def test_small_far_square_has_hidden_corner(self):
        sq = v.square_from(pt(3, 2), pt(Q(1, 4), Q(1, 8)))
        self.assertTrue(v.hidden_corners(sq, Q(1)))

    def test_decision_matches_independent_triangle_certificate(self):
        cases = [(pt(Q(8, 5), Q(9, 5)), pt(Q(-4, 5), 0)), (pt(3, 2), pt(Q(1, 4), Q(1, 8))),
                 (pt(Q(-7, 4), Q(1, 2)), pt(Q(1, 3), Q(1, 2)))]
        for p, u in cases:
            sq = v.square_from(p, u)
            hid = v.hidden_corners(sq, Q(1))
            self.assertTrue(hid)
            for k in hid:
                cert = v.triangle_certificate(sq, k, grid=20)
                self.assertIsNotNone(cert, (p, u, k))
                self.assertTrue(all(Q(x) > 0 for x in cert["barycentric"]))

    def test_big_square_around_disk_has_no_hidden_corner(self):
        self.assertEqual(v.hidden_corners(v.square_from(pt(0, 0), pt(2, 0)), Q(1)), [])

    def test_lemma_is_one_directional(self):
        # a BIG square (side 3*sqrt2) can still hide a corner: the converse is false
        sq = v.square_from(pt(5, 0), pt(3, 0))
        self.assertEqual(v.hidden_corners(sq, Q(1)), [2])


class LemmaCensus(unittest.TestCase):
    def test_exact_grid(self):
        r = v.lemma_grid_census()
        self.assertGreater(r["small_squares"], 10000)
        self.assertEqual(r["small_with_hidden_corner"], r["small_squares"])
        self.assertGreater(r["equality_squares"], 1000)
        self.assertEqual(r["equality_free_off_centre"], 0)
        self.assertEqual(r["equality_without_hidden_corner"], 6)   # six orientations at p = 0

    def test_float_sweep(self):
        r = v.lemma_float_sweep(n=20000, seed=7)
        self.assertEqual(r["without_hidden_corner"], 0)


class ExactConvexCensus(unittest.TestCase):
    def test_every_square_of_convex_polygons_is_big(self):
        out = v.exact_convex_census()
        self.assertGreaterEqual(len(out), 10)
        for name, rec in out.items():
            self.assertTrue(rec["all_side_squared_at_least_2"], name)
        self.assertLess(out["octagon_near_equality"]["min_side_squared_float"], 2.001)

    def test_nonconvex_control_has_small_square(self):
        # Pochuev's polygon contains the unit disk but has a square with side^2 < 2
        P, _, _ = v.read_polygon(v.POCHUEV_JSON)
        self.assertFalse(v.is_strictly_convex_ccw(P))
        c = census_polygon(P)
        self.assertTrue(any(s.side_squared < 2 for s in c.squares))


class Pochuev(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = v.pochuev_certificate()

    def test_small_square_and_case(self):
        self.assertTrue(self.r["side_squared_below_2"])
        self.assertTrue(self.r["both_diagonal_lines_miss_unit_disk"])

    def test_exactly_one_hidden_corner_next_to_reflex_v2(self):
        self.assertEqual(len(self.r["hidden_corner_indices"]), 1)
        self.assertEqual(self.r["hidden_corner_on_edges"], [2])
        self.assertEqual(self.r["reflex_vertices"], [2, 9])
        self.assertLess(self.r["hidden_corner_distance_to_v2_float"], 1e-4)
        b = [Q(x) for x in self.r["barycentric_certificate"]["barycentric"]]
        self.assertTrue(all(x > 0 for x in b) and sum(b) == 1)

    def test_kernel_misses_disk(self):
        self.assertFalse(self.r["star_shaped_about_disk_centre"])
        self.assertTrue(self.r["kernel_misses_unit_disk"])
        self.assertGreater(Q(self.r["kernel_distance_squared"]), Q(14, 10))


@unittest.skipUnless(HAVE_NUMERIC or REQUIRE, "numpy/scipy not installed")
class Numeric(unittest.TestCase):
    def test_float_solver_agrees_with_exact_census(self):
        import numpy as np
        for name in ("lopsided_pentagon", "octagon_near_equality", "triangle"):
            P = v.convex_test_polygons()[name]
            exact = sorted(float(s.side_squared) for s in census_polygon(P).squares)
            polys = v._square_polys(np, [[float(p.x), float(p.y)] for p in P])
            got = sorted(2 * h * h for _, h in polys)
            self.assertEqual(len(got), len(exact), name)
            for a, b in zip(got, exact):
                self.assertAlmostEqual(a, b, places=7)

    def test_random_convex_small_run(self):
        r = v.random_convex_check(count=300, seed=3)
        self.assertGreaterEqual(r["min_half_diagonal_over_inradius"], 1 - 1e-9)

    def test_twin_spike_breaks_strong_form_only(self):
        import numpy as np
        for gap in (0.3, 1.6):
            P = v.twin_spike(gap)
            self.assertGreaterEqual(v._edge_line_dist(np, P).min(), 1 - 1e-12)   # unit disk in kernel
            self.assertGreaterEqual(v._seg_dist(np, P).min(), 1 - 1e-12)         # unit disk inside
        r = v.twin_spike_check()
        for gap, rec in r.items():
            self.assertLess(rec["min_half_diagonal"], 1, gap)
            self.assertGreater(rec["max_half_diagonal"], 1, gap)
            self.assertGreater(rec["small_squares_far_from_disk"], 0, gap)


if __name__ == "__main__":
    unittest.main()
