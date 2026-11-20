"""v266: two-corner lemma, stability, automatic convex cutoff, sectors.

Exact tests use only the standard library. The twin-spike sector test needs
numpy; it is skipped without it unless REQUIRE_NUMERIC_DEPS=1.
"""
import os
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from polygon_orbit_census_v223 import pt, census_polygon  # noqa: E402
import convex_big_peg_v255 as v255  # noqa: E402
import convex_lemmas_v266 as v  # noqa: E402

try:
    import numpy  # noqa: F401
    HAVE_NUMPY = True
except ImportError:
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"


class TwoCornerLemma(unittest.TestCase):
    def test_grid_every_small_square_hides_a_near_corner(self):
        r = v.two_corner_grid_census()
        self.assertGreater(r["small_squares"], 15000)
        self.assertEqual(r["with_hidden_near_corner"], r["small_squares"])

    def test_near_corner_choice(self):
        sq = v255.square_from(pt(3, 2), pt(Q(1, 2), 0))       # p.u > 0, p.Ju > 0
        self.assertEqual(v.near_corner_indices(sq), (2, 3))
        sq = v255.square_from(pt(-3, -2), pt(Q(1, 2), 0))
        self.assertEqual(v.near_corner_indices(sq), (0, 1))

    def test_equality_square_has_no_hidden_near_corner(self):
        for d in v._pythagorean_dirs():
            self.assertEqual(v.near_corners_hidden(v255.square_from(pt(0, 0), d)), [])

    def test_corner_inside_disk_counts_as_hidden(self):
        sq = v255.square_from(pt(Q(1, 2), 0), pt(Q(1, 4), 0))  # whole square inside the disk
        self.assertTrue(v.near_corners_hidden(sq))

    def test_pochuev_hidden_corner_is_a_near_corner(self):
        r = v.pochuev_near_corner()
        self.assertTrue(r["hidden"])
        self.assertTrue(r["hidden_corner_is_near"])


class Stability(unittest.TestCase):
    def test_exact_constants(self):
        K = v.stability_constants(Q(1, 10))
        self.assertEqual(K["eps"], Q(1, 100))
        self.assertEqual(K["margin"], Q(79, 101))
        self.assertEqual(K["T"], Q(101, 3950))
        self.assertLess(K["T"] ** 2, 7 * K["eps"] ** 2)        # used in the convex form
        with self.assertRaises(ValueError):
            v.stability_constants(Q(1, 2))                       # eps = 1/4 > 3 - 2 sqrt 2

    def test_grid_bounds_and_sharpness(self):
        r = v.stability_grid_census(Q(1, 10))
        self.assertGreater(r["unhidden_near_corner_squares"], 1000)
        self.assertEqual(r["satisfy_X_plus_Y_le_T"], r["unhidden_near_corner_squares"])
        self.assertEqual(r["satisfy_convex_centre_bound"], r["fully_unhidden_squares"])
        self.assertEqual(r["max_centre_sq_over_2eps2"], 1.0)

    def test_sharpness_example(self):
        r = v.sharpness_example(Q(1, 100))
        self.assertEqual(r["hidden_corners"], [])
        self.assertTrue(r["corners_outside_open_disk"])
        self.assertEqual(Q(r["centre_distance_squared_over_eps_squared"]), 2)


class MovingDiskCutoff(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fam = v.v264_families()

    def test_bernstein_positivity(self):
        self.assertTrue(v._positive_on([Q(1), Q(-2), Q(2)], Q(0), Q(1), 10))   # 2t^2-2t+1 > 0
        self.assertFalse(v._positive_on([Q(1, 5), Q(-2), Q(2)], Q(0), Q(1), 10))  # roots in (0,1)
        self.assertFalse(v._positive_on([Q(-1), Q(4)], Q(0), Q(1), 10))       # negative at 0

    def test_window_certificate_and_floor(self):
        A, B = self.fam["v254_collision_window"]
        d = v.auto_disk_cutoff(A, B)
        self.assertGreater(d["side_squared_cutoff"], Q(9, 10))
        s = d["slabs"][0]
        c0 = pt(Q(s["centre_at_0"][0]), Q(s["centre_at_0"][1]))
        c1 = pt(Q(s["centre_at_1"][0]), Q(s["centre_at_1"][1]))
        r = Q(s["r"])
        self.assertTrue(v.certify_moving_disk(A, B, c0, c1, r * r))
        self.assertIsNone(v.certify_moving_disk(A, B, c0, c1, (r * Q(11, 10)) ** 2))  # too big fails
        for P in (A, B):                                        # exact census respects the floor
            self.assertTrue(all(sq.side_squared >= d["side_squared_cutoff"] for sq in census_polygon(P).squares))

    def test_beats_v264_floor_on_every_family(self):
        out = v.cutoff_comparison()
        self.assertEqual(len(out), 5)
        for name, rec in out.items():
            self.assertGreater(rec["improvement_factor"], 1.5, name)

    def test_cutoff_unchanged_by_collinear_subdivision(self):
        sq = tuple(pt(x, y) for x, y in ((-1, -1), (1, -1), (1, 1), (-1, 1)))
        ref = tuple(pt(x, y) for x, y in ((-1, -1), (0, -1), (1, -1), (1, 1), (-1, 1)))
        a, b = v.auto_disk_cutoff(sq, sq), v.auto_disk_cutoff(ref, ref)
        self.assertEqual(a["side_squared_cutoff"], b["side_squared_cutoff"])
        self.assertEqual(b["removed_collinear_vertices"], [1])
        A, B = self.fam["v254_collision_window"]
        mid = lambda P: tuple(x for i, p in enumerate(P) for x in (p, pt((p.x + P[(i + 1) % len(P)].x) / 2,
                                                                         (p.y + P[(i + 1) % len(P)].y) / 2)))
        self.assertEqual(v.auto_disk_cutoff(A, B)["side_squared_cutoff"],
                         v.auto_disk_cutoff(mid(A), mid(B))["side_squared_cutoff"])

    def test_nonconvex_motion_fails_closed(self):
        P, _, _ = v.read_polygon(v255.POCHUEV_JSON)
        with self.assertRaises(Exception):
            v.auto_disk_cutoff(P, P)


@unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")
class Sectors(unittest.TestCase):
    def test_twin_spike_graceful_or_off_disk(self):
        r = v.sector_check()
        for gap, rec in r.items():
            self.assertTrue(rec["lemma_holds"], gap)
            self.assertTrue(rec["small_or_ungraceful_have_both_diagonals_off_disk"], gap)
            self.assertTrue(rec["small_squares_all_ungraceful"], gap)
            big = [q for q in rec["squares"] if q["half_diagonal"] >= 1]
            self.assertTrue(big and all(q["graceful"] and q["a_diagonal_crosses_open_disk"] for q in big))


if __name__ == "__main__":
    unittest.main()
