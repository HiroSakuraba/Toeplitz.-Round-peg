"""v338 interior-excursion theorem and cutoff parity. Needs numpy (the v331 square search); set
REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
try:
    from robust_parity_v338 import (bent_base_check, certified_cutoff, excursion_walls, exact_side_ranges,
                                    interior_case, largest_square_positions, side_gaps, square_sides)
    from robust_squares_v336 import perturb, terrain_polygon
    from terrain_gluing_v334 import largest_square, terrain_i
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")


@needs_numpy
class InteriorExcursions(unittest.TestCase):
    def test_walls(self):
        self.assertEqual([excursion_walls([1, 5, 5, 5, 1], a) for a in largest_square_positions([1, 5, 5, 5, 1])], [(1, 4)])
        self.assertEqual(interior_case([1, 3, 3, 1, 4, 2]), 1)
        self.assertIsNone(interior_case([2, 7]))                 # both end walls
        self.assertIsNone(interior_case([1, 8, 8, 8]))           # right end wall
        self.assertIsNone(interior_case([2, 2, 2]))

    def test_bent_base_hypotheses_hold_on_perturbations(self):
        delta, eps = Q(1, 100), Q(1, 1000)
        for hs in ([1, 5, 5, 5, 1], [1, 3, 3, 1, 4, 2], [3, 1, 4, 4, 4, 4, 2, 6]):
            P, i = terrain_polygon(hs), terrain_i(hs)
            for seed in range(4):
                Pp = perturb(P, seed, eps)
                floor, _, _ = bent_base_check(hs, Pp, delta)
                self.assertEqual(floor, i-4*delta)
                self.assertGreaterEqual(largest_square(Pp)[0], floor**2)   # the guaranteed square is there

    def test_end_wall_case_is_rejected(self):
        with self.assertRaises(AssertionError):
            bent_base_check([2, 7], perturb(terrain_polygon([2, 7]), 0, Q(1, 1000)), Q(1, 100))


@needs_numpy
class CutoffParity(unittest.TestCase):
    def test_sharp_example(self):
        P = terrain_polygon([2, 7])
        ranges = sorted(exact_side_ranges(P))
        self.assertIn((Q(2), Q(4)), ranges)                      # the family from the diamond to the corner square
        lists = [square_sides(perturb(P, k, Q(1, 1000))) for k in range(4)]
        self.assertTrue(all(len(s) % 2 for s in lists))          # total count odd
        rho, _ = certified_cutoff(P, lists)
        self.assertAlmostEqual(rho, 2**0.5, places=9)           # certified robust side i/sqrt2

    def test_parity_can_fail_where_robust_squares_exist(self):
        hs = [8, 4, 3, 1, 1, 6]                                  # i = 3, robust side about 5.1
        P = terrain_polygon(hs)
        lists = [square_sides(perturb(P, k, Q(1, 1000))) for k in range(4)]
        rho, table = certified_cutoff(P, lists)
        self.assertAlmostEqual(rho, 2.0, places=9)               # below i/sqrt2 = 2.12
        self.assertTrue(all(not any(par) for (lo, hi), par in table if lo >= 2.0))
        self.assertGreater(min(s[0] for s in lists), 5.0)

    def test_no_gap_around_i_over_sqrt2(self):
        gaps, merged = side_gaps(terrain_polygon([4, 3, 1, 6, 3, 8]))
        self.assertTrue(any(lo <= 3/2**0.5 <= hi for lo, hi in merged))   # exact families cover it


if __name__ == "__main__":
    unittest.main()
