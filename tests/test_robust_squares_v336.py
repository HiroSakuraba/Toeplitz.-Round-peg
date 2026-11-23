"""v336 robust squares of grid terrains. Needs numpy (the v331 square search); set
REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import math
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
try:
    from robust_squares_v336 import corner_square_zero, perturb, robust_estimate, terrain_polygon
    from terrain_gluing_v334 import edges_through, largest_square, terrain_i
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")


@needs_numpy
class SharpExample(unittest.TestCase):
    def test_two_seven_terrain(self):
        P = terrain_polygon([2, 7])
        self.assertEqual(terrain_i([2, 7]), 2)
        self.assertEqual(largest_square(P)[0], 4)                      # exact M = i = 2 (fragile)
        exact, robust, _ = robust_estimate(P, range(8))
        self.assertLess(robust, 2-0.5)                                 # perturbations lose the side-2 square
        self.assertAlmostEqual(robust, math.sqrt(2), delta=0.01)       # the diamond of side i/sqrt2 survives

    def test_perturbations_are_simple_and_small(self):
        P = terrain_polygon([3, 3, 4])
        Pp = perturb(P, 5, Q(1, 1000))
        self.assertTrue(all(abs(a-c) <= Q(1, 1000) and abs(b-d) <= Q(1, 1000) for (a, b), (c, d) in zip(P, Pp)))


@needs_numpy
class CornerZeros(unittest.TestCase):
    def test_corner_square_lies_on_the_boundary(self):
        for hs in ([5, 5, 5], [6, 2], [8, 8, 8, 7, 7, 7, 7, 3], [2, 7], [3, 1, 5]):
            c, P = corner_square_zero(hs), terrain_polygon(hs)
            for v in ((Q(0), Q(0)), (c, Q(0)), (c, c), (Q(0), c)):
                self.assertTrue(edges_through(P, v))


if __name__ == "__main__":
    unittest.main()
