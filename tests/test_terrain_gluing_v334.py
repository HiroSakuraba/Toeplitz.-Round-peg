"""v334 two-terrain gluing. Needs numpy (the v331 square search); set REQUIRE_NUMERIC_DEPS=1 to
make a missing numpy a failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import os
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import PositiveDimensionalSquares, census_polygon, pt
try:
    from terrain_gluing_v334 import (WITNESSES, analyse, column_maps, interior_square, largest_square,
                                     polygon, terrain_i)
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")


@needs_numpy
class SingleTerrains(unittest.TestCase):
    def test_v316_bound_holds_on_random_terrains(self):
        rng = random.Random(7)
        for _ in range(60):
            top = [rng.randint(1, 6) for _ in range(rng.randint(1, 6))]
            P = polygon(0, len(top), {x: h for x, h in enumerate(top)}, {x: 0 for x in range(len(top))})
            side2, _ = largest_square(P)
            self.assertGreaterEqual(side2, terrain_i(top)**2)


@needs_numpy
class Witnesses(unittest.TestCase):
    def test_exact_ratios_and_lost_corners(self):
        a = analyse(*WITNESSES["ratio_2_over_sqrt5"])
        self.assertEqual((a["i"], a["M_squared"]), (5, "20"))
        self.assertEqual(a["base_squares"]["upper"]["corners_lost_to_seam"], [[5, 0]])
        b = analyse(*WITNESSES["ratio_5_over_7"])
        self.assertEqual((b["i"], b["M_squared"], b["i_upper"]), (7, "25", 7))
        self.assertEqual(b["base_squares"]["upper"]["corners_lost_to_seam"], [[7, 0]])
        self.assertEqual(b["polygon"], [[0, 0], [6, 0], [6, -1], [8, -1], [8, 3], [7, 3], [7, 7], [3, 7], [3, 8], [0, 8]])

    def test_exhaustive_census_confirms_five_sevenths(self):
        P = [pt(x, y) for x, y in analyse(*WITNESSES["ratio_5_over_7"])["polygon"]]
        self.assertEqual(len(census_polygon(P, cutoff_squared=Q(25)).squares), 0)  # nothing above side 5
        with self.assertRaises(PositiveDimensionalSquares):                       # a family reaches 25
            census_polygon(P, cutoff_squared=Q(24))


class StraddleReduction(unittest.TestCase):
    def test_max_terrain_square_is_at_least_half(self):
        if not HAVE_NUMPY and not REQUIRE:
            self.skipTest("numpy not installed")
        rng = random.Random(11)
        for _ in range(3000):
            top = [rng.randint(1, 9) for _ in range(rng.randint(1, 9))]
            bot = [rng.randint(1, 9) for _ in range(rng.randint(1, 9))]
            bs = rng.randint(-len(bot)+1, len(top)-1)
            L, R, f, g = column_maps(top, bot, 0, bs)
            i, _ = interior_square(L, R, f, g)
            self.assertGreaterEqual(2*max(terrain_i(top), terrain_i(bot)), i)


if __name__ == "__main__":
    unittest.main()
