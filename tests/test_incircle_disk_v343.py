"""v343 incircle reading and disk version. Needs numpy (the v331 square search and the inradius
search); set REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import math
import os
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
try:
    from incircle_disk_v343 import inradius, named_shapes, random_star, regular_polygon, row, union_polygon
    from robust_squares_v336 import terrain_polygon
    from terrain_gluing_v334 import edges_through, terrain_i
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")


@needs_numpy
class IncircleReading(unittest.TestCase):
    def test_sharp_diamond_is_inscribed_in_the_incircle(self):
        P = terrain_polygon([2, 7])
        diamond = [(Q(1), Q(0)), (Q(2), Q(1)), (Q(1), Q(2)), (Q(0), Q(1))]
        centre, radius2 = (Q(1), Q(1)), Q(1)                      # incircle of the interior square [0,2]^2
        for v in diamond:
            self.assertEqual((v[0]-centre[0])**2+(v[1]-centre[1])**2, radius2)
            self.assertTrue(edges_through(P, v))                  # and it is inscribed in the terrain
        r, c = inradius(P)
        self.assertAlmostEqual(r, 1.0, places=6)                 # the incircle is the largest disk

    def test_constant_chain(self):
        # i/sqrt2 <= sqrt2 r <= i: the disk version implies Conjecture 336.1
        rng = random.Random(3)
        for _ in range(25):
            hs = [rng.randint(1, 7) for _ in range(rng.randint(2, 7))]
            i, (r, _) = terrain_i(hs), inradius(terrain_polygon(hs))
            self.assertLessEqual(i/2, r+1e-7)
            self.assertLessEqual(math.sqrt(2)*r, i+1e-7)


@needs_numpy
class Inradius(unittest.TestCase):
    def test_regular_polygon(self):
        P = regular_polygon(12)
        R = max(math.hypot(float(x), float(y)) for x, y in P)
        self.assertAlmostEqual(inradius(P)[0], R*math.cos(math.pi/12), places=2)   # rational rounding 1e-3

    def test_union_witness(self):
        r, _ = inradius(named_shapes()["v334 glued witness A"])
        self.assertAlmostEqual(r, 3.5, places=6)


@needs_numpy
class DiskVersion(unittest.TestCase):
    def test_equality_cases_stay_at_the_bound(self):
        for name in ("terrain (2,7)", "regular 8-gon"):
            out = row(named_shapes()[name], seeds=range(4))
            self.assertGreater(out["disk_ratio"], 0.99)

    def test_star_polygons_are_simple(self):
        from polygon_orbit_census_v223 import pt, validate_jordan_polygon
        for seed in range(20):
            validate_jordan_polygon([pt(*p) for p in random_star(seed)])


if __name__ == "__main__":
    unittest.main()
