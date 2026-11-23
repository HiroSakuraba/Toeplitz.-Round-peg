"""v335: even-loop counterexample and fragility of v316 squares. Needs numpy; set
REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import math
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
try:
    from even_loop_fragility_v335 import bumped, circle_antipodal_corners, core_loop_winding, perturbed_terrain
    from terrain_gluing_v334 import largest_square
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")


@needs_numpy
class EvenLoopCounterexample(unittest.TestCase):
    def test_circle_band_is_zero_free_with_even_core_loops(self):
        self.assertEqual(circle_antipodal_corners(), (math.sqrt(5), math.sqrt(5)))
        for c in (0.3, 0.5, 0.7):
            w, margin = core_loop_winding(1.0, 1.0, c)
            self.assertEqual(w, 0)
            self.assertGreater(margin, 0.25)
        self.assertAlmostEqual(core_loop_winding(1.0, 1.0, 0.5)[1], math.sqrt(5)-1, places=4)
        self.assertGreater(2*math.sin(0.3*math.pi), math.sqrt(2))     # the claimed certificate exceeds every square

    def test_generic_ellipse_behaves_the_same(self):
        for c in (0.3, 0.5, 0.7):
            w, margin = core_loop_winding(1.05, 1.0, c)
            self.assertEqual(w, 0)
            self.assertGreater(margin, 0.2)


@needs_numpy
class Fragility(unittest.TestCase):
    def test_any_positive_bump_depth_drops_seven_to_five(self):
        self.assertEqual(largest_square(bumped(Q(0)))[0], 49)
        for t in (Q(1, 1000), Q(1, 10), Q(1)):
            self.assertEqual(largest_square(bumped(t))[0], 25)

    def test_bump_must_swallow_the_corner(self):
        self.assertEqual(largest_square(bumped(Q(1), Q(7)))[0], 49)
        self.assertEqual(largest_square(bumped(Q(1), Q(69, 10)))[0], 25)

    def test_tiny_perturbations_of_the_terrain_alone(self):
        sides = [float(largest_square(perturbed_terrain(seed))[0])**0.5 for seed in (0, 1, 2)]
        self.assertLess(sides[0], 5.01)        # seed 0 loses the side-7 square
        self.assertGreater(sides[1], 6.99)     # seed 1 keeps it


if __name__ == "__main__":
    unittest.main()
