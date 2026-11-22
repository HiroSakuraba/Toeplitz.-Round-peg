"""v331 real actions from the preferred capping. Tests that need the residue census need numpy;
set REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import math
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from angular_windows_v328_v330 import area, pi_enclosure, radius_squared
from floer_actions_v331 import PolygonParam, ellipse_major_action, square_actions
from floer_parity_v331 import calibration, square_generators
try:
    import numpy  # noqa: F401
    from floer_actions_v331 import FastParam, carrier_filter
    from residue_census_v331 import certify, float_search
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")

PENTAGON = [(0, 0), (7, 1), (8, 6), (3, 9), (-2, 5)]
TWELVE = [(0, 0), (5, 1), (9, 0), (11, 4), (8, 5), (10, 9), (6, 8), (4, 11), (2, 7), (-2, 9), (-1, 4), (-3, 2)]


class EllipseClosedForm(unittest.TestCase):
    def test_derivative_property_and_endpoints(self):
        a, b = 3.0, 1.0
        for theta in (0.3, 1.0, 1.6, 2.4):
            h = 1e-6
            slope = (ellipse_major_action(a, b, theta+h)-ellipse_major_action(a, b, theta-h))/(2*h)
            q1 = 1/(math.cos(theta/2)**2/a**2+math.sin(theta/2)**2/b**2)
            self.assertAlmostEqual(slope, q1, places=5)
        self.assertAlmostEqual(ellipse_major_action(a, b, 1e-9), 0.0, places=6)
        self.assertAlmostEqual(ellipse_major_action(a, b, math.pi-1e-12), math.pi*a*b, places=6)


@needs_numpy
class ConvexExamples(unittest.TestCase):
    def squares(self, P):
        quads, *_ = float_search(P, qmax2=1e9)
        found, _ = certify(P, quads, qmax2=Q(10**9))
        return [v[0] for v in found.values()]

    def test_pentagon_carrier_and_duality(self):
        P = [(Q(x), Q(y)) for x, y in PENTAGON]
        A, top = area(P), pi_enclosure(30)[1]*radius_squared(P)/2
        (sq,) = self.squares(P)
        acts = square_actions(PolygonParam(P), sq)
        signs = [s for _, s in square_generators(P, sq)]
        even_action = acts[signs.index(calibration())]
        self.assertTrue(A/2 <= even_action <= top)          # the only grading-2 generator is the carrier
        self.assertEqual(sum(acts), A)                      # its odd partner carries l1 = A - l2

    def test_actions_agree_with_residues_and_duality(self):
        for poly in (PENTAGON, TWELVE):
            P = [(Q(x), Q(y)) for x, y in poly]
            A = area(P)
            for loop in (P, P[::-1]):
                par = PolygonParam(loop)
                for sq in self.squares(loop):
                    acts = square_actions(par, sq)
                    self.assertEqual(sum(acts), A)
                    for (r, _), act in zip(square_generators(loop, sq), acts):
                        self.assertEqual((act-r) % A, 0)

    def test_fast_path_matches_exact_path(self):
        P = [(Q(x), Q(y)) for x, y in TWELVE]
        slow, fast = PolygonParam(P), FastParam(P)
        for sq in self.squares(P):
            self.assertEqual(square_actions(slow, sq), fast.square_actions(sq))


@needs_numpy
class SpiralCarriers(unittest.TestCase):
    def test_twenty_four_turns(self):
        row = carrier_filter(24, workers=2)
        self.assertEqual(row["duality_failures"], 0)
        self.assertEqual(row["max_capping_turns"], 47)
        self.assertEqual((row["survivors"], row["no_real_action_in_window"], row["odd_only"],
                          row["vertex_contact"], row["candidates"]), (2617, 1573, 500, 43, 501))


if __name__ == "__main__":
    unittest.main()
