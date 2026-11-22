"""v331 mod-2 Floer grading. The survivor test needs numpy (the residue census); set
REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from angular_windows_v328_v330 import area
from floer_parity_v331 import calibration, ellipse_rectangle_signs, euler_characteristic, square_generators
try:
    import numpy  # noqa: F401
    from floer_parity_v331 import parity_filter
    from residue_census_v331 import certify, float_search
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"


class Calibration(unittest.TestCase):
    def test_ellipse_branches_have_opposite_constant_signs(self):
        for a, b in ((2.0, 1.0), (5.0, 1.0), (1.2, 1.0)):
            signs = {ellipse_rectangle_signs(a, b, t) for t in (0.2, 0.6, 1.0, 1.4)}
            self.assertEqual(signs, {(-1, 1)})
        self.assertEqual(calibration(), -1)


@unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")
class EulerCharacteristic(unittest.TestCase):
    def all_squares(self, P):
        quads, *_ = float_search(P, qmax2=1e9)
        found, _ = certify(P, quads, qmax2=Q(10**9))
        return [v[0] for v in found.values()]

    def test_signs_sum_to_zero_in_both_orientations(self):
        for poly in ([(0, 0), (5, 1), (9, 0), (11, 4), (8, 5), (10, 9), (6, 8), (4, 11), (2, 7), (-2, 9), (-1, 4), (-3, 2)],
                     [(0, 0), (7, 1), (8, 6), (3, 9), (-2, 5)]):
            P = [(Q(x), Q(y)) for x, y in poly]
            for loop in (P, P[::-1]):
                self.assertEqual(euler_characteristic(loop, self.all_squares(loop)), (0, 0))

    def test_square_generators_are_complementary(self):
        P = [(Q(x), Q(y)) for x, y in [(0, 0), (7, 1), (8, 6), (3, 9), (-2, 5)]]
        (sq,) = self.all_squares(P)
        (r, s), (r2, s2) = square_generators(P, sq)
        self.assertEqual(r+r2, area(P))     # actions r and A - r
        self.assertEqual(s, -s2)            # one generator of each parity on a single convex square


@unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")
class SurvivorParity(unittest.TestCase):
    def test_twenty_four_turns(self):
        row = parity_filter(24, workers=2)
        self.assertEqual((row["survivors"], row["even"], row["odd"], row["vertex_contact"], row["ambiguous_window"]),
                         (2617, 1280, 1273, 64, 0))


if __name__ == "__main__":
    unittest.main()
