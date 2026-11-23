from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class ReplacementSquareTests(unittest.TestCase):
    def test_crossing_and_width_orderings(self):
        for H in (Q(31, 4), Q(8), Q(33, 4)):
            c = 4 + (H-2)/2
            for W in (Q(15, 2), H, Q(12)):
                poly = ((0, 0), (W, 0), (W, 2), (4, 2), (c, H), (0, H))
                t = min(W, H)
                v = 4 + (W-2)/2 if W < H else H-2
                u = t-v
                seed = ((u, 0), (t, u), (v, t), (0, v))
                cells = PolygonComponents(poly).trace(tuple(pt(*p) for p in seed))
                floor = min(segment_floor(cell.squares) for cell in cells.values())
                self.assertGreaterEqual(2*floor, t*t)
                self.assertEqual(abs(sum(s for _, s in degree(cells)[1])), 1)
                if W == 12:
                    self.assertEqual(floor, 4+(H-2)**2)

    def test_finite_margin(self):
        e = Q(1, 64)
        self.assertLess(12*e, Q(1, 2))
        # sqrt(80) < 9, so this rational comparison suffices.
        self.assertLess((6+11*e)**2, 40)
        self.assertLess((8+2*e)**2, 72)
