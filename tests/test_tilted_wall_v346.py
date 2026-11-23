from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class TiltedWallTests(unittest.TestCase):
    def check_component(self, polygon, square, floor):
        cells = PolygonComponents(polygon).trace(tuple(pt(*p) for p in square))
        self.assertEqual(min(segment_floor(c.squares) for c in cells.values()), floor)
        self.assertEqual(abs(sum(sign for _, sign in degree(cells)[1])), 1)

    def test_inward(self):
        polygon = ((0, 0), (10, 0), (10, 4), (3, 4), (Q(1, 2), 9), (0, 9))
        self.check_component(polygon, ((4, 0), (6, 4), (2, 6), (0, 2)), 20)

    def test_outward_switch(self):
        for H in (81, 161, 321):
            W, c = H + 9, Q(H + 1, 2)
            polygon = ((0, 0), (W, 0), (W, 1), (1, 1), (c, H), (0, H))
            self.check_component(polygon, ((1, 0), (3, 1), (2, 3), (0, 2)), 5)
            d = 4*c/7
            square = ((d, 2*d-1), (3*d/2, 3*d-1), (d/2, H), (0, 5*d/2-1))
            self.check_component(polygon, square, 5*d*d/4)
            self.assertGreater(2*Q(H, 8)**2, 5)

    def test_corner_transition(self):
        polygon = ((0, 0), (12, 0), (12, 2), (4, 2), (7, 8), (0, 8))
        square = ((4, 2), (6, 6), (2, 8), (0, 4))
        cells = PolygonComponents(polygon).trace(tuple(pt(*p) for p in square))
        self.assertEqual(set(cells), {(2, 3, 4, 5), (3, 3, 4, 5)})
        self.assertEqual(sum(sign for _, sign in degree(cells)[1]), 0)
