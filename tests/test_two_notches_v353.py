from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


def component(polygon, square):
    seed = tuple(pt(*p) for p in square)
    cells = PolygonComponents(polygon).trace(seed)
    return (min(segment_floor(c.squares) for c in cells.values()),
            abs(sum(sign for _, sign in degree(cells)[1])))


def roof(width, height, cuts):
    vertices = [(0, 0), (width, 0), (width, height)]
    for left, right, floor in reversed(cuts):
        vertices.extend(((right,height), (right,floor), (left,floor), (left,height)))
    return (*vertices, (0,height))


class TwoNotchesTests(unittest.TestCase):
    def test_transfer_at_and_across_collision(self):
        for h in (2, 3, Q(7, 2)):
            for H in (Q(19, 2), 10, Q(21, 2), 32):
                poly = roof(11, H, [(3,4,4), (8,10,h)])
                if H < 10:
                    u, v = h, H-h
                elif H == 10:
                    u, v = 5, 5
                else:
                    u, v = 6, 4
                t = u+v
                floor, parity = component(poly, ((u,0),(t,u),(v,t),(0,v)))
                self.assertEqual(parity, 1)
                self.assertEqual(floor, 50 if H == 10 else u*u+v*v)
        poly = roof(11, 10, [(3,4,4), (8,10,3)])
        self.assertEqual(component(poly, ((0,0),(10,0),(10,10),(0,10))), (58,0))

    def test_square_box_multiple_cuts(self):
        for cuts in ([(2,3,2),(5,8,4)], [(1,2,1),(4,5,7),(6,9,3)]):
            poly = roof(10, 10, cuts)
            floor, parity = component(poly, ((0,0),(10,0),(10,10),(0,10)))
            self.assertGreaterEqual(floor, 50)
            self.assertEqual(parity, 1)
            for first, second in zip(cuts, cuts[1:]):
                v = Q(first[1]+second[0], 2)
                u = 10-v
                self.assertEqual(component(poly, ((u,0),(10,u),(v,10),(0,v)))[1], 0)
