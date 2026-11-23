from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from notch_contact_switch_v356 import switch_witnesses
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


def trace(dims, seed):
    W,H,p,q,g,A,B,h = dims
    poly = ((0,0),(W,0),(W,H),(B,H),(B,h),(A,h),(A,H),
            (q,H),(q,g),(p,g),(p,H),(0,H))
    cells = PolygonComponents(poly).trace(tuple(pt(*z) for z in seed))
    return min(segment_floor(c.squares) for c in cells.values()), abs(sum(s for _,s in degree(cells)[1]))


class SwitchTests(unittest.TestCase):
    def test_roof_contact(self):
        for H in (Q(15,2), Q(799,100), 8, Q(801,100), 14):
            dims = (14,H,3,5,4,7,12,5)
            method, seed, floor = switch_witnesses(*dims)[0]
            self.assertEqual(method, "roof switch")
            self.assertEqual(trace(dims, seed), (floor, 1))

    def test_far_wall_collision_and_takeover(self):
        for g in (Q(17,2), Q(899,100), 9, Q(901,100), Q(19,2)):
            dims = (15,14,2,5,g,8,13,4)
            (method, seed, floor), = switch_witnesses(*dims)
            self.assertEqual(method, "floor square")
            self.assertEqual(trace(dims, seed), (floor, 1))
            if g <= 9:
                for t in {g, 9}:
                    self.assertEqual(trace(dims, ((4,0),(t,4),(t-4,t),(0,t-4))),
                                     (16+(t-4)**2, int(g<9)))

    def test_two_events_at_once(self):
        for g in (Q(89,10), 9, Q(91,10)):
            dims = (16,15,2,5,g,8,14,4)
            (_, seed, floor), = switch_witnesses(*dims)
            self.assertEqual(floor, max(g,9)**2)
            self.assertEqual(trace(dims, seed), (floor, 1))
