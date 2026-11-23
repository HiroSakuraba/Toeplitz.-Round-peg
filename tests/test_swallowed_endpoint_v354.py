from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from swallowed_endpoint_v354 import endpoint_witness
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class EndpointTests(unittest.TestCase):
    def test_component_at_all_contact_transitions(self):
        for A,B,d,h,p in ((8,10,1,3,3), (9,12,1,4,2), (10,13,1,5,2)):
            for q in (B-h, Q(B-h+A, 2)):
                for g in (Q(1,8), A-Q(1,100), A, A+Q(1,100), B-Q(1,8)):
                    with self.subTest(A=A, q=q, g=g):
                        seed, floor = endpoint_witness(A,B,d,h,p,q,g)
                        poly = ((0,0),(B+d,0),(B+d,B),(B,B),(B,h),(A,h),
                                (A,B),(q,B),(q,g),(p,g),(p,B),(0,B))
                        cells = PolygonComponents(poly).trace(tuple(pt(*z) for z in seed))
                        self.assertEqual(min(segment_floor(c.squares) for c in cells.values()),
                                         floor)
                        self.assertEqual(abs(sum(s for _,s in degree(cells)[1])), 1)
                        self.assertGreaterEqual(2*floor, A*A)

    def test_same_square_before_at_and_after_cut(self):
        seed = tuple(pt(*z) for z in ((5,0),(8,5),(3,8),(0,3)))
        for q in (Q(13,2), 7, Q(15,2)):
            poly = ((0,0),(11,0),(11,10),(10,10),(10,3),(8,3),(8,10),
                    (q,10),(q,4),(3,4),(3,10),(0,10))
            cells = PolygonComponents(poly).trace(seed)
            self.assertEqual(min(segment_floor(c.squares) for c in cells.values()), 34)
            self.assertEqual(abs(sum(s for _,s in degree(cells)[1])), 1)
