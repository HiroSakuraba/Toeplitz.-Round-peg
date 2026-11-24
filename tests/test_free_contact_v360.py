from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from free_contact_v360 import free_contact_certificate
from polygon_component_v341 import PolygonComponents, degree
from polygon_orbit_census_v223 import pt


def fixture(height, folded=False):
    left, right = Q(49,10), Q(51,10)
    bottom = [(Q(11,2),0),(Q(13,2),0)]
    red = [(11,Q(13,2)),(11,Q(11,2))]
    blue = [(0,Q(11,2)),(0,Q(9,2))]
    if folded:
        for arc, axis in ((bottom,0),(red,1),(blue,1)):
            a,b = arc
            for i,(along,offset) in enumerate(((3,1),(3,2),(2,2),(2,3),(4,3))):
                p = list(a)
                p[axis] += (b[axis]-a[axis])*along/5
                p[1-axis] += Q(offset,1024)
                arc.insert(i+1,tuple(p))
    poly = [(0,0),*bottom,(15,0),(15,12),(11,12),*red,(11,4),(8,4),
            (8,12),(right,12),(right,10),(5,10),(5,height),(left,height),
            (left,9),(2,9),(2,12),(0,12),*blue]
    return poly, [(poly.index(arc[0]),poly.index(arc[-1])) for arc in (bottom,red,blue)]


class FreeContactTests(unittest.TestCase):
    def test_free_contact_transfer(self):
        for height in (12, 11, Q(21,2)):
            for folded in (False,True):
                poly, windows = fixture(height,folded)
                a,b,floor = free_contact_certificate(poly, windows, 11, Q(1,256))
                self.assertLess(a,b)
                self.assertGreater(floor,Q(777,100)**2)
            poly,_ = fixture(height)
            engine = PolygonComponents(poly)
            expected = {Q(51,10): 1}
            if height > 11:
                expected.update({Q(49,10):1, Q(5):1})
            if height == 11:
                expected[Q(99,20)] = 0
            for w, parity in expected.items():
                seed = ((11-w,0),(11,11-w),(w,11),(0,w))
                cells = engine.trace(tuple(pt(*v) for v in seed))
                self.assertEqual(abs(sum(s for _,s in degree(cells)[1])),parity)

    def test_clearance_and_traversal_failures(self):
        poly,windows = fixture(12,True)
        self.assertIsNone(free_contact_certificate(poly,windows,11,0))
        self.assertIsNone(free_contact_certificate(poly,windows,11,Q(1,8)))
        with self.assertRaises(ValueError):
            free_contact_certificate(poly,[windows[0]]*3,11,Q(1,256))
