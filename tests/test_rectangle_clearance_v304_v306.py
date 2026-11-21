from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from rectangle_clearance_v305 import BOX, T, fixture, local_graph, matrix, pt, report, search
from radial_and_rough_v304_v306 import ETA, XS, radial_bounds, rough_shift, rough_witness, wave


class CrossingTests(unittest.TestCase):
    def test_radial_margins(self):
        bounds = radial_bounds()
        self.assertLess(bounds['radial_error'],bounds['angular_wall'])
        self.assertLess(5*ETA,Q(1,100))
        self.assertGreater(bounds['side_squared'],Q(69,100)**2)

    def test_tracked_edges_and_tie(self):
        r = report()
        self.assertGreater(Q(r['normalized_clearance']),Q(1,1024))
        self.assertLess(Q(r['endpoints'][0][0]),0)
        self.assertLess(Q(r['endpoints'][1][3]),0)

    def test_unmarked_reversed_and_reindexed(self):
        p = fixture()[::-1]
        p = p[4:]+p[:4]
        self.assertTrue(any(c['origin']==pt(0,0) and c['axis']==pt(1,1)
                            for c in search(p)))

    def test_extra_boundary_arc_rejected(self):
        edges = [(pt(-1,0),pt(1,0)),(pt(-1,Q(1,8)),pt(1,Q(1,8)))]
        self.assertIsNone(local_graph(edges,pt(0,0),0))

    def test_twelve_rough_rectangles_in_degree_box(self):
        A = matrix([Q(0),Q(1,4),Q(-1,2),Q(1,2)],1)
        for n in range(1,5):
            for k in (T,Q(1),1/T):
                v = rough_witness(k,n)
                self.assertLess(max(abs(sum(a*b for a,b in zip(row,v))) for row in A),BOX)
                self.assertGreater(min(1,k)**2*((1+v[2])**2+v[3]**2),Q(94,100)**2)

    def test_infinite_tail_and_retired_contacts(self):
        for i in range(4):
            self.assertNotEqual(rough_shift(i,Q(XS[i]),4),0)
            for x in (Q(1,13),Q(1,5),Q(1,3)):
                for n in range(5):
                    self.assertLessEqual(abs(rough_shift(i,XS[i]+x,9)-rough_shift(i,XS[i]+x,n)),
                                         Q(1,2048*2**n))
        for n in range(2,6):
            l = Q(1,4**n)
            self.assertEqual(wave(3*l/2,n+3)-(wave(l,n+3)+wave(2*l,n+3))/2,Q(1,2**n))

    def test_shrinking_bump_physical_scale(self):
        # v272's four bump contacts: a unit-square normalisation hides their scale.
        for n in (4,16,64):
            e = Q(1,n)
            h = Q(1,3)
            corners = (pt(h,1),pt(h+e,1),pt(h+e,1+e),pt(h,1+e))
            a,b,c,_ = corners
            self.assertEqual(((c.x-a.x)**2+(c.y-a.y)**2)/4,e**2/2)
            self.assertLess(b.x-a.x,Q(94,100))
