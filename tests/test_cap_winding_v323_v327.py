import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cap_winding_v323_v327 import (
    gadget, fixture, integral, simple, winding, joint_formula, scaled_payment,
    payment, endpoint_model
)


class CapWindingTests(unittest.TestCase):
    def test_exact_unscaled_sum(self):
        for n in (1,2,3,8,17,32):
            p=gadget(n)
            i=p.index((0,0)); j=p.index((0,1))
            self.assertEqual(integral(p[:i+1])+integral(p[j:]),joint_formula(n))

    def test_registered_polygons(self):
        for k in (2,3,4):
            row=payment(k)
            self.assertEqual(row["side"],Q(1,k))
            self.assertGreater(row["action_gap"],Q(1,2))
            self.assertLessEqual(abs(row["area"]-16),Q(4,k*k))

    def test_exact_limit(self):
        last=Q(0)
        for k in (2,4,8,16,32,64,1024):
            v=scaled_payment(k)
            self.assertGreater(v,last)
            self.assertEqual(v,Q(19,24)-Q(1,2*k*k)-Q(1,384*k**4))
            self.assertLess(v,Q(19,24))
            last=v

    def test_joint_winding_large(self):
        for n in (1,2,7,16):
            p=gadget(n)
            i=p.index((0,0))
            self.assertEqual(winding(p[:i+1]+[p[0]],(Q(-1,8),Q(1,8))),n)

    def test_full_boundary_winding_small(self):
        for k in (2,3,4):
            p=fixture(k)
            self.assertEqual(winding(p+[p[0]],(Q(2),Q(2))),-1)
            self.assertEqual(winding(p+[p[0]],(Q(5),Q(5))),0)

    def test_crossing_and_backtrack_negative_controls(self):
        self.assertFalse(simple([(0,0),(2,2),(0,2),(2,0)]))
        self.assertFalse(simple([(0,0),(2,0),(1,0),(1,1)]))
        self.assertTrue(simple([(0,0),(1,0),(2,0),(2,1),(0,1)]))

    def test_dual_endpoint_models(self):
        for b in (Q(1,100),Q(1,4),Q(1,2)):
            for m in (Q(0),Q(1,100),Q(1)):
                f,g=endpoint_model(b,Q(1,2),m)
                self.assertEqual(f,b-m*b/2)
                self.assertEqual(g,m*b/2)
                for i in range(33):
                    x=Q(i,32)
                    fx,gx=endpoint_model(b,x,m)
                    self.assertLessEqual(gx,fx)
                    self.assertGreaterEqual(gx,0)
                    self.assertLessEqual(fx,b)
                    self.assertEqual(fx+endpoint_model(b,1-x,m)[1],b)
                    for j in range(33-i):
                        y=Q(j,32)
                        self.assertGreaterEqual(
                            fx+endpoint_model(b,y,m)[0],
                            endpoint_model(b,x+y,m)[0])

    def test_domain_checks(self):
        for n in (0,-1,1.5):
            with self.assertRaises(ValueError):
                gadget(n)
        with self.assertRaises(ValueError):
            fixture(1)
        with self.assertRaises(ValueError):
            scaled_payment(1)
        with self.assertRaises(ValueError):
            winding([(0,0),(1,1)],(0,0))
        with self.assertRaises(ValueError):
            endpoint_model(Q(3,4),Q(1,2))


if __name__=="__main__":
    unittest.main()
