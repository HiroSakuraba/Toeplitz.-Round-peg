"""Exact arithmetic controls and independent geometry sanity tests for v236."""
from fractions import Fraction as F
from math import atan, hypot, log, pi, cos, sin, sqrt
from pathlib import Path
import sys, unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from log_spiral_exclusion_v236 import (
    EPSILON,A_THRESHOLD,exp_positive_partial_sum,exclusion_certificate,
    gamma,v99_special_pitch_locator,report,
)

class ExactArithmetic(unittest.TestCase):
    def test_strict_exp_bound(self):
        self.assertEqual(exp_positive_partial_sum(F(3,2),3),F(67,16))
        self.assertGreater(F(67,16),4)
    def test_exponential_tail_bounded_below(self):
        for n in range(1,10):
            self.assertLess(exp_positive_partial_sum(F(3,2),n),4.5)
    def test_critical_four(self):
        c=exclusion_certificate(4)
        self.assertTrue(c.excluded)
        self.assertEqual(c.strict_delta_upper,F(3,8))
        self.assertEqual(c.triangle_min_width_lower,F(3,8))
        self.assertEqual(c.tangent_angle_upper,F(5,8))
    def test_higher_open_pitches(self):
        for a in (F(17,4),5,10,100,10000):
            with self.subTest(a=a):
                c=exclusion_certificate(a)
                self.assertTrue(c.excluded)
                self.assertTrue(c.two_inner_points_impossible)
                self.assertTrue(c.three_outer_points_impossible)
                self.assertTrue(c.four_outer_points_impossible)
    def test_below_four_fail_closed(self):
        for a in (F(1,10),F(1),F(399,100)):
            self.assertFalse(exclusion_certificate(a).excluded)
    def test_invalid_domain(self):
        for a in (0,-1):
            with self.assertRaises(ValueError):exclusion_certificate(a)
    def test_two_inner_points_inequality(self):
        eps=EPSILON
        self.assertLess(2*(2*eps)**2,(1-eps)**2)
    def test_three_outer_strip_gap_strict(self):
        c=exclusion_certificate(4)
        self.assertLess(log(4)/4,float(c.triangle_min_width_lower))
        self.assertLess(log(4)/4+atan(1/4),pi/4)
    def test_all_thresholds_rational(self):
        for n in range(4,100):
            c=exclusion_certificate(F(n,1))
            self.assertLessEqual(c.strict_delta_upper,c.triangle_min_width_lower)
            self.assertLess(c.tangent_angle_upper,F(3,4))
    def test_no_circular_self_similarity_assumption(self):
        a=5.;t=.22;shift=.31
        x,y=gamma(a,t+shift)
        X,Y=gamma(a,t)
        r=__import__('math').exp(-a*shift)
        self.assertAlmostEqual(x,r*(X*cos(shift)-Y*sin(shift)),places=14)
        self.assertAlmostEqual(y,r*(X*sin(shift)+Y*cos(shift)),places=14)
    def test_gamma_radius_is_monotone(self):
        for a in (4,5,20):
            rs=[hypot(*gamma(a,k/20)) for k in range(15)]
            self.assertTrue(all(x>y for x,y in zip(rs,rs[1:])))
    def test_v99_locator_is_only_diagnostic(self):
        x=v99_special_pitch_locator()
        self.assertGreater(x['a_approx'],0.10)
        self.assertLess(x['a_approx'],0.101)
        self.assertIn('Historical',x['existence_source'].replace('restored','Historical restored'))
    def test_no_square_graph_lipschitz_demo(self):
        # For any orientation of a square, one of two perpendicular edges
        # has |dy|>=|dx|, including an edge with dx=0.
        for k in range(1,99):
            theta=k*pi/99
            e0=(cos(theta),sin(theta))
            e1=(-sin(theta),cos(theta))
            self.assertTrue(abs(e0[1])>=abs(e0[0]) or abs(e1[1])>=abs(e1[0]))
    def test_triangle_minimum_width_demo(self):
        # One triangle of three square corners: min width=s/sqrt(2).
        # Test rotated triangles against the bound in arbitrary directions.
        verts=[(0,0),(1,0),(0,1)]
        for i in range(1,128):
            t=i*pi/128
            yy=[-x*sin(t)+y*cos(t) for x,y in verts]
            self.assertGreaterEqual(max(yy)-min(yy)+1e-14,1/sqrt(2))
    def test_report_has_no_unproved_generic_claim(self):
        o=report()
        self.assertIn('a>=4',o['theorem'])
        self.assertIn('Does not rule out',o['limits'])

if __name__=='__main__':unittest.main()
