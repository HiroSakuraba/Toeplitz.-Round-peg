"""v240 exact-rational regression checks and diagnostic geometric controls."""
from fractions import Fraction as F
from math import atan, cos, exp, log, pi, sin, sqrt
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from log_spiral_exclusion_v240 import (
    THRESHOLD, CUTOFF, WIDTH_FLOOR, Certificate, exact_margins,
    exp_lower, exclusion_certificate, y_over_R, report,
)
from log_spiral_exclusion_v239 import exclusion_certificate as old_exclusion


class V240CertificateTests(unittest.TestCase):
    def test_cutoff_constants(self):
        self.assertEqual(THRESHOLD, F(17,20))
        self.assertEqual(CUTOFF, F(4,15))
        self.assertEqual(WIDTH_FLOOR, F(11,30))
        self.assertLess(3*CUTOFF, 1)

    def test_pi_bound_24_gon_exact(self):
        m = exact_margins()
        u = m["nested_root_strict_upper"]
        self.assertGreater(m["sqrt3_strict_upper"]**2, 3)
        self.assertGreater(u**2, 2+m["sqrt3_strict_upper"])
        self.assertEqual(144*(2-u), F(25,8)**2)
        self.assertEqual(m["pi_strict_lower"], F(25,8))

    def test_logarithmic_inequality_exact(self):
        m = exact_margins()
        self.assertGreater(exp_lower(F(53,40),5), F(15,4))
        self.assertEqual(m["tangent_span_strict_upper"], F(53,34))
        self.assertLess(m["tangent_span_strict_upper"], F(25,16))

    def test_sine_cosine_alternating_bounds_exact(self):
        m = exact_margins()
        self.assertGreater(m["cos_taylor_strict_lower"], 0)
        self.assertLess(17*m["sin_taylor_strict_upper"],
                        20*m["cos_taylor_strict_lower"])
        self.assertEqual(m["atan_strict_lower"], F(108,125))

    def test_sqrt_denominator_exact(self):
        m=exact_margins()
        self.assertLess(m["sqrt_denominator_strict_lower"]**2,
                        1+THRESHOLD**2)

    def test_transverse_envelope_margin_exact(self):
        m=exact_margins()
        self.assertGreater(m["exp_transverse_partial_sum"],
                           F(3000,1441))
        self.assertEqual(1/(m["sqrt_denominator_strict_lower"]*WIDTH_FLOOR),
                         F(3000,1441))

    def test_exact_endpoint(self):
        c=exclusion_certificate(F(17,20))
        self.assertTrue(c.excluded)
        self.assertEqual(c.cutoff, F(4,15))
        self.assertEqual(c.triangle_width_floor_over_R, F(11,30))

    def test_monotonic_extension(self):
        for a in (F(850001,1000000), F(86,100), F(9,10), 1, 4, 10000):
            with self.subTest(a=a):
                self.assertTrue(exclusion_certificate(a).excluded)

    def test_not_established_below(self):
        for a in (F(849999,1000000), F(4,5), F(3,4), F(1,10)):
            with self.subTest(a=a):
                self.assertEqual(exclusion_certificate(a).status,
                                 "not-established")

    def test_reject_rounded_floats_and_bools(self):
        for a in (0.85, 1.0, True, False):
            with self.subTest(a=a):
                with self.assertRaises(TypeError):
                    exclusion_certificate(a)
        self.assertTrue(exclusion_certificate("0.85").excluded)
        self.assertFalse(exclusion_certificate("0.84999999999999").excluded)

    def test_reject_nonpositive_domain(self):
        for a in (0,-1,"-1/3"):
            with self.assertRaises(ValueError):
                exclusion_certificate(a)

    def test_exponential_partial_sum_validation(self):
        for x,n in ((0,2),(-1,2),(1,-1)):
            with self.assertRaises(ValueError):
                exp_lower(x,n)

    def test_monotonic_y_envelope_numerical_only(self):
        for a in (0.85,0.86,0.9,1,4):
            delta = log(15/4)/a
            for j in range(201):
                t = delta*j/200
                y = y_over_R(a,t)
                self.assertGreaterEqual(y,-1e-14)
                self.assertLess(y,11/30)

    def test_tangent_axis_numerical_only(self):
        for a in (0.85,0.86,0.9,1,4):
            delta = log(15/4)/a
            axis = delta/2-atan(1/a)
            self.assertLess(delta, pi/2)
            for j in range(201):
                theta=delta*j/200
                self.assertLess(abs(theta-atan(1/a)-axis),pi/4)

    def test_regression_stronger_than_v239(self):
        self.assertTrue(exclusion_certificate(F(17,20)).excluded)
        self.assertFalse(old_exclusion(F(17,20)).excluded)
        for a in (F(9,10),1,4):
            self.assertTrue(exclusion_certificate(a).excluded)
            self.assertTrue(old_exclusion(a).excluded)

    def test_square_edge_orthogonality_sanity(self):
        for j in range(1,99):
            theta=j*pi/99
            u=(cos(theta),sin(theta))
            v=(-sin(theta),cos(theta))
            self.assertTrue(abs(u[1]) >= abs(u[0]) or
                            abs(v[1]) >= abs(v[0]))

    def test_report_limits(self):
        r=report()
        self.assertIn("17/20",r["theorem"])
        self.assertIn("NOT a closed Jordan curve",r["scope"])
        self.assertIn("not established",r["smaller_a"])


if __name__=="__main__":
    unittest.main()
