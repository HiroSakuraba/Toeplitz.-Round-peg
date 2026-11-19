"""v239 regression suite: exact rational margins, numerical sanity and scope."""
from fractions import Fraction as F
from math import atan, exp, log, pi, sin, sqrt
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from log_spiral_exclusion_v239 import (
    THRESHOLD, INNER_CUTOFF, THREE_VERTEX_WIDTH, exact_margins,
    exp_lower, exclusion_certificate, report, transverse_coordinate,
)
from log_spiral_exclusion_v238 import exclusion_certificate as v238


class V239ExactTests(unittest.TestCase):
    def test_cutoff_and_side_floor(self):
        self.assertEqual(INNER_CUTOFF, F(13, 50))
        self.assertEqual(THREE_VERTEX_WIDTH, F(37,100))
        self.assertEqual(THRESHOLD, F(9,10))

    def test_two_inner_separation_exact(self):
        m = exact_margins()
        self.assertEqual(m["two_inner_left_squared"], F(338,625))
        self.assertEqual(m["two_inner_right_squared"], F(1369,2500))
        self.assertLess(m["two_inner_left_squared"], m["two_inner_right_squared"])

    def test_angular_exponential_exact(self):
        m = exact_margins()
        self.assertGreater(m["exp_27_20_lower"], F(50,13))
        self.assertEqual(m["tangent_span_upper_at_threshold"], F(3,2))

    def test_arctangent_lower_bound_exact(self):
        m = exact_margins()
        self.assertLess(m["tan_4_5_sine_upper"],
                        F(10,9)*m["tan_4_5_cosine_lower"])
        self.assertGreater(m["tan_4_5_cosine_lower"], 0)

    def test_root_lower_bound_exact(self):
        self.assertGreater(F(181,100), F(16,9))

    def test_transverse_exponential_exact(self):
        m = exact_margins()
        self.assertGreater(m["exp_18_25_lower"], m["exp_18_25_target"])
        self.assertEqual(m["exp_18_25_target"], F(75,37))

    def test_threshold_endpoint(self):
        c = exclusion_certificate(F(9,10))
        self.assertTrue(c.excluded)
        self.assertEqual(c.minimum_triangle_width_over_R, F(37,100))

    def test_monotonic_extension(self):
        for a in (F(9000001,10000000), F(91,100), F(1), F(2), F(4), F(20)):
            with self.subTest(a=a):
                self.assertTrue(exclusion_certificate(a).excluded)

    def test_not_established_below(self):
        for a in (F(8999999,10000000), F(4,5), F(1,2), F(1,10)):
            with self.subTest(a=a):
                c = exclusion_certificate(a)
                self.assertFalse(c.excluded)
                self.assertEqual(c.status,"not-established")

    def test_input_precision_guard(self):
        for a in (0.9, 1.0, True, False):
            with self.subTest(a=a):
                with self.assertRaises(TypeError):
                    exclusion_certificate(a)
        self.assertTrue(exclusion_certificate("0.9").excluded)
        self.assertFalse(exclusion_certificate("0.899999999999999").excluded)

    def test_invalid_domain(self):
        for a in (0, -1, "-1/2"):
            with self.assertRaises(ValueError):
                exclusion_certificate(a)

    def test_invalid_exponential_inputs(self):
        for x,n in ((0,2),(-1,2),(1,-1)):
            with self.assertRaises(ValueError):
                exp_lower(x,n)

    def test_transverse_strip_numerical_only(self):
        for a in (0.9, 0.91, 1, 2, 4):
            delta = log(50/13)/a
            for j in range(201):
                theta = delta*j/200
                y = transverse_coordinate(a,theta)
                self.assertGreaterEqual(y,-1e-14)
                self.assertLess(y,37/100)

    def test_rotated_tangent_graph_numerical_only(self):
        for a in (0.9, 0.91, 1, 4):
            delta = log(50/13)/a
            axis = delta/2 - atan(1/a)
            self.assertLess(delta,pi/2)
            for j in range(201):
                theta = delta*j/200
                self.assertLess(abs(theta-atan(1/a)-axis),pi/4)

    def test_old_threshold_extension(self):
        self.assertFalse(v238(F(9,10)).excluded)
        self.assertTrue(exclusion_certificate(F(9,10)).excluded)
        for a in (1, 4, 10):
            self.assertTrue(v238(a).excluded)
            self.assertTrue(exclusion_certificate(a).excluded)

    def test_report_scope(self):
        r=report()
        self.assertIn("9/10",r["result"])
        self.assertIn("not a Jordan-curve",r["status"])
        self.assertIn("do not infer",r["below_threshold"])


if __name__=="__main__":
    unittest.main()
