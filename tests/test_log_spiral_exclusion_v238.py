"""v238: exact inequalities, fail-closed API, and geometry sanity checks."""
from fractions import Fraction as F
from math import atan, cos, exp, log, pi, sin, sqrt
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from log_spiral_exclusion_v238 import (
    INNER_RADIUS, PROVED_THRESHOLD, exclusion_certificate,
    positive_exponential_partial_sum, report, transverse_coordinate,
)


class ExactThresholdTests(unittest.TestCase):
    def test_e_lower_bound_strict(self):
        self.assertEqual(positive_exponential_partial_sum(1, 3), F(8, 3))
        self.assertEqual(positive_exponential_partial_sum(F(3, 2), 3), F(67, 16))
        self.assertGreater(F(67, 16), 4)

    def test_exact_boundary(self):
        c = exclusion_certificate(1)
        self.assertTrue(c.excluded)
        self.assertEqual(c.openness, 1)
        self.assertEqual(c.inner_cutoff, F(1, 4))
        self.assertEqual(c.triangle_width_floor_over_R, F(3, 8))
        self.assertEqual(c.strip_width_strict_upper_over_R, F(3, 8))
        self.assertEqual(c.tangent_half_span_strict_upper, F(3, 4))

    def test_upward_monotonicity(self):
        for a in (F(1000001, 1000000), F(17, 10), F(4), F(100), F(10**8)):
            with self.subTest(a=a):
                c = exclusion_certificate(a)
                self.assertTrue(c.excluded)
                self.assertLessEqual(c.strip_width_strict_upper_over_R,
                                     c.triangle_width_floor_over_R)
                self.assertLessEqual(c.tangent_half_span_strict_upper, F(3, 4))

    def test_below_one_fails_closed(self):
        for a in ("0.999999999999999999999", F(999, 1000), F(1, 2), F(1, 10)):
            with self.subTest(a=a):
                c = exclusion_certificate(a)
                self.assertFalse(c.excluded)
                self.assertEqual(c.status, "not-established")

    def test_float_and_bool_rejected(self):
        for a in (1.0, 0.999999999999, True, False):
            with self.subTest(a=a):
                with self.assertRaises(TypeError):
                    exclusion_certificate(a)

    def test_invalid_positive_domain(self):
        for a in (0, -1, "-1/2"):
            with self.subTest(a=a):
                with self.assertRaises(ValueError):
                    exclusion_certificate(a)

    def test_exponential_input_validation(self):
        for x, n in ((0, 2), (-1, 2), (1, -1)):
            with self.assertRaises(ValueError):
                positive_exponential_partial_sum(x, n)

    def test_two_inner_exact_separation(self):
        eps = INNER_RADIUS
        self.assertLess(2 * (2*eps)**2, (1-eps)**2)
        self.assertLess(F(1, 2), F(9, 16))

    def test_strict_strip_gap_at_one_numerical_sanity(self):
        # The proof uses e>8/3, NOT this floating-point comparison.
        self.assertLess(1/exp(1), 3/8)
        self.assertLess(log(4), 3/2)
        self.assertLess(log(4)/2, pi/4)

    def test_curve_points_fit_analytic_strip(self):
        for a in (1, 1.01, 1.5, 4):
            delta = log(4)/a
            for j in range(201):
                theta = delta*j/200
                y = transverse_coordinate(a, theta)
                self.assertGreaterEqual(y, -1e-15)
                self.assertLess(y, 3/(8*a))

    def test_graph_midpoint_axis_works_at_threshold(self):
        for a in (1, 1.01, 2, 4):
            delta = log(4)/a
            axis = delta/2 - atan(1/a)
            for j in range(201):
                theta = delta*j/200
                tangent = theta - atan(1/a)
                self.assertLess(abs(tangent-axis), pi/4)

    def test_perpendicular_adjacent_edges_obstruct_graph(self):
        for j in range(1, 80):
            angle = j*pi/80
            ux, uy = cos(angle), sin(angle)
            vx, vy = -uy, ux
            self.assertTrue(abs(uy) >= abs(ux) or abs(vy) >= abs(vx))

    def test_right_isosceles_triangle_min_width(self):
        verts = ((0,0),(1,0),(0,1))
        for j in range(1, 160):
            phi = pi*j/160
            ys = [-x*sin(phi)+y*cos(phi) for x,y in verts]
            self.assertGreaterEqual(max(ys)-min(ys)+1e-14, 1/sqrt(2))

    def test_new_theorem_extends_v236(self):
        from log_spiral_exclusion_v236 import exclusion_certificate as old
        for a in (4, 5, 10):
            self.assertTrue(old(a).excluded)
            self.assertTrue(exclusion_certificate(a).excluded)
        self.assertFalse(old(1).excluded)
        self.assertTrue(exclusion_certificate(1).excluded)

    def test_result_is_an_arc_claim_only(self):
        r = report()
        self.assertIn("a>=1", r["theorem"])
        self.assertIn("open arc", r["scope"])
        self.assertIn("Jordan closures", r["limits"])

    def test_threshold_constant(self):
        self.assertEqual(PROVED_THRESHOLD, 1)


if __name__ == "__main__":
    unittest.main()
