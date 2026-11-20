"""Exact fixtures for the scaled degree, finite-coil envelope and tangent control."""

from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from polygon_orbit_census_v223 import _solve
from protected_rectangle_band_v284 import norm_inf
import rough_flat_contacts_v296_v299 as rough


class RoughFlatContactTests(unittest.TestCase):
    def test_replay(self):
        self.assertEqual(rough.report(), json.loads(rough.ARTIFACT.read_text(encoding="utf-8")))

    def test_scaled_side_estimate_with_quadratic_position_errors(self):
        for end, jet in rough.normalized_jets(Q(16)):
            with self.subTest(end=end):
                a, _, c, _ = jet
                b, d = Q(2, 3), -Q(4, 5)
                certificate = rough.certify_scaled_box(a, b, c, d)
                h = certificate["parameter_box_halfwidth"]
                e = certificate["quadratic_error_tolerance"]
                left = lambda t: (a*t*t+e*t*abs(t), t+b*t*t-e*t*abs(t))
                right = lambda t: (1+c*t*t-e*t*abs(t), t+d*t*t+e*t*abs(t))
                J = rough.limiting_jacobian(a, c)
                r = certificate["half_ratio_upper_before_source_domain_cutoff"]
                for signs in product((-1, 1), repeat=4):
                    z = tuple(v+sgn*h for v, sgn in zip(rough.FLAT_TAGS, signs))
                    residual = rough.rescaled_residual(left, right, z, r)
                    processed = _solve(J, residual)
                    error = tuple(v-sgn*h for v, sgn in zip(processed, signs))
                    self.assertLessEqual(norm_inf(error), 3*h/8)

    def test_parabolic_rectangle_at_a_known_positive_ratio(self):
        a, c, s = Q(24, 5), -Q(12, 5), Q(1, 1000)
        left = lambda t: (a*t*t, t)
        right = lambda t: (1+c*t*t, t)
        long = 1+(c-a)*s*s
        r = s/long
        z = tuple(v*long for v in rough.FLAT_TAGS)
        self.assertEqual(rough.rescaled_residual(left, right, z, r), (0,)*4)

    def test_coil_budget_after_shrinking_support(self):
        budget = Q(1, 100000)
        for length in (Q(1, 128), Q(1, 8192)):
            n = rough.first_coil_index(length, budget)
            self.assertLessEqual(rough.coil_quadratic_error(n, length), budget)
            self.assertGreater(rough.coil_quadratic_error(n-1, length), budget)
            self.assertLess(rough.coil_quadratic_error(n+1, length),
                            rough.coil_quadratic_error(n, length))

    def test_singular_jets_and_cubic_control(self):
        with self.assertRaisesRegex(ArithmeticError, "degenerate"):
            rough.certify_scaled_box(0, 0, 0, 0)
        self.assertTrue(rough.cubic_control()["same_flat_chord_and_tangent_directions"])
        self.assertEqual(rough.cubic_control()["genuine_rectangles_in_the_tagged_flat_collar"], 0)
