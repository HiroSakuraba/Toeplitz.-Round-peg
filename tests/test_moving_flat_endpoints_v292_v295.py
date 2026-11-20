"""Exact numerical fixtures and replay for the moving-flat theorems."""

from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

import moving_flat_endpoints_v292_v295 as flats


def parabola(offset, quadratic, center):
    return (lambda t: offset+quadratic*(t-center)**2,
            lambda t: 2*quadratic*(t-center))


class MovingFlatEndpointTests(unittest.TestCase):
    def test_replay(self):
        self.assertEqual(flats.report(), json.loads(flats.ARTIFACT.read_text(encoding="utf-8")))

    def test_both_flat_jet_signs_and_transition_limits(self):
        result = flats.reference_endpoint_audit()
        self.assertEqual([result[str(x)]["normalized_transition_limit_determinant"] for x in (0, 1)],
                         ["-288", "4/5"])
        self.assertEqual([result[str(x)]["hessian_determinant_sign_for_curvature_parameter_at_least_one"]
                          for x in (0, 1)], [1, -1])

    def test_chord_hessian_and_angular_center_jacobian_agree(self):
        for length, left, right in ((Q(5), Q(48, 5), -Q(24, 5)), (Q(2), Q(2), Q(4, 5))):
            hessian = flats.chord_hessian(length, left, right)
            centers = flats.binormal_jacobian(length, left, right)
            self.assertEqual(flats.rational_det(hessian), -length*flats.rational_det(centers))
            self.assertNotEqual(flats.rational_det(hessian), 0)

    def test_moving_parabolic_binormal_and_join_rectangles(self):
        shift = Q(1, 100)
        f = parabola(Q(1, 200), Q(24, 5), shift)
        g = parabola(5+Q(1, 300), -Q(12, 5), shift)
        self.assertNotEqual(flats.fixed_aspect_collar_residual(f, g, Q(0), Q(0), Q(0), Q(0)), (0, 0, 0))
        for r in (Q(0), Q(1, 1000), Q(1, 100)):
            long = g[0](shift+r) - f[0](shift+r)
            ratio = 2*r/long
            self.assertGreater(long, 4)
            self.assertEqual(flats.fixed_aspect_collar_residual(f, g, shift, shift, r, ratio),
                             (0, 0, 0))

    def test_tilt_removes_tagged_flat_and_fold_has_degree_zero(self):
        self.assertEqual(flats.parallel_tilt_control()["binormals_between_the_tilted_straight_arcs"], 0)
        self.assertEqual(flats.fold_control()["positive_time_local_degrees"], [-1, 1])
        self.assertEqual(flats.fold_control()["transverse_fold_degree"], 0)
