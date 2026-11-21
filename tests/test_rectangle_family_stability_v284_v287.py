"""Exact regression checks for the near-square and moved-support certificates."""

from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

import moving_spiral_support_v285_v286 as moving
import protected_rectangle_band_v284 as band
from polygon_orbit_census_v223 import _solve
from relative_witness_refinement_v282_v283 import diagonal_dot, difference, dot


class RectangleFamilyStabilityTests(unittest.TestCase):
    def test_frozen_certificates(self):
        for module in (band, moving):
            with self.subTest(module=module.__name__):
                self.assertEqual(module.report(), json.loads(module.ARTIFACT.read_text(encoding="utf-8")))

    def test_degree_guard_rejects_excessive_perturbation(self):
        with self.assertRaisesRegex(ArithmeticError, "homotopy may hit a side"):
            band.certify_band(band.diamond_edges(), Q(1, 1000), Q(1, 20000), Q(1, 100), Q(9, 10))

    def test_independently_bent_arcs_have_rectangles_at_three_angles(self):
        record = band.report()["diamond_spiral"]
        windows = []
        halfwidth = Q(record["contact_box_halfwidth"])
        bends = ((1, 1), (1, -1), (-1, 1), (-1, -1))
        for edge, center, bend in zip(band.diamond_edges(), record["square_edge_parameters"], bends):
            center = Q(center)
            lo = band.lerp(edge, center-halfwidth)
            hi = band.lerp(edge, center+halfwidth)
            kink = tuple(v+Q(s, 100000) for v, s in zip(band.lerp(edge, center), bend))
            windows.append(((lo, kink), (kink, hi)))
        low, high = map(Q, record["angle_coordinate_interval"])
        for x, sign in ((low, 1), (Q(1, 2), 0), (high, -1)):
            with self.subTest(x=x):
                found = []
                for edges in product(*windows):
                    matrix, constant = band.affine_matrix(edges, x)
                    params = _solve(matrix, tuple(-v for v in constant))
                    if params is None or any(t < 0 or t > 1 for t in params):
                        continue
                    corners = tuple(band.lerp(edge, t) for edge, t in zip(edges, params))
                    self.assertEqual(band.rectangle_residual(corners, x), (0,)*4)
                    diag = difference(corners[2], corners[0])
                    self.assertGreater(dot(diag, diag)/4, Q(81, 100))
                    value = diagonal_dot(corners)
                    self.assertEqual(int(value > 0)-int(value < 0), sign)
                    found.append(corners)
                self.assertTrue(found)

    def test_flat_collar_derivatives(self):
        record = moving.rounding_audit()
        self.assertEqual(record["horizontal_flat"]["inward_midpoint_derivatives"], ["15", "1"])
        self.assertEqual(record["vertical_flat"]["inward_midpoint_derivatives"], ["7/8", "7/8"])

    def test_finite_coil_placement(self):
        record = band.report()["finite_coil_in_concave_polygon"]
        self.assertEqual(record["second_edge_height_gap"], "1/3")
        self.assertGreater(Q(record["gap_to_first_tagged_window"]), 1)
