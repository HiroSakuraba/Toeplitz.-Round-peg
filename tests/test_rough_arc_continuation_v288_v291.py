"""Replay the rough-arc certificates and resolve the independent tent fixtures."""

from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from moving_polygon_events_v224 import evaluate
from polygon_orbit_census_v223 import _solve, pt, validate_jordan_polygon
from protected_rectangle_band_v284 import affine_matrix, concave_edges, lerp, norm_inf
from relative_witness_refinement_v282_v283 import diagonal_dot, difference, dot
import rough_arc_continuation_v288_v291 as rough


def tent_solutions(x):
    """Independently solve every tagged piecewise-linear cell at fixed x."""
    reference = concave_edges()
    patch = rough.tent_edges()
    mid = rough.baseline_parameters(Q(1, 2))
    cells = []
    for i, (a, b) in enumerate(reference):
        if i not in patch:
            cells.append(((Q(0), Q(1), (a, b)),))
            continue
        l, m, h = patch[i]
        lo, hi = mid[i]-rough.HALFWIDTH, mid[i]+rough.HALFWIDTH
        cells.append(((Q(0), lo, (a, l)), (lo, mid[i], (l, m)),
                      (mid[i], hi, (m, h)), (hi, Q(1), (h, b))))
    baseline = rough.baseline_parameters(x)
    result = set()
    for group in product(*cells):
        segment_edges = tuple(e for lo, hi, e in group)
        A, constant = affine_matrix(segment_edges, x)
        p = _solve(A, tuple(-v for v in constant))
        if p is None or not all(0 <= v <= 1 for v in p):
            continue
        actual = tuple(lo+(hi-lo)*t for (lo, hi, e), t in zip(group, p))
        if norm_inf(difference(actual, baseline)) <= rough.HALFWIDTH:
            result.add(tuple(lerp(e, t) for e, t in zip(segment_edges, p)))
    return result


class RoughArcContinuationTests(unittest.TestCase):
    def test_exact_ledger(self):
        self.assertEqual(rough.report(), json.loads(rough.ARTIFACT.read_text(encoding="utf-8")))

    def test_two_lifts_match_baseline_in_their_windows(self):
        for edge, x in ((0, Q(1, 8)), (2, Q(7, 9))):
            with self.subTest(edge=edge):
                preimage = rough.baseline_parameters(x)[edge]
                lo, hi = rough.BOTTOM_WINDOW if edge == 0 else rough.OPPOSITE_WINDOW
                self.assertTrue(lo < preimage < hi)
                corners = rough.baseline_corners(x)
                root = corners[2][0] if edge == 0 else corners[0][0]
                self.assertEqual(evaluate(rough.contact_polynomial(edge, corners[edge]), root), 0)
                self.assertEqual(rough.lifted_corners(edge, corners[edge], root), corners)

    def test_tents_are_a_simple_polygon(self):
        e = concave_edges()
        p = rough.tent_edges()
        polygon = (e[0][0], *p[0], e[1][0], e[2][0], *p[2], e[3][0])
        self.assertEqual(len(validate_jordan_polygon(tuple(pt(*q) for q in polygon))), 10)

    def test_tent_tube_meets_exact_join_rectangles_and_opposite_signs(self):
        lo, hi = rough.ANGLE_WINDOW
        for x, direction in ((lo, 1), (Q(9, 20), 1), (Q(11, 20), -1), (hi, -1)):
            with self.subTest(x=x):
                found = tent_solutions(x)
                self.assertEqual(len(found), 1)
                corners = next(iter(found))
                if x in (lo, hi):
                    self.assertEqual(corners, rough.baseline_corners(x))
                self.assertGreater(direction*diagonal_dot(corners), 0)
                diag = difference(corners[2], corners[0])
                self.assertGreaterEqual(dot(diag, diag)/4, (Q(1, 2)-rough.EPSILON)**2)

    def test_anchored_coil_error_forces_the_reference_square(self):
        audit = rough.simultaneous_coil_audit()
        self.assertEqual(audit["physical_pointwise_error_over_source_distance_upper"], "8/255")
        self.assertLess(Q(audit["preconditioned_square_fiber_error_factor_upper"]), 1)
        self.assertTrue(audit["same_square_meets_both_coil_accumulation_points"])

    def test_coupled_square_touches_both_perturbed_segments(self):
        found = tent_solutions(Q(1, 2))
        self.assertEqual(len(found), 1)
        corners = next(iter(found))
        self.assertEqual(diagonal_dot(corners), 0)
        self.assertLess(corners[0][1], 0)
        self.assertGreater(2*corners[2][0]+corners[2][1], 3)
