"""Independent algebraic checks at both patch joins and throughout their enclosures."""

import importlib.util
from itertools import product
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The rough-arc algebraic lift checks use SageMath")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from sage.all import AA, QQ, PolynomialRing
from protected_rectangle_band_v284 import concave_edges, lerp
from relative_witness_refinement_v282_v283 import difference, dot, on_segment
import rough_arc_continuation_v288_v291 as rough


R = PolynomialRing(QQ, "z")


def patch_root(edge, q):
    polynomial = R(list(rough.contact_polynomial(edge, q)))
    lo, hi = (0, QQ(1)/2) if edge == 0 else (QQ(3)/20, QQ(1)/2)
    roots = tuple(r for r in polynomial.roots(AA, multiplicities=False) if lo < r < hi)
    if len(roots) != 1:
        raise AssertionError("patch has an unexpected algebraic root count")
    return roots[0]


class RoughArcAlgebraTests(unittest.TestCase):
    def test_all_four_joins_agree_with_the_reference_path(self):
        edges = tuple(tuple(tuple(QQ(v) for v in q) for q in e) for e in concave_edges())
        for edge, window in ((0, rough.BOTTOM_WINDOW), (2, rough.OPPOSITE_WINDOW)):
            for t in map(QQ, window):
                with self.subTest(edge=edge, t=t):
                    H = R(list(rough.DEN))
                    equation = R(list(rough.NUMS[edge]))-t*H
                    roots = [x for x in equation.roots(AA, multiplicities=False) if 0 < x < 1]
                    self.assertEqual(len(roots), 1)
                    x = roots[0]
                    parameters = tuple(R(list(n))(x)/H(x) for n in rough.NUMS)
                    expected = tuple(lerp(e, u) for e, u in zip(edges, parameters))
                    q = lerp(edges[edge], t)
                    self.assertEqual(rough.lifted_corners(edge, q, patch_root(edge, q)), expected)

    def test_exact_box_vertices_and_center_are_ccw_and_have_clearance(self):
        edges = tuple(tuple(tuple(QQ(v) for v in q) for q in e) for e in concave_edges())
        for edge, box in ((0, rough.COIL_BOX), (2, rough.OPPOSITE_BOX)):
            corners_of_box = list(product(*(tuple(map(QQ, v)) for v in box)))
            corners_of_box.append(tuple(sum(point[i] for point in corners_of_box)/4 for i in range(2)))
            for q in corners_of_box:
                with self.subTest(edge=edge, q=q):
                    square = rough.lifted_corners(edge, q, patch_root(edge, q))
                    self.assertEqual(square[edge], q)
                    self.assertTrue(all(on_segment(p, *e) for i, (p, e) in enumerate(zip(square, edges))
                                        if i != edge))
                    u, v = difference(square[1], square[0]), difference(square[3], square[0])
                    self.assertEqual(dot(u, v), 0)
                    self.assertGreater(u[0]*v[1]-u[1]*v[0], 0)
                    diag = difference(square[2], square[0])
                    bound = QQ(1) if edge == 0 else QQ(117)/200
                    self.assertGreaterEqual(dot(diag, diag)/4, bound**2)
