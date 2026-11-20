"""Independent symbolic and algebraic checks for the v285 support motion."""

import importlib.util
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The v285 symbolic identity checks use SageMath")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from sage.all import AA, QQ, PolynomialRing
from moving_spiral_support_v285_v286 import (
    ARM_HEIGHTS, CAP_DEPTHS, coefficients, comparisons, polynomial_value, walking_vertices,
)
from relative_witness_refinement_v282_v283 import check_contacts, diagonal_dot, diamond_vertices, difference, dot


class MovingSupportIdentitiesTests(unittest.TestCase):
    def test_all_comparisons_are_polynomial_identities(self):
        ring = PolynomialRing(QQ, names=("h", "k", "a", "b")).fraction_field()
        h, k, a, b = ring.gens()
        G = coefficients(h, k, (a, b))
        self.assertTrue(all(polynomial_value(G, r) == value for r, value in comparisons(h, k, (a, b))))
        L = 7+7*a+8*b
        B = 1+b
        self.assertEqual(tuple(8*c for c in coefficients(ring(1), ring(QQ(7)/8), (a, b))),
                         (L*(B-a), 16*a-14*B, ring(-2)))

    def test_algebraic_rectangles_on_moved_supports(self):
        depths = tuple(QQ(v) for v in CAP_DEPTHS)
        heights = tuple(QQ(v) for v in ARM_HEIGHTS)
        for h, k in ((depths[0], heights[1]), (depths[1], heights[0]), (depths[1], heights[1]),
                     (sum(depths)/2, sum(heights)/2)):
            fixed = (((-1, 0), (0, -k)), ((-1, 0), (0, -h)), ((0, -h), (1, 0)))
            vertices = tuple(tuple(QQ(v) for v in p) for p in walking_vertices(k))
            sample = set(vertices)
            sample.update(tuple(QQ(v) for v in p) for p in diamond_vertices())
            sample.update(tuple((a+b)/2 for a, b in zip(p, q)) for p, q in zip(vertices, vertices[1:]))
            for q in sample:
                a, b = q
                c0, c1, c2 = coefficients(h, k, q)
                r = (-c1-AA(c1*c1-4*c0*c2).sqrt())/(2*c2)
                ell = ((h+k)*r-(k*(1+a)+b))/(h-k)
                corners = ((a-ell-r, b+h*(ell-r)), (-ell, h*(ell-1)), (r, h*(r-1)), q)
                check_contacts(fixed, corners)
                v = difference(corners[2], corners[1])
                w = difference(corners[3], corners[2])
                self.assertEqual(dot(v, w), 0)
                if q == vertices[0]:
                    self.assertEqual(diagonal_dot(corners), 4)
                elif q == vertices[-1]:
                    self.assertEqual(diagonal_dot(corners), -(h-k)**2)
                else:
                    self.assertGreater(v[0]*w[1]-v[1]*w[0], 0)
                d = difference(corners[2], corners[0])
                self.assertGreater(dot(d, d)/4, QQ(1)/625)
