"""Derive the rescaled limit and cubic obstruction with SageMath."""

import importlib.util
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The rescaled endpoint algebra checks use SageMath")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from sage.all import QQ, PolynomialRing, matrix
from moving_flat_endpoints_v292_v295 import chord_hessian
from rough_flat_contacts_v296_v299 import FLAT_TAGS, limiting_jacobian, limiting_residual, rescaled_residual


class RoughFlatAlgebraTests(unittest.TestCase):
    def test_rescaled_parabolic_limit_and_jacobian(self):
        p = PolynomialRing(QQ, names=("a", "b", "c", "d", "r", "z0", "z1", "z2", "z3"))
        field = p.fraction_field()
        a, b, c, d, r, *z = map(field, p.gens())
        left = lambda t: (a*t*t, t+b*t*t)
        right = lambda t: (1+c*t*t, t+d*t*t)
        rhs = rescaled_residual(left, right, z, r)
        limit = []
        for item in rhs:
            self.assertEqual(item.denominator(), 1)
            limit.append(item.numerator().subs({p.gen(4): 0}))
        self.assertEqual(tuple(limit), limiting_residual(a, c, z))
        tag_point = dict(zip(p.gens()[5:], FLAT_TAGS))
        jacobian = matrix(field, [[v.derivative(w).subs(tag_point) for w in p.gens()[5:]]
                                  for v in limit])
        self.assertEqual(jacobian, matrix(field, limiting_jacobian(a, c)))
        self.assertEqual(jacobian.det(), 2*matrix(field, chord_hessian(1, 2*a, 2*c)).det())

    def test_cubic_secant_expression_is_strict_for_positive_length(self):
        p = PolynomialRing(QQ, names=("v", "r"))
        v, r = p.fraction_field().gens()
        self.assertEqual(((1+(v+r)**3)-(1+(v-r)**3))/(2*r), 3*v*v+r*r)
