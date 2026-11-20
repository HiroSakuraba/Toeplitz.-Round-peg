"""Independent symbolic checks for the moving flat endpoint algebra."""

import importlib.util
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The flat-jet identity checks use SageMath")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))

from sage.all import QQ, PolynomialRing, matrix, vector
from moving_flat_endpoints_v292_v295 import chord_hessian, normalized_contact_jacobian


class FlatJetsAlgebraTests(unittest.TestCase):
    def test_normalization_matches_direct_rectangle_linearization(self):
        ring = PolynomialRing(QQ, names=("s", "k", "a", "b", "c", "d", "w", "l")).fraction_field()
        s, k, a, b, c, d, w, l = ring.gens()
        A = vector(ring, (a, b))
        B = vector(ring, (c, d))
        T = vector(ring, (s*k*w, 1))
        velocities = (k*A-T, k*A+T, k*B+T, k*B-T)
        columns = []
        for j in range(4):
            mA, dA, mB, dB = [int(i == j) for i in range(4)]
            weights = (mA-dA, mA+dA, mB+dB, mB-dB)
            q = tuple(weight*v for weight, v in zip(weights, velocities))
            U, V = (q[1]-q[0])/2, (q[2]-q[3])/2
            W = (q[2]+q[3]-q[0]-q[1])/2
            short = U-V
            long = s*W-l*vector(ring, (U[1], -U[0]))
            columns.append((short[0]/k, short[1], long[0], long[1]/k))
        direct = matrix(ring, columns).transpose()
        expected = matrix(ring, normalized_contact_jacobian(s, k, s/k, a, b, c, d, w, l))
        self.assertEqual(expected, direct)
        limiting = matrix(ring, normalized_contact_jacobian(0, k, 0, a, b, c, d, w, l))
        self.assertEqual(limiting.det(), l*l*a*c)

    def test_hessian_is_the_exact_chord_length_hessian(self):
        ring = PolynomialRing(QQ, names=("u", "v", "a", "c", "l", "j"))
        u, v, a, c, l, j = ring.gens()
        energy = ((l+j*c*v*v-j*a*u*u)**2+(v-u)**2)/2
        direct = matrix(ring, [[energy.derivative(x).derivative(y)(0, 0, a, c, l, j)
                                for y in (u, v)] for x in (u, v)])
        expected = matrix(ring, chord_hessian(l, 2*j*a, 2*j*c))
        self.assertEqual(direct, expected)
