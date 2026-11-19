"""v240: classification of inscribed squares on open log spirals, 3/50 <= a <= 9/10.

Needs python-flint (Arb). Set REQUIRE_SPIRAL_AUDIT_DEPS=1 to make a missing
dependency a failure instead of a skip; the dedicated CI job does this.
"""
import cmath
import math
import os
import random
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
try:
    import flint  # noqa: F401
    import log_spiral_classification_v240 as m
except ImportError:  # pragma: no cover - environment dependent
    m = None
try:
    import mpmath
except ImportError:  # pragma: no cover
    mpmath = None

REQUIRE = os.environ.get('REQUIRE_SPIRAL_AUDIT_DEPS') == '1'


class ElementaryLemmas(unittest.TestCase):
    """Standard-library sanity checks of the reduction lemmas (not proofs)."""

    def test_adjacent_vertices_have_radius_at_least_sqrt2_minus_1(self):
        rng = random.Random(240)
        floor = math.sqrt(2) - 1
        for _ in range(20000):
            r, th = math.sqrt(rng.random()), rng.uniform(0, 2 * math.pi)
            v1 = cmath.rect(r, th)
            for sg in (1, -1):
                v3 = 1 + sg * 1j * (v1 - 1)
                self.assertGreaterEqual(abs(v3), floor - 1e-12)

    def test_log_spiral_tangent_turns_exactly_by_parameter(self):
        a = 0.37
        arg = lambda t: cmath.phase((-a + 1j) * cmath.exp((-a + 1j) * t))
        for t in (0.1, 0.5, 1.0):
            d = (arg(t) - arg(0.0)) % (2 * math.pi)
            self.assertAlmostEqual(d, t, places=12)


@unittest.skipIf(m is None and not REQUIRE, 'python-flint not installed')
class Classification(unittest.TestCase):
    def setUp(self):
        if m is None:
            self.fail('REQUIRE_SPIRAL_AUDIT_DEPS=1 but python-flint is missing')
        m.ctx()

    def test_constants_are_certified(self):
        arb = m.arb
        self.assertTrue((1 + arb(2).sqrt()).log() < m.ball(m.T1C))
        self.assertTrue(m.ball(m.TAU) < arb.pi() / 2)
        self.assertGreater(m.SQRT2_UP ** 2, 2)

    def test_full_classification(self):
        res = m.classify()
        self.assertEqual(len(res['bands']), 168)
        self.assertEqual(res['bands'][0][0], F(3, 50))
        self.assertEqual(res['bands'][-1][1], F(9, 10))
        for (a0, b0, *_), (a1, *__) in zip(res['bands'], res['bands'][1:]):
            self.assertEqual(b0, a1)
        self.assertEqual(set(res['kboxes']), {1, -1})
        lo, hi = res['kboxes'][1][2], res['kboxes'][1][3]
        self.assertTrue(F('0.10027465') < lo < hi < F('0.10027466'))

    def test_band_with_the_square_fails_without_krawczyk_box(self):
        r = m.certify_band(F(1, 10), F(21, 200))
        self.assertTrue(r['stuck'])
        for sg, t, a in r['stuck']:
            self.assertAlmostEqual(a, 0.1002746561, places=8)

    def test_krawczyk_controls(self):
        t, a = m.RESONANCE[1]
        self.assertTrue(m.krawczyk(1, F(t), F(a)))
        self.assertFalse(m.krawczyk(1, F(t) + 20 * m.K_RADIUS_T, F(a)))
        self.assertFalse(m.krawczyk(1, F(t), F(a) + 50 * m.K_RADIUS_A))
        self.assertFalse(m.krawczyk(-1, F(t), F(a)))

    def test_independent_overlap_with_v239(self):
        for x in range(90, 150, 10):
            self.assertFalse(m.certify_band(F(x, 100), F(x + 10, 100))['stuck'])

    def test_range_stops_where_search_stalls(self):
        r = m.certify_band(F(55, 1000), F(56, 1000), depth_max=40)
        self.assertTrue(r['stuck'])

    def test_lemma_L_bottom_band(self):
        eps, M, _ = m.lemma_L(F(3, 50), F(61, 1000))
        self.assertGreater(eps, F(1, 10))
        self.assertLess(m.SQRT2_UP * (1 + F(61, 1000) ** 2 / 2) * eps, M)


@unittest.skipIf((m is None or mpmath is None) and not REQUIRE, 'flint/mpmath missing')
class SquareAtAStar(unittest.TestCase):
    def test_the_square_is_genuine_and_both_labels_agree(self):
        if m is None or mpmath is None:
            self.fail('REQUIRE_SPIRAL_AUDIT_DEPS=1 but a dependency is missing')
        mp = mpmath
        mp.mp.dps = 40
        t1, a = (mp.mpf(x) for x in m.RESONANCE[1])
        z = lambda t: mp.exp((-a + 1j) * t)
        v0, v1 = mp.mpc(1), z(t1)
        v3 = 1 + 1j * (v1 - 1)
        v2 = v1 + 1j * (v1 - 1)
        params = []
        for w in (v2, v3):
            t = -mp.log(abs(w)) / a
            self.assertLess(abs(z(t) - w), mp.mpf('1e-18'))
            params.append(t)
        sides = [abs(v1 - v0), abs(v2 - v1), abs(v3 - v2), abs(v0 - v3)]
        self.assertLess(max(sides) - min(sides), mp.mpf('1e-30'))
        self.assertTrue(all(abs(w) < 1 for w in (v1, v2, v3)))
        # sigma=-1 labelling uses the other neighbour of v0: its t1 is v3's parameter
        self.assertLess(abs(params[1] - mp.mpf(m.RESONANCE[-1][0])), mp.mpf('1e-18'))


if __name__ == '__main__':
    unittest.main()
