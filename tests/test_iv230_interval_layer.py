"""v230 audit: soundness of the corrected Decimal interval layer iv230.

iv230 is the recovered iv101 with two fixes (see its header). These tests
check enclosures against references that do not reuse its series code:
exact squaring, a hard-coded 100-decimal value of pi, Decimal's correctly
rounded ln at higher precision, and (when installed) Arb ball arithmetic.

Set REQUIRE_SPIRAL_AUDIT_DEPS=1 to turn missing python-flint into a failure
instead of a skip (the dedicated CI job does this).
"""
import os
import random
import sys
import unittest
from decimal import Decimal as D, localcontext, ROUND_CEILING
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'related' / 'power-spiral-v99-v102'))
import iv101  # recovered original, byte-identical
import iv230  # corrected copy

try:
    import flint
except ImportError:  # pragma: no cover - environment dependent
    flint = None

REQUIRE = os.environ.get('REQUIRE_SPIRAL_AUDIT_DEPS') == '1'
# pi to 100 decimals as literals (no Decimal arithmetic: importing iv101/iv230
# sets the global Decimal precision to 40):  PI_TRUNC < pi < PI_UPPER.
PI_TRUNC = D('3.1415926535897932384626433832795028841971693993751058209749445923'
             '078164062862089986280348253421170679')
PI_UPPER = D('3.1415926535897932384626433832795028841971693993751058209749445923'
             '078164062862089986280348253421170680')


def exact_le(a, b):
    with localcontext() as c:
        c.prec = 400
        return a <= b


def sq(s):
    with localcontext() as c:
        c.prec = 400
        return s * s


def sample_inputs(seed, n):
    rng = random.Random(seed)
    xs = [D(2), D(3), D(5), D(7), D(10), D('0.3'), D(1), D(4), D(100),
          D('1e-30'), D('12345678901234567890123456789012345678.9')]
    for _ in range(n):
        xs.append(D(rng.randint(1, 10 ** 40)) / D(10) ** rng.randint(0, 45))
    return xs


class SqrtSoundness(unittest.TestCase):
    def test_iv230_sqrt_endpoints_verified_by_exact_squaring(self):
        for x in sample_inputs(230, 400):
            r = iv230.IV(x).sqrt()
            self.assertTrue(r.lo >= 0)
            self.assertTrue(exact_le(sq(r.lo), x), (x, r))
            self.assertTrue(exact_le(x, sq(r.hi)), (x, r))
            # tight: at most a few units in the last of 40 digits
            self.assertLess((r.hi - r.lo) / r.hi, D('1e-37'))

    def test_iv230_sqrt_of_interval_and_zero(self):
        r = iv230.IV(D(2), D(3)).sqrt()
        self.assertTrue(exact_le(sq(r.lo), D(2)) and exact_le(D(3), sq(r.hi)))
        z = iv230.IV(0).sqrt()
        self.assertEqual((z.lo, z.hi), (D(0), D(0)))
        with self.assertRaises(iv230.DomainError):
            iv230.IV(D(-1), D(1)).sqrt()

    def test_recovered_iv101_sqrt_bug_reproduces(self):
        """Documents WHY iv230 exists: Decimal.sqrt ignores directed rounding."""
        low_too_high = high_too_low = 0
        for x in sample_inputs(101, 200):
            r = iv101.IV(x).sqrt()
            low_too_high += not exact_le(sq(r.lo), x)
            high_too_low += not exact_le(x, sq(r.hi))
        self.assertGreater(low_too_high, 0)
        self.assertGreater(high_too_low, 0)


class TranscendentalEnclosures(unittest.TestCase):
    def test_machin_pi_encloses_pi(self):
        self.assertTrue(iv230.PI.lo <= PI_TRUNC)
        self.assertTrue(PI_UPPER <= iv230.PI.hi)
        self.assertLess(iv230.PI.hi - iv230.PI.lo, D('1e-35'))

    def test_atan_one_encloses_quarter_pi(self):
        a = iv230.atan_pt(D(1))
        with localcontext() as c:
            c.prec = 200
            self.assertTrue(4 * a.lo <= PI_TRUNC)
            self.assertTrue(PI_UPPER <= 4 * a.hi)

    def test_atanh_negative_arguments_against_correctly_rounded_ln(self):
        rng = random.Random(7)
        for _ in range(40):
            t = -D(rng.randint(1, 67 * 10 ** 6)) / D(10 ** 8)   # t in (-0.67, 0)
            enc = iv230.atanh_pt(t)
            with localcontext() as c:
                c.prec = 90
                ref = ((1 + t) / (1 - t)).ln() / 2      # error < 1e-85
                slack = D('1e-80')
                self.assertTrue(enc.lo <= ref - slack and ref + slack <= enc.hi, (t, enc, ref))

    def test_atanh_tail_rounding_hole_in_iv101(self):
        """iv101 bounded term*t*t under ROUND_CEILING; for t<0 that can undershoot."""
        rng = random.Random(2)
        undershoots = 0
        for _ in range(200):
            a = D(rng.randint(1, 10 ** 39)) / D(10) ** 60
            t = -D(rng.randint(1, 68 * 10 ** 37)) / D(10) ** 39
            with localcontext() as c:
                c.prec = 40
                c.rounding = ROUND_CEILING
                old = a * t * t                 # iv101 expression
                u = abs(t)
                new = a * u * u                 # iv230 expression
            with localcontext() as c:
                c.prec = 300
                exact = a * t * t
            undershoots += old < exact
            self.assertTrue(new >= exact)
        self.assertGreater(undershoots, 0)


@unittest.skipIf(flint is None and not REQUIRE, 'python-flint not installed')
class ArbCrossCheck(unittest.TestCase):
    """Every iv230 enclosure must contain Arb's independent ball."""

    def setUp(self):
        if flint is None:
            self.fail('REQUIRE_SPIRAL_AUDIT_DEPS=1 but python-flint is missing')
        flint.ctx.prec = 400
        self.arb = flint.arb

    def ball(self, x):
        f = Fraction(x)
        return self.arb(f.numerator) / self.arb(f.denominator)

    def contains(self, enc, ref):
        return ref >= self.ball(enc.lo) and self.ball(enc.hi) >= ref

    def test_reference_pi_digits(self):
        pi = self.arb.pi()
        self.assertTrue(pi >= self.ball(PI_TRUNC))
        self.assertTrue(self.ball(PI_UPPER) >= pi)

    def test_pi_and_sqrt(self):
        self.assertTrue(self.contains(iv230.PI, self.arb.pi()))
        for x in sample_inputs(3, 60):
            self.assertTrue(self.contains(iv230.IV(x).sqrt(), self.ball(x).sqrt()), x)

    def test_series_evaluators(self):
        rng = random.Random(230)
        def pts(lo, hi, n=25):
            return [D(rng.randint(int(lo * 10 ** 6), int(hi * 10 ** 6))) / D(10 ** 6)
                    for _ in range(n)]
        A = self.arb
        cases = [
            (iv230.exp_pt, pts(0, 1.6), lambda x: x.exp()),
            (iv230.sin_pt, pts(0, 0.8), lambda x: x.sin()),
            (iv230.cos_pt, pts(0, 0.8), lambda x: x.cos()),
            (iv230.sinc_pt, pts(0.001, 0.8), lambda x: x.sin() / x),
            (iv230.coshC_pt, pts(0.001, 1.6), lambda x: (x.cosh() - 1) / (x * x)),
            (iv230.atanh_pt, pts(-0.68, 0.68), lambda x: x.atanh()),
            (iv230.L1_pt, [p for p in pts(-0.8, 2.5) if p != 0], lambda x: x.log1p() / x),
            (iv230.Efun_pt, [p for p in pts(-1, 3) if p != 0], lambda x: -((-x).expm1()) / x),
            (iv230.atan_pt, pts(0, 5), lambda x: x.atan()),
        ]
        for fn, xs, ref in cases:
            for x in xs:
                self.assertTrue(self.contains(fn(x), ref(self.ball(x))), (fn.__name__, x))


if __name__ == '__main__':
    unittest.main()
