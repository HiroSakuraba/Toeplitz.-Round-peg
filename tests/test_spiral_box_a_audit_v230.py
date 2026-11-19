"""v230 audit: checks that can actually fail, behind the v230 double-spiral note.

- Standard library: complex-step check of the exact y dx identity (with a
  negative control proving the check has teeth) and the p<2 / p>=2 behaviour
  of the non-oscillating primitive.
- python-flint: Arb certification of the three Theorem 2 inequalities on ALL
  of Box A (every p>=2), plus a corner check that the hand-proof constants are
  not vacuous.
- mpmath: numerical continuation of the actual Box-A square (not a proof),
  including the Asano-Ike theta>=e truncation caveat.

Set REQUIRE_SPIRAL_AUDIT_DEPS=1 to make missing python-flint/mpmath a failure
instead of a skip (the dedicated CI job does this).
"""
import math
import os
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
import spiral_box_a_audit_v230 as audit

try:
    import flint
except ImportError:  # pragma: no cover
    flint = None
try:
    import mpmath
except ImportError:  # pragma: no cover
    mpmath = None

REQUIRE = os.environ.get('REQUIRE_SPIRAL_AUDIT_DEPS') == '1'
GRID = [(t, p) for p in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0)
        for t in (0.05, 0.1, 0.3, 0.5, 0.9)]


class YdxIdentity(unittest.TestCase):
    def test_identity_holds_by_complex_step(self):
        # Residual is limited only by float argument reduction of cos(t^-p)
        # for large t^-p (about 1e-10 here), far below the control below.
        for t, p in GRID:
            self.assertLess(audit.y_dx_identity_residual(t, p), 1e-8, (t, p))

    def test_negative_control_detects_wrong_coefficient(self):
        import cmath
        def wrong(t, p):   # (p/3) instead of (p/2): must be caught
            x = lambda s: s * cmath.cos(s ** (-p))
            Fn = lambda s: s * s / 4 * cmath.sin(2 * s ** (-p))
            lhs = t * math.sin(t ** (-p)) * audit._cstep(x, t)
            rhs = audit._cstep(Fn, t) + p / 3 * t ** (1 - p)
            return abs(lhs - rhs) / max(1.0, abs(lhs), abs(rhs), p * t ** (1 - p))
        self.assertGreater(max(wrong(t, p) for t, p in GRID), 1e-3)

    def test_primitive_bounded_below_two_unbounded_from_two(self):
        tails = lambda p: [audit.liouville_primitive_tail(p, 10.0 ** -k) for k in (4, 8, 16)]
        for p in (0.5, 1.0, 1.5, 1.9):
            a = tails(p)
            self.assertLess(a[-1] - a[0], 1e-3 * max(1, a[-1]) + p / (2 * (2 - p)) * 1e-4 ** (2 - p))
            self.assertLessEqual(a[-1], p / (2 * (2 - p)) + 1e-12)
        for p in (2.0, 2.5, 3.0):
            a = tails(p)
            self.assertTrue(a[0] < a[1] < a[2])
            self.assertGreater(a[2], 30)


@unittest.skipIf(flint is None and not REQUIRE, 'python-flint not installed')
class BoxACertificate(unittest.TestCase):
    def setUp(self):
        if flint is None:
            self.fail('REQUIRE_SPIRAL_AUDIT_DEPS=1 but python-flint is missing')

    def test_three_inequalities_certified_on_whole_box(self):
        out = audit.certify_box_a()
        self.assertEqual(set(out), {'I1', 'I2', 'I3'})
        for name, (cells, margin) in out.items():
            self.assertGreaterEqual(cells, 1)
            self.assertGreater(margin, 0, name)

    def test_inequalities_are_not_vacuous_at_worst_corner(self):
        """The I1 threshold is close to tight: a 1/19 threshold must FAIL somewhere."""
        flint.ctx.prec = 200
        arb = flint.arb
        c, d, q = arb(1), arb(9) / 20, arb(1) / 3
        val = (1 - (-(c * (1 - q))).exp()) / (3 * arb.pi() - d * q)
        self.assertTrue(val > arb(1) / 25)
        self.assertTrue(val < arb(1) / 19)


@unittest.skipIf(mpmath is None and not REQUIRE, 'mpmath not installed')
class BoxASquareTracking(unittest.TestCase):
    """NUMERICAL regression, not a proof: the square the theorem talks about."""

    @classmethod
    def setUpClass(cls):
        if mpmath is None:
            raise AssertionError('REQUIRE_SPIRAL_AUDIT_DEPS=1 but mpmath is missing')
        cls.rows = audit.track_family([2, 2.5, 3, 5, 10, 30, 100])

    def test_square_satisfies_theorem_bounds(self):
        for r in self.rows:
            p = float(r['p'])
            self.assertLess(r['residual'], 1e-40, p)
            self.assertLess(r['side_spread'], 1e-40, p)
            self.assertTrue(r['all_on_first_arm'], p)
            self.assertTrue(r['in_box_a'], p)
            self.assertGreater(r['min_radius'], 0.2, p)
            self.assertGreater(r['min_radius'], 25 ** (-1 / p), p)
            self.assertLess(r['max_radius'], 1, p)
            self.assertGreater(r['side'], 0.25, p)
            self.assertGreater(r['side'], r['min_radius'], p)

    def test_asano_ike_theta_e_cut_caveat(self):
        """Asano-Ike cut their spiral at theta>=e; the Box-A square leaves that piece for large p."""
        by_p = {float(r['p']): r for r in self.rows}
        for p in (2, 2.5, 3, 5, 10):
            self.assertGreater(by_p[p]['outer_theta'], math.e, p)
        for p in (30, 100):
            self.assertLess(by_p[p]['outer_theta'], math.e, p)
        cross = float(audit.asano_ike_cut_crossover())
        self.assertTrue(23 < cross < 25, cross)


if __name__ == '__main__':
    unittest.main()
