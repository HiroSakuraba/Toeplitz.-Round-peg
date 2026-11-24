from fractions import Fraction as Q
from math import sqrt
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from asymmetric_spiral_v374 import (  # noqa: E402
    assert_jordan, asymmetric_drawing, comb_radius, default_phi, large_graceful_parity, radius_over_side)
from polygon_orbit_census_v223 import pt, validate_jordan_polygon  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
S = Q(1, 20)


def max_center_radius(squares):
    return max(radius_over_side(e) * sqrt(float(e[0])) for e in squares.values())


class GapPath(unittest.TestCase):
    def test_drawings_are_jordan_and_arms_keep_their_radii(self):
        for lam in (Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1)):
            V = asymmetric_drawing(4, lam)
            validate_jordan_polygon([pt(*p) for p in V])
            n = (len(V) - 3) // 2
            A, B = V[:n], V[n:2 * n]
            # arm two is arm one turned vertex by vertex: the radii agree exactly
            self.assertEqual([x * x + y * y for x, y in A], [x * x + y * y for x, y in reversed(B)])

    def test_deep_drawings_are_rejected_when_chords_cross(self):
        V = asymmetric_drawing(75)
        with self.assertRaises(ValueError):
            assert_jordan(V)
        with self.assertRaises(ValueError):
            large_graceful_parity(V, S)
        assert_jordan(asymmetric_drawing(75, Q(0)))

    def test_parity_along_the_gap_path(self):
        expected = {Q(0): (1, 29, 29), Q(1, 2): (0, 32, 27), Q(1): (1, 35, 29)}
        for lam, row in expected.items():
            par, counts, _ = large_graceful_parity(asymmetric_drawing(4, lam), S)
            self.assertEqual((par, counts['aligned'], counts['reversed']), row)
            self.assertEqual(counts['crossed'], 0)

    def test_first_event_is_an_aligned_square_on_the_outer_turn(self):
        lo, hi = Q(2567, 8192), Q(5135, 16384)
        p_lo, c_lo, _ = large_graceful_parity(asymmetric_drawing(4, lo), S)
        p_hi, c_hi, sq_hi = large_graceful_parity(asymmetric_drawing(4, hi), S)
        self.assertNotEqual(p_lo, p_hi)
        self.assertEqual(c_hi['aligned'], c_lo['aligned'] + 1)
        self.assertEqual(c_hi['reversed'], c_lo['reversed'])
        top = max(sq_hi.values(), key=lambda e: e[0])
        self.assertGreater(top[0], (S * Q(99999, 100000)) ** 2)
        self.assertAlmostEqual(radius_over_side(top), 6.4, delta=0.1)


class CombRadius(unittest.TestCase):
    def test_small_squares_stop_at_the_comb_radius(self):
        for lam in (Q(0), Q(1)):
            _, _, sq = large_graceful_parity(asymmetric_drawing(4, lam), S)
            r = comb_radius(0.05, phi=default_phi, lam=float(lam))
            self.assertLess(abs(max_center_radius(sq) / r - 1), 0.02)

    def test_recorded_scan(self):
        rows = json.loads((ROOT / 'results/v374/comb_scan.json').read_text())['comb_scan']
        self.assertEqual(len(rows), 8)
        for r in rows:
            self.assertLess(abs(r['max_center_radius'] / r['comb_radius'] - 1), 0.02)
            self.assertLess(r['inner'], r['comb_radius'])
            self.assertEqual(r['crossed'], 0)
            if r['lam'] == '0':
                self.assertEqual(r['aligned'], r['reversed'])
        ratios = [r['max_center_over_side'] for r in rows if r['lam'] == '0']
        self.assertEqual(ratios, sorted(ratios))


if __name__ == '__main__':
    unittest.main()
