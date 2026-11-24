from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from path_sweep_v373 import large_graceful_parity, locate_event, nearest_to_cutoff, twisted_drawing  # noqa: E402
from polygon_orbit_census_v223 import pt, validate_jordan_polygon  # noqa: E402
from wild_point_parity_v370 import orientation, sector, sector_counts  # noqa: E402

S = Q(1, 20)


class TwistPath(unittest.TestCase):
    def test_endpoints(self):
        for u, parity in ((Q(0), 1), (Q(1, 200), 1), (Q(1, 20), 0)):
            V = twisted_drawing(3, u)
            validate_jordan_polygon([pt(*p) for p in V])
            self.assertEqual(large_graceful_parity(V, S)[0], parity)

    def test_the_flip_is_a_square_crossing_the_cutoff(self):
        (lo, hi), _, small_hi = locate_event(3, S, Q(1, 200), Q(1, 20), steps=14)
        self.assertEqual((lo, hi), (Q(38029, 3276800), Q(19019, 1638400)))
        entering = nearest_to_cutoff(small_hi, S)
        self.assertLess(entering[0], S * S)
        self.assertGreater(entering[0], (S * Q(9999, 10000)) ** 2)
        V_lo, V_hi = twisted_drawing(3, lo), twisted_drawing(3, hi)
        _, degenerate, wider = sector_counts(V_lo, Q(51, 1000))
        self.assertEqual(degenerate, [])
        above = [e for e in wider.values() if e[0] >= S * S]
        shift = lambda e: max(abs(p[0] - q[0]) + abs(p[1] - q[1]) for p, q in zip(sorted(e[1]), sorted(entering[1])))
        match = min(above, key=shift)
        self.assertLess(shift(match), Q(1, 10**5))
        self.assertLess(match[0], (S * Q(10001, 10000)) ** 2)
        self.assertEqual(sector(match, len(V_lo), orientation(V_lo)), sector(entering, len(V_hi), orientation(V_hi)))
        self.assertEqual(sector(entering, len(V_hi), orientation(V_hi)), 'reversed')


if __name__ == '__main__':
    unittest.main()
