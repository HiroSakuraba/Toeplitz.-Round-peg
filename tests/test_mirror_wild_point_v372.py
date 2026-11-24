from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from mirror_wild_point_v372 import in_cone, mirror_arm, mirror_drawing, reflect, symmetric_chords  # noqa: E402
from polygon_orbit_census_v223 import pt, validate_jordan_polygon  # noqa: E402
from wild_point_parity_v370 import orientation, sector, sector_counts  # noqa: E402


def by_sector(V, s):
    counts, degenerate, sq = sector_counts(V, s)
    o = orientation(V)
    out = {'aligned': set(), 'reversed': set(), 'crossed': set()}
    for e in sq.values():
        out[sector(e, len(V), o)].add(frozenset(e[1]))
    return out, degenerate


def control(zig):
    arm = [(Q(0), Q(0))] + [tuple(map(Q, z)) for z in zig] + [(Q(2), Q(0))]
    return arm, [(-x, y) for x, y in reversed(arm)][:-1] + arm + [(Q(2), Q(1)), (Q(-2), Q(1))]


class ConeArms(unittest.TestCase):
    def test_curled_arms_pair_without_fixed_squares(self):
        for depth, s in ((1, Q(1, 16)), (2, Q(1, 16))):
            arm = mirror_arm(depth, turns=1, per_turn=6)
            self.assertTrue(in_cone(arm))
            self.assertEqual(symmetric_chords(arm), [])
            V = mirror_drawing(depth, turns=1, per_turn=6)
            validate_jordan_polygon([pt(*p) for p in V])
            sec, degenerate = by_sector(V, s)
            self.assertEqual(degenerate, [])
            for name in ('aligned', 'reversed'):
                self.assertEqual({reflect(q) for q in sec[name]}, sec[name])
                self.assertFalse(any(reflect(q) == q for q in sec[name]))
                self.assertEqual(len(sec[name]) % 2, 0)

    def test_large_aligned_count_is_odd(self):
        V = mirror_drawing(1, turns=1, per_turn=4)
        total, degenerate, _ = sector_counts(V, Q(100))
        self.assertEqual(degenerate, [])
        self.assertEqual((total['aligned'] % 2, total['reversed'] % 2), (1, 0))
        small, _ = by_sector(V, Q(1, 16))
        self.assertEqual((total['aligned'] - len(small['aligned'])) % 2, 1)


class DegenerateChords(unittest.TestCase):
    def test_overlapping_edges_raise_only_inside_the_cutoff(self):
        q = lambda pts: [tuple(map(Q, p)) for p in pts]
        with self.assertRaises(ArithmeticError):
            symmetric_chords(q([(1, 0), (1, 2)]))
        with self.assertRaises(ArithmeticError):
            symmetric_chords(q([(0, 0), (1, 0), (1, 3)]))
        self.assertEqual(symmetric_chords(q([(0, 0), (1, 0), (1, 3)]), max_u=Q(1, 2)), [])
        with self.assertRaises(ArithmeticError):
            symmetric_chords(q([(0, 0), (1, 0), (1, 3)]), max_u=Q(2))


class OutsideTheCone(unittest.TestCase):
    def test_one_symmetric_chord_makes_one_small_sector_odd(self):
        cases = [([('23/100', '-1/20'), ('19/50', '-9/50'), ('19/100', '3/50'), ('2/25', '7/25')], 'aligned'),
                 ([('17/100', '17/50'), ('1/4', '19/50'), ('2/25', '-11/100'), ('11/100', '-1/20')], 'reversed')]
        for zig, odd in cases:
            arm, V = control(zig)
            self.assertFalse(in_cone(arm))
            validate_jordan_polygon([pt(*p) for p in V])
            chords = symmetric_chords(arm, max_u=Q(1, 8))
            self.assertEqual(len(chords), 1)
            sec, degenerate = by_sector(V, Q(1, 4))
            self.assertEqual(degenerate, [])
            self.assertEqual(len(sec[odd]), 1)
            (square,) = sec[odd]
            self.assertEqual(reflect(square), square)          # the odd square is the symmetric one
            u, P = chords[0]
            self.assertIn(P, square)
            self.assertIn((P[0], P[1] + 2 * u), square)


if __name__ == '__main__':
    unittest.main()
