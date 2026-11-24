from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import census_polygon, pt, validate_jordan_polygon  # noqa: E402
from wild_point_parity_v370 import drawing, orientation, sector, sector_counts  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BIG = Q(100)


def neg(points):
    return frozenset((-x, -y) for x, y in points)


class PolygonSectorDegrees(unittest.TestCase):
    """Lemma 370.1: every Jordan polygon has aligned degree 1 and reversed degree 0 (mod 2)."""

    def test_small_polygons_agree_with_the_census(self):
        shapes = [[(0, 0), (7, 1), (3, 6)], [(0, 0), (9, 1), (2, 2)],
                  [(0, 0), (5, 1), (6, 5), (1, 4)], [(0, 0), (6, 1), (2, 2), (1, 6)]]
        for V in shapes:
            V = [tuple(map(Q, p)) for p in V]
            counts, degenerate, squares = sector_counts(V, BIG)
            self.assertEqual(degenerate, [])
            census = census_polygon([pt(*p) for p in V])
            self.assertEqual({frozenset(e[1]) for e in squares.values()},
                             {frozenset((q.x, q.y) for q in s.vertices) for s in census.squares})
            self.assertEqual((counts['aligned'] % 2, counts['reversed'] % 2), (1, 0))

    def test_spiral_drawings(self):
        for h, per in ((2, 8), (2, 12), (3, 12)):
            V, _ = drawing(0.5, 3.0, h, per)
            validate_jordan_polygon([pt(*p) for p in V])
            counts, degenerate, _ = sector_counts(V, BIG)
            self.assertEqual(degenerate, [])
            self.assertEqual((counts['aligned'] % 2, counts['reversed'] % 2), (1, 0))


class HalfTurnPairing(unittest.TestCase):
    """The half-turn pairs small aligned squares with small reversed squares (Theorem 370.2)."""

    def test_exact_pairing_and_odd_large_graceful_count(self):
        for h, per, s in ((2, 12, Q(1, 10)), (3, 12, Q(1, 10)), (4, 24, Q(1, 20))):
            V, _ = drawing(0.5, 3.0, h, per)
            o = orientation(V)
            _, _, small = sector_counts(V, s)
            by_sector = {'aligned': set(), 'reversed': set(), 'crossed': set()}
            for e in small.values():
                by_sector[sector(e, len(V), o)].add(frozenset(e[1]))
            self.assertEqual({neg(q) for q in by_sector['aligned']}, by_sector['reversed'])
            self.assertEqual({neg(q) for q in by_sector['crossed']}, by_sector['crossed'])
            if h <= 3:
                total, degenerate, _ = sector_counts(V, BIG)
                self.assertEqual(degenerate, [])
                large = total['aligned'] + total['reversed'] - 2 * len(by_sector['aligned'])
                self.assertEqual(large % 2, 1)

    def test_an_asymmetric_connector_breaks_the_pairing(self):
        s = Q(1, 20)
        V, _ = drawing(0.5, 3.0, 6, 24, 'bent')
        counts, degenerate, _ = sector_counts(V, s)
        self.assertEqual(degenerate, [])
        self.assertEqual((counts['aligned'], counts['reversed']), (77, 78))
        V, _ = drawing(0.5, 3.0, 6, 24, 'straight')
        counts, _, _ = sector_counts(V, s)
        self.assertEqual((counts['aligned'], counts['reversed']), (78, 78))

    def test_recorded_scan_is_symmetric(self):
        data = json.loads((ROOT / 'results/v370/spiral_small_square_scan.json').read_text())
        self.assertGreaterEqual(len(data['rows']), 8)
        for row in data['rows']:
            self.assertEqual(row['degenerate'], 0)
            for s in ('1/10', '1/20', '1/40'):
                self.assertEqual(row[s]['aligned'], row[s]['reversed'])


if __name__ == '__main__':
    unittest.main()
