from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from angular_windows_v328_v330 import area, radius_squared
from polygon_orbit_census_v223 import census_polygon, pt
from residue_census_v331 import census, certify, float_search, midside_fixture, residue_distance


def squares_by_fast_search(P):
    quads, _, _, _ = float_search(P, qmax2=1e9)
    found, _ = certify(P, quads, qmax2=Q(10**9))
    return found


class AgreesWithExhaustiveCensus(unittest.TestCase):
    def test_irregular_twelve_gon(self):
        poly = [(0, 0), (5, 1), (9, 0), (11, 4), (8, 5), (10, 9), (6, 8), (4, 11), (2, 7), (-2, 9), (-1, 4), (-3, 2)]
        P = [(Q(x), Q(y)) for x, y in poly]
        reference = {frozenset((v.x, v.y) for v in s.vertices) for s in census_polygon([pt(*p) for p in P]).squares}
        self.assertEqual(set(squares_by_fast_search(P)), reference)
        self.assertEqual(len(reference), 7)


class SpiralCensus(unittest.TestCase):
    def test_fixture_keeps_the_square_diameter(self):
        for N in (4, 20, 40):
            self.assertEqual(radius_squared(midside_fixture(N)), 8)
        self.assertGreater(area(midside_fixture(20)), 16)

    def test_no_survivors_below_the_threshold(self):
        for N in (4, 16):
            row = census(N, workers=1)
            self.assertEqual(row["survivors"], [])
            self.assertLess(row["max_float_residue_distance"], float(row["spectral_gap_lower"]))

    def test_twenty_turns(self):
        row = census(20, workers=2)
        self.assertEqual(len(row["survivors"]), 6)
        self.assertLess(row["marked_square_residue"], row["spectral_gap_lower"])   # marked square excluded
        top_distance, top = row["survivors"][0]
        self.assertEqual(set(top), {(Q(-1), Q(2)), (Q(-1, 2), Q(2)), (Q(-1, 2), Q(5, 2)), (Q(-1), Q(5, 2))})
        P = midside_fixture(20)
        self.assertEqual(residue_distance(P, top), top_distance)       # independent residue routine
        self.assertTrue(all(d >= row["spectral_gap_lower"] for d, _ in row["survivors"]))


if __name__ == "__main__":
    unittest.main()
