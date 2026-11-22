"""v331 residue census. Needs numpy; set REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a
failure instead of a skip (the dedicated CI job does this)."""
from fractions import Fraction as Q
from pathlib import Path
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from angular_windows_v328_v330 import area, pi_enclosure, radius_squared
from polygon_orbit_census_v223 import census_polygon, pt
try:
    from residue_census_v331 import census, certify, degenerate_audit, float_search, midside_fixture, residue_distance
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")


def squares_by_fast_search(P):
    quads, _, _, _ = float_search(P, qmax2=1e9)
    found, _ = certify(P, quads, qmax2=Q(10**9))
    return found


@needs_numpy
class AgreesWithExhaustiveCensus(unittest.TestCase):
    def test_irregular_twelve_gon(self):
        poly = [(0, 0), (5, 1), (9, 0), (11, 4), (8, 5), (10, 9), (6, 8), (4, 11), (2, 7), (-2, 9), (-1, 4), (-3, 2)]
        P = [(Q(x), Q(y)) for x, y in poly]
        reference = {frozenset((v.x, v.y) for v in s.vertices) for s in census_polygon([pt(*p) for p in P]).squares}
        self.assertEqual(set(squares_by_fast_search(P)), reference)
        self.assertEqual(len(reference), 7)


@needs_numpy
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

    def test_skipped_collinear_cases(self):
        P = midside_fixture(20)
        A = area(P)
        gap = min(A/2, A-pi_enclosure(30)[1]*radius_squared(P)/2)
        counts, families, point_survivors = degenerate_audit(P, gap, workers=2)
        self.assertEqual(counts["other"], 0)
        self.assertEqual(point_survivors, [])
        self.assertEqual(len(families), 1)                       # the frame-box family of tilted squares
        residues = {d for d in families[0]["residue_distances"]}
        self.assertEqual(len(residues), 1)                       # action is constant along the family
        self.assertLess(residues.pop(), gap)

    def test_family_minimum_can_be_interior(self):
        # [0, 1.99]^2 with the bottom-left corner cut: the family of tilted squares on the four
        # original sides has its smallest square strictly inside the parameter interval, below
        # every sampled endpoint and the midpoint (review finding on PR #85).
        L = Q(199, 100)
        P = [(Q(199, 2000), Q(0)), (L, Q(0)), (L, L), (Q(0), L), (Q(0), Q(199, 400))]
        _, families, _ = degenerate_audit(P, Q(0), qmax2=Q(1))
        self.assertEqual(len(families), 1)
        self.assertTrue(all(s*s/2 > 1 for s in families[0]["sides"]))   # every sample misses it


if __name__ == "__main__":
    unittest.main()
