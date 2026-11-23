from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from inward_tilt_v350 import anchor_interval, witness
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class InwardTests(unittest.TestCase):
    def check(self, polygon, seed, floor):
        cells = PolygonComponents(polygon).trace(tuple(pt(*p) for p in seed))
        self.assertEqual(min(segment_floor(c.squares) for c in cells.values()), floor)
        self.assertEqual(abs(sum(s for _, s in degree(cells)[1])), 1)

    def test_inward_components(self):
        for W in range(4, 8):
            for H in (3, 5, 7):
                for a in range(2, W):
                    for c in range(1, a):
                        for b in (1, H-1):
                            dims = W, H, a, b, c
                            with self.subTest(dims=dims):
                                self.check(((0, 0), (W, 0), (W, b), (a, b),
                                            (c, H), (0, H)), *witness(*dims))
        dims = (9, 8, Q(13, 2), Q(1, 3), Q(3, 2))
        W, H, a, b, c = dims
        self.check(((0, 0), (W, 0), (W, b), (a, b), (c, H), (0, H)), *witness(*dims))

    def test_multiple_steps(self):
        frontier = ((9, 0), (9, 1), (6, 1), (6, 3), (4, 3), (4, 6), (2, 8), (0, 8))
        t, lo, hi = anchor_interval(frontier)
        u = max(lo, min(t/2, hi))
        v = t-u
        self.check(((0, 0), *frontier), ((u, 0), (t, u), (v, t), (0, v)), u*u+v*v)
