from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from outward_tilt_v349 import witnesses
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class OutwardTests(unittest.TestCase):
    def test_component_floors_and_parity(self):
        # Width, height, elbow x, ledge y, top x. Include exact cap ties.
        for dims in ((12, 12, 2, 2, 7), (12, 6, 2, 2, 4),
                     (6, 8, 2, 2, 5), (6, 6, 2, 2, 4),
                     (10, 7, 2, 3, 6), (12, 5, 6, 1, 8),
                     (12, 6, 6, 2, 8), (90, 81, 1, 1, 41),
                     (20, 16, 1, 8, 5), (8, 9, 3, 2, Q(15, 2))):
            W, H, a, b, c = dims
            engine = PolygonComponents(((0, 0), (W, 0), (W, b),
                                        (a, b), (c, H), (0, H)))
            for kind, square, floor in witnesses(*dims):
                with self.subTest(dims=dims, kind=kind):
                    cells = engine.trace(tuple(pt(*p) for p in square))
                    self.assertEqual(min(segment_floor(v.squares) for v in cells.values()), floor)
                    self.assertEqual(abs(sum(s for _, s in degree(cells)[1])), 1)

    def test_disk_comparisons_in_two_square_regime(self):
        count = 0
        for H in range(5, 17):
            for c in range(2, H):
                for b in range(1, H):
                    rows = witnesses(H+1, H, 1, b, c)
                    if len(rows) == 2:
                        count += 1
                        k = Q(c-1, H-b)
                        elbow = rows[0][2]
                        self.assertGreaterEqual(elbow, max(1, b)**2)
                        self.assertGreaterEqual(elbow, 2/(2-k)**2)
                        # Upper triangle inradius = c/(1+k+sqrt(1+k*k)).
                        # Lower-bound the radical by 1 to test a stronger estimate.
                        self.assertGreaterEqual(rows[1][2], 2*c*c/(2+k)**2)
        self.assertGreater(count, 100)
