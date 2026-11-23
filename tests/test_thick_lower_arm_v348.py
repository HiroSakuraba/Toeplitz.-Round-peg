from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class ThickLowerArmTests(unittest.TestCase):
    def test_component_floor_and_degree(self):
        W = 6
        for b in (6, Q(601, 100), 8):
            for a in (2, 3, 4):
                for interior in ((), ((5, b+1), (5, b+3), (2, b+2), (1, b+4))):
                    with self.subTest(b=b, a=a, folded=bool(interior)):
                        H = b+5
                        polygon = ((0, 0), (W, 0), (W, b), (a, b),
                                   *interior, (3, H), (0, H))
                        seed = ((0, 0), (W, 0), (W, W), (0, W))
                        cells = PolygonComponents(polygon).trace(tuple(pt(*p) for p in seed))
                        d = min(W-a, Q(W, 2)) if b == W else 0
                        self.assertEqual(min(segment_floor(c.squares) for c in cells.values()),
                                         d*d + (W-d)**2)
                        self.assertEqual(abs(sum(s for _, s in degree(cells)[1])), 1)
