from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from rectangular_notch_v351 import notch_certificate
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class NotchTests(unittest.TestCase):
    def test_wide_notch_disk_comparison(self):
        for left in range(1, 7):
            for right in range(1, 7):
                for h in range(1, 6):
                    W = left+right+h+1
                    for H in (h+1, W, W+1):
                        model, _, floor = notch_certificate(W, H, left, right, h)
                        self.assertTrue(model.reaches_disk_target(floor))

    def test_components_on_actual_notch(self):
        for W, H, a, n, h, seed, floor, parity in (
            (14, 10, 7, 5, 3, ((3,0),(10,3),(7,10),(0,7)), 58, 1),
            (8, 7, 3, 2, 2, ((2,0),(5,2),(3,5),(0,3)), 13, 0),
            (8, 7, 3, 2, 2, ((0,0),(7,0),(7,7),(0,7)), 49, 1),
            (13, 7, 6, 1, 3, ((0,0),(7,0),(7,7),(0,7)), Q(49, 2), 1),
            (12, 7, 3, 6, 2, ((2,0),(5,2),(3,5),(0,3)), 13, 1),
        ):
            polygon = ((0,0),(W,0),(W,H),(a+n,H),(a+n,h),(a,h),(a,H),(0,H))
            cells = PolygonComponents(polygon).trace(tuple(pt(*p) for p in seed))
            self.assertEqual(min(segment_floor(c.squares) for c in cells.values()), floor)
            self.assertEqual(abs(sum(s for _, s in degree(cells)[1])), parity)
