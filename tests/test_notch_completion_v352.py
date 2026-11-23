from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt
from rectangular_notch_v352 import complete_notch_certificate


class NotchCompletionTests(unittest.TestCase):
    def test_regular_and_exceptional_coverage(self):
        for A in range(1, 7):
            for d in range(1, A+1):
                for n in range(1, 4):
                    W = A+n+d
                    for H in range(2, W+3):
                        for h in range(1, H):
                            model, _, floor = complete_notch_certificate(W, H, d, A, h)
                            self.assertTrue(model.reaches_disk_target(floor))

    def test_exceptional_marked_components(self):
        for A, n, d, h in ((6,1,4,3), (6,1,6,3), (5,2,3,2), (3,2,3,2),
                            (3,2,2,3), (Q(7,2),1,2,Q(5,2))):
            W = A+n+d
            for H in (A+n, W):
                poly = ((0,0),(W,0),(W,H),(A+n,H),(A+n,h),(A,h),(A,H),(0,H))
                seed = ((0,0),(H,0),(H,H),(0,H))
                cells = PolygonComponents(poly).trace(tuple(pt(*p) for p in seed))
                floor = min(segment_floor(c.squares) for c in cells.values())
                self.assertGreaterEqual(2*floor, H*H)
                self.assertEqual(abs(sum(s for _,s in degree(cells)[1])), 1)
