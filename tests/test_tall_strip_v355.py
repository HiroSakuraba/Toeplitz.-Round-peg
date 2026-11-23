from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from tall_strip_v355 import roof_disk_fits, tall_strip_witness
from convex_l_disk_v344_v345 import OrthogonalL
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class TallStripTests(unittest.TestCase):
    def test_square_components_and_disk_comparison(self):
        W,H,A,B = 15,14,8,13
        for p in (1,2,3):
            for h in (4,5,6,7,8,9):
                for g in (6,8,9,10):
                    q = Q(15,2)
                    result = tall_strip_witness(W,H,p,q,g,A,B,h)
                    self.assertIsNotNone(result)
                    seed, floor = result
                    poly = ((0,0),(W,0),(W,H),(B,H),(B,h),(A,h),(A,H),
                            (q,H),(q,g),(p,g),(p,H),(0,H))
                    cells = PolygonComponents(poly).trace(tuple(pt(*z) for z in seed))
                    self.assertEqual(min(segment_floor(c.squares) for c in cells.values()), floor)
                    self.assertEqual(abs(sum(s for _,s in degree(cells)[1])), 1)
                    self.assertTrue(OrthogonalL(W,H,A,h).reaches_disk_target(floor))

    def test_disk_can_exceed_old_bound(self):
        cuts = [(1,5,9),(8,10,6)]
        self.assertTrue(roof_disk_fits(Q(21,2),10,cuts,Q(21,5)))
        self.assertFalse(roof_disk_fits(Q(21,2),10,cuts,Q(421,100)))

    def test_disk_tangencies_and_single_cut_oracle(self):
        # Two cuts force the diameter-four disk between their vertical walls.
        cuts = [(1,2,1),(6,7,1)]
        self.assertTrue(roof_disk_fits(8,6,cuts,2))
        self.assertFalse(roof_disk_fits(8,6,cuts,Q(201,100)))
        for h in range(1,7):
            model = OrthogonalL(12,8,6,h)
            for k in range(1,18):
                r = Q(k,4)
                self.assertEqual(roof_disk_fits(12,8,[(4,6,h)],r),model.fits_disk(r))
