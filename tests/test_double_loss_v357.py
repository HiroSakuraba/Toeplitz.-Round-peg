from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from double_loss_v357 import double_loss_witness
from convex_l_disk_v344_v345 import OrthogonalL
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


def component(engine, seed):
    cells = engine.trace(tuple(pt(*p) for p in seed))
    return min(segment_floor(cell.squares) for cell in cells.values()), abs(sum(s for _, s in degree(cells)[1]))


def engine_for(W, H, p, q, g, A, B, h):
    return PolygonComponents(((0,0),(W,0),(W,H),(B,H),(B,h),(A,h),(A,H),
                              (q,H),(q,g),(p,g),(p,H),(0,H)))


class DoubleLossTests(unittest.TestCase):
    def test_replacements_and_size(self):
        for A in (6,7,8):
            for H in (Q(19,2), 10, 11, 12, 16):
                for g in (Q(89,10), 9, Q(91,10)):
                    dims = (15,H,2,5,g,A,11,4)
                    seed, floor = double_loss_witness(*dims)
                    self.assertEqual(component(engine_for(*dims),seed), (floor,1))
                    self.assertTrue(OrthogonalL(15,H,A,4).reaches_disk_target(floor))

    def test_two_pairs_at_once(self):
        for g in (Q(89,10), 9):
            engine = engine_for(15,12,2,5,g,8,11,4)
            for t in {g, 9}:
                tilted = ((4,0),(t,4),(t-4,t),(0,t-4))
                axis = ((11-t,0),(11,0),(11,t),(11-t,t))
                self.assertEqual(component(engine, tilted), (16+(t-4)**2,int(g<9)))
                self.assertEqual(component(engine, axis), (t*t,int(g<9)))
