from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from square_relabel_path_v362 import replay_path


class SquarePathTests(unittest.TestCase):
    polygon = [(1,0),(Q(3,5),Q(4,5)),(0,1),(-Q(4,5),Q(3,5)),
               (-1,0),(-Q(3,5),-Q(4,5)),(0,-1),(Q(4,5),-Q(3,5))]
    knots = [Q(i,8) for i in range(9)]
    start = tuple(Q(i,4) for i in range(4))

    def test_all_nonidentity_relabelings(self):
        for k in (1,2,3):
            end = tuple(t+Q(k,4) for t in self.start)
            result = replay_path(self.polygon,self.knots,[self.start,end])
            self.assertEqual(result.relabeling,k)
            self.assertEqual(result.max_residual_squared,0)
            self.assertEqual(result.parameter_distance,Q(1,2))
            self.assertGreater(result.chord_squared,0)

    def test_interior_knot_error_on_asymmetric_redraw(self):
        end = tuple(t+Q(1,4) for t in self.start)
        for e in (Q(1,16),Q(1,64),Q(1,256)):
            polygon = self.polygon.copy()
            polygon[1] = (Q(3,5)+e,Q(4,5))
            result = replay_path(polygon,self.knots,[self.start,end])
            self.assertEqual(result.max_residual_squared,2*e**2)

    def test_tight_cluster_must_spread_and_incurs_residual(self):
        row = (Q(0),Q(1,100),Q(2,100),Q(3,100))
        end = row[1:]+(Q(1),)
        r = replay_path(self.polygon,self.knots,[row,end])
        self.assertGreaterEqual(r.parameter_distance,Q(1,3))
        self.assertGreater(r.max_residual_squared,1)
        with self.assertRaises(ValueError):
            replay_path(self.polygon,self.knots,[row,row[1:]+(Q(0),)])
        with self.assertRaises(ValueError):
            replay_path(self.polygon,self.knots,[row,row])

    def test_blind_grid_finds_and_replays_path(self):
        from square_relabel_path_v362 import find_vertex_path
        for e in (Q(0),Q(1,64)):
            polygon = self.polygon.copy()
            polygon[1] = (Q(3,5)+e,Q(4,5))
            rows = find_vertex_path(polygon,2*e**2)
            self.assertIsNotNone(rows)
            r = replay_path(polygon,self.knots,rows)
            self.assertLessEqual(r.max_residual_squared,2*e**2)
            self.assertGreaterEqual(r.parameter_distance,Q(1,3))
        self.assertIsNone(find_vertex_path([(0,0),(4,0),(2,3)],0))
