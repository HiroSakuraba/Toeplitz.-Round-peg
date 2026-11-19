"""v243 complete exact Big Peg polygon checks (standard-library only)."""
from fractions import Fraction as Q
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from polygon_orbit_census_v223 import pt, census_polygon
from big_peg_exact_v243 import (
    maximum_square, contained_disk, certify_big_peg,
    read_polygon, _candidate_vertices, _point_inside,
)

SQUARE=tuple(pt(*xy) for xy in ((-1,-1),(1,-1),(1,1),(-1,1)))
DIAMOND=tuple(pt(*xy) for xy in ((1,0),(0,1),(-1,0),(0,-1)))
TRIANGLE=tuple(pt(*xy) for xy in ((0,0),(3,0),(0,3)))


class ExactFeasibility(unittest.TestCase):
    def test_full_rank_unique_hit(self):
        points,rank=_candidate_vertices(
            ((Q(-1),Q(1),Q(0)),(Q(-2),Q(0),Q(1))),
            ((Q(0),Q(1),Q(0)),(Q(3),Q(-1),Q(0)),
             (Q(0),Q(0),Q(1)),(Q(3),Q(0),Q(-1))))
        self.assertEqual(rank,2)
        self.assertEqual(points,((Q(1),Q(2)),))

    def test_rank_one_segment(self):
        bounds=((Q(0),Q(1),Q(0)),(Q(1),Q(-1),Q(0)),
                (Q(0),Q(0),Q(1)),(Q(1),Q(0),Q(-1)))
        points,rank=_candidate_vertices(
            ((Q(-1),Q(1),Q(1)),(Q(-2),Q(2),Q(2))),bounds)
        self.assertEqual(rank,1)
        self.assertEqual(set(points),{(Q(0),Q(1)),(Q(1),Q(0))})

    def test_inconsistent_parallel_lines(self):
        points,rank=_candidate_vertices(
            ((Q(-1),Q(1),Q(1)),(Q(-2),Q(1),Q(1))),
            ((Q(0),Q(1),Q(0)),(Q(1),Q(-1),Q(0)),
             (Q(0),Q(0),Q(1)),(Q(1),Q(0),Q(-1))))
        self.assertEqual(points,())
        self.assertEqual(rank,-1)

    def test_rank_zero_2d_polytope(self):
        bounds=((Q(0),Q(1),Q(0)),(Q(1),Q(-1),Q(0)),
                (Q(0),Q(0),Q(1)),(Q(1),Q(0),Q(-1)))
        points,rank=_candidate_vertices(
            ((Q(0),Q(0),Q(0)),(Q(0),Q(0),Q(0))),bounds)
        self.assertEqual(rank,0)
        self.assertEqual(len(points),4)


class PolygonCensus(unittest.TestCase):
    def test_square_complete_maximum_exact(self):
        result=maximum_square(SQUARE)
        self.assertEqual(result.side_squared,Q(4))
        self.assertEqual(result.charts,512)
        self.assertTrue(result.positive)
        self.assertGreater(result.singular_feasible_charts,0)

    def test_diamond_complete_maximum_exact(self):
        result=maximum_square(DIAMOND)
        self.assertEqual(result.side_squared,Q(2))
        self.assertEqual(result.charts,512)

    def test_triangle_vs_existing_v223(self):
        result=maximum_square(TRIANGLE)
        previous=census_polygon(TRIANGLE)
        self.assertEqual(result.side_squared,
                         max(s.side_squared for s in previous.squares))
        self.assertGreater(result.side_squared,0)

    def test_reverse_orientation_invariance(self):
        a=maximum_square(TRIANGLE)
        b=maximum_square(tuple(reversed(TRIANGLE)))
        self.assertEqual(a.side_squared,b.side_squared)

    def test_scaled_square(self):
        S=tuple(pt(x*3,y*3) for x,y in
                ((-1,-1),(1,-1),(1,1),(-1,1)))
        self.assertEqual(maximum_square(S).side_squared,Q(36))

    def test_no_singular_charts_silently_dropped(self):
        result=maximum_square(SQUARE)
        self.assertGreater(result.singular_feasible_charts,0)
        self.assertGreater(result.feasible_charts,0)


class InradiusCertification(unittest.TestCase):
    def test_square_certifies_unit_disk(self):
        self.assertTrue(contained_disk(SQUARE,pt(0,0),Q(1)))
        self.assertFalse(contained_disk(SQUARE,pt(0,0),Q(10001,10000)))

    def test_diamond_centre_clearance(self):
        self.assertTrue(contained_disk(DIAMOND,pt(0,0),Q(1,2)))
        self.assertFalse(contained_disk(DIAMOND,pt(0,0),Q(5001,10000)))

    def test_outside_centre_rejected(self):
        self.assertFalse(_point_inside(SQUARE,pt(3,0)))
        self.assertFalse(contained_disk(SQUARE,pt(3,0),Q(1,10)))

    def test_zero_radius_rejected(self):
        with self.assertRaises(ValueError):
            contained_disk(SQUARE,pt(0,0),Q(0))

    def test_false_counterexample_not_misreported(self):
        result=certify_big_peg(SQUARE,pt(0,0),Q(1))
        self.assertFalse(result["counterexample_certified"])
        self.assertEqual(result["maximum_side_squared"],"4")
        self.assertEqual(result["J_squared_upper_bound"],"4")

    def test_claimed_disk_is_verified_not_assumed(self):
        with self.assertRaises(ValueError):
            certify_big_peg(SQUARE,pt(0,0),Q(2))


class ExactInputTests(unittest.TestCase):
    def test_rational_and_decimal_input(self):
        with TemporaryDirectory() as d:
            f=Path(d)/"polygon.json"
            f.write_text(json.dumps({"vertices":[
                ["-1","-1"],["1","-1"],["1","1"],["-1","1"]],
                "disk":{"center":["0.0","0"],"radius_squared":"1"}}))
            P,c,r2=read_polygon(f)
            self.assertEqual(P,SQUARE)
            self.assertEqual(c,pt(0,0))
            self.assertEqual(r2,Q(1))

    def test_float_input_rejected(self):
        with TemporaryDirectory() as d:
            f=Path(d)/"polygon.json"
            f.write_text(json.dumps({"vertices":[
                [-1,-1],[1.0,-1],[1,1],[-1,1]]}))
            with self.assertRaises(ValueError):
                read_polygon(f)

    def test_self_intersection_rejected(self):
        bow=tuple(pt(*v) for v in ((0,0),(2,2),(0,2),(2,0)))
        with self.assertRaises(ValueError):
            maximum_square(bow)


if __name__=="__main__":
    unittest.main()
