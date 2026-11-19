"""v245: rational candidate certification against a pinned primary source."""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from big_peg_exact_v243 import read_polygon, contained_disk, _distance_sq_to_segment
from big_peg_improvement_v245 import (
    verify, BASE, CANDIDATE, EXPECTED_BOUNDARY_D2, EXPECTED_MAX_S2,
    EXPECTED_RATIO_SQ, BASELINE_CERTIFIED_RATIO_SQ, ORIGINAL_L2,
    ORIGINAL_ORIGIN_DISK_D2,
)


class ExactImprovementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=verify()
        cls.baseline,_,_=read_polygon(BASE)
        cls.poly,cls.centre,cls.r2=read_polygon(CANDIDATE)
        cls.claim=json.loads(CANDIDATE.read_text())

    def test_real_exact_grid_perturbation(self):
        deltas=self.claim["provenance"]["integer_deltas_from_original_times_10"]
        self.assertEqual(len(deltas),10)
        for p0,p1,(dx,dy) in zip(self.baseline,self.poly,deltas):
            self.assertEqual(p1,pt(p0.x+Q(dx,12000),p0.y+Q(dy,12000)))

    def test_coordinates_reject_float(self):
        self.assertTrue(all(type(x) is str for pair in self.claim["vertices"] for x in pair))

    def test_jordan_polygon_and_disk_containment(self):
        self.assertEqual(validate_jordan_polygon(self.poly),self.poly)
        self.assertEqual(self.centre,pt(Q(3,10000),0))
        self.assertTrue(contained_disk(self.poly,self.centre,self.r2))

    def test_minimum_rational_boundary_distance(self):
        d2=min(_distance_sq_to_segment(self.centre,self.poly[i],self.poly[(i+1)%10])
               for i in range(10))
        self.assertEqual(d2,EXPECTED_BOUNDARY_D2)
        self.assertEqual(d2,self.r2)
        self.assertGreater(d2,1)

    def test_exact_complete_maximum(self):
        self.assertEqual(Q(self.report["exact_maximum_square_side_squared"]),
                         EXPECTED_MAX_S2)
        self.assertEqual(self.report["edge_chart_cases"],20000)
        self.assertEqual(len(self.report["maximizing_square_vertices"]),4)

    def test_candidate_counterexample_strict(self):
        self.assertLess(EXPECTED_MAX_S2,Q(2)*EXPECTED_BOUNDARY_D2)
        self.assertTrue(self.report["sharp_sqrt2_big_peg_counterexample"])

    def test_ratio_strictly_improved_over_original_explicit_disk_bound(self):
        self.assertEqual(ORIGINAL_L2, Q(56381275521625791352241234309,
                                       28343589374868261752462760000))
        self.assertEqual(ORIGINAL_ORIGIN_DISK_D2,Q(2403854890969,
                                                 2400326560000))
        self.assertEqual(EXPECTED_RATIO_SQ,EXPECTED_MAX_S2/EXPECTED_BOUNDARY_D2)
        self.assertLess(EXPECTED_RATIO_SQ, BASELINE_CERTIFIED_RATIO_SQ)
        self.assertTrue(self.report["strict_upper_bound_improved"])

    def test_claimed_json_numerators_are_audited(self):
        claim=self.claim["claimed_exact_invariants"]
        self.assertEqual(Q(claim["maximum_square_side_squared"]),EXPECTED_MAX_S2)
        self.assertEqual(Q(claim["minimum_boundary_distance_squared"]),EXPECTED_BOUNDARY_D2)
        self.assertEqual(Q(claim["candidate_squared_ratio_bound"]),EXPECTED_RATIO_SQ)

    def test_no_false_optimality_claim(self):
        self.assertIn("does not prove a globally minimal ratio",self.report["scope"])
        self.assertIn("Toeplitz",self.report["scope"])


if __name__=="__main__":
    unittest.main()
