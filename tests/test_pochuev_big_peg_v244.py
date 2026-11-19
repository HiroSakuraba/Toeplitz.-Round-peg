"""v244 Pochuev witness: exact, replayable and non-heuristic checks."""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from big_peg_exact_v243 import (
    read_polygon, contained_disk, maximum_square, _distance_sq_to_segment,
)
from pochuev_big_peg_v244 import (
    audit, WITNESS, SOURCE_COMMIT, SOURCE_BLOB,
    EXPECTED_VERTICES_GRID, EXPECTED_MIN_BOUNDARY_D2,
    EXPECTED_MAX_SIDE2, EXPECTED_GAP,
)


class ExactPochuevWitness(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads(WITNESS.read_text(encoding="utf-8"))
        cls.P,cls.center,cls.r2=read_polygon(WITNESS)
        cls.report=audit()  # Full 20,000-chart exact census, once for class.

    def test_original_provenance_pinned(self):
        self.assertEqual(self.data["provenance"]["pinned_commit"],SOURCE_COMMIT)
        self.assertEqual(self.data["provenance"]["source_blob_sha"],SOURCE_BLOB)
        self.assertIn("Grisha-Pochuev",self.data["provenance"]["repository"])
        self.assertEqual(self.data["provenance"]["coordinate_denominator"],1200)

    def test_ten_exact_source_vertices(self):
        self.assertEqual(len(self.P),10)
        self.assertEqual(tuple((int(p.x*1200),int(p.y*1200)) for p in self.P),
                         EXPECTED_VERTICES_GRID)
        self.assertTrue(all((p.x*1200).denominator==1 and
                            (p.y*1200).denominator==1 for p in self.P))

    def test_jordan_simplicity(self):
        self.assertEqual(validate_jordan_polygon(self.P),self.P)

    def test_unit_disk_inside(self):
        self.assertEqual(self.center,pt(0,0))
        self.assertEqual(self.r2,1)
        self.assertTrue(contained_disk(self.P,self.center,self.r2))

    def test_independent_exact_boundary_distance(self):
        md2=min(_distance_sq_to_segment(
            self.center,self.P[i],self.P[(i+1)%10]) for i in range(10))
        self.assertEqual(md2,EXPECTED_MIN_BOUNDARY_D2)
        self.assertGreater(md2,Q(1))
        self.assertEqual(md2,self.report["minimum_boundary_distance_squared"] and
                         Q(self.report["minimum_boundary_distance_squared"]))

    def test_exact_maximum_matches_original_source(self):
        self.assertEqual(Q(self.report["maximum_inscribed_square_side_squared"]),
                         EXPECTED_MAX_SIDE2)
        self.assertEqual(self.report["exact_edge_charts"],20000)

    def test_strict_inequality_is_exact(self):
        self.assertEqual(2-EXPECTED_MAX_SIDE2,EXPECTED_GAP)
        self.assertGreater(EXPECTED_GAP,0)
        self.assertLess(EXPECTED_MAX_SIDE2,Q(2))
        self.assertTrue(self.report["sharp_big_peg_refuted_for_polygon"])

    def test_maximum_square_has_attained_witness(self):
        self.assertEqual(len(self.report["maximizing_square_vertices"]),4)
        self.assertEqual(len({tuple(x) for x in self.report["maximizing_square_vertices"]}),4)
        self.assertEqual(len(self.report["maximizing_square_edge_assignment"]),4)

    def test_source_claims_equal_recomputed_rationals(self):
        expected=self.data["claimed_exact_invariants"]
        self.assertEqual(Q(expected["maximum_square_side_squared"]),EXPECTED_MAX_SIDE2)
        self.assertEqual(Q(expected["minimum_boundary_distance_squared"]),
                         EXPECTED_MIN_BOUNDARY_D2)
        self.assertEqual(Q(expected["gap_to_two"]),EXPECTED_GAP)

    def test_no_general_big_peg_or_toeplitz_claim(self):
        self.assertIn("does not refute",self.report["limits"])
        self.assertIn("Toeplitz",self.report["limits"])

    def test_no_float_coordinates(self):
        self.assertTrue(all(isinstance(x,str) for row in self.data["vertices"] for x in row))


if __name__=="__main__":
    unittest.main()
