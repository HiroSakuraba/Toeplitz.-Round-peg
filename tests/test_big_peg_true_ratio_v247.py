"""v247: exact rational witness, global inradius, ratio ordering and controls."""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from big_peg_exact_v243 import read_polygon, contained_disk, maximum_square
from big_peg_inradius_v246 import (
    certify_global_radius_upper, min_boundary_distance_squared,
)
from big_peg_true_ratio_v247 import (
    audit, NEW_PATH, BASE_PATH, GRID_DENOMINATOR,
    NEW_DISK_CENTER, NEW_DISK_RADIUS2, NEW_GLOBAL_RADIUS_UPPER,
    NEW_MAX_SQUARE_SIDE2, OLD_GLOBAL_RADIUS_UPPER, OLD_MAX_SQUARE_SIDE2,
    EXPECTED_NEW_COVER_BOXES, EXPECTED_OLD_COVER_BOXES,
)


class RationalTrueRatioImprovement(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.new,cls.center,cls.lower2=read_polygon(NEW_PATH)
        cls.old,_,_=read_polygon(BASE_PATH)
        cls.data=json.loads(NEW_PATH.read_text(encoding="utf-8"))
        cls.certificate=audit()

    def test_source_fidelity_and_exact_integer_grid(self):
        deltas=self.data["provenance"]["vertex_deltas_from_v245_in_grid_units"]
        self.assertEqual(len(deltas),10)
        for original,modified,(dx,dy) in zip(self.old,self.new,deltas):
            self.assertEqual(modified,pt(
                original.x+Q(dx,GRID_DENOMINATOR),
                original.y+Q(dy,GRID_DENOMINATOR)))
            self.assertEqual((modified.x*GRID_DENOMINATOR).denominator,1)
            self.assertEqual((modified.y*GRID_DENOMINATOR).denominator,1)

    def test_exact_json_has_no_rounded_float_coordinates(self):
        self.assertTrue(all(type(x) is str for pair in self.data["vertices"]
                            for x in pair))
        self.assertTrue(all(type(x) is int
                            for row in self.data["provenance"]["vertex_deltas_from_v245_in_grid_units"]
                            for x in row))

    def test_both_polygons_are_jordan(self):
        self.assertEqual(validate_jordan_polygon(self.new),self.new)
        self.assertEqual(validate_jordan_polygon(self.old),self.old)

    def test_positive_new_disk_inside(self):
        self.assertEqual(self.center,NEW_DISK_CENTER)
        self.assertEqual(self.lower2,NEW_DISK_RADIUS2)
        self.assertTrue(contained_disk(self.new,self.center,self.lower2))
        self.assertGreater(self.lower2,1)

    def test_exact_disk_clearance_recomputed(self):
        self.assertEqual(min_boundary_distance_squared(self.new,self.center),
                         NEW_DISK_RADIUS2)

    def test_new_radius_cover_complete_and_deterministic(self):
        proof=certify_global_radius_upper(self.new,NEW_GLOBAL_RADIUS_UPPER)
        self.assertTrue(proof.complete)
        self.assertEqual(proof.visited_boxes,EXPECTED_NEW_COVER_BOXES)
        self.assertEqual(proof.visited_boxes,4*proof.split_boxes+1)
        self.assertGreater(proof.covered_boxes,0)
        self.assertGreater(proof.exterior_boxes,0)
        self.assertLessEqual(proof.maximum_depth,48)

    def test_previous_global_radius_reverified(self):
        proof=certify_global_radius_upper(self.old,OLD_GLOBAL_RADIUS_UPPER)
        self.assertTrue(proof.complete)
        self.assertEqual(proof.visited_boxes,EXPECTED_OLD_COVER_BOXES)

    def test_completed_polygonal_square_census_exact(self):
        found=maximum_square(self.new)
        self.assertEqual(found.charts,20000)
        self.assertEqual(found.side_squared,NEW_MAX_SQUARE_SIDE2)
        self.assertTrue(found.positive)
        self.assertEqual(len(found.witness),4)

    def test_previous_square_recomputed(self):
        self.assertEqual(maximum_square(self.old).side_squared,
                         OLD_MAX_SQUARE_SIDE2)

    def test_strict_true_J_inequality_exact_fraction(self):
        new_j2_upper=NEW_MAX_SQUARE_SIDE2/NEW_DISK_RADIUS2
        old_j2_lower=OLD_MAX_SQUARE_SIDE2/OLD_GLOBAL_RADIUS_UPPER**2
        self.assertLess(new_j2_upper,old_j2_lower)
        self.assertGreater(old_j2_lower-new_j2_upper,Q(0))
        self.assertGreater(old_j2_lower-new_j2_upper,Q(1,10000))

    def test_new_global_radius_interval_very_narrow(self):
        self.assertLess(
            NEW_GLOBAL_RADIUS_UPPER**2-NEW_DISK_RADIUS2,Q(1,10**9))
        self.assertGreater(
            NEW_GLOBAL_RADIUS_UPPER**2-NEW_DISK_RADIUS2,Q(0))

    def test_report_declares_scope_not_global_optimality(self):
        r=self.certificate
        self.assertIn("not ten-gon optimality",r["scope"])
        self.assertIn("Toeplitz remains open",r["scope"])
        self.assertIn("PASS",r["status"])

    def test_false_global_radius_fails_closed(self):
        with self.assertRaises(RuntimeError):
            certify_global_radius_upper(self.new,Q(1),max_depth=12)

    def test_unreasonably_low_node_limit_fails_closed(self):
        with self.assertRaises(RuntimeError):
            certify_global_radius_upper(self.new,NEW_GLOBAL_RADIUS_UPPER,
                                        max_boxes=1)

    def test_float_radius_rejected(self):
        with self.assertRaises(ValueError):
            certify_global_radius_upper(self.new,1.00133200316)

    def test_claimed_box_count_independent_of_source_json(self):
        self.assertEqual(self.certificate["new_global_cover_boxes"],1257)
        self.assertEqual(self.certificate["previous_global_cover_boxes"],917)


if __name__=="__main__":
    unittest.main()
