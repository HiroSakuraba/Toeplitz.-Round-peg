"""v246: test complete rational global inradius upper bounds and ratio separation."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from polygon_orbit_census_v223 import pt
from big_peg_exact_v243 import read_polygon
from big_peg_inradius_v246 import (
    ORIGINAL,IMPROVED,ORIGINAL_CENTER,IMPROVED_CENTER,
    ORIGINAL_RADIUS_UPPER,IMPROVED_RADIUS_UPPER,
    ORIGINAL_RADIUS2_LOWER,IMPROVED_RADIUS2_LOWER,
    ORIGINAL_SIDE2,IMPROVED_SIDE2,
    min_boundary_distance_squared,certify_global_radius_upper,
    audit,
)


class GlobalRationalCoverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original,_,_=read_polygon(ORIGINAL)
        cls.improved,_,_=read_polygon(IMPROVED)
        cls.report=audit()

    def test_original_exact_lower_disk(self):
        self.assertEqual(
            min_boundary_distance_squared(self.original,ORIGINAL_CENTER),
            ORIGINAL_RADIUS2_LOWER)
        self.assertGreater(ORIGINAL_RADIUS2_LOWER,1)
        self.assertLess(ORIGINAL_RADIUS2_LOWER,ORIGINAL_RADIUS_UPPER**2)

    def test_candidate_exact_lower_disk(self):
        self.assertEqual(
            min_boundary_distance_squared(self.improved,IMPROVED_CENTER),
            IMPROVED_RADIUS2_LOWER)
        self.assertGreater(IMPROVED_RADIUS2_LOWER,1)
        self.assertLess(IMPROVED_RADIUS2_LOWER,IMPROVED_RADIUS_UPPER**2)

    def test_both_exhaustive_rectangle_covers(self):
        for poly,R in ((self.original,ORIGINAL_RADIUS_UPPER),
                       (self.improved,IMPROVED_RADIUS_UPPER)):
            proof=certify_global_radius_upper(poly,R)
            self.assertTrue(proof.complete)
            self.assertEqual(proof.visited_boxes,4*proof.split_boxes+1)
            self.assertGreater(proof.covered_boxes,0)
            self.assertGreater(proof.exterior_boxes,0)
            self.assertLessEqual(proof.maximum_depth,48)

    def test_expected_audited_box_counts_are_reproducible(self):
        self.assertEqual(self.report["original"]["quadtree_boxes_checked"],885)
        self.assertEqual(self.report["improved"]["quadtree_boxes_checked"],917)

    def test_true_J_ranges_exact_and_disjoint(self):
        # New J <= sqrt(s_new^2/r_new_lower^2);
        # old J >= sqrt(s_old^2/r_old_upper^2).
        new_upper_sq=IMPROVED_SIDE2/IMPROVED_RADIUS2_LOWER
        old_lower_sq=ORIGINAL_SIDE2/ORIGINAL_RADIUS_UPPER**2
        self.assertLess(new_upper_sq,old_lower_sq)
        self.assertGreater(old_lower_sq-new_upper_sq,0)

    def test_ratio_bounds_each_have_correct_order(self):
        self.assertLessEqual(
            ORIGINAL_SIDE2/ORIGINAL_RADIUS_UPPER**2,
            ORIGINAL_SIDE2/ORIGINAL_RADIUS2_LOWER)
        self.assertLessEqual(
            IMPROVED_SIDE2/IMPROVED_RADIUS_UPPER**2,
            IMPROVED_SIDE2/IMPROVED_RADIUS2_LOWER)

    def test_original_interval_is_tight(self):
        self.assertLess(
            ORIGINAL_RADIUS_UPPER**2-ORIGINAL_RADIUS2_LOWER,Q(1,10**7))

    def test_improved_interval_is_tight(self):
        self.assertLess(
            IMPROVED_RADIUS_UPPER**2-IMPROVED_RADIUS2_LOWER,Q(1,10**7))

    def test_irrational_sqrt_not_rounded_into_decision(self):
        r=self.report
        self.assertIn("TRUE v245 Big Peg ratio",r["theorem"])
        self.assertIn("global rational-covering",r["status"])
        self.assertGreater(Q(r["strict_j_squared_gap_lower_exact"]),0)
        self.assertTrue(all("/" in r[key]["J_squared_upper_exact"]
                            for key in ("original","improved")))

    def test_smaller_than_known_interior_disk_fails_closed(self):
        # A radius upper bound of 1 is impossible: original has a
        # certified disk radius strictly >1. Max-depth causes an error,
        # not a false success.
        with self.assertRaises(RuntimeError):
            certify_global_radius_upper(self.original,Q(1),max_depth=12)

    def test_invalid_radii_fail_closed(self):
        for invalid in (Q(0),Q(-1),1,1.0,"1.2"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    certify_global_radius_upper(self.original,invalid)

    def test_resource_limit_fails_closed(self):
        with self.assertRaises(RuntimeError):
            certify_global_radius_upper(self.original,ORIGINAL_RADIUS_UPPER,
                                        max_boxes=1)

    def test_depth_limit_fails_closed(self):
        with self.assertRaises(RuntimeError):
            certify_global_radius_upper(self.original,ORIGINAL_RADIUS_UPPER,
                                        max_depth=0)

    def test_easy_exact_square_control(self):
        square=tuple(pt(*xy) for xy in (
            (-1,-1),(1,-1),(1,1),(-1,1)))
        proof=certify_global_radius_upper(square,Q(11,10))
        self.assertTrue(proof.complete)
        self.assertLess(proof.maximum_depth,20)

    def test_bounding_boxes_include_every_vertex(self):
        for poly,R in ((self.original,ORIGINAL_RADIUS_UPPER),
                       (self.improved,IMPROVED_RADIUS_UPPER)):
            b=certify_global_radius_upper(poly,R).bounding_box
            self.assertTrue(all(b[0]<=p.x<=b[1] and b[2]<=p.y<=b[3]
                                for p in poly))


if __name__=="__main__":
    unittest.main()
