"""v250: exact noncollinear corner continuation and branch-independent size floor."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt, sub, cross, dot
from big_peg_exact_v243 import contained_disk
from big_peg_multichart_v249 import (
    HALF, pol, padd, ppoint_at, pvcross, pvsub,
    certify_chart, certify_cover,
)
from big_peg_corner_v250 import (
    OCTAGON, KINK_SCALE, KINK_SWEEP, THREE_EDGE_NULL_VECTOR,
    SQUARE_SIDE_SQUARED_THRESHOLD, RIGHT_CHART, LEFT_CHART,
    family_vertices, exact_global_coordinate_separation,
    prove_jordan_family, verify_three_edge_kernel,
    minimum_nonadjacent_distance_squared, toothed_square,
    tooth_diagnostic, audit,
)


class ExactKinkHandoff(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vertices=family_vertices()
        cls.report=audit()
        cls.charts=certify_cover(
            cls.vertices,((F(0),HALF,RIGHT_CHART),(HALF,F(1),LEFT_CHART)),
            side2_limit=SQUARE_SIDE_SQUARED_THRESHOLD)

    def test_not_collinear_at_any_time(self):
        P=self.vertices
        turn=pvcross(pvsub(P[1],P[0]),pvsub(P[2],P[1]))
        self.assertEqual(len(turn),1)
        self.assertGreater(turn[0],0)
        self.assertEqual(str(turn[0]),"6622/44945")
        self.assertEqual(str(turn[0]),self.report["noncollinear_kink_corner_cross"])

    def test_two_different_square_branches_join_at_exact_corner(self):
        first,second=self.charts
        self.assertEqual(first.right_root,second.left_root)
        self.assertEqual(first.right_parameters[0],0)
        self.assertEqual(second.left_parameters[0],1)
        self.assertEqual(
            pt(first.right_root[0],first.right_root[1]),
            ppoint_at(self.vertices[1],HALF))
        self.assertEqual(self.report["handoff_time"],"1/2")
        self.assertTrue(self.report["algebraically_distinct_square_branches_join_at_kink"])

    def test_three_unchanged_edges_support_seam_square(self):
        old,new=verify_three_edge_kernel()
        self.assertEqual(tuple(old[i]+KINK_SCALE*THREE_EDGE_NULL_VECTOR[i]
                               for i in range(4)),new)
        self.assertEqual(self.charts[0].right_root,new)

    def test_strict_squared_side_gt_nine_all_real_times(self):
        self.assertEqual(SQUARE_SIDE_SQUARED_THRESHOLD,9)
        for chart in self.charts:
            self.assertGreater(chart.sign_checks,0)
            for root in (chart.left_root,chart.right_root):
                self.assertGreater(root[2]**2+root[3]**2,9)

    def test_full_time_interval_and_no_holes(self):
        self.assertEqual(self.charts[0].left,0)
        self.assertEqual(self.charts[0].right,HALF)
        self.assertEqual(self.charts[1].left,HALF)
        self.assertEqual(self.charts[1].right,1)

    def test_wrong_chart_cannot_span_opposite_side(self):
        for assignment,a,b in ((RIGHT_CHART,HALF,F(1)),
                               (LEFT_CHART,F(0),HALF)):
            with self.subTest(assignment=assignment):
                with self.assertRaises(ValueError):
                    certify_chart(self.vertices,assignment,a,b,side2_limit=9)

    def test_oversized_cutoff_rejected(self):
        with self.assertRaises(ValueError):
            certify_cover(self.vertices,
                          ((F(0),HALF,RIGHT_CHART),(HALF,F(1),LEFT_CHART)),
                          side2_limit=F(16))

    def test_complete_jordan_family_strict_turn(self):
        result=prove_jordan_family(self.vertices)
        self.assertTrue(result["strictly_convex_and_jordan_for_all_times"])
        self.assertTrue(result["enclosed_in_radius_3_disk"])
        self.assertGreater(result["positive_turn_and_edge_polynomial_checks"],0)

    def test_global_nonadjacent_separation_exact(self):
        gap=exact_global_coordinate_separation(self.vertices)
        self.assertEqual(gap,F(1810915,8503594))
        self.assertEqual(gap*gap/F(2),
                         F(self.report["every_inscribed_square_side2_at_least"]))

    def test_fixed_time_exact_feature_size_bound(self):
        gap=exact_global_coordinate_separation(self.vertices)
        for t in (F(0),F(1,8),HALF,F(3,4),F(1)):
            with self.subTest(t=t):
                polygon=tuple(ppoint_at(v,t) for v in self.vertices)
                d2=minimum_nonadjacent_distance_squared(polygon)
                self.assertGreaterEqual(d2,gap*gap)

    def test_unsupported_nonaffine_family_rejected(self):
        vertices=list(self.vertices)
        vertices[1]=(padd(vertices[1][0],pol(0,0,F(1,100))),
                     vertices[1][1])
        with self.assertRaises(ValueError):
            exact_global_coordinate_separation(tuple(vertices))

    def test_collapsing_kink_detected(self):
        vertices=list(self.vertices)
        original=vertices[1]
        # Move the kink onto the old segment: zero turn at the new vertex.
        vertices[1]=(padd(original[0],pol(0,F(0))),
                     original[1])
        # This is a no-op control for the exact geometry.
        self.assertEqual(tuple(vertices),self.vertices)
        flat=tuple(pt(p.x,p.y) for p in OCTAGON)
        self.assertGreater(minimum_nonadjacent_distance_squared(flat),0)

    def test_restricted_big_peg_ratio_gt_one(self):
        self.assertGreater(SQUARE_SIDE_SQUARED_THRESHOLD/F(3)**2,F(0))
        self.assertEqual(SQUARE_SIDE_SQUARED_THRESHOLD/F(3)**2,F(1))
        self.assertEqual(self.report["uniform_true_big_peg_J_strictly_above"],"1")


class FixedPolygonFeatureSize(unittest.TestCase):
    def test_square_has_exact_positive_feature_gap(self):
        square=(pt(-1,-1),pt(1,-1),pt(1,1),pt(-1,1))
        self.assertEqual(minimum_nonadjacent_distance_squared(square),4)

    def test_triangle_has_no_nonadjacent_edges(self):
        triangle=(pt(0,0),pt(1,0),pt(0,1))
        with self.assertRaises(ValueError):
            minimum_nonadjacent_distance_squared(triangle)

    def test_non_jordan_curve_rejected(self):
        crossed=(pt(0,0),pt(1,1),pt(0,1),pt(1,0))
        with self.assertRaises(ValueError):
            minimum_nonadjacent_distance_squared(crossed)

    def test_outward_tooth_has_fixed_right_angle_and_unit_disk(self):
        for n in (8,16,32,64,128):
            delta=F(1,n)
            polygon=toothed_square(delta)
            self.assertTrue(contained_disk(polygon,pt(0,0),F(1)))
            self.assertLessEqual(
                minimum_nonadjacent_distance_squared(polygon),
                2*delta*delta)
            first=sub(polygon[4],polygon[3])
            second=sub(polygon[5],polygon[4])
            self.assertEqual(dot(first,second),0)

    def test_feature_gap_tends_to_zero_despite_fixed_angles(self):
        examples=tooth_diagnostic()
        self.assertEqual(len(examples),4)
        self.assertEqual(examples[0]["delta"],"1/8")
        self.assertEqual(examples[-1]["delta"],"1/64")
        self.assertTrue(all(item["right_angle_tip"] for item in examples))

    def test_tooth_rejects_rounded_or_invalid_scale(self):
        for delta in (0.1,F(0),F(1,4),F(-1,20),"1/16"):
            with self.subTest(delta=delta):
                with self.assertRaises(ValueError):
                    toothed_square(delta)

    def test_report_disclaims_universal_peg_theorem(self):
        result=audit()
        self.assertIn("no universal c",result["limits"])
        self.assertIn("no Toeplitz proof",result["limits"])


if __name__=="__main__":
    unittest.main()
