"""v258: conditional near-square -> exact square, and singular-geometry controls."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import (
    pt,add,sub,dot,cross,validate_jordan_polygon,census_polygon,
)
from big_peg_persistence_v248 import (
    square_corner,certify_persistence,
)
from central_symmetry_square_v256 import (
    NOTCHED_12,certify_near_square,certify_centered_square,
    disk_certificate,is_centrally_symmetric,
)
from near_square_transversality_v258 import (
    BASE_TEN,BASE_ASSIGNMENT,R2,DELTA,GRAPH_EPSILON,
    asymmetric_ten,det4,certify_polygon_square_lift,
    certify_v257_near_square_lift,audit,
)


class TransverseNearSquareLifting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=asymmetric_ten()
        cls.near,cls.lift=certify_v257_near_square_lift(cls.P,R2,DELTA)
        cls.report=audit()

    def test_reference_is_genuinely_nonconvex_centrally_symmetric_jordan(self):
        P=BASE_TEN
        self.assertEqual(validate_jordan_polygon(P),P)
        self.assertTrue(is_centrally_symmetric(P))
        self.assertTrue(any(cross(sub(P[i],P[i-1]),
                                  sub(P[(i+1)%len(P)],P[i]))<0
                            for i in range(len(P))))
        self.assertGreater(disk_certificate(P,R2),R2)

    def test_target_asymmetric_nonconvex_jordan_and_disk(self):
        P=self.P
        self.assertEqual(validate_jordan_polygon(P),P)
        self.assertFalse(is_centrally_symmetric(P))
        self.assertGreater(disk_certificate(P,R2),R2)
        self.assertTrue(any(cross(sub(P[i],P[i-1]),
                                  sub(P[(i+1)%len(P)],P[i]))<0
                            for i in range(len(P))))

    def test_v257_witness_is_approximate_not_claimed_exact(self):
        self.assertFalse(self.near.exact_square_claimed)
        self.assertGreater(self.near.paired_vertex_defect_squared,0)
        self.assertLess(self.near.paired_vertex_defect_squared,DELTA**2)
        self.assertLessEqual(self.near.max_corner_error_squared,
                             self.near.paired_vertex_defect_squared)
        self.assertGreaterEqual(self.near.ideal_side_squared,2*R2)

    def test_transverse_incidence_four_by_four_is_nonsingular(self):
        self.assertEqual(self.lift.assignment,(3,5,8,0))
        self.assertEqual(self.lift.determinant,Q(-64))
        self.assertEqual(self.report["four_edge_incidence_determinant"],"-64")

    def test_conditioning_margin_quantifies_exact_lift(self):
        c=self.lift
        self.assertEqual(c.condition_response_bound,Q(6))
        self.assertEqual(c.guaranteed_parameter_displacement,Q(1,5))
        self.assertLessEqual(c.actual_parameter_displacement,
                             c.guaranteed_parameter_displacement)
        self.assertTrue(all(b<m for b,m in zip(
            c.maximum_relative_edge_parameter_bounds,
            c.verified_interior_contact_margins)))

    def test_every_corrected_corner_is_exactly_on_actual_polygon(self):
        c=self.lift
        for k,i in enumerate(c.assignment):
            a=self.P[i];b=self.P[(i+1)%len(self.P)]
            q=c.exact_square_corners[k]
            self.assertEqual(c.exact_square_corners[k],
                             square_corner(c.exact_square_parameters,k))
            self.assertEqual(cross(sub(b,a),sub(q,a)),0)
            self.assertGreater(c.exact_edge_parameters[k],Q(0))
            self.assertLess(c.exact_edge_parameters[k],Q(1))

    def test_corrected_square_is_nonzero_and_big(self):
        c=self.lift
        self.assertGreater(c.ideal_side_squared,2*R2)
        self.assertGreater(c.guaranteed_side_squared_lower,Q(28))
        self.assertGreater(c.exact_side_squared,c.guaranteed_side_squared_lower)
        self.assertEqual(c.exact_side_squared,Q(94179510917,2934392832))

    def test_exact_root_is_not_just_the_original_approximate_square(self):
        c=self.lift
        self.assertNotEqual(c.exact_square_parameters,
                            c.ideal_square_parameters)
        self.assertGreater(c.actual_parameter_displacement,0)
        self.assertGreater(c.max_corner_error_squared if hasattr(c,'max_corner_error_squared')
                           else c.maximum_observed_corner_error_squared,Q(0))

    def test_derived_rational_square_root_matches_independent_census(self):
        # Independent v223 complete ORBIT census, not the same 4-line solver.
        result=census_polygon(self.P)
        expected=frozenset(c.exact_square_corners for c in [self.lift])
        self.assertTrue(any(frozenset(s.vertices) in expected
                            for s in result.squares))
        self.assertGreaterEqual(result.count,1)

    def test_repeatable_exact_lift_from_raw_v257_corners(self):
        c=certify_polygon_square_lift(
            self.P,self.near.ideal_centered_square,
            self.near.exact_boundary_points,DELTA,disk_r2=R2)
        self.assertEqual(c,self.lift)

    def test_general_delta_bound_may_be_loose_but_must_hold(self):
        c=self.lift
        self.assertLessEqual(c.maximum_observed_corner_error_squared,DELTA**2)
        self.assertLessEqual(c.actual_parameter_displacement,
                             c.guaranteed_parameter_displacement)

    def test_reference_square_rough_persistence_for_arbitrary_continuous_graphs(self):
        # Strong v248 theorem for ALL continuous normal graphs <=1/100,
        # not just polygonal perturbations or smooth arcs.
        cert=certify_persistence(BASE_TEN,BASE_ASSIGNMENT,GRAPH_EPSILON)
        self.assertGreater(cert.minimum_square_side_squared,Q(0))
        self.assertGreater(cert.original_square_side_squared,Q(30))
        self.assertEqual(self.report["reference_normal_graph_rough_epsilon"],"1/100")

    def test_singular_axis_aligned_near_square_does_not_meet_transversality(self):
        # v256 notched 12-gon has centered diamond points on four
        # axis-aligned edges; its four-line incidence determinant is 0.
        P=NOTCHED_12
        ideal=(pt(4,0),pt(0,4),pt(-4,0),pt(0,-4))
        cert=certify_centered_square(P,Q(4))
        self.assertGreaterEqual(cert.side_squared,Q(8))
        with self.assertRaisesRegex(ValueError,"singular"):
            certify_polygon_square_lift(
                P,ideal,ideal,Q(0),assignment=(5,10,11,4),disk_r2=Q(4))

    def test_insufficient_delta_rejects_unproved_near_square_proximity(self):
        with self.assertRaises(ValueError):
            certify_polygon_square_lift(
                self.P,self.near.ideal_centered_square,
                self.near.exact_boundary_points,Q(1,1000),disk_r2=R2)

    def test_excessively_loose_delta_fails_margin_check(self):
        with self.assertRaisesRegex(ValueError,"conditioning"):
            certify_polygon_square_lift(
                self.P,self.near.ideal_centered_square,
                self.near.exact_boundary_points,Q(1,2),disk_r2=R2)

    def test_bad_edge_assignment_rejected(self):
        with self.assertRaises(ValueError):
            certify_polygon_square_lift(
                self.P,self.near.ideal_centered_square,
                self.near.exact_boundary_points,DELTA,
                assignment=(3,5,9,0),disk_r2=R2)

    def test_repeated_edge_assignment_rejected(self):
        with self.assertRaises(ValueError):
            certify_polygon_square_lift(
                self.P,self.near.ideal_centered_square,
                self.near.exact_boundary_points,DELTA,
                assignment=(3,5,3,0),disk_r2=R2)

    def test_wrong_disk_hypothesis_rejected(self):
        with self.assertRaises(ValueError):
            certify_polygon_square_lift(
                self.P,self.near.ideal_centered_square,
                self.near.exact_boundary_points,DELTA,disk_r2=Q(10))

    def test_non_square_ideal_is_rejected(self):
        ideal=list(self.near.ideal_centered_square)
        ideal[2]=add(ideal[2],pt(Q(1,10),Q(0)))
        with self.assertRaises(ValueError):
            certify_polygon_square_lift(
                self.P,ideal,self.near.exact_boundary_points,DELTA)

    def test_invalid_float_delta_rejected(self):
        with self.assertRaises(ValueError):
            certify_polygon_square_lift(
                self.P,self.near.ideal_centered_square,
                self.near.exact_boundary_points,1/30)

    def test_invalid_polygon_rejected_before_lift(self):
        p=list(self.P)
        p[2]=p[4]  # crosses or duplicates vertices: not Jordan
        with self.assertRaises(ValueError):
            certify_polygon_square_lift(
                tuple(p),self.near.ideal_centered_square,
                self.near.exact_boundary_points,DELTA)

    def test_no_general_toeplitz_or_unconditional_symmetry_claim(self):
        self.assertIn("Conditional",self.report["limitation"])
        self.assertIn("ALONE",self.report["limitation"])
        self.assertIn("No universal Big Peg",self.report["limitation"])


if __name__=="__main__":
    unittest.main()
