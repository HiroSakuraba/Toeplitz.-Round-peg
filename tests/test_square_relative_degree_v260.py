"""v260 exact degree gluing at positive-square vertex pair annihilation."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt,sub,dot,cross,census_polygon
from big_peg_pair_death_v254 import (
    family_vertices,certify_collision,certify_jordan_family,
    COLLIDING_LEFT,COLLIDING_RIGHT,QPOLY,ROOT_LEFT,ROOT_RIGHT,
)
from big_peg_multichart_v249 import peval,ppoint_at
from moving_polygon_events_v224 import divide,iszero
from square_relative_degree_v260 import (
    START,END,BOX,INTERIOR_GAP,ENTRY_BOUND,DETERMINANT_FLOOR,
    UNIFORM_INVERSE_NORM_BOUND,BOUNDARY_RESIDUAL_FLOOR,
    C0_EPSILON,PERTURBATION_RESIDUAL_RESPONSE,
    SIDE_X_FLOOR,SIDE_Y_FLOOR,
    _snapshot,_square_points,_map,_matrix,_chart_at,
    _local_numerators,_normalized_event_polynomial,
    certify_cell,certify_interface_gluing,
    certify_sides_and_robust_boundary,_verify_independent_snapshots,audit,
)
from near_square_transversality_v258 import det4


class SignedGluedSquareDegree(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.left=certify_cell(cls.P,COLLIDING_LEFT)
        cls.right=certify_cell(cls.P,COLLIDING_RIGHT)
        cls.report=audit()

    def test_all_time_convex_jordan_pentagon(self):
        certified=certify_jordan_family(self.P)
        self.assertTrue(certified["convex_Jordan_for_all_real_times"])
        self.assertGreater(certified["global_nonadjacent_edge_rect_gap"],0)

    def test_noncollinear_vertex_hinge_is_continuous(self):
        for t in (START,(START+END)/2,END):
            p=_snapshot(self.P,t)
            pre=_square_points(self.P,t,(Q(91,100),Q(0),Q(44,100),Q(75,100)))
            self.assertEqual(pre[1],p[2])
            self.assertNotEqual(sub(p[2],p[1]),sub(p[3],p[2]))
            self.assertGreater(cross(sub(p[2],p[1]),sub(p[3],p[2])),0)

    def test_each_cell_det_has_opposite_nonzero_orientation(self):
        self.assertEqual(self.left.orientation_sign,-1)
        self.assertEqual(self.right.orientation_sign,1)
        self.assertEqual((self.left.determinant_floor,self.right.determinant_floor),
                         (Q(1,10),Q(1,10)))

    def test_independent_matrix_det_agrees_with_polynomial_at_both_ends(self):
        for ids,c in ((COLLIDING_LEFT,self.left),(COLLIDING_RIGHT,self.right)):
            chart=_chart_at(self.P,ids)
            for tt,tau in ((START,Q(0)),(END,Q(1))):
                self.assertEqual(det4(_matrix(self.P,tt,ids)),
                                 peval(chart.determinant,tau))
                self.assertGreater(abs(det4(_matrix(self.P,tt,ids))),
                                   DETERMINANT_FLOOR)

    def test_entry_bound_verified_for_both_moving_cells(self):
        for ids in (COLLIDING_LEFT,COLLIDING_RIGHT):
            for tt in (START,END):
                mat=_matrix(self.P,tt,ids)
                self.assertLessEqual(max(abs(v) for row in mat for v in row),
                                     ENTRY_BOUND)

    def test_strict_exact_outer_root_margin_for_all_real_time(self):
        self.assertGreater(self.left.parameter_coordinate_sign_checks,8)
        self.assertGreater(self.right.parameter_coordinate_sign_checks,8)
        for c in (self.left,self.right):
            self.assertEqual(c.gap_from_outer_faces,INTERIOR_GAP)
            for root in (c.signed_roots_at_start,c.signed_roots_at_end):
                for j,(lo,hi) in enumerate(BOX):
                    self.assertTrue(lo+INTERIOR_GAP<root[j]<hi-INTERIOR_GAP)

    def test_algebraic_time_event_quadratic_is_reused_exactly(self):
        q=_normalized_event_polynomial()
        self.assertGreater(peval(q,Q(0)),0)
        self.assertLess(peval(q,Q(1)),0)
        self.assertEqual(len(q),3)
        self.assertEqual(tuple(QPOLY),
                         tuple(Q(x) for x in self.report["algebraic_collision_polynomial"]))

    def test_both_critical_parameters_have_unique_common_event_factor(self):
        q=_normalized_event_polynomial()
        for ids in (COLLIDING_LEFT,COLLIDING_RIGHT):
            chart=_chart_at(self.P,ids)
            nums=_local_numerators(chart,ids)
            k,r=divide(nums[1],q)
            self.assertTrue(iszero(r))
            self.assertNotEqual(k,(Q(0),))
            self.assertGreater(peval(k,Q(0)),0)
            self.assertGreater(peval(k,Q(1)),0)

    def test_full_four_parameter_roots_glue_at_vertex(self):
        quotients=certify_interface_gluing(self.P)
        self.assertEqual(len(quotients),4)
        self.assertTrue(any(any(v!=0 for v in q) for q in quotients))
        self.assertTrue(certify_collision(self.P)["shared_second_corner_at_pentagon_vertex_2"])

    def test_before_event_there_are_two_distinct_actual_local_squares(self):
        a=self.left.signed_roots_at_start
        b=self.right.signed_roots_at_start
        self.assertLess(a[1],0)
        self.assertGreater(b[1],0)
        self.assertNotEqual(a,b)
        self.assertEqual(_map(self.P,START,a),(Q(0),)*4)
        self.assertEqual(_map(self.P,START,b),(Q(0),)*4)
        for root in (a,b):
            self.assertTrue(all(lo<z<hi for z,(lo,hi) in zip(root,BOX)))
            points=_square_points(self.P,START,root)
            self.assertGreater(dot(sub(points[1],points[0]),
                                   sub(points[1],points[0])),Q(1))

    def test_after_event_both_extended_solutions_have_wrong_hinge_sign(self):
        a=self.left.signed_roots_at_end
        b=self.right.signed_roots_at_end
        self.assertGreater(a[1],0)
        self.assertLess(b[1],0)
        self.assertTrue(any(_map(self.P,END,a)))
        self.assertTrue(any(_map(self.P,END,b)))

    def test_algebraic_collision_inner_seam_not_outer_boundary(self):
        self.assertTrue(self.report["algebraic_collision_strictly_inside_time_bracket"])
        self.assertTrue(self.report["inner_interface_undefined_cell_degrees_at_event"])
        self.assertEqual(self.report["glued_degree_at_nonsmooth_algebraic_event"],0)

    def test_pre_event_orientation_indices_cancel(self):
        self.assertEqual((self.left.orientation_sign,self.right.orientation_sign),(-1,1))
        self.assertEqual(self.report["signed_pre_collision_sum"],0)
        self.assertEqual(self.report["signed_post_collision_sum"],0)

    def test_no_degree_can_be_assigned_to_individual_halfbox_at_event(self):
        self.assertEqual(self.report["inner_interface_undefined_cell_degrees_at_event"],True)
        self.assertIn("degree_zero",self.report["scope"].replace("Degree zero","degree_zero"))

    def test_inverse_norm_bound_and_uniform_outer_residual_floor(self):
        self.assertEqual(UNIFORM_INVERSE_NORM_BOUND,Q(1920))
        self.assertEqual(BOUNDARY_RESIDUAL_FLOOR,Q(1,19200000))
        self.assertGreater(BOUNDARY_RESIDUAL_FLOOR,0)

    def test_small_arbitrary_continuous_perturbation_preserves_zero_degree(self):
        robust=certify_sides_and_robust_boundary(self.P,C0_EPSILON)
        self.assertEqual(robust["safe_continuous_coordinate_epsilon"],Q(1,10**9))
        self.assertLess(robust["max_continuous_perturbed_defect"],
                        robust["rigorous_unperturbed_boundary_defect_floor"])
        self.assertEqual(robust["max_continuous_perturbed_defect"],
                         PERTURBATION_RESIDUAL_RESPONSE*C0_EPSILON)

    def test_critical_epsilon_at_boundary_margin_fails_closed(self):
        critical=BOUNDARY_RESIDUAL_FLOOR/PERTURBATION_RESIDUAL_RESPONSE
        with self.assertRaises(ValueError):
            certify_sides_and_robust_boundary(self.P,critical)
        with self.assertRaises(ValueError):
            certify_sides_and_robust_boundary(self.P,critical+Q(1,10**12))

    def test_float_or_negative_epsilon_fails_closed(self):
        for eps in (1e-9,Q(-1,10**9),None,"1/1000000000"):
            with self.subTest(eps=eps):
                with self.assertRaises(ValueError):
                    certify_sides_and_robust_boundary(self.P,eps)

    def test_multiaffine_side_noncollapse_all_real_times(self):
        floor=certify_sides_and_robust_boundary(self.P,Q(0))
        self.assertEqual(floor["checked_multiaffine_side_vertices"],12)
        self.assertEqual(SIDE_X_FLOOR,Q(74,100))
        self.assertEqual(SIDE_Y_FLOOR,Q(79,100))
        self.assertGreater(floor["certified_side_squared_floor"],1)

    def test_independent_unfiltered_whole_polygon_census_is_three_then_one(self):
        counts=_verify_independent_snapshots(self.P)
        self.assertEqual(counts,[(START,3),(END,1)])
        self.assertEqual(
            [s["total_geometric_square_orbits"] for s in
             self.report["independent_all_square_census_regular_times"]],[3,1])

    def test_degree_zero_does_not_force_a_square_to_exist(self):
        self.assertTrue(self.report["positive_square_collision_but_degree_zero"])
        self.assertIn("does NOT force a square",self.report["scope"])
        self.assertIn("not square existence",self.report["scope"])

    def test_box_or_margin_tampering_rejected(self):
        wrong=list(BOX)
        wrong[0]=(Q(91,100),Q(93,100))
        with self.assertRaises(ValueError):
            certify_cell(self.P,COLLIDING_LEFT,box=tuple(wrong))
        with self.assertRaises(ValueError):
            certify_cell(self.P,COLLIDING_RIGHT,gap=Q(1,1000))

    def test_noncanonical_geometry_rejected_before_degree_claim(self):
        modified=list(self.P)
        x,y=modified[1]
        modified[1]=(tuple(a+Q(1,10000) if i==0 else a for i,a in enumerate(x)),y)
        for fn,args in (
                (certify_cell,(tuple(modified),COLLIDING_LEFT)),
                (certify_interface_gluing,(tuple(modified),)),
                (certify_sides_and_robust_boundary,(tuple(modified),))):
            with self.subTest(fn=fn.__name__):
                with self.assertRaises(ValueError):
                    fn(*args)

    def test_no_global_all_time_total_orbit_count_claimed(self):
        self.assertIn("No general event atlas",self.report["scope"])
        self.assertIn("Toeplitz",self.report["scope"])


if __name__=="__main__":
    unittest.main()
