"""v261 exact oriented atlas and relabeling obstruction regressions."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from big_peg_pair_death_v254 import family_vertices
from polygon_orbit_census_v223 import add,sub,mul
from square_degree_cube_v259 import square_defect
from square_relative_degree_v260 import (
    START,END,BOX as PAIR_BOX,_square_points,
    BOUNDARY_RESIDUAL_FLOOR as PAIR_RESIDUAL_FLOOR,
    C0_EPSILON,
)
from square_oriented_atlas_v261 import (
    PERSISTENT_BOX,PERSISTENT_GAP,PERSISTENT_DET_FLOOR,
    PERSISTENT_ENTRY_BOUND,PERSISTENT_INVERSE_BOUND,
    PERSISTENT_RESIDUAL_FLOOR,
    CYCLIC_RESIDUAL_TRANSFORM,CYCLIC_SOURCE_PERMUTATION,
    SIDE_X_FLOOR,SIDE_Y_FLOOR,
    persistent_points,_unwrapped_q_persistent,_unwrapped_q_pair,
    global_boundary_point,certify_persistent_box,
    certify_cyclic_relabeling,audit,
)
from near_square_transversality_v258 import det4


class OrientedMultiRegionAtlas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.certificate=certify_persistent_box(cls.P)
        cls.report=audit()

    def test_all_real_times_persistent_index_is_positive(self):
        self.assertEqual(self.certificate.orientation_sign,+1)
        self.assertGreaterEqual(self.certificate.parameter_sign_checks,9)
        self.assertEqual(self.report["persistent_index"],+1)

    def test_exact_rational_persistent_box_margin(self):
        self.assertEqual(PERSISTENT_GAP,Q(1,1000))
        for root in self.certificate.exact_endpoint_roots:
            for x,(lo,hi) in zip(root,PERSISTENT_BOX):
                self.assertLess(lo+PERSISTENT_GAP,x)
                self.assertLess(x,hi-PERSISTENT_GAP)

    def test_persistent_exact_square_at_both_endpoints(self):
        for t,root in zip((START,END),
                          self.certificate.exact_endpoint_roots):
            p=persistent_points(self.P,t,root)
            self.assertEqual(square_defect(p),(Q(0),)*4)

    def test_persistent_square_uses_actual_ordered_boundary_arcs(self):
        for t,root in zip((START,END),
                          self.certificate.exact_endpoint_roots):
            q=_unwrapped_q_persistent(root)
            self.assertTrue(Q(0)<q[0]<q[1]<q[2]<q[3]<Q(5))
            direct=tuple(global_boundary_point(self.P,t,v) for v in q)
            self.assertEqual(direct,persistent_points(self.P,t,root))

    def test_disappearing_pair_is_a_single_genuine_vertex_glued_chart(self):
        for t in (START,END):
            for h in (-Q(1,1000),Q(0),Q(1,1000)):
                s=(Q(91,100),h,Q(44,100),Q(75,100))
                q=_unwrapped_q_pair(s)
                self.assertTrue(Q(0)<q[0]<q[1]<q[2]<q[3]<Q(5))
                self.assertEqual(
                    tuple(global_boundary_point(self.P,t,v) for v in q),
                    _square_points(self.P,t,s))

    def test_persistent_and_pair_boxes_disjoint_in_global_coordinates(self):
        self.assertLess(PERSISTENT_BOX[0][1],PAIR_BOX[0][0])
        self.assertLess(Q(2)+PERSISTENT_BOX[2][1],
                        Q(3)+PAIR_BOX[2][0])

    def test_same_orientation_on_all_local_edge_parameters(self):
        # All q charts are translations of their four local parameters.
        # The pair's h translates to q1=2+h on BOTH sides.
        for u in self.certificate.exact_endpoint_roots:
            q=_unwrapped_q_persistent(u)
            self.assertEqual(tuple(q[i]-u[i] for i in range(4)),
                             (Q(0),Q(1),Q(2),Q(4)))
        for h in (-Q(1,2000),Q(1,2000)):
            s=(Q(91,100),h,Q(44,100),Q(75,100))
            q=_unwrapped_q_pair(s)
            self.assertEqual(tuple(q[i]-s[i] for i in range(4)),
                             (Q(0),Q(2),Q(3),Q(4)))

    def test_exact_det_and_entry_floor(self):
        self.assertEqual(PERSISTENT_DET_FLOOR,Q(1,5))
        self.assertEqual(PERSISTENT_ENTRY_BOUND,Q(2))
        self.assertEqual(PERSISTENT_INVERSE_BOUND,Q(960))
        self.assertEqual(PERSISTENT_RESIDUAL_FLOOR,Q(1,960000))
        self.assertEqual(self.certificate.inverse_bound,Q(960))

    def test_arbitrary_small_continuous_perturbation_keeps_plus_one(self):
        self.assertEqual(C0_EPSILON,Q(1,10**9))
        self.assertLess(Q(4)*C0_EPSILON,
                        min(PERSISTENT_RESIDUAL_FLOOR,
                            PAIR_RESIDUAL_FLOOR))
        self.assertEqual(self.report["oriented_disconnected_union_degree"],1)
        self.assertGreater(self.certificate.side_squared_floor,Q(1))

    def test_size_floor_on_persistent_box(self):
        self.assertEqual(SIDE_X_FLOOR,Q(1))
        self.assertEqual(SIDE_Y_FLOOR,Q(9,100))
        self.assertGreater(
            (SIDE_X_FLOOR-2*C0_EPSILON)**2+
            (SIDE_Y_FLOOR-2*C0_EPSILON)**2,Q(1))

    def test_two_disappearing_indices_cancel_in_bigger_atlas(self):
        self.assertEqual(self.report["colliding_left_index"],-1)
        self.assertEqual(self.report["colliding_right_index"],+1)
        self.assertEqual(self.report["pair_glued_index"],0)
        self.assertEqual(self.report["oriented_disconnected_union_degree"],1)

    def test_independent_complete_census_at_rational_slices_only(self):
        counts=self.report[
            "independent_whole_polygon_census_at_two_rational_slices"]
        self.assertEqual([r["total_geometric_square_orbits"] for r in counts],
                         [3,1])
        self.assertFalse(
            self.report["complete_all_real_time_global_square_census_proved"])

    def test_cyclic_relabeling_linear_identity_and_orientation(self):
        c=certify_cyclic_relabeling()
        self.assertEqual(det4(CYCLIC_SOURCE_PERMUTATION),-1)
        self.assertEqual(det4(CYCLIC_RESIDUAL_TRANSFORM),1)
        self.assertEqual(c["local_degree_relabeling_factor"],-1)
        self.assertEqual(c["number_of_exact_basis_checks"],8)

    def test_unlabeled_global_integer_degree_not_asserted(self):
        self.assertFalse(
            self.report["canonical_unlabeled_integer_degree_proved"])
        self.assertIn("NOT a canonical global integer degree",
                      self.report["scope"])
        self.assertIn("reverses sign",self.report["scope"])

    def test_wrong_pentagon_fails_closed(self):
        P=list(self.P)
        x,y=P[1]
        P[1]=(tuple(a+Q(1,10000) if j==0 else a
                    for j,a in enumerate(x)),y)
        with self.assertRaises(ValueError):
            certify_persistent_box(tuple(P))

    def test_wrong_region_or_gap_fails_closed(self):
        other=list(PERSISTENT_BOX)
        other[0]=(Q(15,100),Q(17,100))
        with self.assertRaises(ValueError):
            certify_persistent_box(self.P,box=tuple(other))
        with self.assertRaises(ValueError):
            certify_persistent_box(self.P,gap=Q(1,100))

    def test_critical_perturbation_fails_closed(self):
        critical=PERSISTENT_RESIDUAL_FLOOR/Q(4)
        with self.assertRaises(ValueError):
            certify_persistent_box(self.P,epsilon=critical)
        with self.assertRaises(ValueError):
            certify_persistent_box(self.P,epsilon=critical+Q(1,10**9))

    def test_invalid_numeric_perturbation_fails_closed(self):
        for bad in (-Q(1,10**9),1e-9,None,"1/1000000000"):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    certify_persistent_box(self.P,epsilon=bad)

    def test_outside_global_unwrapped_parameter_is_rejected(self):
        for q in (-Q(1,10),Q(5),Q(51,10)):
            with self.subTest(q=q):
                with self.assertRaises(ValueError):
                    global_boundary_point(self.P,START,q)


if __name__=="__main__":
    unittest.main()
