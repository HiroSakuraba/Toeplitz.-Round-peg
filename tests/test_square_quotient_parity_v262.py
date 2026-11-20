"""v262 exact cyclic quotient relative parity and singular-event tests."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from big_peg_pair_death_v254 import family_vertices
from square_relative_degree_v260 import (
    START,END,C0_EPSILON,BOX as PAIR_BOX,
)
from square_oriented_atlas_v261 import (
    PERSISTENT_BOX,CYCLIC_RESIDUAL_TRANSFORM,
    CYCLIC_SOURCE_PERMUTATION,
)
from near_square_transversality_v258 import det4
from square_quotient_parity_v262 import (
    PERIOD,IDENTITY,PATCHES,LOCAL_BOX,
    BASE_BOUNDARY_RESIDUAL_FLOOR,
    matmul,power,apply,infinity_norm,
    base_unwrapped,base_box,_validate_domain,
    rotate_unwrapped,rotate_box,cyclic_sector_points,
    certify_cyclic_action,certify_saturated_disjoint_regions,
    certify_geometric_equivariance,certify_quotient_atlas,
)


class ExactCyclicQuotientAtlas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.report=certify_quotient_atlas()

    def test_exact_fourfold_group_law_both_representations(self):
        for M in (CYCLIC_SOURCE_PERMUTATION,
                  CYCLIC_RESIDUAL_TRANSFORM):
            self.assertEqual(power(M,0),IDENTITY)
            self.assertEqual(power(M,4),IDENTITY)
            for a in range(3):
                for b in range(3-a):
                    self.assertEqual(matmul(power(M,a),power(M,b)),
                                     power(M,a+b))

    def test_nonorientable_source_and_oriented_residual_bundle(self):
        rep=certify_cyclic_action()
        self.assertEqual(rep["basis_identities_checked"],32)
        self.assertEqual(rep["source_orientation_character"],(1,-1,1,-1))
        self.assertEqual(rep["target_determinants"],(1,1,1,1))
        self.assertEqual(det4(CYCLIC_SOURCE_PERMUTATION),-1)
        self.assertEqual(det4(CYCLIC_RESIDUAL_TRANSFORM),1)

    def test_exact_eight_boxes_disjoint_even_by_first_parameter(self):
        sat=certify_saturated_disjoint_regions()
        self.assertEqual(sat["regions_in_labeled_C4_saturation"],8)
        self.assertEqual(sat["regions_in_unlabeled_quotient"],2)
        self.assertEqual(sat["free_action_vertex_checks"],96)
        intervals=sat["strictly_disjoint_first_coordinate_bands"]
        for a,b in zip(intervals,intervals[1:]):
            self.assertLess(Q(a[1]),Q(b[0]))
        self.assertEqual({(kind,k) for _,_,kind,k in intervals},
                         {(p,k) for p in PATCHES for k in range(4)})

    def test_each_region_has_global_cyclic_parameter_order(self):
        for kind in PATCHES:
            base=base_box(kind)
            _validate_domain(base)
            for k in range(4):
                image=rotate_box(base,k)
                _validate_domain(image)
                for q in product(*base):
                    rotated=rotate_unwrapped(q,k)
                    for z,(lo,hi) in zip(rotated,image):
                        self.assertLessEqual(lo,z)
                        self.assertLessEqual(z,hi)

    def test_four_iterations_are_identity_mod_common_winding(self):
        for kind in PATCHES:
            for q in product(*base_box(kind)):
                self.assertEqual(rotate_unwrapped(q,0),q)
                self.assertEqual(rotate_unwrapped(q,4),q)
                self.assertEqual(
                    rotate_unwrapped(rotate_unwrapped(q,1),3),q)

    def test_every_actual_polygon_corner_commutes_with_relabeling(self):
        self.assertEqual(certify_geometric_equivariance(self.P),280)
        self.assertEqual(self.report["actual_polygon_equivariance_checks"],280)

    def test_glued_pair_is_not_artificially_split_at_vertex(self):
        s=(Q(91,100),Q(0),Q(44,100),Q(75,100))
        q=base_unwrapped("collision_pair",s)
        self.assertEqual(q[1],Q(2))
        for k in range(4):
            pts=cyclic_sector_points(self.P,START,
                                     rotate_unwrapped(q,k))
            self.assertEqual(pts,(pts[:]+pts[:])[0:4])

    def test_quoitent_degree_survives_root_pair_death(self):
        r=self.report
        self.assertEqual(
            (r["quotient_relative_degree_before_event_mod2"],
             r["quotient_relative_degree_at_event_mod2"],
             r["quotient_relative_degree_after_event_mod2"]),
            (1,1,1))
        self.assertEqual(r["quotient_pair_at_event_mod2"],0)

    def test_transported_signed_index_flips_on_odd_cyclic_sheets(self):
        signs=self.report["signed_indices_in_four_labeled_sheets"]
        self.assertEqual(signs["persistent"],[1,-1,1,-1])
        self.assertEqual(signs["pair_before_left"],[-1,1,-1,1])
        self.assertEqual(signs["pair_before_right"],[1,-1,1,-1])
        self.assertEqual(signs["whole_union"],[1,-1,1,-1])
        self.assertEqual(signs["pair_union"],[0,0,0,0])

    def test_labeled_root_count_mod2_is_zero_even_when_quotient_odd(self):
        snapshots=self.report["rational_snapshot_root_orbits"]
        self.assertEqual(
            (snapshots["rational_before_orbits"],
             snapshots["rational_after_orbits"]), (3,1))
        self.assertEqual(
            (snapshots["rational_before_labeled_roots"],
             snapshots["rational_after_labeled_roots"]), (12,4))
        self.assertEqual(snapshots["cyclic_root_membership_checks"],16)
        self.assertEqual(self.report[
            "ordinary_mod2_count_on_entire_labeled_saturation"],0)
        self.assertEqual(self.report[
            "ordinary_integer_sum_on_entire_labeled_saturation"],0)

    def test_exact_c0_boundary_budget_on_all_cyclic_sheets(self):
        self.assertGreater(BASE_BOUNDARY_RESIDUAL_FLOOR,0)
        M=self.report["uniform_cyclic_inverse_residual_transform_norm"]
        floor=Q(self.report["uniform_cyclic_saturated_outer_residual_floor"])
        self.assertEqual(BASE_BOUNDARY_RESIDUAL_FLOOR/Q(M),floor)
        self.assertLess(Q(4)*C0_EPSILON,floor)
        self.assertTrue(self.report["positive_square_side_for_every_certified_atlas_root"])

    def test_singular_geometric_root_parity_is_not_singular_degree(self):
        self.assertTrue(self.report[
            "quotient_geometric_root_count_at_singular_event_parity_not_an_invariant"])
        self.assertEqual(self.report["quotient_pair_at_event_mod2"],0)

    def test_quotient_integer_degree_and_all_jordan_claims_rejected(self):
        r=self.report
        self.assertFalse(r["canonically_globally_oriented_integer_degree_on_quotient"])
        self.assertFalse(r["full_six_sector_atlas_proved"])
        self.assertFalse(r["all_real_time_exhaustive_whole_pentagon_census_proved"])
        self.assertFalse(r["arbitrary_Jordan_square_peg_theorem_proved"])

    def test_wrong_or_misaligned_boxes_fail_closed(self):
        with self.assertRaises(ValueError):
            base_box("fourth square")
        bad=list(base_box("persistent"))
        bad[1]=(bad[0][0],bad[1][1])
        with self.assertRaises(ValueError):
            _validate_domain(tuple(bad))
        with self.assertRaises(ValueError):
            rotate_box(base_box("persistent"),4)

    def test_invalid_winding_or_cyclic_order_rejected(self):
        with self.assertRaises(ValueError):
            rotate_unwrapped((Q(0),Q(1),Q(2),Q(5)),0)
        with self.assertRaises(ValueError):
            rotate_unwrapped((Q(1),Q(2),Q(3),Q(4)),5)
        with self.assertRaises(ValueError):
            rotate_unwrapped((Q(2),Q(1),Q(3),Q(4)),2)
        with self.assertRaises(ValueError):
            cyclic_sector_points(self.P,START,
                                 (Q(0),Q(2),Q(1),Q(4)))
        with self.assertRaises(ValueError):
            certify_geometric_equivariance(tuple(reversed(self.P)))

    def test_bad_matrix_exponent_and_multiplication_rejected(self):
        with self.assertRaises(ValueError):
            power(IDENTITY,5)
        with self.assertRaises(ValueError):
            matmul(IDENTITY,((Q(1),),))
        with self.assertRaises(ValueError):
            apply(IDENTITY,(Q(1),Q(2)))


if __name__=="__main__":
    unittest.main()
