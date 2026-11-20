"""v267: exact secant-spin and diagonal square forcing / no-go checks."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt,add,mul
from collision_quotient_v219 import notch_square_side_squared
from configuration_symmetry_v218 import (
    CYCLIC_CHART_JACOBIAN,residual_cyclic_matrix,determinant,
)
from secant_spin_rectangle_v267 import (
    SecantAuditInconclusive,secant,exact_secant_winding,
    _breakpoints,_quadratic_unisolvent_points,_source_sign,
    _square_pairing_residual,all_pairings,is_positive_square_diagonals,
    ROT,REFLECT,ROT_TARGET,REFLECT_TARGET,
    certify_diagonal_D4_equivariance,
    rectangle_family,certify_rectangle_orthogonality_crossing,
    notch_shrinkout_comparison,audit,
)


SQUARE=(pt(0,0),pt(2,0),pt(2,2),pt(0,2))
TRIANGLE=(pt(0,0),pt(3,0),pt(1,2))
NONCONVEX=(pt(0,0),pt(4,0),pt(4,4),pt(3,4),
           pt(3,1),pt(1,1),pt(1,4),pt(0,4))


class SecantSpinAudit(unittest.TestCase):
    def test_antipodal_reversal_is_exact_at_all_segment_boundaries(self):
        for P in (SQUARE,TRIANGLE,NONCONVEX):
            for t in _breakpoints(P):
                self.assertEqual(secant(P,t+Q(1,2)),
                                 mul(secant(P,t),Q(-1)))

    def test_exact_square_secant_winding_degree_one(self):
        result=exact_secant_winding(SQUARE)
        self.assertEqual(result["secant_winding"],1)
        self.assertTrue(result["odd_degree"])
        self.assertEqual(
            result["minimum_antipodal_secant_length_squared"],Q(4))

    def test_clockwise_square_secant_winding_degree_minus_one(self):
        self.assertEqual(exact_secant_winding(tuple(reversed(SQUARE)))[
            "secant_winding"],-1)

    def test_triangle_winding_odd(self):
        self.assertEqual(exact_secant_winding(TRIANGLE)["secant_winding"],1)

    def test_nonconvex_embedded_curve_winding_is_still_odd(self):
        self.assertEqual(exact_secant_winding(NONCONVEX)[
            "secant_winding"],1)

    def test_nonconvex_reversed_curve_flips_winding(self):
        self.assertEqual(exact_secant_winding(
            tuple(reversed(NONCONVEX)))["secant_winding"],-1)

    def test_secant_nonzero_length_is_explicit_and_curve_dependent(self):
        for P in (SQUARE,TRIANGLE,NONCONVEX):
            a=exact_secant_winding(P)
            self.assertGreater(
                a["minimum_antipodal_secant_length_squared"],0)
            self.assertEqual(
                a["uniform_coordinate_displacement_epsilon_squared_budget"],
                a["minimum_antipodal_secant_length_squared"]/8)

    def test_invalid_crossed_polygon_is_not_silently_counted(self):
        crossed=(pt(0,0),pt(2,2),pt(2,0),pt(0,2))
        with self.assertRaises(ValueError):
            exact_secant_winding(crossed)

    def test_mobius_chord_line_has_odd_orientation_monodromy(self):
        a=exact_secant_winding(SQUARE)
        self.assertEqual(
            a["w1_of_unoriented_chord_orientation_line_on_core"],1)
        self.assertEqual(a["secant_winding"]%2,1)

    def test_complete_three_pairings_do_not_assume_graceful_sector(self):
        pairing=all_pairings(SQUARE)
        self.assertEqual(len(pairing),3)
        self.assertEqual(
            sum(f==(Q(0),)*4 for f in pairing.values()),1)
        self.assertTrue(is_positive_square_diagonals(SQUARE))
        for reorder in ((0,2,1,3),(1,3,0,2),(3,0,2,1)):
            P=tuple(SQUARE[i] for i in reorder)
            self.assertTrue(is_positive_square_diagonals(P))

    def test_collapsed_four_point_solution_is_rejected_as_square(self):
        P=(pt(0,0),)*4
        self.assertEqual(
            _square_pairing_residual(P),(Q(0),)*4)
        self.assertFalse(is_positive_square_diagonals(P))

    def test_rectangle_not_square_has_exact_orthogonality_obstruction(self):
        P=rectangle_family(Q(3,2))
        f=_square_pairing_residual(P)
        self.assertEqual(f[:3],(Q(0),)*3)
        self.assertNotEqual(f[3],0)
        self.assertFalse(is_positive_square_diagonals(P))

    def test_rectangle_orthogonality_sign_crossing_forces_positive_square(self):
        r=certify_rectangle_orthogonality_crossing(Q(1,2),Q(3,2))
        self.assertEqual(r["orthogonality_at_ends"],(Q(3),Q(-5)))
        self.assertEqual(r["exact_square_at_parameter"],Q(1))
        self.assertTrue(r["positive_size_protected_on_chosen_family"])

    def test_rectangle_path_without_sign_crossing_does_not_force_square(self):
        with self.assertRaises(SecantAuditInconclusive):
            certify_rectangle_orthogonality_crossing(Q(3,2),Q(2))
        with self.assertRaises(ValueError):
            certify_rectangle_orthogonality_crossing(Q(0),Q(2))

    def test_universal_exact_quadratic_D4_action(self):
        r=certify_diagonal_D4_equivariance()
        self.assertEqual(r["total_exact_polynomial_identity_probes"],90)
        self.assertEqual(len(tuple(_quadratic_unisolvent_points())),45)
        self.assertEqual(ROT_TARGET,(-1,-1,-1,-1))
        self.assertEqual(REFLECT_TARGET,(1,1,1,-1))

    def test_cyclic_rotation_still_flips_integer_local_index(self):
        d=certify_diagonal_D4_equivariance()
        self.assertEqual(d["cyclic_rotation"],{
            "source_orientation":-1,
            "target_orientation":1,
            "local_index_transport":-1,
        })
        self.assertEqual(
            determinant(CYCLIC_CHART_JACOBIAN),-1)
        self.assertEqual(
            determinant(residual_cyclic_matrix(1)),1)

    def test_chord_endpoint_reversal_has_matching_orientation_twists(self):
        d=certify_diagonal_D4_equivariance()
        self.assertEqual(d["diagonal_endpoint_reversal"],{
            "source_orientation":-1,
            "target_orientation":-1,
            "local_index_transport":1,
        })
        self.assertEqual(_source_sign(REFLECT),-1)

    def test_notch_shrinking_square_survives_odd_winding(self):
        for n in (2,3,4,8,64,128):
            r=notch_shrinkout_comparison(n)
            self.assertEqual(r["secant_winding"],1)
            self.assertEqual(r["chosen_square_side_squared"],Q(1,n*n))
            self.assertEqual(r["chosen_square_side_squared"],
                             notch_square_side_squared(n))
            self.assertTrue(r[
                "no_uniform_square_size_floor_deduced_from_winding"])

    def test_chosen_square_can_be_far_smaller_than_secant_length(self):
        r=notch_shrinkout_comparison(128)
        self.assertGreater(r["minimum_antipodal_chord_squared"],0)
        self.assertLess(r["chosen_square_side_squared"],
                        r["minimum_antipodal_chord_squared"]/100)

    def test_all_conclusions_explicitly_stop_before_global_toeplitz(self):
        r=audit()
        self.assertTrue(r["ccw_square_secant"]["odd_degree"])
        self.assertFalse(
            r["unconditional_all_Jordan_square_existence_inferred"])
        self.assertIn("nondegenerate inscribed rectangles",
                      r["actionable_conditional_theorem"])


if __name__=="__main__":
    unittest.main()
