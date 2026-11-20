"""v268: exact cyclic/dihedral monodromy and two closed rectangle-loop tests."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from moving_polygon_events_v224 import (
    iszero,divide,evaluate,isolate_roots,polynomial,
)
from rectangle_monodromy_v268 import (
    P,Rat,ZERO,ONE,U,rectangle,rotate,unit_quarter,
    rectangle_defect,diagonal_lengths_squared,corner_radius_squared,
    diagonal_dot_character,even_moebius_loop,odd_twisted_loop,
    _monodromy,_pointset,_assert_rectangle_identically,
    _positive_rational_function,_nonnegative_rational_function,
    monodromy_no_go,audit,certify_circle_rectangle_loop,
)


class ExactRectangleMonodromy(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.even=even_moebius_loop()
        cls.odd=odd_twisted_loop()
        cls.all=audit()

    def test_exact_rational_function_cancellation(self):
        self.assertEqual((ONE+U)*(ONE-U),ONE-U*U)
        self.assertEqual((ONE+U)/(ONE+U),ONE)
        self.assertEqual((3-U)/(6),Rat((Q(1,2),Q(-1,6))))
        self.assertEqual((ONE+U).at(Q(1,3)),Q(4,3))

    def test_unacceptable_poles_rejected(self):
        with self.assertRaises(ZeroDivisionError):
            Rat.of(1)/ZERO
        with self.assertRaises(ZeroDivisionError):
            (ONE/(ONE-U)).at(1)

    def test_full_rational_quarterturn_is_exact_unit_circle(self):
        c,s=unit_quarter(U)
        self.assertEqual(c*c+s*s,ONE)
        self.assertEqual((c.at(0),s.at(0)),(Q(1),Q(0)))
        self.assertEqual((c.at(1),s.at(1)),(Q(0),Q(1)))

    def test_wrong_matrix_rejected_by_exact_rotation_identity(self):
        with self.assertRaises(ValueError):
            rotate(rectangle(3,4),ONE,ONE)

    def test_even_loop_endpoints_match_under_halfturn_T_squared(self):
        r=self.even
        self.assertEqual(r["cyclic_monodromy_T_power"],2)
        self.assertEqual(r["tensor_w1"],0)
        self.assertEqual(r["individual_chord_w1_values"],[1,1])

    def test_even_loop_has_no_square_ANYWHERE_on_both_continuous_arcs(self):
        r=self.even
        self.assertEqual(r["diagonal_dot_constant"],"28")
        self.assertEqual(r["nondegenerate_square_count_on_loop"],0)
        self.assertEqual(r["two_rational_circle_arcs"],2)
        self.assertTrue(r["all_real_parameter_identities_checked"])

    def test_even_loop_diagonals_fixed_length_and_corners_on_circle(self):
        a,b=Rat.of(3),Rat.of(4)
        for c,s in (unit_quarter(U),
                    (-unit_quarter(U)[1],unit_quarter(U)[0])):
            corners=rotate(rectangle(a,b),c,s)
            self.assertTrue(all(d==Rat.of(100)
                    for d in diagonal_lengths_squared(corners)))
            self.assertTrue(all(q==Rat.of(25)
                    for q in corner_radius_squared(corners)))
            self.assertEqual(
                rectangle_defect(corners),
                (ZERO,ZERO,ZERO,Rat.of(28)))

    def test_even_loop_composition_matches_midpoint(self):
        base=rectangle(3,4)
        c,s=unit_quarter(U)
        path0=rotate(base,c,s)
        path1=rotate(base,-s,c)
        self.assertEqual(_pointset(path0,1),_pointset(path1,0))
        self.assertEqual(
            _monodromy(_pointset(path0,0),_pointset(path1,1)),2)

    def test_odd_loop_has_genuine_quarterturn_T_monodromy(self):
        self.assertEqual(self.odd["cyclic_monodromy_T_power"],1)
        self.assertEqual(self.odd["tensor_w1"],1)
        self.assertFalse(self.odd["individual_diagonal_lines_globally_labeled"])
        self.assertEqual(
            self.odd["diagonal_dot_endpoints"],["28","-28"])

    def test_odd_rectangle_loop_is_full_time_noncollapsed(self):
        r=self.odd
        self.assertEqual(r["semiwidth_endpoints"],["3","4"])
        self.assertEqual(r["semiheight_endpoints"],["4","3"])
        self.assertEqual(r["minimum_rectangle_side_squared"],"36")
        self.assertEqual(r["fixed_circumcircle_radius_squared"],"25")
        self.assertEqual(r["square_side_squared_at_event"],"50")

    def test_odd_event_exact_minimal_polynomial(self):
        self.assertEqual(self.odd["event_minimal_polynomial"],
                         ["9","-18","1"])
        q=polynomial((9,-18,1))
        self.assertGreater(evaluate(q,Q(0)),0)
        self.assertLess(evaluate(q,Q(1)),0)
        roots=isolate_roots(q,Q(0),Q(1),width=Q(1,2**22))
        self.assertEqual(len(roots),1)
        self.assertFalse(roots[0].rational)

    def test_odd_event_has_certified_unique_irrational_zero(self):
        lo,hi=map(Q,self.odd["event_isolation_bracket"])
        self.assertLess(Q(0),lo)
        self.assertLess(lo,hi)
        self.assertLess(hi,Q(1))
        self.assertLess(hi-lo,Q(1,2**20))
        q=polynomial((9,-18,1))
        self.assertLess(evaluate(q,lo)*evaluate(q,hi),0)
        self.assertEqual(self.odd["number_of_distinct_square_crossings"],1)

    def test_odd_loop_exact_analytic_continuous_family(self):
        p=(3-U)/6
        a=5*(1-p*p)/(1+p*p)
        b=10*p/(1+p*p)
        self.assertEqual(a*a+b*b,Rat.of(25))
        for r in (a,b):
            self.assertTrue(_positive_rational_function(r))
            self.assertTrue(_nonnegative_rational_function(r-3))
            self.assertTrue(_nonnegative_rational_function(4-r))
        c,s=unit_quarter(U)
        corners=rotate(rectangle(a,b),c,s)
        dot=_assert_rectangle_identically(corners)
        self.assertEqual(dot,4*(b*b-a*a))
        self.assertEqual(_monodromy(
            _pointset(corners,0),_pointset(corners,1)),1)

    def test_even_vs_odd_monodromy_character(self):
        self.assertEqual(
            [diagonal_dot_character(k) for k in range(4)],
            [1,-1,1,-1])
        self.assertEqual(self.even["tensor_w1"],0)
        self.assertEqual(self.odd["tensor_w1"],1)

    def test_noninteger_monodromy_rejected(self):
        for x in (Q(1,2),1.0,"1",True):
            with self.subTest(x=x):
                with self.assertRaises(ValueError):
                    diagonal_dot_character(x)

    def test_all_four_corner_points_must_be_distinct(self):
        s=((Q(0),Q(0)),(Q(1),Q(0)),
           (Q(1),Q(0)),(Q(0),Q(1)))
        with self.assertRaises(ValueError):
            _monodromy(s,s)

    def test_wrong_endpoint_identification_rejected(self):
        s=((Q(3),Q(4)),(-Q(3),Q(4)),
           (-Q(3),-Q(4)),(Q(3),-Q(4)))
        t=list(s);t[0]=(Q(4),Q(4))
        with self.assertRaises(ValueError):
            _monodromy(s,tuple(t))

    def test_generic_checker_accepts_static_non_square_loop(self):
        c=certify_circle_rectangle_loop(
            (rectangle(3,4),),Q(25),Q(30))
        self.assertEqual(c["geometric_C4_monodromy_T_power"],0)
        self.assertEqual(c["distinct_isolated_square_event_count"],0)
        self.assertEqual(c["all_side_inequalities_certified"],4)

    def test_generic_checker_validates_both_example_loops(self):
        e=self.even["reusable_piecewise_path_audit"]
        o=self.odd["reusable_piecewise_path_audit"]
        self.assertEqual(e["number_of_rational_patches"],2)
        self.assertEqual(e["all_side_inequalities_certified"],8)
        self.assertEqual(e["distinct_isolated_square_event_count"],0)
        self.assertEqual(o["number_of_rational_patches"],1)
        self.assertEqual(o["all_side_inequalities_certified"],4)
        self.assertEqual(o["distinct_isolated_square_event_count"],1)
        self.assertTrue(o["nondegenerate_at_every_real_parameter"])

    def test_wrong_circle_radius_rejected(self):
        with self.assertRaises(ValueError):
            certify_circle_rectangle_loop(
                (rectangle(3,4),),Q(24),Q(30))

    def test_false_size_floor_rejected(self):
        with self.assertRaises(ValueError):
            certify_circle_rectangle_loop(
                (rectangle(3,4),),Q(25),Q(40))

    def test_square_continuum_is_not_a_finite_event_ledger(self):
        with self.assertRaises(ValueError):
            certify_circle_rectangle_loop(
                (rectangle(3,3),),Q(18),Q(10))

    def test_degenerate_flat_endpoint_rectangle_is_rejected(self):
        p=ONE-U/2
        a=5*(ONE-p*p)/(ONE+p*p)
        b=10*p/(ONE+p*p)
        with self.assertRaises(ValueError):
            certify_circle_rectangle_loop(
                (rectangle(a,b),),Q(25),Q(1))

    def test_unclosed_rectangle_path_is_rejected(self):
        c,s=unit_quarter(U)
        with self.assertRaises(ValueError):
            certify_circle_rectangle_loop(
                (rotate(rectangle(3,4),c,s),),Q(25),Q(30))

    def test_piecewise_discontinuous_rectangle_path_is_rejected(self):
        c,s=unit_quarter(U)
        with self.assertRaises(ValueError):
            certify_circle_rectangle_loop(
                (rectangle(3,4),rotate(rectangle(3,4),-s,c)),
                Q(25),Q(30))

    def test_wrong_cyclic_orientation_is_rejected(self):
        with self.assertRaises(ValueError):
            certify_circle_rectangle_loop(
                (rectangle(-3,4),),Q(25),Q(30))

    def test_odd_component_no_square_is_topologically_impossible(self):
        r=monodromy_no_go()
        self.assertTrue(r[
            "odd_monodromy_forces_nonzero_size_square"])
        self.assertTrue(r[
            "nonsquare_closed_loop_must_have_even_monodromy"])
        self.assertFalse(r[
            "existence_of_odd_rectangle_loop_for_arbitrary_Jordan_curve_proved"])

    def test_automatic_report_rejects_universal_toeplitz_claim(self):
        self.assertTrue(self.all["not_an_unrestricted_Toeplitz_proof"])
        self.assertEqual(
            self.all["even_no_square_rectangle_loop"][
                "nondegenerate_square_count_on_loop"],0)
        self.assertEqual(
            self.all["odd_algebraic_square_forcing_rectangle_loop"][
                "number_of_distinct_square_crossings"],1)


if __name__=="__main__":
    unittest.main()
