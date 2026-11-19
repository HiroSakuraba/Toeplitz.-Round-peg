"""v253: exact relative cutoff parity and global-max bottleneck regression tests."""
import sys
from fractions import Fraction as Q
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import census_polygon
from big_peg_multichart_v249 import (
    pol,pmul,padd,psub,pneg,peval,certify_sign,ppoint_at,
)
from big_peg_corner_v250 import family_vertices
from big_peg_size_atlas_v251 import (
    BRANCH_A,BRANCH_B,WINDOW_START,WINDOW_END,HIGH_S2,
    side_comparison_polynomial,
)
from big_peg_cutoff_ledger_v252 import (
    B_ENTER_BRACKET,A_EXIT_BRACKET,cutoff_numerator,
)
from big_peg_relative_parity_v253 import (
    CROSS_LO,CROSS_HI,BOTTLE_S2_LO,BOTTLE_S2_HI,
    _squared_side,_side_derivative_sign_polynomial,
    certify_global_envelope,rational_time_state,audit,
)


class RelativeSquareCountAndBottleneck(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.report=audit()
        cls.envelope=certify_global_envelope(cls.P)

    def test_full_126_chart_prior_ledger_rechecked(self):
        self.assertEqual(self.report["complete_high_square_charts_categorized"],126)
        self.assertEqual(self.report["number_of_possible_high_square_branches"],2)
        self.assertIn("EVERY real",self.report["global_maximum_square_side_squared_identity"])

    def test_global_size_envelope_is_exact_max_of_two_squares(self):
        for time in (Q(1,2),Q(500545,10**6),
                     Q(10011,20000),Q(501,1000)):
            with self.subTest(time=time):
                a=_squared_side(self.P,BRANCH_A,time)
                b=_squared_side(self.P,BRANCH_B,time)
                s=rational_time_state(self.P,time)
                self.assertEqual(s.largest_square_side_squared,max(a,b))
                self.assertGreater(s.largest_square_side_squared,HIGH_S2)

    def test_both_branch_derivatives_strictly_opposite_for_all_times(self):
        self.assertGreater(self.envelope.derivative_sign_checks_A,0)
        self.assertGreater(self.envelope.derivative_sign_checks_B,0)
        A=_side_derivative_sign_polynomial(self.P,BRANCH_A)
        B=_side_derivative_sign_polynomial(self.P,BRANCH_B)
        self.assertGreater(certify_sign(pneg(A),WINDOW_START,WINDOW_END,
                                        strict=True),0)
        self.assertGreater(certify_sign(B,WINDOW_START,WINDOW_END,
                                        strict=True),0)

    def test_bottleneck_time_is_unique_equal_size_root(self):
        self.assertEqual(CROSS_LO,Q(50054929969,10**11))
        self.assertEqual(CROSS_HI,Q(50054929970,10**11))
        self.assertEqual(CROSS_HI-CROSS_LO,Q(1,10**11))
        G=side_comparison_polynomial(self.P,BRANCH_A,BRANCH_B)
        self.assertGreater(peval(G,CROSS_LO),0)
        self.assertLess(peval(G,CROSS_HI),0)

    def test_tight_interval_contains_true_global_bottleneck(self):
        self.assertEqual(BOTTLE_S2_LO,Q(94257027585,10**10))
        self.assertEqual(BOTTLE_S2_HI,Q(94257027586,10**10))
        self.assertEqual(BOTTLE_S2_HI-BOTTLE_S2_LO,Q(1,10**10))
        self.assertGreater(BOTTLE_S2_LO,HIGH_S2)
        A_hi=_squared_side(self.P,BRANCH_A,CROSS_HI)
        B_lo=_squared_side(self.P,BRANCH_B,CROSS_LO)
        A_lo=_squared_side(self.P,BRANCH_A,CROSS_LO)
        B_hi=_squared_side(self.P,BRANCH_B,CROSS_HI)
        self.assertLess(BOTTLE_S2_LO,max(A_hi,B_lo))
        self.assertLessEqual(max(A_hi,B_lo),min(A_lo,B_hi))
        self.assertLess(min(A_lo,B_hi),BOTTLE_S2_HI)

    def test_global_bottleneck_above_previous_cutoff(self):
        self.assertGreater(BOTTLE_S2_LO-HIGH_S2,Q(1,10**6))
        # Exact lower bound on J: r_in<=3 for ALL family members,
        # so J² >= s_max²/9 > BOTTLE_S2_LO/9.
        self.assertGreater(BOTTLE_S2_LO,Q(102337,10**5)**2 * 9)

    def test_four_representative_time_states(self):
        expected=[
            (Q(1,2),1,0,0,0,1),
            (Q(500545,10**6),2,1,0,1,1),
            (Q(10011,20000),2,1,0,1,1),
            (Q(501,1000),1,1,1,0,1),
        ]
        for time,count,entered,exited,flux,index in expected:
            with self.subTest(time=time):
                state=rational_time_state(self.P,time)
                self.assertEqual(
                    (state.raw_large_square_count,state.entered_cutoff_count,
                     state.exited_cutoff_count,state.net_boundary_flux,
                     state.corrected_integer_index),
                    (count,entered,exited,flux,index))
                self.assertEqual(state.corrected_mod2_index,1)

    def test_every_rational_query_between_event_brackets_uses_exact_signs(self):
        for numerator in range(500510,500601,3):
            time=Q(numerator,10**6)
            state=rational_time_state(self.P,time)
            self.assertEqual(state.corrected_integer_index,1)
            self.assertEqual(
                state.raw_large_square_count-state.net_boundary_flux,1)

    def test_signed_event_count_matches_two_exact_cutoff_polynomials(self):
        for time in (Q(1,2),Q(500545,10**6),
                     Q(10011,20000),Q(501,1000)):
            state=rational_time_state(self.P,time)
            signs=[]
            for assignment in (BRANCH_A,BRANCH_B):
                from big_peg_multichart_v249 import chart_polynomials
                D,N=chart_polynomials(self.P,assignment)
                signs.append(peval(cutoff_numerator(D,N,HIGH_S2),time))
            self.assertEqual(state.exited_cutoff_count,int(signs[0]<=0))
            self.assertEqual(state.entered_cutoff_count,int(signs[1]>0))

    def test_event_time_conventions_in_strict_filter(self):
        # The algebraic roots are not assumed rational. Monotonicity
        # determines how equality belongs to the filtered count.
        # At B entry F_B=0 -> not counted; at A exit F_A=0 -> excluded.
        self.assertEqual(self.report["integer_signed_boundary_event_coefficients"],
                         {"B_enters_strict_cutoff":1,"A_exits_strict_cutoff":-1})
        self.assertIn("entry_before",self.report["corrected_integer_count_for_ALL_real_times"])
        self.assertIn("exit_by",self.report["corrected_integer_count_for_ALL_real_times"])

    def test_cutoff_events_precede_and_follow_envelope_bottleneck(self):
        self.assertTrue(B_ENTER_BRACKET[1]<CROSS_LO<CROSS_HI<A_EXIT_BRACKET[0])
        self.assertEqual(self.report["global_minimizer_time_bracket"],
                         [str(CROSS_LO),str(CROSS_HI)])

    def test_size_equality_event_does_not_change_relative_count(self):
        before=rational_time_state(self.P,CROSS_LO)
        after=rational_time_state(self.P,CROSS_HI)
        self.assertEqual(before.raw_large_square_count,2)
        self.assertEqual(after.raw_large_square_count,2)
        self.assertEqual(before.corrected_integer_index,1)
        self.assertEqual(after.corrected_integer_index,1)
        self.assertEqual(before.maximizing_chart,BRANCH_A)
        self.assertEqual(after.maximizing_chart,BRANCH_B)

    def test_existing_complete_orbit_census_matches_global_max_at_switch(self):
        # Independent v223 all-square enumeration tests the new global
        # envelope identity at a rational time near the bottleneck.
        time=Q(10011,20000)
        polygon=tuple(ppoint_at(p,time) for p in self.P)
        census=census_polygon(polygon)
        actual=max(s.side_squared for s in census.squares)
        self.assertEqual(actual,rational_time_state(self.P,time).largest_square_side_squared)

    def test_left_endpoint_independent_all_square_max(self):
        polygon=tuple(ppoint_at(p,WINDOW_START) for p in self.P)
        census=census_polygon(polygon)
        actual=max(s.side_squared for s in census.squares)
        self.assertEqual(actual,rational_time_state(self.P,WINDOW_START).largest_square_side_squared)

    def test_arbitrary_family_or_wrong_cutoff_cannot_claim_invariant(self):
        modified=list(self.P)
        x,y=modified[2]
        modified[2]=(padd(x,pol(Q(1,100000))),y)
        for time in (Q(1,2),Q(10011,20000)):
            with self.assertRaises(ValueError):
                rational_time_state(tuple(modified),time)
        with self.assertRaises(ValueError):
            rational_time_state(self.P,WINDOW_START,Q(9))

    def test_unverified_time_and_floating_point_are_rejected(self):
        for invalid in (0.5,0.50055,Q(499,1000),Q(503,1000),
                        "0.50055",None):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    rational_time_state(self.P,invalid)

    def test_global_envelope_rejects_uncertified_geometry(self):
        bad=list(self.P)
        x,y=bad[4]
        bad[4]=(padd(x,pol(Q(1,100))),y)
        with self.assertRaises(ValueError):
            certify_global_envelope(tuple(bad))

    def test_report_is_proven_and_explicitly_restricted(self):
        self.assertEqual(self.report["cutoff_square_side_squared"],str(HIGH_S2))
        self.assertEqual(self.report["minimum_of_global_max_square_side_squared_rigorous_interval"],
                         [str(BOTTLE_S2_LO),str(BOTTLE_S2_HI)])
        self.assertIn("no invariant across arbitrary orbit births/deaths",self.report["scope"])
        self.assertIn("no Toeplitz proof",self.report["scope"])
        self.assertIn("1.02337",self.report["restricted_Big_Peg_consequence"])


if __name__=="__main__":
    unittest.main()
