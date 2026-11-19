"""v251: every-time square SIZE cutoff with discontinuous selected witnesses."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from big_peg_multichart_v249 import (
    HALF, pol, peval, certify_chart, chart_polynomials, root_at,
    certify_sign,
)
from big_peg_corner_v250 import (
    family_vertices, prove_jordan_family, LEFT_CHART, RIGHT_CHART,
)
from big_peg_size_atlas_v251 import (
    BRANCH_A,BRANCH_B,WINDOW_START,WINDOW_END,
    ATLAS_SWITCH,HIGH_S2,FULL_S2,CROSS_LEFT,CROSS_RIGHT,
    certify_existence_atlas,side_comparison_polynomial,derivative,
    exact_switch_proof,square_side2,audit,
)


class ExactSizeAwareBranchCover(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.report=audit()
        cls.switch=exact_switch_proof(cls.P)

    def test_entire_family_remains_jordan_and_bounded(self):
        proof=prove_jordan_family(self.P)
        self.assertTrue(proof["strictly_convex_and_jordan_for_all_times"])
        self.assertTrue(proof["enclosed_in_radius_3_disk"])

    def test_window_two_charts_cover_every_real_time(self):
        atlas=self.switch["atlas"]
        self.assertTrue(atlas.fully_covered)
        self.assertEqual(atlas.start,HALF)
        self.assertEqual(atlas.end,Q(501,1000))
        self.assertEqual(atlas.side_squared_strictly_above,Q(94257,10000))
        self.assertEqual(len(atlas.certificates),2)

    def test_two_distinct_squares_at_atlas_switch(self):
        self.assertEqual(self.switch["atlas"].distinct_seam_squares,1)
        self.assertTrue(self.switch["A_and_B_distinct_at_seam"])
        A,B=self.switch["atlas"].certificates
        self.assertNotEqual(A.right_root,B.left_root)
        self.assertGreater(square_side2(A.right_root),HIGH_S2)
        self.assertGreater(square_side2(B.left_root),HIGH_S2)

    def test_different_candidates_reverse_size_ordering(self):
        s=self.switch
        self.assertGreater(s["candidate_A_side2_at_start"],
                           s["candidate_B_side2_at_start"])
        self.assertLess(s["candidate_A_side2_at_end"],
                        s["candidate_B_side2_at_end"])

    def test_neither_candidate_alone_proves_high_cutoff(self):
        A=certify_chart(self.P,BRANCH_A,WINDOW_START,WINDOW_END,side2_limit=FULL_S2)
        B=certify_chart(self.P,BRANCH_B,WINDOW_START,WINDOW_END,side2_limit=FULL_S2)
        self.assertLess(square_side2(A.right_root),HIGH_S2)
        self.assertLess(square_side2(B.left_root),HIGH_S2)
        with self.assertRaises(ValueError):
            certify_chart(self.P,BRANCH_A,WINDOW_START,WINDOW_END,
                          side2_limit=HIGH_S2)
        with self.assertRaises(ValueError):
            certify_chart(self.P,BRANCH_B,WINDOW_START,WINDOW_END,
                          side2_limit=HIGH_S2)

    def test_polynomial_crossing_has_both_exact_signs(self):
        G=side_comparison_polynomial(self.P,BRANCH_A,BRANCH_B)
        self.assertGreater(peval(G,WINDOW_START),0)
        self.assertLess(peval(G,WINDOW_END),0)
        self.assertGreater(peval(G,CROSS_LEFT),0)
        self.assertLess(peval(G,CROSS_RIGHT),0)

    def test_crossover_is_unique_by_derivative(self):
        G=side_comparison_polynomial(self.P,BRANCH_A,BRANCH_B)
        deriv=derivative(G)
        self.assertGreater(
            certify_sign(tuple(-c for c in deriv),
                         WINDOW_START,WINDOW_END,strict=True),0)
        self.assertGreater(self.switch["strict_derivative_sign_checks"],0)

    def test_crossover_bracket_is_tighter_than_one_millionth(self):
        self.assertEqual(CROSS_LEFT,Q(500549,1000000))
        self.assertEqual(CROSS_RIGHT,Q(10011,20000))
        self.assertGreater(CROSS_RIGHT-CROSS_LEFT,0)
        self.assertLess(CROSS_RIGHT-CROSS_LEFT,Q(1,100000))

    def test_whole_family_branch_neutral_cover(self):
        selected=(
            (Q(0),HALF,RIGHT_CHART),
            (HALF,ATLAS_SWITCH,BRANCH_A),
            (ATLAS_SWITCH,WINDOW_END,BRANCH_B),
            (WINDOW_END,Q(1),BRANCH_A),
        )
        result=certify_existence_atlas(self.P,selected,Q(9))
        self.assertTrue(result.fully_covered)
        self.assertEqual(result.distinct_seam_squares,2)
        self.assertEqual(len(result.certificates),4)
        self.assertEqual(result.start,Q(0))
        self.assertEqual(result.end,Q(1))

    def test_overlapping_charts_can_be_accepted(self):
        atlas=certify_existence_atlas(
            self.P,((HALF,ATLAS_SWITCH,BRANCH_A),
                    (HALF,WINDOW_END,BRANCH_B)),
            FULL_S2,start=HALF,end=WINDOW_END)
        self.assertTrue(atlas.fully_covered)

    def test_empty_or_untouched_atlas_rejected(self):
        with self.assertRaises(ValueError):
            certify_existence_atlas(self.P,(),HIGH_S2,start=HALF,end=WINDOW_END)

    def test_open_gap_rejected(self):
        with self.assertRaises(ValueError):
            certify_existence_atlas(
                self.P,((HALF,Q(10001,20000),BRANCH_A),
                        (ATLAS_SWITCH,WINDOW_END,BRANCH_B)),
                HIGH_S2,start=HALF,end=WINDOW_END)

    def test_missing_final_endpoint_rejected(self):
        with self.assertRaises(ValueError):
            certify_existence_atlas(
                self.P,((HALF,ATLAS_SWITCH,BRANCH_A),),
                HIGH_S2,start=HALF,end=WINDOW_END)

    def test_invalid_float_cutoffs_and_times_rejected(self):
        for invalid in (9.4257,"94257/10000",Q(0),Q(-1)):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    certify_existence_atlas(
                        self.P,((HALF,ATLAS_SWITCH,BRANCH_A),),
                        invalid,start=HALF,end=ATLAS_SWITCH)
        with self.assertRaises(ValueError):
            certify_existence_atlas(
                self.P,((0.5,ATLAS_SWITCH,BRANCH_A),),
                HIGH_S2,start=HALF,end=ATLAS_SWITCH)

    def test_invalid_single_chart_fails_closed(self):
        with self.assertRaises(ValueError):
            certify_existence_atlas(
                self.P,((HALF,WINDOW_END,BRANCH_A),),
                HIGH_S2,start=HALF,end=WINDOW_END)

    def test_endpoint_roots_are_real_squares_in_selected_charts(self):
        for c in self.switch["atlas"].certificates:
            for root,pars in ((c.left_root,c.left_parameters),
                              (c.right_root,c.right_parameters)):
                self.assertTrue(all(Q(0)<=x<=Q(1) for x in pars))
                self.assertGreater(square_side2(root),HIGH_S2)

    def test_replay_report_scoped_and_reproducible(self):
        self.assertEqual(self.report["target_short_interval"],
                         ["1/2","501/1000"])
        self.assertEqual(self.report["short_interval_side_squared_strictly_above"],
                         "94257/10000")
        self.assertTrue(self.report["candidates_reverse_relative_size"])
        self.assertIn("no proof of global maximal-orbit switching",self.report["scope"])
        self.assertIn("Toeplitz",self.report["scope"])


if __name__=="__main__":
    unittest.main()
