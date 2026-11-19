"""v252: exhaustive exact cutoff square-orbit ledger and parity failure."""
import sys
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import census_polygon, pt
from big_peg_multichart_v249 import (
    pol, padd, psub, pneg, pmul, peval, ppoint_at, bcoeff,
)
from big_peg_corner_v250 import family_vertices, prove_jordan_family
from big_peg_size_atlas_v251 import (
    BRANCH_A, BRANCH_B, WINDOW_START, WINDOW_END,
    HIGH_S2, ATLAS_SWITCH,
)
from big_peg_cutoff_ledger_v252 import (
    EXPECTED_CLASSIFICATION, EXPECTED_ORBITS,
    B_ENTER_BRACKET, A_EXIT_BRACKET,
    _strict_sign_by_bernstein, _positive_except_possible_left,
    _edge_box, _max_pair_distance_squared,
    cutoff_numerator, poly_derivative, exhaustive_cutoff_classification,
    exact_endpoint_census, unique_cutoff_event, audit,
)


class CompleteHighSquareLedger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.report=audit()

    def test_complete_chart_partition_of_all_126_possibilities(self):
        parts=self.report["classifications"]
        self.assertEqual(parts,EXPECTED_CLASSIFICATION)
        self.assertEqual(sum(parts.values()),126)
        self.assertEqual(self.report["all_distinct_edge_charts_checked"],126)

    def test_two_and_only_two_surviving_large_square_branches(self):
        actual={tuple(a) for a in self.report["only_two_possible_orbit_branches"]}
        self.assertEqual(actual,set(EXPECTED_ORBITS))
        self.assertEqual(len(actual),2)

    def test_all_time_convex_jordan_geometry_is_certified(self):
        jordan=prove_jordan_family(self.P)
        self.assertTrue(jordan["strictly_convex_and_jordan_for_all_times"])
        self.assertTrue(jordan["enclosed_in_radius_3_disk"])

    def test_87_geometry_rejections_are_actual_global_max_distance_bounds(self):
        from big_peg_cutoff_ledger_v252 import _edge_box, _max_pair_distance_squared
        boxes=[_edge_box(self.P,i,WINDOW_START,WINDOW_END)
               for i in range(len(self.P))]
        count=0
        for ids in combinations(range(9),4):
            if any(_max_pair_distance_squared(boxes[i],boxes[j])<HIGH_S2
                   for i,j in combinations(ids,2)):
                count+=1
        self.assertEqual(count,87)

    def test_edge_length_exclusion_repeated_edge_is_exact(self):
        from big_peg_multichart_v249 import pvsub, pvdot, certify_sign
        for i in range(9):
            e=pvsub(self.P[(i+1)%9],self.P[i])
            result=certify_sign(
                psub(pol(HIGH_S2),pvdot(e,e)),
                WINDOW_START,WINDOW_END,strict=True)
            self.assertGreater(result,0)

    def test_36_off_edge_and_one_below_cutoff_rejected(self):
        self.assertEqual(self.report["classifications"]["off_edge"],36)
        self.assertEqual(self.report["classifications"]["below_cutoff"],1)
        self.assertGreater(self.report["exact_auxiliary_sign_checks"],10)

    def test_endpoint_geometric_census_deduplicates_vertex_chart(self):
        first,last=self.report["independent_complete_v223_endpoint_censuses"]
        self.assertEqual(first["geometric_high_square_orbits"],1)
        self.assertEqual(last["geometric_high_square_orbits"],1)
        self.assertTrue(first["positive_square_on_polygon_vertex"])
        self.assertGreater(first["repeated_labelings_deduplicated"],0)
        self.assertEqual(first["time"],"1/2")
        self.assertEqual(last["time"],"501/1000")

    def test_exact_interior_census_counts_two_geometric_orbits(self):
        t=ATLAS_SWITCH
        P=tuple(ppoint_at(p,t) for p in self.P)
        census=census_polygon(P,cutoff_squared=HIGH_S2)
        self.assertEqual(census.count,2)
        self.assertTrue(all(s.side_squared>HIGH_S2 for s in census.squares))

    def test_entry_and_exit_brackets_are_strictly_separated(self):
        self.assertLess(B_ENTER_BRACKET[1],A_EXIT_BRACKET[0])
        self.assertGreater(B_ENTER_BRACKET[0],WINDOW_START)
        self.assertLess(A_EXIT_BRACKET[1],WINDOW_END)

    def test_B_enters_by_one_unique_exact_polynomial_root(self):
        r=unique_cutoff_event(
            self.P,BRANCH_B,HIGH_S2,
            WINDOW_START,WINDOW_END,B_ENTER_BRACKET,direction=1)
        self.assertEqual(r["direction"],"enters above cutoff")
        self.assertTrue(r["unique_algebraic_crossing"])
        self.assertGreater(r["derivative_sign_checks"],0)
        self.assertLess(peval(r["cutoff_polynomial"],B_ENTER_BRACKET[0]),0)
        self.assertGreater(peval(r["cutoff_polynomial"],B_ENTER_BRACKET[1]),0)

    def test_A_exits_by_one_unique_exact_polynomial_root(self):
        r=unique_cutoff_event(
            self.P,BRANCH_A,HIGH_S2,
            WINDOW_START,WINDOW_END,A_EXIT_BRACKET,direction=-1)
        self.assertEqual(r["direction"],"exits above cutoff")
        self.assertTrue(r["unique_algebraic_crossing"])
        self.assertGreater(r["derivative_sign_checks"],0)
        self.assertGreater(peval(r["cutoff_polynomial"],A_EXIT_BRACKET[0]),0)
        self.assertLess(peval(r["cutoff_polynomial"],A_EXIT_BRACKET[1]),0)

    def test_size_filtered_parity_changes_odd_even_odd(self):
        self.assertEqual(self.report["geometric_orbit_counts_strictly_above_cutoff_in_time_order"],
                         [1,2,1])
        self.assertEqual(self.report["filtered_orbit_parities_in_time_order"],
                         [1,0,1])
        self.assertTrue(self.report["filter_parity_is_not_invariant"])

    def test_equal_size_crossing_is_distinct_from_two_size_cutoff_events(self):
        center=self.report["separate_same_size_crossing_bracket"]
        xlo,xhi=map(Q,center)
        self.assertTrue(B_ENTER_BRACKET[1]<xlo<xhi<A_EXIT_BRACKET[0])
        self.assertTrue(self.report["same_size_crossing_does_not_change_filtered_count"])

    def test_closed_right_exclusion_allows_only_special_left_endpoint(self):
        # All Bernstein basis terms are strictly positive in the
        # OPEN interval, so a nonzero sum of nonnegative terms is >0.
        a,b=WINDOW_START,WINDOW_END
        self.assertTrue(_positive_except_possible_left(pol(-a,1),a,b))
        self.assertFalse(_positive_except_possible_left(pol(0),a,b))
        self.assertFalse(_positive_except_possible_left(pol(a,-1),a,b))

    def test_bernstein_sign_rejects_zero_or_changing_determinant(self):
        self.assertEqual(_strict_sign_by_bernstein(pol(1),WINDOW_START,WINDOW_END),1)
        self.assertEqual(_strict_sign_by_bernstein(pol(-1),WINDOW_START,WINDOW_END),-1)
        with self.assertRaises(ValueError):
            _strict_sign_by_bernstein(pol(0),WINDOW_START,WINDOW_END)
        with self.assertRaises(ValueError):
            _strict_sign_by_bernstein(pol(-ATLAS_SWITCH,1),WINDOW_START,WINDOW_END)

    def test_scope_rejects_claims_about_all_square_orbits_or_toeplitz(self):
        scope=self.report["scope"]
        self.assertIn("EXHAUSTIVE only",scope)
        self.assertIn("no assertion on smaller square orbits",scope)
        self.assertIn("Toeplitz",scope)

    def test_different_cutoff_fails_closed(self):
        with self.assertRaises(ValueError):
            exhaustive_cutoff_classification(
                self.P,Q(9),WINDOW_START,WINDOW_END)

    def test_changed_family_fails_closed(self):
        tampered=list(self.P)
        x,y=tampered[2]
        tampered[2]=(padd(x,pol(Q(1,100000))),y)
        with self.assertRaises(ValueError):
            exhaustive_cutoff_classification(
                tuple(tampered),HIGH_S2,WINDOW_START,WINDOW_END)

    def test_wrong_interval_and_float_fail_closed(self):
        for a,b,k in ((Q(0),WINDOW_END,HIGH_S2),
                      (WINDOW_START,Q(1),HIGH_S2),
                      (0.5,WINDOW_END,HIGH_S2),
                      (WINDOW_START,WINDOW_END,9.4257)):
            with self.subTest(a=a,b=b,k=k):
                with self.assertRaises(ValueError):
                    exhaustive_cutoff_classification(self.P,k,a,b)

    def test_invalid_event_brackets_rejected(self):
        with self.assertRaises(ValueError):
            unique_cutoff_event(
                self.P,BRANCH_B,HIGH_S2,WINDOW_START,WINDOW_END,
                (A_EXIT_BRACKET[0],A_EXIT_BRACKET[1]),direction=1)
        with self.assertRaises(ValueError):
            unique_cutoff_event(
                self.P,BRANCH_A,HIGH_S2,WINDOW_START,WINDOW_END,
                B_ENTER_BRACKET,direction=-1)


if __name__=="__main__":
    unittest.main()
