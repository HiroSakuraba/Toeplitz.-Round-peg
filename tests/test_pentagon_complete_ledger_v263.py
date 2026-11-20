"""v263 all-real-time, whole-pentagon complete positive-square census tests."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from big_peg_pair_death_v254 import (
    family_vertices, QPOLY, ROOT_LEFT,ROOT_RIGHT,
    PERSISTENT,COLLIDING_LEFT,COLLIDING_RIGHT,
    certify_collision,certify_jordan_family
)
from big_peg_multichart_v249 import (
    ppoint_at,certify_sign,pneg,pscale
)
from moving_polygon_events_v224 import (
    chart_polynomials as build_chart
)
from polygon_orbit_census_v223 import census_polygon
from pentagon_complete_ledger_v263 import (
    START,END,ALL_EDGE_CHARTS,TRACKED,adjacent_edges,
    sign_constant,edge_parameter_sign_polynomials,
    singular_inconsistency_witness,certified_side_cutoff,
    classify_chart,diagnostic,audit,
)


class ExhaustivePentagonSquareLedger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.report=audit()
        cls.eta,cls.gap=certified_side_cutoff(cls.P)

    def test_all_weak_increasing_cyclic_edge_quadruples_covered(self):
        self.assertEqual(len(ALL_EDGE_CHARTS),70)
        self.assertEqual(ALL_EDGE_CHARTS,
                         tuple(combinations_with_replacement(range(5),4)))
        self.assertEqual(len(set(ALL_EDGE_CHARTS)),70)

    def test_no_unresolved_regular_or_singular_chart(self):
        self.assertEqual(self.report["unresolved_charts"],{})
        self.assertTrue(self.report["all_real_time_exhaustiveness_proved"])
        self.assertEqual(sum(self.report["classes"].values()),70)

    def test_all_four_exclusion_categories_accounted(self):
        self.assertEqual(self.report["classes"],{
            "three_collinear":25,
            "adjacent_two_edges":5,
            "off_edge":37,
            "tracked":3,
        })
        self.assertNotIn("below_certified_cutoff",self.report["classes"])
        self.assertNotIn("singular_inconsistent",self.report["classes"])
        self.assertNotIn("unresolved_determinant_event",self.report["classes"])

    def test_every_nontracked_chart_is_excluded(self):
        report=self.report["full_classification"]
        tracked={tuple(map(int,k.split(","))) for k,r in report.items()
                 if r["class"]=="tracked"}
        self.assertEqual(tracked,set(TRACKED))
        self.assertEqual(tracked,{
            PERSISTENT,COLLIDING_LEFT,COLLIDING_RIGHT})

    def test_rational_bernstein_witness_for_every_off_edge_chart(self):
        P0=tuple(ppoint_at(v,START) for v in self.P)
        P1=tuple(ppoint_at(v,END) for v in self.P)
        total=0
        for key,entry in self.report["full_classification"].items():
            if entry["class"]!="off_edge":
                continue
            ids=tuple(map(int,key.split(",")))
            chart=build_chart(P0,P1,ids,orientation=1,
                              cutoff_squared=self.eta)
            detsign=sign_constant(chart.determinant)
            self.assertIn(detsign,(-1,1))
            witness=edge_parameter_sign_polynomials(
                chart,detsign)[entry["witness"]]
            # A strict negative bound for all real times excludes root,
            # including every algebraic time inside the interval.
            self.assertGreater(
                certify_sign(pneg(witness),Q(0),Q(1),strict=True),0)
            self.assertLess(witness[0],0)
            self.assertLess(sum(witness),0)
            total+=1
        self.assertEqual(total,37)

    def test_25_three_collinear_charts_cannot_support_positive_square(self):
        total=0
        for key,entry in self.report["full_classification"].items():
            if entry["class"]=="three_collinear":
                ids=tuple(map(int,key.split(",")))
                self.assertGreaterEqual(
                    max(ids.count(i) for i in set(ids)),3)
                total+=1
        self.assertEqual(total,25)

    def test_five_adjacent_two_edge_charts_excluded_by_convex_support(self):
        total=0
        for key,entry in self.report["full_classification"].items():
            if entry["class"]=="adjacent_two_edges":
                ids=tuple(map(int,key.split(",")))
                self.assertEqual(len(set(ids)),2)
                i,j=sorted(set(ids))
                self.assertTrue(adjacent_edges(i,j))
                self.assertEqual(ids.count(i),2)
                self.assertEqual(ids.count(j),2)
                total+=1
        self.assertEqual(total,5)

    def test_convex_jordan_property_is_all_real_time(self):
        j=certify_jordan_family(self.P)
        self.assertTrue(j["convex_Jordan_for_all_real_times"])
        self.assertEqual(j["global_nonadjacent_edge_rect_gap"],self.gap)
        self.assertEqual(self.gap,Q(3901,10000))

    def test_strict_size_lower_bound_for_any_square(self):
        self.assertEqual(self.eta,Q(15217801,400000000))
        self.assertEqual(self.eta,self.gap*self.gap/Q(4))
        self.assertGreater(self.gap*self.gap/Q(2),self.eta)
        self.assertLess(self.eta,1)
        # There are NO classified below-cutoff active charts.
        self.assertEqual(self.report["classes"].get(
            "below_certified_cutoff",0),0)

    def test_exact_algebraic_single_event_replayed(self):
        collision=certify_collision(self.P)
        self.assertTrue(collision[
            "root_identity_divisible_by_event_in_four_components"])
        self.assertEqual(
            self.report["unique_square_event_polynomial"],
            [str(x) for x in QPOLY])

    def test_all_real_time_geometric_orbit_counts(self):
        self.assertEqual(
            self.report["all_real_time_positive_geometric_square_count"],
            {"before_algebraic_collision":3,
             "at_algebraic_collision":2,
             "after_algebraic_collision":1})

    def test_independent_all_square_census_at_rational_endpoints(self):
        expected=(3,1)
        for t,n in zip((START,END),expected):
            polygon=tuple(ppoint_at(v,t) for v in self.P)
            census=census_polygon(polygon)
            self.assertEqual(census.count,n)
            self.assertTrue(all(s.side_squared>1 for s in census.squares))
        self.assertEqual(
            [x["geometric_squares"] for x in self.report[
             "independent_complete_rational_slices"]],[3,1])

    def test_noncanonical_geometry_fails_closed(self):
        changed=list(self.P)
        x,y=changed[1]
        changed[1]=((x[0]+Q(1,100),),y)
        with self.assertRaises(ValueError):
            certified_side_cutoff(tuple(changed))
        with self.assertRaises(ValueError):
            classify_chart(tuple(changed),PERSISTENT,
                           self.eta,START,END)

    def test_invalid_cutoff_time_or_assignment_fails_closed(self):
        for wrong in (Q(0),Q(-1),1e-3,None):
            with self.subTest(wrong=wrong):
                with self.assertRaises(ValueError):
                    classify_chart(self.P,PERSISTENT,wrong,START,END)
        with self.assertRaises(ValueError):
            classify_chart(self.P,PERSISTENT,self.eta,START,END+Q(1,1000))
        with self.assertRaises(ValueError):
            classify_chart(self.P,(4,3,2,1),self.eta,START,END)

    def test_unresolved_chart_causes_audit_to_fail_closed(self):
        # Force the arithmetic certificate primitive to refuse proof;
        # the global checker MUST NOT silently count only known roots.
        with patch("pentagon_complete_ledger_v263.sign_constant",
                   return_value=None):
            with self.assertRaises(AssertionError):
                audit()

    def test_nobody_can_promote_to_general_jordan_theorem(self):
        self.assertIn("ONLY for v254",self.report["scope"])
        self.assertIn("Does NOT prove",self.report["scope"])


if __name__=="__main__":
    unittest.main()
