"""v227 exact event-germ provenance and geometric-orbit incidence tests."""
import sys
import unittest
from dataclasses import replace
from pathlib import Path
from fractions import Fraction as Q
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))

from polygon_orbit_census_v223 import pt, convex_vertex_alignment_polygon, census_polygon
from moving_polygon_events_v224 import polynomial, RootBracket, times, minus, T, ONE
from algebraic_orbit_events_v226 import EventCoverInconclusive
from persistent_strata_v225 import classify_polygon_motion
from local_orbit_germs_v227 import (
    AlgebraicTime, sign_at_time, germ_signs, _same_geometric_orbit,
    _classify, GermInconclusive, audit_local_event_germs,
    classify_rational_or_algebraic_event, report,
)


class AlgebraicSignTests(unittest.TestCase):
    def test_rational_simple_crossing(self):
        a=AlgebraicTime((-Q(8,9),1),Q(8,9),Q(8,9))
        self.assertEqual(germ_signs((-Q(8,9),1),a),(-1,0,1))
        self.assertEqual(germ_signs((Q(8,9),-1),a),(1,0,-1))

    def test_rational_even_tangency(self):
        a=AlgebraicTime((-Q(1,3),1),Q(1,3),Q(1,3))
        f=times(polynomial((-Q(1,3),1)),polynomial((-Q(1,3),1)))
        self.assertEqual(germ_signs(f,a),(1,0,1))
        self.assertEqual(germ_signs(tuple(-x for x in f),a),(-1,0,-1))

    def test_rational_triple_crossing(self):
        a=AlgebraicTime((-Q(1,3),1),Q(1,3),Q(1,3))
        f=polynomial((-Q(1,3),1))
        cube=times(times(f,f),f)
        self.assertEqual(germ_signs(cube,a),(-1,0,1))

    def test_irrational_root_sign_and_tangent(self):
        a=AlgebraicTime((-2,0,1),Q(7,5),Q(3,2))
        self.assertEqual(sign_at_time((-2,0,1),a),0)
        self.assertEqual(germ_signs((-2,0,1),a),(-1,0,1))
        self.assertEqual(sign_at_time((-1,1),a),1)
        self.assertEqual(sign_at_time((Q(3,2),-1),a),1)
        self.assertEqual(sign_at_time((2,-2),a),-1)

    def test_irrational_quadratic_double_root(self):
        a=AlgebraicTime((-2,0,1),Q(7,5),Q(3,2))
        f=times(polynomial((-2,0,1)),polynomial((-2,0,1)))
        self.assertEqual(germ_signs(f,a),(1,0,1))

    def test_zero_polynomial_has_zero_germs(self):
        a=AlgebraicTime((-2,0,1),Q(7,5),Q(3,2))
        self.assertEqual(germ_signs((0,),a),(0,0,0))

    def test_nearby_distinct_irrational_roots_not_equal(self):
        a=AlgebraicTime((-2,0,1),Q(7,5),Q(3,2))
        shifted=polynomial((Q(1,10000)-2,-Q(1,50),1))
        self.assertNotEqual(sign_at_time(shifted,a),0)

    def test_invalid_algebraic_brackets_rejected(self):
        with self.assertRaises(ValueError):AlgebraicTime((-2,0,1),1,Q(11,10))
        with self.assertRaises(ValueError):AlgebraicTime((-2,0,1),2,1)
        with self.assertRaises(ValueError):AlgebraicTime((Q(2,9),-1,1),0,1)

    def test_algebraic_sign_and_germ_are_pure_repeatable(self):
        a=AlgebraicTime((-2,0,1),Q(7,5),Q(3,2))
        f=(-1,1)
        self.assertEqual(sign_at_time(f,a),sign_at_time(f,a))
        self.assertEqual(germ_signs(f,a),germ_signs(f,a))


class ConvexVertexGermTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P0=convex_vertex_alignment_polygon(0)
        cls.P1=convex_vertex_alignment_polygon(Q(9,8))
        cls.base,cls.events=audit_local_event_germs(cls.P0,cls.P1,Q(1,100),Q(3,4),1)
        cls.e=cls.events[0]

    def test_single_rational_event_at_eight_ninths(self):
        self.assertEqual(len(self.events),1)
        self.assertTrue(self.e.algebraic_time.rational)
        self.assertEqual(self.e.algebraic_time.lo,Q(8,9))
        self.assertEqual(self.e.raw_candidate_entries,224)

    def test_two_event_geometric_squares_exactly(self):
        self.assertEqual(len(self.e.groups),2)
        keys={g.event_key for g in self.e.groups}
        self.assertEqual(len(keys),2)
        self.assertEqual(len(census_polygon(convex_vertex_alignment_polygon(1),Q(1,100)).squares),2)

    def test_genuine_birth_and_continuation_not_chart_duplicates(self):
        bykind={g.classification:g for g in self.e.groups}
        self.assertEqual(set(bykind),{'pair-birth','continuation'})
        self.assertEqual((bykind['pair-birth'].left_count,bykind['pair-birth'].right_count),(0,2))
        self.assertEqual((bykind['continuation'].left_count,bykind['continuation'].right_count),(1,1))
        self.assertEqual(bykind['pair-birth'].chart_presentations,16)
        self.assertEqual(bykind['continuation'].chart_presentations,8)
        self.assertEqual((self.e.left_chamber_count,self.e.right_chamber_count),(1,3))
        self.assertEqual(self.e.event_slice_count,2)
        self.assertTrue(self.e.balances())

    def test_birth_event_is_at_a_polygon_vertex(self):
        birth=next(g for g in self.e.groups if g.classification=='pair-birth')
        self.assertTrue(birth.vertex_incidence)
        self.assertEqual(set(birth.event_key),{(Q(1),Q(0)),(Q(2),Q(0)),(Q(2),Q(1)),(Q(1),Q(1))})

    def test_no_parity_change_when_two_squares_born(self):
        self.assertEqual(self.e.net_change,2)
        self.assertEqual(self.e.left_chamber_count%2,self.e.right_chamber_count%2)

    def test_all_group_germs_account_for_chamber_orbits(self):
        self.assertEqual(sum(g.left_count for g in self.e.groups),1)
        self.assertEqual(sum(g.right_count for g in self.e.groups),3)
        self.assertEqual(set(k for g in self.e.groups for k in g.left_keys),set(self.base.chambers[0].orbit_keys))
        self.assertEqual(set(k for g in self.e.groups for k in g.right_keys),set(self.base.chambers[1].orbit_keys))

    def test_bad_index_rejected(self):
        with self.assertRaises(ValueError):
            classify_rational_or_algebraic_event(self.P0,self.P1,Q(1,100),self.base,2)

    def test_report(self):
        r=report()
        self.assertEqual(r['convex_event_time'],'8/9')
        self.assertEqual(r['net_change'],2)


class IrrationalCutoffGermTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P0=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        cls.P1=tuple(pt(2*p.x,2*p.y) for p in cls.P0)
        cls.base,cls.events=audit_local_event_germs(cls.P0,cls.P1,Q(144,25),0,1)

    def test_six_clusters_only_one_true_crossing(self):
        self.assertEqual(len(self.events),6)
        self.assertEqual(self.base.raw_candidates,88)
        self.assertEqual([e.net_change for e in self.events],[0,0,0,1,0,0])
        self.assertTrue(all(e.balances() for e in self.events))

    def test_irrational_size_entry_is_unique_orbit_germ(self):
        ev=self.events[3]
        self.assertFalse(ev.algebraic_time.rational)
        self.assertEqual(len(ev.groups),1)
        g=ev.groups[0]
        self.assertEqual(g.classification,'size-entry')
        self.assertEqual(g.size_status,'on-cutoff')
        self.assertEqual((g.left_count,g.right_count),(0,1))
        self.assertEqual(g.chart_presentations,8)
        self.assertEqual((ev.left_chamber_count,ev.right_chamber_count),(0,1))
        self.assertIsNone(ev.event_slice_count)

    def test_size_entry_root_is_sqrt2_minus_one(self):
        ev=self.events[3]
        lo,hi=ev.algebraic_time.lo,ev.algebraic_time.hi
        poly=(Q(-1),Q(2),Q(1))
        self.assertLess(lo,hi)
        self.assertLess((lo+1)**2-2,0)
        self.assertGreater((hi+1)**2-2,0)
        self.assertEqual(sign_at_time(poly,ev.algebraic_time),0)
        self.assertEqual(germ_signs(poly,ev.algebraic_time),(-1,0,1))

    def test_three_false_early_clusters_have_no_geometric_germ(self):
        self.assertTrue(all(len(self.events[i].groups)==0 for i in (0,1,2)))

    def test_two_late_clusters_only_continue_unchanged_square(self):
        for i in (4,5):
            ev=self.events[i]
            self.assertEqual(len(ev.groups),1)
            self.assertEqual(ev.groups[0].classification,'continuation')
            self.assertEqual((ev.groups[0].left_count,ev.groups[0].right_count),(1,1))

    def test_irrational_event_has_no_fabricated_rational_key(self):
        self.assertIsNone(self.events[3].groups[0].event_key)
        self.assertEqual(self.events[3].groups[0].chart_presentations,8)

    def test_total_orbits_not_labeled_charts(self):
        self.assertEqual(sum(g.right_count for g in self.events[3].groups),1)
        self.assertEqual(sum(g.chart_presentations for g in self.events[3].groups),8)


class FailClosedTests(unittest.TestCase):
    def test_positive_dimensional_rectangle_rejected(self):
        P=tuple(pt(*xy) for xy in ((0,0),(2,0),(2,1),(0,1)))
        with self.assertRaises((EventCoverInconclusive,ValueError)):
            audit_local_event_germs(P,P,Q(1,100))

    def test_rational_event_fails_if_chart_evidence_is_deleted(self):
        P0=convex_vertex_alignment_polygon(0)
        P1=convex_vertex_alignment_polygon(Q(9,8))
        base,_=audit_local_event_germs(P0,P1,Q(1,100),Q(3,4),1)
        strata=classify_polygon_motion(P0,P1,Q(1,100),Q(3,4),1)
        empty=replace(strata,certificates=())
        with self.assertRaises(GermInconclusive):
            classify_rational_or_algebraic_event(P0,P1,Q(1,100),base,0,empty)

    def test_irrational_event_fails_on_possible_singular_exception(self):
        P0=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        P1=tuple(pt(2*p.x,2*p.y) for p in P0)
        base,_=audit_local_event_germs(P0,P1,Q(144,25),0,1)
        strata=classify_polygon_motion(P0,P1,Q(144,25),0,1)
        artificial=replace(strata.certificates[0],status='exceptional-time-only',
                           obstruction=(Q(-1),Q(2),Q(1)))
        altered=replace(strata,certificates=(artificial,)+strata.certificates[1:])
        with self.assertRaisesRegex(GermInconclusive,'possibly feasible exceptional'):
            classify_rational_or_algebraic_event(P0,P1,Q(144,25),base,3,altered)

    def test_time_reversal_converts_pair_birth_to_death(self):
        P0=convex_vertex_alignment_polygon(Q(9,8))
        P1=convex_vertex_alignment_polygon(0)
        _,events=audit_local_event_germs(P0,P1,Q(1,100),0,Q(1,4))
        self.assertEqual(len(events),1)
        self.assertEqual(events[0].algebraic_time.lo,Q(1,9))
        self.assertEqual(set(g.classification for g in events[0].groups),
                         {'pair-death','continuation'})
        self.assertEqual(events[0].net_change,-2)

    def test_cutoff_must_be_positive(self):
        P=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        with self.assertRaises(ValueError):audit_local_event_germs(P,P,0)


if __name__=='__main__':unittest.main()
