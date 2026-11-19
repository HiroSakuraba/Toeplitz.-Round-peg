"""v226 exact algebraic root equivalence + independent polygon orbit census."""
import random
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))

from moving_polygon_events_v224 import (
    RootBracket,polynomial,times,minus,T,ONE,isolate_roots,evaluate,
    interpolated_polygon,
)
from polygon_orbit_census_v223 import (
    pt,convex_vertex_alignment_polygon,census_polygon,
)
from persistent_strata_v225 import classify_polygon_motion
from algebraic_orbit_events_v226 import (
    AlgebraicCandidate,EventCoverInconclusive,
    _has_root_in_closed,_same_root,cluster_candidate_roots,
    candidate_atoms_from_strata,audit_moving_polygon_orbits,report,
)


def atom(p,lo,hi,mark):
    return AlgebraicCandidate(polynomial(p),RootBracket(Q(lo),Q(hi)),
                              (((0,1,2,3),1,mark),))


class ExactClusterTests(unittest.TestCase):
    def test_rational_identical_root_unequal_polynomials(self):
        # Both vanish at 1/3, one also at 2/3 (outside chosen isolator).
        p=polynomial((-Q(1,3),1))
        q=times(p,polynomial((-Q(2,3),1)))
        x=atom(p,Q(1,4),Q(2,5),'a')
        y=atom(q,Q(1,4),Q(2,5),'b')
        self.assertTrue(_same_root(p,x.bracket,q,y.bracket))
        grouped=cluster_candidate_roots((x,y))
        self.assertEqual(len(grouped),1)
        self.assertEqual(grouped[0].exact_time,Q(1,3))
        self.assertEqual(grouped[0].raw_entries,2)

    def test_overlapping_intervals_different_close_rational_roots(self):
        x=atom((-Q(1,3),1),Q(3,10),Q(9,25),'a')
        y=atom((-Q(17,50),1),Q(3,10),Q(9,25),'b')
        grouped=cluster_candidate_roots((x,y))
        self.assertEqual([c.exact_time for c in grouped],[Q(1,3),Q(17,50)])
        self.assertLess(grouped[0].bracket.upper,grouped[1].bracket.lower)

    def test_common_irrational_root_proved_using_gcd(self):
        p=polynomial((-2,0,1));q=times(p,polynomial((1,1)))
        x=atom(p,Q(5,4),Q(3,2),'p')
        y=atom(q,Q(1),Q(2),'q')
        grouped=cluster_candidate_roots((x,y),Q(1,2**15))
        self.assertEqual(len(grouped),1)
        self.assertIsNone(grouped[0].exact_time)
        self.assertLess(grouped[0].bracket.upper-grouped[0].bracket.lower,Q(1,2**14))
        self.assertLess(evaluate(p,grouped[0].bracket.lower),0)
        self.assertGreater(evaluate(p,grouped[0].bracket.upper),0)

    def test_overlapping_intervals_distinct_irrational_roots(self):
        p=polynomial((-2,0,1))
        q=polynomial((Q(1,10000)-2,-Q(1,50),1))
        x=atom(p,1,Q(3,2),'p')
        y=atom(q,1,Q(3,2),'q')
        grouped=cluster_candidate_roots((x,y),Q(1,2**12))
        self.assertEqual(len(grouped),2)
        self.assertLess(grouped[0].bracket.upper,grouped[1].bracket.lower)
        self.assertTrue(all(c.exact_time is None for c in grouped))

    def test_polynomial_with_two_roots_not_collapsed(self):
        p=polynomial((Q(2,9),-1,1))  # (t-1/3)(t-2/3)
        a=atom(p,Q(1,4),Q(2,5),'first')
        b=atom(p,Q(3,5),Q(3,4),'second')
        grouped=cluster_candidate_roots((a,b))
        self.assertEqual(len(grouped),2)
        self.assertTrue(grouped[0].bracket.lower<Q(1,3)<grouped[0].bracket.upper)
        self.assertTrue(grouped[1].bracket.lower<Q(2,3)<grouped[1].bracket.upper)

    def test_repeated_root_squarefree_cluster(self):
        p=times(polynomial((-Q(1,3),1)),polynomial((-Q(1,3),1)))
        grouped=cluster_candidate_roots((atom(p,Q(1,4),Q(2,5),'double'),))
        self.assertEqual(grouped[0].exact_time,Q(1,3))

    def test_pairwise_root_equality_symmetric(self):
        rng=random.Random(226)
        for _ in range(20):
            v=Q(rng.randint(1,15),20)
            p=polynomial((-v,1))
            q=times(p,polynomial((1,Q(rng.randint(1,5)))))
            b=RootBracket(v-Q(1,100),v+Q(1,100))
            self.assertTrue(_same_root(p,b,q,b))
            self.assertTrue(_same_root(q,b,p,b))

    def test_input_requires_true_isolator(self):
        with self.assertRaises(ValueError):
            cluster_candidate_roots((atom((1,1),0,Q(1,2),'no-root'),))
        with self.assertRaises(ValueError):
            cluster_candidate_roots((atom((Q(2,9),-1,1),0,1,'two-roots'),))
        with self.assertRaises(ValueError):cluster_candidate_roots((),width=0)

    def test_root_on_exact_time(self):
        a=atom((-Q(2,3),1),Q(2,3),Q(2,3),'exact')
        b=atom(times(polynomial((-Q(2,3),1)),polynomial((1,1))),Q(1,2),Q(3,4),'wide')
        grouped=cluster_candidate_roots((a,b))
        self.assertEqual(len(grouped),1)
        self.assertEqual(grouped[0].exact_time,Q(2,3))

    def test_empty_root_list(self):
        self.assertEqual(cluster_candidate_roots(()),())

    def test_root_in_closed_interval_including_endpoints(self):
        p=polynomial((-Q(8,9),1))
        self.assertTrue(_has_root_in_closed(p,Q(8,9),Q(8,9)))
        self.assertTrue(_has_root_in_closed(p,Q(2,3),Q(8,9)))
        self.assertFalse(_has_root_in_closed(p,Q(8,9)+Q(1,1000),1))


class WholePolygonOrbitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P0=convex_vertex_alignment_polygon(Q(0))
        cls.P1=convex_vertex_alignment_polygon(Q(9,8))
        cls.audit=audit_moving_polygon_orbits(cls.P0,cls.P1,Q(1,100),Q(3,4),1,Q(1,2**11))

    def test_previous_224_raw_entries_preserved(self):
        self.assertEqual(self.audit.charts_total,408)
        self.assertEqual(self.audit.raw_candidates,224)
        self.assertEqual(sum(c.raw_entries for c in self.audit.clusters),224)

    def test_one_real_algebraic_time_cluster(self):
        self.assertEqual(len(self.audit.clusters),1)
        c=self.audit.clusters[0]
        self.assertEqual(c.exact_time,Q(8,9))
        self.assertEqual(c.minimal_factor,(-Q(8,9),Q(1)))

    def test_true_two_square_birth(self):
        self.assertEqual(tuple(c.count for c in self.audit.chambers),(1,3))
        ev=self.audit.events[0]
        self.assertEqual(ev.at_event_count,2)
        self.assertEqual(ev.net_change,2)
        self.assertEqual((ev.left_parity,ev.right_parity),(1,1))
        self.assertFalse(ev.has_parity_change)

    def test_crosscheck_v223_on_both_chambers(self):
        for c in self.audit.chambers:
            P=interpolated_polygon(self.P0,self.P1,c.sample)
            result=census_polygon(P,Q(1,100))
            self.assertEqual(result.count,c.count)
            self.assertEqual(tuple(s.side_squared for s in result.squares),c.square_sizes_squared)
            self.assertTrue(all(q>Q(1,100) for q in c.square_sizes_squared))

    def test_corner_event_square_exact(self):
        P=interpolated_polygon(self.P0,self.P1,Q(8,9))
        result=census_polygon(P,Q(1,100))
        self.assertEqual(result.count,2)
        expected={pt(1,0),pt(2,0),pt(2,1),pt(1,1)}
        self.assertTrue(any(set(s.vertices)==expected for s in result.squares))
        self.assertTrue(any(s.has_polygon_vertex for s in result.squares))

    def test_sectors_all_geometric_orbits_not_labelings(self):
        self.assertEqual(sum(len(c.orbit_keys) for c in self.audit.chambers),4)
        for c in self.audit.chambers:
            self.assertEqual(len(set(c.orbit_keys)),c.count)
            self.assertEqual(len(c.sectors),c.count)

    def test_full_report_consistency(self):
        x=report()
        self.assertEqual(x['raw_candidates'],224)
        self.assertEqual(x['distinct_algebraic_times'],1)
        self.assertEqual(x['event_time'],'8/9')
        self.assertEqual(x['event_chamber_orbit_counts'],(1,3))
        self.assertFalse(x['parity_changed'])

    def test_truly_singular_rectangle_fails_closed(self):
        P=tuple(pt(*x) for x in ((0,0),(2,0),(2,1),(0,1)))
        with self.assertRaises(EventCoverInconclusive):
            audit_moving_polygon_orbits(P,P,Q(1,100),0,1,Q(1,2**10))

    def test_invalid_cutoff_refused(self):
        for k in (0,-Q(1,8)):
            with self.assertRaises(ValueError):
                audit_moving_polygon_orbits(self.P0,self.P1,k,Q(3,4),1)


class ScalingPolygonControl(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P0=tuple(pt(*x) for x in ((0,0),(2,0),(3,3),(0,2)))
        cls.P1=tuple(pt(2*p.x,2*p.y) for p in cls.P0)
        cls.audit=audit_moving_polygon_orbits(cls.P0,cls.P1,Q(144,25),0,1,Q(1,2**10))

    def test_false_positive_roots_distinguished(self):
        self.assertEqual(self.audit.raw_candidates,88)
        self.assertEqual(len(self.audit.clusters),6)
        self.assertEqual(tuple(c.count for c in self.audit.chambers),(0,0,0,0,1,1,1))
        self.assertEqual([e.net_change for e in self.audit.events],[0,0,0,1,0,0])

    def test_irrational_true_cutoff_time(self):
        c=self.audit.clusters[3]
        self.assertIsNone(c.exact_time)
        lo,hi=c.bracket.lower,c.bracket.upper
        f=lambda t:(t+1)**2-2
        self.assertLess(f(lo),0)
        self.assertGreater(f(hi),0)
        self.assertEqual(self.audit.events[3].net_change,1)
        self.assertTrue(self.audit.events[3].has_parity_change)

    def test_scaling_has_no_geometric_birth(self):
        all_orbits=census_polygon(self.P0)
        self.assertEqual(all_orbits.count,1)
        self.assertEqual(all_orbits.squares[0].side_squared,Q(72,25))
        self.assertEqual(self.audit.chambers[-1].count,1)

    def test_no_event_root_is_replaced_by_decimal_guess(self):
        for e in self.audit.events:
            if e.cluster.exact_time is None:
                self.assertIsNone(e.at_event_count)
                self.assertLess(e.cluster.bracket.lower,e.cluster.bracket.upper)
        self.assertEqual(self.audit.unresolved_chart_count,0)


if __name__=='__main__':unittest.main()
