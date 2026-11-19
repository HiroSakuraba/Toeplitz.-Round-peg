"""v225 exact persistent-face normalization and moving singular-rank certificates."""
import random
import sys
import unittest
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations,product
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from persistent_strata_v225 import (
    adjacent_two_pairs,affine_system_polynomials,generic_polynomial_rank,
    singular_inconsistency_certificate,classify_chart,owned_square_at_rational_time,
    singular_chart_snapshot,classify_polygon_motion,generic_snapshot_check,report,
)
from moving_polygon_events_v224 import (
    polynomial,ZERO,ONE,T, evaluate,chart_polynomials,polygon_event_cover,
    interpolated_polygon,certify_jordan_isotopy,
)
from polygon_orbit_census_v223 import (
    pt,_affine_system,_rank,_box_vertices,_evaluate,_side_sq,
    validate_jordan_polygon,census_polygon,PositiveDimensionalSquares,
    convex_vertex_alignment_polygon,
)


class SingularAlgebraTests(unittest.TestCase):
    def test_polynomial_rank_does_not_use_one_bad_time(self):
        A=((T,ZERO),(ZERO,ONE))
        self.assertEqual(generic_polynomial_rank(A),2)
        self.assertEqual(generic_polynomial_rank(((ZERO,ZERO),(ZERO,ZERO))),0)

    def test_singular_inconsistency_whole_interval(self):
        A=((ONE,ZERO,ZERO,ZERO),(ZERO,ONE,ZERO,ZERO),
           (ZERO,ZERO,ONE,ZERO),(ZERO,ZERO,ZERO,ZERO))
        b=(ZERO,ZERO,ZERO,ONE)
        r,s,w,roots=singular_inconsistency_certificate(A,b)
        self.assertEqual((r,s),(3,4))
        self.assertEqual(roots,())
        self.assertNotEqual(evaluate(w,Q(1,2)),0)

    def test_exceptional_rational_time_has_no_false_existence_claim(self):
        A=((ONE,ZERO,ZERO,ZERO),(ZERO,ONE,ZERO,ZERO),
           (ZERO,ZERO,ONE,ZERO),(ZERO,ZERO,ZERO,ZERO))
        b=(ZERO,ZERO,ZERO,polynomial((-Q(1,3),1)))
        r,s,w,roots=singular_inconsistency_certificate(A,b)
        self.assertEqual((r,s),(3,4))
        self.assertEqual(len(roots),1)
        self.assertLess(roots[0].lower,Q(1,3))
        self.assertGreater(roots[0].upper,Q(1,3))
        # At the root the system actually has an affine family, and we
        # deliberately report only 'possible' rather than an orbit parity.
        self.assertEqual(evaluate(w,Q(1,3)),0)
        self.assertNotEqual(evaluate(w,Q(1,2)),0)

    def test_exceptional_irrational_time_is_isolated_exactly(self):
        A=((ONE,ZERO,ZERO,ZERO),(ZERO,ONE,ZERO,ZERO),
           (ZERO,ZERO,ONE,ZERO),(ZERO,ZERO,ZERO,ZERO))
        b=(ZERO,ZERO,ZERO,polynomial((-Q(1,2),0,1)))
        r,s,w,roots=singular_inconsistency_certificate(A,b,width=Q(1,2**18))
        self.assertEqual((r,s),(3,4))
        self.assertEqual(len(roots),1)
        z=roots[0]
        self.assertFalse(z.rational)
        self.assertLess(z.upper-z.lower,Q(1,2**17))
        self.assertLess(evaluate(w,z.lower)*evaluate(w,z.upper),0)

    def test_consistent_generic_singular_system_not_dismissed(self):
        A=((ONE,ZERO,ZERO,ZERO),(ZERO,ONE,ZERO,ZERO),
           (ZERO,ZERO,ONE,ZERO),(ZERO,ZERO,ZERO,ZERO))
        b=(ZERO,ZERO,ZERO,ZERO)
        r,s,w,roots=singular_inconsistency_certificate(A,b)
        self.assertEqual((r,s),(3,3))
        self.assertIsNone(w)
        self.assertEqual(roots,())

    def test_non_polynomial_4x4_data_rejected(self):
        with self.assertRaises(ValueError):
            singular_inconsistency_certificate(((ONE,),),(ONE,))
        A=((ONE,ZERO,ZERO,ZERO),)*4
        with self.assertRaises(ValueError):
            singular_inconsistency_certificate(A,(ZERO,)*4,1,0)


class AdjacentTwoPairTests(unittest.TestCase):
    def test_all_pair_labelings_excluded(self):
        for n in (4,5,8):
            for a in range(n):
                b=(a+1)%n
                for ids in set(permutations((a,a,b,b))):
                    self.assertTrue(adjacent_two_pairs(ids,n))
        self.assertFalse(adjacent_two_pairs((0,2,0,2),4))
        self.assertFalse(adjacent_two_pairs((0,0,0,1),4))

    def test_geometry_exclusion_independently_by_static_feasibility(self):
        polygons=[tuple(pt(*v) for v in points) for points in (
            ((0,0),(3,1),(2,3),(-1,2)),
            ((0,0),(2,0),(3,3),(0,2)),
            ((0,0),(2,0),(3,1),(2,2),(0,1)),
        )]
        for P in polygons:
            validate_jordan_polygon(P)
            edges=tuple((P[i],P[(i+1)%len(P)]) for i in range(len(P)))
            for i in range(len(P)):
                j=(i+1)%len(P)
                for ids in set(permutations((i,i,j,j))):
                    for sigma in (1,-1):
                        A,b=_affine_system(tuple(edges[k] for k in ids),sigma)
                        roots=list(_box_vertices(A,b))
                        self.assertFalse(any(_side_sq(_evaluate(tuple(edges[k] for k in ids),u))>0
                                             for u in roots), (P,ids,sigma))

    def test_adjacent_exclusion_also_used_for_moving_deformation(self):
        P=convex_vertex_alignment_polygon(Q(0))
        R=convex_vertex_alignment_polygon(Q(9,8))
        for ids in ((0,0,1,1),(0,3,0,3),(3,3,0,0)):
            for sigma in (-1,1):
                c=classify_chart(P,R,ids,sigma,Q(1,100))
                self.assertEqual(c.status,'structurally-empty')
                self.assertTrue(c.resolved)


class PersistentFacetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=convex_vertex_alignment_polygon(Q(0))
        cls.R=convex_vertex_alignment_polygon(Q(9,8))
        cls.global_audit=classify_polygon_motion(cls.P,cls.R,Q(1,100),Q(3,4),1)

    def test_all_408_charts_resolved_in_convex_model(self):
        c=self.global_audit
        self.assertEqual(len(c.certificates),408)
        self.assertEqual(dict(c.status_counts),{
            'regular':280,'vertex-owned':40,
            'incoming-edge-duplicate':40,'structurally-empty':48})
        self.assertTrue(c.completely_stratified)
        self.assertEqual(c.unresolved_charts,())

    def test_v224_128_unresolved_reclassified_not_silently_discarded(self):
        old,unresolved=polygon_event_cover(self.P,self.R,Q(1,100),
                             Q(3,4),1,Q(1,2**12))
        self.assertEqual(len(old),408)
        self.assertEqual(len(unresolved),128)
        now={(x.edges,x.orientation):x for x in self.global_audit.certificates}
        for ids,sigma,reason in unresolved:
            self.assertTrue(now[(ids,sigma)].resolved)
            self.assertNotEqual(now[(ids,sigma)].status,'regular')

    def test_regular_and_persistent_facets_still_cover_physical_square(self):
        for t,expected in ((Q(3,4),1),(Q(17,20),1),
                           (Q(8,9),2),(Q(11,12),3),(Q(1),3)):
            v=generic_snapshot_check(self.P,self.R,t,Q(1,100))
            self.assertTrue(v['matches'],t)
            self.assertEqual(len(v['known']),expected)
            self.assertEqual(len(v['found']),expected)
            self.assertEqual(v['unresolved'],())

    def test_persistent_start_endpoint_assigns_square_once_per_d4_orbit(self):
        P=tuple(pt(*v) for v in ((0,0),(2,0),(2,Q(3,2)),(0,Q(1,2))))
        self.assertEqual(census_polygon(P).count,2)
        out=generic_snapshot_check(P,P,Q(1,2),Q(1,100))
        self.assertTrue(out['matches']);self.assertEqual(len(out['found']),2)
        ids=(0,1,1,2)
        c=classify_chart(P,P,ids,1,Q(1,100))
        self.assertEqual(c.status,'vertex-owned')
        self.assertEqual(c.forced_start,(1,))
        result=owned_square_at_rational_time(P,P,ids,1,Q(1,2),Q(1,100))
        self.assertIsNotNone(result)
        self.assertEqual(result[1],(Q(1,2),Q(0),Q(2,3),Q(1,2)))
        self.assertEqual(result[2],Q(1))

    def test_persistent_incoming_endpoint_is_not_counted_twice(self):
        P=tuple(pt(*v) for v in ((0,0),(2,0),(2,Q(3,2)),(0,Q(1,2))))
        # The incoming edge describes a TRUE square at u_1=1.
        ids=(0,0,1,2)
        c=classify_chart(P,P,ids,1,Q(1,100))
        self.assertEqual(c.status,'incoming-edge-duplicate')
        self.assertEqual(c.forced_end,(1,))
        algebraic=chart_polynomials(P,P,ids,1,Q(1,100))
        self.assertEqual(algebraic.algebraic_parameters(Q(1,2)),
                         (Q(1,2),Q(1),Q(2,3),Q(1,2)))
        self.assertIsNone(owned_square_at_rational_time(P,P,ids,1,Q(1,2),Q(1,100)))
        # The same geometric square is genuinely present on the outgoing
        # edge chart (0,1,1,2), with that parameter now zero.
        outgoing=owned_square_at_rational_time(P,P,(0,1,1,2),1,Q(1,2),Q(1,100))
        self.assertIsNotNone(outgoing)
        self.assertEqual(outgoing[1][1],Q(0))
        self.assertTrue(generic_snapshot_check(P,P,Q(1,2),Q(1,100))['matches'])

    def test_root_cover_is_finite_for_canonical_regular_charts(self):
        c=self.global_audit
        self.assertGreater(c.raw_candidate_events,0)
        self.assertTrue(all(not x.exceptional_roots for x in c.certificates
                            if x.status in ('incoming-edge-duplicate','structurally-empty')))
        self.assertTrue(all(all(z!=ZERO for name,z in x.active_polynomials)
                            for x in c.certificates if x.has_regular_branch))

    def test_exact_vertex_alignment_is_preserved(self):
        target=(pt(1,0),pt(2,0),pt(2,1),pt(1,1))
        t=Q(8,9)
        actual=generic_snapshot_check(self.P,self.R,t,Q(1,100))
        self.assertIn(tuple(sorted(target)),actual['found'])
        self.assertEqual(len(actual['found']),2)

    def test_cutoff_filter_can_remove_all_small_squares(self):
        for t in (Q(0),Q(1,3),Q(8,9),Q(1)):
            x=generic_snapshot_check(self.P,self.R,t,Q(1000))
            self.assertTrue(x['matches'])
            self.assertEqual(x['found'],set())

    def test_whole_interval_embedded(self):
        self.assertEqual(certify_jordan_isotopy(self.P,self.R).triples_checked,4)

    def test_half_open_ownership_agrees_with_reversed_polygon(self):
        t=Q(3,4)
        a=generic_snapshot_check(self.P,self.R,t,Q(1,100))
        b=generic_snapshot_check(tuple(reversed(self.P)),tuple(reversed(self.R)),t,Q(1,100))
        self.assertEqual(a['known'],b['known'])
        self.assertEqual(a['found'],b['found'])
        self.assertTrue(a['matches'] and b['matches'])

    def test_rational_symmetric_perturbations_keep_exact_snapshot_census(self):
        rng=random.Random(225)
        for _ in range(1):
            P=tuple(pt(v.x+Q(rng.randint(-1,1),80),
                       v.y+Q(rng.randint(-1,1),80)) for v in self.P)
            R=tuple(pt(v.x+Q(rng.randint(-1,1),80),
                       v.y+Q(rng.randint(-1,1),80)) for v in self.R)
            certify_jordan_isotopy(P,R)
            for t in (Q(3,4),):
                result=generic_snapshot_check(P,R,t,Q(1,100))
                self.assertTrue(result['matches'])
                self.assertEqual(result['unresolved'],())

    def test_bad_polygon_isotopy_rejected(self):
        P=self.P
        Qbad=(P[0],P[2],P[1],P[3])
        with self.assertRaises(ValueError):
            classify_polygon_motion(P,Qbad,Q(1,100))

    def test_resource_budget_is_enforced(self):
        with self.assertRaises(RuntimeError):
            classify_polygon_motion(self.P,self.R,Q(1,100),max_charts=10)


class RankDeficientExamples(unittest.TestCase):
    def test_rectangle_parallel_diagonals_generically_inconsistent(self):
        P=tuple(pt(*v) for v in ((0,0),(2,0),(2,1),(0,1)))
        c=classify_chart(P,P,(0,2,0,2),1,Q(1,64))
        self.assertEqual(c.status,'always-empty')
        self.assertEqual((c.rank_A,c.rank_augmented),(3,4))
        self.assertFalse(c.exceptional_roots)

    def test_rectangle_true_square_continuum_not_discarded(self):
        P=tuple(pt(*v) for v in ((0,0),(2,0),(2,1),(0,1)))
        c=classify_chart(P,P,(0,0,2,2),1,Q(1,64))
        self.assertEqual(c.status,'genuinely-unresolved')
        self.assertFalse(c.resolved)
        s=singular_chart_snapshot(P,(0,0,2,2),1,Q(1,64))
        self.assertEqual(s.status,'positive-dimensional-family')
        self.assertGreater(s.max_side_squared,Q(1,64))
        with self.assertRaises(PositiveDimensionalSquares):census_polygon(P)

    def test_singular_chart_below_huge_cutoff_is_empty_at_snapshot(self):
        P=tuple(pt(*v) for v in ((0,0),(2,0),(2,1),(0,1)))
        s=singular_chart_snapshot(P,(0,0,2,2),1,Q(10000))
        self.assertEqual(s.status,'no-positive-size-square')
        self.assertLessEqual(s.max_side_squared,Q(10000))

    def test_moving_affine_system_matches_static_rank_and_feasibility(self):
        P=tuple(pt(*v) for v in ((0,0),(2,0),(3,3),(0,2)))
        R=tuple(pt(*v) for v in ((0,0),(2,0),(3,2),(0,2)))
        ids=(0,2,0,2)
        for sigma in (-1,1):
            A,b=affine_system_polynomials(P,R,ids,sigma)
            for t in (Q(0),Q(1,4),Q(1,2),Q(1)):
                moved=interpolated_polygon(P,R,t)
                edges=tuple((moved[i],moved[(i+1)%4]) for i in ids)
                As,bs=_affine_system(edges,sigma)
                self.assertEqual(tuple(tuple(evaluate(z,t) for z in row) for row in A),
                                 tuple(tuple(row) for row in As))
                self.assertEqual(tuple(evaluate(z,t) for z in b),bs)

    def test_report_contains_scientific_scope_and_no_unresolved_for_control(self):
        r=report()
        self.assertEqual(r['unresolved_charts'],0)
        self.assertEqual(r['total_charts'],408)
        self.assertTrue(r['sample_t_3_4'])
        self.assertTrue(r['sample_t_11_12'])
        self.assertIn('no global parity',r['scope'])

if __name__=='__main__':unittest.main()
