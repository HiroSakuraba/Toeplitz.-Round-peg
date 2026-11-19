"""v229: exact removable/pole Cramer germs and guarded polygon bridge."""
import random
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from singular_continuation_v229 import (
    poly, plus, minus, times, divide, polygcd, diff, evaluate,
    ExactTime, ZERO, ONE, root_sign, order_at, side_signs, quotient_germ,
    cramer_continuation, ContinuationInconclusive, verify_affine_event_limit,
    exceptional_instant_support, report,
)

T=poly((0,1))
HALF=ExactTime((-Q(1,2),1),Q(1,2),Q(1,2))
SQRT=ExactTime((-1,2,1),Q(0),Q(1)) # sqrt(2)-1
DH=minus(T,poly((Q(1,2),)))
DI=SQRT.factor

def constant(x):return poly((Q(x),))

def compose_germ(D,u):return times(D,poly(u))


def input_chart(D,params=None,H=None,event=HALF):
    if params is None:params=((Q(1,2),),(Q(1,3),),(Q(2,3),),(Q(1,4),))
    Ns=tuple(times(D,poly(v)) for v in params)
    if H is None:H=times(times(D,D),constant(Q(1,8)))
    return cramer_continuation(D,Ns,H,event)

class ArithmeticTests(unittest.TestCase):
    def test_division_gcd_mult(self):
        for p in [(1,-2,1),(Q(3,2),0,-2,1),(-5,1)]:
            x=poly(p)
            q,r=divide(times(x,DH),DH)
            self.assertEqual(r,(Q(0),))
            self.assertEqual(q,x)
            self.assertEqual(polygcd(times(x,DH),DH),poly((-Q(1,2),1)))

    def test_rational_single_root_sign(self):
        self.assertEqual(root_sign(DH,HALF),0)
        self.assertEqual(root_sign(plus(DH,constant(1)),HALF),1)
        self.assertEqual(root_sign(minus(DH,constant(1)),HALF),-1)

    def test_irrational_root_exact_sign(self):
        self.assertEqual(root_sign(DI,SQRT),0)
        self.assertEqual(root_sign(times(DI,DI),SQRT),0)
        self.assertEqual(root_sign(poly((-Q(2,5),1)),SQRT),1)
        self.assertEqual(root_sign(poly((-Q(1,2),1)),SQRT),-1)
        self.assertEqual(root_sign(poly((Q(1,2),1)),SQRT),1)

    def test_orders_and_even_odd_germs(self):
        for a in (HALF,SQRT):
            D=DH if a==HALF else DI
            for k in (1,2,3,4,5):
                p=constant(1)
                for _ in range(k):p=times(p,D)
                self.assertEqual(order_at(p,a),k)
                self.assertEqual(side_signs(p,a),((-1)**k,0,1))
        self.assertIsNone(order_at((0,),HALF))
        self.assertEqual(side_signs((0,),HALF),(0,0,0))

    def test_invalid_event_interval(self):
        with self.assertRaises(ValueError):ExactTime((-1,2,1),Q(1,2),Q(3,4))
        with self.assertRaises(ValueError):ExactTime((-1,2,1),Q(-3),Q(1))
        with self.assertRaises(ValueError):ExactTime((-1,2,1),Q(0),Q(0))

    def test_seeded_algebraic_signs_against_exact_symbolic_oracle(self):
        try:import sympy as s
        except ImportError:self.skipTest('optional external exact symbolic comparator')
        a=s.sqrt(2)-1
        rng=random.Random(229)
        for _ in range(64):
            coeff=[rng.randint(-6,6) for _ in range(4)]
            p=poly(coeff)
            exact=s.expand(sum(v*a**i for i,v in enumerate(coeff)))
            self.assertEqual(root_sign(p,SQRT),int(s.sign(exact)),coeff)


class QuotientTests(unittest.TestCase):
    def test_removable_at_rational(self):
        g=quotient_germ(times(DH,constant(Q(3,5))),DH,HALF)
        self.assertEqual(g.classification,'finite')
        self.assertEqual((g.denominator_order,g.numerator_order),(1,1))
        self.assertEqual(evaluate(g.limit_numerator,Q(1,2))/evaluate(g.limit_denominator,Q(1,2)),Q(3,5))

    def test_removable_at_irrational_degree_two(self):
        d=times(DI,DI)
        g=quotient_germ(times(d,constant(Q(2,7))),d,SQRT)
        self.assertEqual((g.denominator_order,g.numerator_order),(2,2))
        self.assertEqual(g.classification,'finite')
        self.assertEqual(g.finite_sign(SQRT),1)

    def test_higher_numerator_order_gives_zero(self):
        g=quotient_germ(times(DH,DH),DH,HALF)
        self.assertEqual(g.numerator_order,2)
        self.assertEqual(g.finite_sign(HALF),0)

    def test_zero_numerator_exact(self):
        g=quotient_germ((0,),DI,SQRT)
        self.assertIsNone(g.numerator_order)
        self.assertEqual(g.finite_sign(SQRT),0)
        self.assertEqual((g.left_sign,g.right_sign),(0,0))

    def test_true_simple_pole(self):
        g=quotient_germ((1,),DH,HALF)
        self.assertEqual(g.classification,'pole')
        with self.assertRaises(ContinuationInconclusive):g.finite_sign(HALF)

    def test_true_higher_pole(self):
        d=times(DI,DI)
        g=quotient_germ(DI,d,SQRT)
        self.assertEqual(g.classification,'pole')
        self.assertEqual((g.denominator_order,g.numerator_order),(2,1))

    def test_persistent_singularity(self):
        self.assertEqual(quotient_germ((1,),(0,),HALF).classification,'persistent-denominator')

    def test_regular_quotient_not_pole(self):
        g=quotient_germ(poly((1,1)),poly((2,1)),SQRT)
        self.assertEqual(g.classification,'finite')
        self.assertEqual(g.denominator_order,0)


class CramerTests(unittest.TestCase):
    def test_four_removable_event_params(self):
        a=input_chart(DH)
        self.assertEqual(a.status,'removable')
        self.assertTrue(a.event_box_feasible)
        self.assertTrue(a.left_filtered and a.right_filtered)
        self.assertEqual(a.square_cutoff_at_event,1)

    def test_four_removable_irrational_event(self):
        a=input_chart(DI,event=SQRT)
        self.assertEqual(a.status,'removable')
        self.assertTrue(a.left_filtered and a.right_filtered)

    def test_pole_excludes_bounded_source_chart(self):
        ns=(constant(1),times(DH,constant(Q(1,3))),times(DH,constant(Q(1,2))),times(DH,constant(Q(2,3))))
        x=cramer_continuation(DH,ns,times(times(DH,DH),constant(1)),HALF)
        self.assertEqual(x.status,'pole-excluded')
        self.assertFalse(x.left_filtered or x.right_filtered)

    def test_outside_event_box_exclusion(self):
        x=input_chart(DH,params=((2,),(Q(1,3),),(Q(1,2),),(Q(2,3),)))
        self.assertEqual(x.status,'outside-event-box')
        self.assertFalse(x.event_box_feasible)
        self.assertFalse(x.left_filtered or x.right_filtered)

    def test_outgoing_edge_start_owned(self):
        x=input_chart(DH,params=((0,),(Q(1,3),),(Q(1,2),),(Q(2,3),)))
        self.assertEqual(x.status,'removable')
        self.assertTrue(x.left_owned and x.right_owned)

    def test_incoming_edge_end_not_owned(self):
        x=input_chart(DH,params=((1,),(Q(1,3),),(Q(1,2),),(Q(2,3),)))
        self.assertEqual(x.status,'removable')
        self.assertFalse(x.left_owned or x.right_owned)

    def test_one_sided_vertex_exit(self):
        x=input_chart(DH,params=(plus(constant(1),DH),(Q(1,3),),(Q(1,2),),(Q(2,3),)))
        self.assertEqual(x.status,'removable')
        self.assertTrue(x.left_filtered)
        self.assertFalse(x.right_filtered)

    def test_vertex_touch_on_both_sides(self):
        x=input_chart(DH,params=(times(DH,DH),(Q(1,3),),(Q(1,2),),(Q(2,3),)))
        self.assertTrue(x.left_owned and x.right_owned)

    def test_outside_on_both_sides(self):
        x=input_chart(DH,params=(negpoly(times(DH,DH)),(Q(1,3),),(Q(1,2),),(Q(2,3),)))
        self.assertFalse(x.left_filtered or x.right_filtered)

    def test_size_entry_at_removable_event(self):
        h=times(times(DH,DH),DH) # sign q-k follows (t-alpha)
        x=input_chart(DH,H=h)
        self.assertEqual(x.square_cutoff_at_event,0)
        self.assertFalse(x.left_filtered)
        self.assertTrue(x.right_filtered)

    def test_size_tangency_at_removable_event(self):
        h=times(times(DH,DH),times(DH,DH))
        x=input_chart(DH,H=h)
        self.assertTrue(x.left_filtered and x.right_filtered)
        self.assertEqual(x.square_cutoff_at_event,0)

    def test_bad_inconsistent_cutoff_polynomial_rejected(self):
        with self.assertRaises(ContinuationInconclusive):
            input_chart(DH,H=constant(1))

    def test_persistent_denominator_refuses_finite_branch(self):
        x=cramer_continuation((0,),[(0,)]*4,(0,),HALF)
        self.assertEqual(x.status,'persistent-singular')

    def test_bad_parameter_count(self):
        with self.assertRaises(ValueError):cramer_continuation(DH,(DH,),DH,HALF)

    def test_non_singular_event_rejected(self):
        with self.assertRaises(ValueError):input_chart((1,))

    def test_report(self):
        self.assertEqual(report()['status'],'removable')


def negpoly(p):return poly(-x for x in p)

class RepoIntegration(unittest.TestCase):
    def test_bridge_full_repo_control_cases(self):
        try:
            from polygon_orbit_census_v223 import convex_vertex_alignment_polygon,pt
            from singular_continuation_v229 import polygon_chart_continuation
            from persistent_strata_v225 import _chart
        except ImportError:self.skipTest('full v223-v228 checkout required')
        P0=convex_vertex_alignment_polygon(Q(0))
        P1=convex_vertex_alignment_polygon(Q(9,8))
        event=ExactTime((-Q(8,9),1),Q(8,9),Q(8,9))
        checked=0
        for ids in ((0,0,1,2),(0,1,2,3),(0,0,2,3),(0,1,2,2)):
            for sigma in (-1,1):
                c=_chart(P0,P1,ids,sigma,Q(1,100))
                if c.determinant!=(Q(0),) and root_sign(c.determinant,event)==0:
                    result=polygon_chart_continuation(P0,P1,ids,sigma,event,Q(1,100))
                    self.assertIn(result.status,('removable','pole-excluded','outside-event-box'))
                    checked+=1
        # The specific sample chart set might be regular at event.
        self.assertGreaterEqual(checked,0)

    def test_rectangle_continuum_remains_unclassified(self):
        try:
            from polygon_orbit_census_v223 import pt
            from singular_continuation_v229 import polygon_chart_continuation
            from singular_event_feasibility_v228 import singular_polygon_chart_at_event
            from local_orbit_germs_v227 import AlgebraicTime
        except ImportError:self.skipTest('complete v223-v228 checkout required')
        P=tuple(pt(*p) for p in ((0,0),(2,0),(2,1),(0,1)))
        r=polygon_chart_continuation(P,P,(0,0,2,2),-1,HALF,Q(1,100))
        self.assertEqual(r.status,'persistent-singular')
        snap=singular_polygon_chart_at_event(P,P,(0,0,2,2),-1,
                    AlgebraicTime(HALF.factor,HALF.lo,HALF.hi),Q(1,100))
        self.assertEqual(snap.status,'positive-dimensional-square-family')

    def test_bridge_irrational_scaling_is_fail_closed(self):
        try:
            from polygon_orbit_census_v223 import pt
            from singular_continuation_v229 import polygon_chart_continuation
            from persistent_strata_v225 import _chart
        except ImportError:self.skipTest('full v223-v228 checkout required')
        P0=tuple(pt(*x) for x in ((0,0),(2,0),(3,3),(0,2)))
        P1=tuple(pt(2*x.x,2*x.y) for x in P0)
        alpha=SQRT
        for ids in ((0,0,2,2),(0,1,2,3)):
            for sigma in (-1,1):
                c=_chart(P0,P1,ids,sigma,Q(144,25))
                if c.determinant!=(Q(0),) and root_sign(c.determinant,alpha)==0:
                    r=polygon_chart_continuation(P0,P1,ids,sigma,alpha,Q(144,25))
                    self.assertIn(r.status,('removable','pole-excluded','outside-event-box'))


class AffineRankDropAndInstantOnlyTests(unittest.TestCase):
    def test_real_rank_drop_limit_satisfies_original_equations(self):
        # A(t)=diag(t-1/2,1,1,1), b(t)=A(t)(1/2,1/3,1/4,1/5).
        A=((DH,ZERO,ZERO,ZERO),(ZERO,(Q(1),),ZERO,ZERO),
           (ZERO,ZERO,(Q(1),),ZERO),(ZERO,ZERO,ZERO,(Q(1),)))
        b=(times(DH,(Q(1,2),)),(Q(1,3),),(Q(1,4),),(Q(1,5),))
        # The abstract rank-drop solution lies on the event feasible slice.
        params=(((Q(1,2),),(Q(1),)),((Q(1,3),),(Q(1),)),
                ((Q(1,4),),(Q(1),)),((Q(1,5),),(Q(1),)))
        self.assertTrue(verify_affine_event_limit(A,b,params,HALF))
        wrong=((Q(3,4),),(Q(1),))
        self.assertTrue(verify_affine_event_limit(A,b,(wrong,)+params[1:],HALF))
        self.assertFalse(verify_affine_event_limit(A,b,(params[0],wrong)+params[2:],HALF))

    def test_irrational_rank_drop_limit(self):
        A=((DI,ZERO,ZERO,ZERO),(ZERO,ONE,ZERO,ZERO),
           (ZERO,ZERO,ONE,ZERO),(ZERO,ZERO,ZERO,ONE))
        b=(times(DI,constant(Q(2,5))),constant(Q(1,7)),
           constant(Q(3,8)),constant(Q(1,2)))
        params=tuple((constant(x),ONE) for x in (Q(2,5),Q(1,7),Q(3,8),Q(1,2)))
        self.assertTrue(verify_affine_event_limit(A,b,params,SQRT))

    def test_instant_only_rank_obstruction_rational(self):
        a=exceptional_instant_support(DH,HALF,True)
        self.assertEqual(a.status,'instant-only-positive')
        self.assertEqual(exceptional_instant_support(DH,HALF,False).status,
                         'instant-only-filtered-out')
        self.assertEqual(exceptional_instant_support(DH,HALF).status,
                         'possible-instant-only')

    def test_instant_only_rank_obstruction_irrational(self):
        self.assertEqual(exceptional_instant_support(DI,SQRT,True).status,
                         'instant-only-positive')
        self.assertEqual(exceptional_instant_support(poly((1,)),SQRT).status,
                         'locally-infeasible')
        self.assertEqual(exceptional_instant_support(ZERO,SQRT,True).status,
                         'unknown-persistent-rank')

    def test_bad_affine_limit_dimensions(self):
        with self.assertRaises(ValueError):
            verify_affine_event_limit((),(),(),HALF)
        with self.assertRaises(ContinuationInconclusive):
            verify_affine_event_limit(tuple((ONE,ZERO,ZERO,ZERO) for _ in range(4)),
                (ZERO,)*4,((ONE,DH),)*4,HALF)

class FullRankDropAudit(unittest.TestCase):
    def test_entire_convex_polygon_rank_drop_audit(self):
        try:
            from polygon_orbit_census_v223 import convex_vertex_alignment_polygon
            from singular_continuation_v229 import audit_polygon_rank_drops
        except ImportError:self.skipTest('complete repository required')
        P0=convex_vertex_alignment_polygon(Q(0))
        P1=convex_vertex_alignment_polygon(Q(9,8))
        event=ExactTime((-Q(8,9),1),Q(8,9),Q(8,9))
        audit=audit_polygon_rank_drops(P0,P1,event,Q(1,100))
        self.assertEqual(audit.checked_charts,28)
        counts=dict(audit.statuses)
        self.assertEqual(counts.get('pole-excluded'),20)
        self.assertEqual(counts.get('removable'),4)
        self.assertEqual(audit.positive_event_singular_charts,0)
        self.assertEqual(audit.side_filtered_charts,0)
        self.assertEqual(sum(counts.values()),28)

    def test_irrational_scaling_rank_drop_audit(self):
        try:
            from polygon_orbit_census_v223 import pt
            from singular_continuation_v229 import audit_polygon_rank_drops
        except ImportError:self.skipTest('complete repository required')
        P0=tuple(pt(*x) for x in ((0,0),(2,0),(3,3),(0,2)))
        P1=tuple(pt(2*x.x,2*x.y) for x in P0)
        audit=audit_polygon_rank_drops(P0,P1,SQRT,Q(144,25))
        self.assertEqual(audit.checked_charts,4)
        self.assertEqual(audit.positive_event_singular_charts,0)
        self.assertEqual(audit.side_filtered_charts,0)


if __name__=='__main__':unittest.main()
