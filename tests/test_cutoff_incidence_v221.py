"""v221 exact cutoff/parameter-boundary incidence and adversarial controls."""
import sys
import unittest
import random
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from cutoff_incidence_v221 import (
    pt,sub,add,scale,rot,dot,residual,MovingSegment,NonGenericIncidence,
    LinearSquareBranch,planted_square_family,linear_square_branch,
    incidence_ledger,shrinking_notch_tagged_count,polynomial_value,
    _crossings_of_quadratic,report,
)


class IncidenceTests(unittest.TestCase):
    def test_linear_square_branch_all_rational_times(self):
        f=planted_square_family(pt(Q(1,8),Q(1,8)),pt(Q(1,9),-Q(1,12)))
        self.assertEqual(f.params0,(Q(1,2),)*4)
        self.assertEqual(f.params1,(Q(0),)*4)
        for t in (Q(0),Q(1,13),Q(1,2),Q(5,6),Q(1)):
            p=f.square(t)
            self.assertEqual(residual(p,1),(0,0,0,0))
            self.assertEqual(dot(sub(p[1],p[0]),sub(p[1],p[0])),f.side_squared(t))
            self.assertGreater(f.side_squared(t),0)

    def test_one_exact_cutoff_crossing(self):
        f=planted_square_family(pt(0,0),pt(Q(1,8),Q(1,8)))
        d=incidence_ledger(f,Q(1,4),Q(3,4),Q(1,128))
        self.assertEqual(f.cutoff_polynomial(Q(1,128)),(Q(1,32),Q(0),-Q(1,128)))
        self.assertEqual(f.side_squared(Q(1,2)),Q(1,128))
        self.assertEqual((d.start_count,d.end_count),(0,1))
        self.assertEqual((d.source_events,d.cutoff_crossings),(0,1))
        self.assertTrue(d.balances)

    def test_two_cutoff_crossings_have_even_parity(self):
        f=planted_square_family(pt(Q(1,4),0),pt(-Q(1,2),0))
        d=incidence_ledger(f,Q(1,4),Q(3,4),Q(1,256))
        self.assertEqual((d.start_count,d.end_count),(1,1))
        self.assertEqual((d.source_events,d.cutoff_crossings),(0,2))
        self.assertEqual(f.side_squared(Q(3,8)),Q(1,256))
        self.assertEqual(f.side_squared(Q(5,8)),Q(1,256))
        self.assertEqual(f.side_squared(Q(1,2)),0)
        self.assertTrue(d.balances)

    def test_irrational_cutoff_crossing_decided_with_rationals(self):
        f=planted_square_family(pt(Q(1,8),0),pt(0,Q(1,8)))
        k=Q(3,128)
        d=incidence_ledger(f,Q(0),Q(1),k)
        self.assertEqual(f.cutoff_polynomial(k),(Q(1,64),Q(0),-Q(1,128)))
        self.assertEqual((d.start_count,d.end_count),(0,1))
        self.assertEqual(d.cutoff_crossings,1)
        self.assertTrue(d.balances)
        # The time is sqrt(1/2), which is irrational; no float approximation used.
        self.assertEqual(polynomial_value(f.cutoff_polynomial(k),Q(1,2)),-Q(1,256))
        self.assertEqual(polynomial_value(f.cutoff_polynomial(k),Q(3,4)),Q(1,1024))

    def test_source_wall_exit_is_not_a_size_crossing(self):
        f=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0))
        x=list(f.segments)
        x[1]=MovingSegment(x[1].start,x[1].end,pt(0,1))
        g=linear_square_branch(tuple(x))
        d=incidence_ledger(g,Q(0),Q(1,5),Q(1,1024))
        self.assertEqual((d.start_count,d.end_count),(1,0))
        self.assertEqual(d.source_wall_times,(Q(55,19968),))
        self.assertEqual(d.cutoff_crossings,0)
        self.assertTrue(d.balances)
        self.assertGreater(g.side_squared(Q(55,19968)),Q(1,1024))
        self.assertTrue(any(u in (0,1) for u in g.parameters(Q(55,19968))))

    def test_fixed_side_translation_has_no_events(self):
        f=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0))
        d=incidence_ledger(f,0,1,Q(1,256))
        self.assertEqual((d.start_count,d.end_count,d.source_events,d.cutoff_crossings),
                         (1,1,0,0))

    def test_tangency_cutoff_refused_not_false_transverse_event(self):
        f=planted_square_family(pt(Q(1,8),-Q(1,16)),pt(0,Q(1,8)))
        self.assertEqual(f.side_squared(Q(1,2)),Q(1,64))
        with self.assertRaisesRegex(NonGenericIncidence,'tangential'):
            incidence_ledger(f,Q(0),Q(1),Q(1,64))

    def test_boundary_cutoff_coincidence_refused(self):
        f=planted_square_family(pt(0,0),pt(Q(1,8),Q(1,8)))
        with self.assertRaisesRegex(NonGenericIncidence,'time slice'):
            incidence_ledger(f,Q(1,2),Q(3,4),Q(1,128))

    def test_bad_cutoffs_time_intervals_rejected(self):
        f=planted_square_family(pt(0,0),pt(Q(1,8),Q(1,8)))
        for args in ((0,1,0),(1,0,Q(1,128)),(Q(1,2),Q(1,2),Q(1,128))):
            with self.assertRaises(ValueError):incidence_ledger(f,*args)
        for k in (0,-Q(1,4)):
            with self.assertRaises(ValueError):shrinking_notch_tagged_count(k)

    def test_singular_chart_refused(self):
        p=pt(0,0)
        segs=tuple(MovingSegment(p,pt(1,0),pt(0,0)) for _ in range(4))
        with self.assertRaisesRegex(NonGenericIncidence,'singular'):
            linear_square_branch(segs)

    def test_constant_zero_cutoff_polynomial_detected(self):
        with self.assertRaisesRegex(NonGenericIncidence,'entire branch'):
            _crossings_of_quadratic((Q(0),Q(0),Q(0)),Q(0),Q(1))

    def test_quadratic_internal_double_root_refused(self):
        with self.assertRaisesRegex(NonGenericIncidence,'tangential'):
            _crossings_of_quadratic((Q(1),-Q(1),Q(1,4)),Q(0),Q(1))

    def test_two_irregular_quadratic_crossings_without_sqrt(self):
        # t^2-t+1/8=0 has two irrational roots in (0,1).
        self.assertEqual(_crossings_of_quadratic((Q(1),-Q(1),Q(1,8)),0,1),2)

    def test_exact_random_cutoff_crosscheck_by_sampling_signs(self):
        # Dense exact rational sample comparison for generic monotone cases.
        # Sampling is merely a regression, NOT the mathematical proof.
        rng=random.Random(221)
        for _ in range(80):
            b=Q(rng.randint(-20,20),10)
            c=Q(rng.randint(-20,20),10)
            A=Q(rng.randint(1,10),10)
            f=(A,b,c)
            if polynomial_value(f,0)==0 or polynomial_value(f,1)==0:continue
            vertex=-b/(2*A)
            if 0<vertex<1 and polynomial_value(f,vertex)==0:continue
            count=_crossings_of_quadratic(f,0,1)
            grid=[polynomial_value(f,Q(j,1000)) for j in range(1001)]
            # An exact sample can land on a rational root, so count changes
            # from sign transitions across zero samples too.
            signs=[(z>0)-(z<0) for z in grid]
            nonzero=[s for s in signs if s]
            changes=sum(x!=y for x,y in zip(nonzero,nonzero[1:]))
            self.assertEqual(changes,count)

    def test_time_reversal_exchanges_source_exit_and_entry(self):
        f=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0))
        segs=list(f.segments)
        segs[1]=MovingSegment(segs[1].start,segs[1].end,pt(0,1))
        forward=linear_square_branch(tuple(segs))
        original=incidence_ledger(forward,0,Q(1,5),Q(1,1024))
        reversed_segments=tuple(MovingSegment(add(s.start,s.velocity),
                                               add(s.end,s.velocity),
                                               scale(-1,s.velocity)) for s in segs)
        reverse=linear_square_branch(reversed_segments)
        backwards=incidence_ledger(reverse,Q(4,5),1,Q(1,1024))
        self.assertEqual((original.start_count,original.end_count),(1,0))
        self.assertEqual((backwards.start_count,backwards.end_count),(0,1))
        self.assertEqual(backwards.source_wall_times,
                         (1-original.source_wall_times[0],))
        self.assertEqual(backwards.cutoff_crossings,0)
        self.assertTrue(backwards.balances)

    def test_rational_rotation_and_translation_preserve_filter_ledger(self):
        f=planted_square_family(pt(Q(1,8),0),pt(0,Q(1,8)))
        base=incidence_ledger(f,0,1,Q(3,128))
        def turn(p):return add(pt(-p.y,p.x),pt(Q(113,7),-Q(53,11)))
        turned=tuple(MovingSegment(turn(s.start),turn(s.end),
                                   pt(-s.velocity.y,s.velocity.x))
                     for s in f.segments)
        rotated=linear_square_branch(turned)
        changed=incidence_ledger(rotated,0,1,Q(3,128))
        self.assertEqual(base,changed)
        self.assertEqual(rotated.side_squared(Q(1,3)),f.side_squared(Q(1,3)))

    def test_tagged_notches_cutoff_counts_exact(self):
        for n in list(range(1,30))+[50,80,110]:
            squared=Q(1,16**(n+1))
            self.assertEqual(shrinking_notch_tagged_count(squared),n)
            prev=Q(1,16**n)
            self.assertEqual(shrinking_notch_tagged_count((squared+prev)/2),n-1)
            self.assertEqual((shrinking_notch_tagged_count(squared)-
                              shrinking_notch_tagged_count((squared+prev)/2))%2,1)

    def test_tagged_notches_not_all_squares_on_curve(self):
        # n=0 is not counted even though the ambient rectangle is itself
        # known to have macroscopic inscribed squares: named counts are lower
        # bounds for marked configurations, not an exhaustive global parity.
        self.assertEqual(shrinking_notch_tagged_count(Q(1)),0)
        self.assertEqual(shrinking_notch_tagged_count(Q(1,256)),1)

    def test_full_repo_v220_notch_square_crosscheck(self):
        try:
            from accumulating_notches_v220 import notch,side_squared,square_residual
        except ImportError:
            self.skipTest('v220 available in full repository checkout')
        for n in range(1,25):
            p=notch(n)
            self.assertEqual(side_squared(p),Q(1,16**(n+1)))
            self.assertEqual(square_residual(p),(0,0,0,0))

    def test_full_repo_v216_affine_square_crosscheck(self):
        try:
            from boundary_degree_v216 import pt as vpt,residual as vres
        except ImportError:
            self.skipTest('v216 available in full repository checkout')
        f=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,Q(1,8)))
        for t in (Q(0),Q(1,3),Q(1)):
            pts=tuple(vpt(p.x,p.y) for p in f.square(t))
            self.assertEqual(vres(pts,1),(0,0,0,0))

    def test_report(self):
        r=report()
        self.assertEqual(r['scaling_cutoff_crossings'],1)
        self.assertEqual(r['source_exit_times'],('55/19968',))
        self.assertEqual(r['named_notches_cutoff_at_n5'],5)
        self.assertTrue(r['scaling_balances'] and r['source_balances'])


if __name__=='__main__':unittest.main()
