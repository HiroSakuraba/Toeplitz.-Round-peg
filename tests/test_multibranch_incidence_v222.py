"""v222 exact multi-chart, internal-knot, cutoff and orbit tests."""
import unittest
import random
import sys
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from cutoff_incidence_v221 import (
    pt, add, sub, scale, MovingSegment, planted_square_family,
    linear_square_branch, incidence_ledger, NonGenericIncidence,
)
from multibranch_incidence_v222 import (
    Point, MovingArc, moving_arc, moving_segments, _rank, _box_feasible,
    _singular_chart_can_hold_square, _active_closed_time,
    multibranch_ledger, polygon_example_two_roots,
    canonical_square_orbit, dihedral_labels, fold_model_root_count, report,
)


def from_segments(segments):
    return tuple(MovingArc((s.start,s.end),s.velocity) for s in segments)


class MultiBranchIncidenceTests(unittest.TestCase):
    def test_tent_two_squares_disappear_through_knot(self):
        d=multibranch_ledger(polygon_example_two_roots(),0,1,Q(1,1024))
        self.assertEqual((d.start_count,d.end_count),(2,0))
        self.assertEqual((len(d.internal_events),len(d.source_events),d.cutoff_crossings),(1,0,0))
        self.assertEqual(d.internal_events[0].time,Q(1,2))
        self.assertEqual(set(d.internal_events[0].chart_indices),
                         {(0,0,0,0),(0,0,1,0)})
        self.assertEqual(len(d.chart_branches),2)
        self.assertTrue(d.balances)

    def test_tent_reversal_two_squares_appear(self):
        arcs=polygon_example_two_roots()
        rev=tuple(MovingArc(tuple(add(v,arc.velocity) for v in arc.vertices),
                            scale(-1,arc.velocity)) for arc in arcs)
        d=multibranch_ledger(rev,0,1,Q(1,1024))
        self.assertEqual((d.start_count,d.end_count),(0,2))
        self.assertEqual(len(d.internal_events),1)
        self.assertEqual(d.internal_events[0].time,Q(1,2))
        self.assertTrue(d.balances)

    def test_tent_mixed_size_and_knot_boundary(self):
        arcs=polygon_example_two_roots()
        d=multibranch_ledger(arcs,0,1,Q(1,40))
        self.assertEqual((d.start_count,d.end_count),(1,0))
        self.assertEqual(d.cutoff_crossings,1)
        self.assertEqual(len(d.internal_events),1)
        self.assertEqual(len(d.source_events),0)
        self.assertTrue(d.balances)

    def test_tent_size_threshold_at_knot_rejected(self):
        with self.assertRaisesRegex(NonGenericIncidence,'cutoff coincides'):
            multibranch_ledger(polygon_example_two_roots(),0,1,Q(1,32))

    def test_internal_knot_transit_not_birth(self):
        base=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0),half_width=Q(1,32))
        segs=list(base.segments)
        segs[1]=MovingSegment(segs[1].start,segs[1].end,pt(0,Q(1,20)))
        full=linear_square_branch(segs)
        knot=full.parameters(Q(1,2))[1]
        midpoint=add(segs[1].start,scale(knot,sub(segs[1].end,segs[1].start)))
        arcs=list(from_segments(segs))
        arcs[1]=MovingArc((segs[1].start,midpoint,segs[1].end),segs[1].velocity)
        d=multibranch_ledger(tuple(arcs),0,1,Q(1,1024))
        self.assertEqual((d.start_count,d.end_count),(1,1))
        self.assertEqual((len(d.internal_events),len(d.source_events),d.cutoff_crossings),(1,0,0))
        self.assertEqual(d.internal_events[0].time,Q(1,2))
        self.assertTrue(d.balances)

    def test_redundant_collinear_knot_on_stationary_square(self):
        base=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0),half_width=Q(1,32))
        arc=list(from_segments(base.segments))
        s=base.segments[1]
        mid=add(s.start,scale(Q(1,2),sub(s.end,s.start)))
        arc[1]=MovingArc((s.start,mid,s.end),s.velocity)
        with self.assertRaisesRegex(NonGenericIncidence,'persistent internal/source-face'):
            multibranch_ledger(arc,0,1,Q(1,1024))

    def test_one_branch_v221_cutoff_crosscheck(self):
        f=planted_square_family(pt(0,0),pt(Q(1,8),Q(1,8)))
        d=multibranch_ledger(from_segments(f.segments),Q(1,4),Q(3,4),Q(1,128))
        old=incidence_ledger(f,Q(1,4),Q(3,4),Q(1,128))
        self.assertEqual((d.start_count,d.end_count), (old.start_count,old.end_count))
        self.assertEqual(d.cutoff_crossings,old.cutoff_crossings)
        self.assertEqual(len(d.source_events),old.source_events)
        self.assertTrue(d.balances)

    def test_one_branch_v221_irrational_cutoff(self):
        f=planted_square_family(pt(Q(1,8),0),pt(0,Q(1,8)))
        d=multibranch_ledger(from_segments(f.segments),0,1,Q(3,128))
        self.assertEqual((d.start_count,d.end_count,d.cutoff_crossings),(0,1,1))
        self.assertTrue(d.balances)

    def test_one_branch_v221_source_exit(self):
        f=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0))
        x=list(f.segments)
        x[1]=MovingSegment(x[1].start,x[1].end,pt(0,1))
        prev=incidence_ledger(linear_square_branch(x),0,Q(1,5),Q(1,1024))
        d=multibranch_ledger(from_segments(x),0,Q(1,5),Q(1,1024))
        self.assertEqual((d.start_count,d.end_count),(1,0))
        self.assertEqual(len(d.source_events),1)
        self.assertEqual(d.source_events[0].time,Q(55,19968))
        self.assertEqual(d.source_events[0].time,prev.source_wall_times[0])
        self.assertEqual(d.cutoff_crossings,0)
        self.assertTrue(d.balances)

    def test_one_branch_two_cutoff_crossings(self):
        f=planted_square_family(pt(Q(1,4),0),pt(-Q(1,2),0))
        d=multibranch_ledger(from_segments(f.segments),Q(1,4),Q(3,4),Q(1,256))
        self.assertEqual((d.start_count,d.end_count,d.cutoff_crossings),(1,1,2))
        self.assertTrue(d.balances)

    def test_time_reversal_of_source_exit(self):
        f=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0))
        x=list(f.segments);x[1]=MovingSegment(x[1].start,x[1].end,pt(0,1))
        rv=tuple(MovingSegment(add(s.start,s.velocity),add(s.end,s.velocity),
                               scale(-1,s.velocity)) for s in x)
        d=multibranch_ledger(from_segments(rv),Q(4,5),1,Q(1,1024))
        self.assertEqual((d.start_count,d.end_count),(0,1))
        self.assertEqual(d.source_events[0].time,1-Q(55,19968))
        self.assertTrue(d.balances)

    def test_dihedral_orbit_eight_distinct_labelings_one_geometry(self):
        p=(pt(0,0),pt(2,0),pt(2,2),pt(0,2))
        labels=dihedral_labels(p)
        self.assertEqual(len(labels),8)
        self.assertEqual(len(set(labels)),8)
        self.assertEqual(len({canonical_square_orbit(r) for r in labels}),1)
        from cutoff_incidence_v221 import residual
        self.assertEqual(sum(residual(x,+1)==(0,0,0,0) for x in labels),4)
        self.assertEqual(sum(residual(x,-1)==(0,0,0,0) for x in labels),4)
        self.assertEqual(len(labels)%2,0)

    def test_geometric_orbit_rejects_nonsquare_and_collapse(self):
        with self.assertRaises(ValueError):canonical_square_orbit((pt(0,0),)*4)
        with self.assertRaises(ValueError):canonical_square_orbit((pt(0,0),pt(1,0),pt(2,0),pt(3,0)))
        with self.assertRaises(ValueError):dihedral_labels((pt(0,0),)*4)

    def test_abstract_fold_birth_death_even(self):
        self.assertEqual((fold_model_root_count(-Q(1,7)),fold_model_root_count(0),fold_model_root_count(Q(1,7))),(0,1,2))
        self.assertEqual((fold_model_root_count(-1)+fold_model_root_count(1))%2,0)

    def test_singular_empty_chart_skipped_exactly(self):
        a=(moving_arc(((0,1),(1,1))),moving_arc(((0,0),(1,0))),
           moving_arc(((0,-1),(1,-1))),moving_arc(((0,0),(1,0))))
        d=multibranch_ledger(a,0,1,Q(1,10))
        self.assertEqual((d.start_count,d.end_count,len(d.chart_branches)),(0,0,0))
        self.assertTrue(d.balances)

    def test_singular_feasible_chart_fails_closed(self):
        p=(moving_arc(((0,0),(1,0))),)*4
        with self.assertRaisesRegex(NonGenericIncidence,'singular'):
            multibranch_ledger(p,0,1,Q(1,100))

    def test_exact_feasibility_rank_deficient_inconsistent(self):
        self.assertEqual(_rank(((Q(1),Q(0)),(Q(2),Q(0)))),1)
        self.assertTrue(_box_feasible(((Q(1),Q(0)),(Q(2),Q(0))), (Q(1,2),Q(1)),2))
        self.assertFalse(_box_feasible(((Q(1),Q(0)),(Q(2),Q(0))), (Q(1,2),Q(2)),2))

    def test_invalid_inputs(self):
        f=polygon_example_two_roots()
        for args in ((0,0,Q(1,1024)),(0,1,0)):
            with self.assertRaises(ValueError):multibranch_ledger(f,*args)
        with self.assertRaises(ValueError):multibranch_ledger(f,0,1,Q(1,1024),orientation=0)
        with self.assertRaises(ValueError):multibranch_ledger(f,0,1,Q(1,1024),shifts=(1,0,0,0))
        with self.assertRaises(ValueError):moving_arc(((0,0),(0,0)))

    def test_rotated_tent_same_multibranch_ledger(self):
        arcs=polygon_example_two_roots()
        turn=lambda p:pt(-p.y+Q(11,5),p.x-Q(7,9))
        modified=tuple(MovingArc(tuple(turn(v) for v in a.vertices),
                                 pt(-a.velocity.y,a.velocity.x)) for a in arcs)
        x=multibranch_ledger(modified,0,1,Q(1,1024))
        self.assertEqual((x.start_count,x.end_count,len(x.internal_events)),(2,0,1))
        self.assertTrue(x.balances)

    def test_seeded_independent_one_chart_crosschecks(self):
        rng=random.Random(222)
        for _ in range(36):
            x=pt(Q(rng.randint(-4,4),128),Q(rng.randint(-4,4),128))
            e=pt(Q(rng.randint(-4,4),64),Q(rng.randint(-4,4),64))
            f=planted_square_family(pt(Q(1,8),Q(1,8)),e)
            old=incidence_ledger(f,0,1,Q(1,256))
            d=multibranch_ledger(from_segments(f.segments),0,1,Q(1,256))
            self.assertEqual((old.start_count,old.end_count,old.source_events,old.cutoff_crossings),
                             (d.start_count,d.end_count,len(d.source_events),d.cutoff_crossings))
            self.assertTrue(d.balances)

    def test_fixed_physical_lifts_not_quotient_only(self):
        from cutoff_incidence_v221 import rot
        p0=pt(Q(1,2),0);e=pt(-1,-1);j=rot(e,1)
        pts=(p0,add(p0,e),add(add(p0,e),j),add(p0,j))
        shifts=(0,-1,0,1)
        slopes=(-Q(2),Q(3),Q(7,2),-Q(1,3))
        half=Q(1,32)
        arcs=tuple(moving_arc(((v.x-k-half,v.y-m*half),
                               (v.x-k+half,v.y+m*half)))
                   for v,k,m in zip(pts,shifts,slopes))
        correct=multibranch_ledger(arcs,0,1,Q(1,2),shifts=shifts)
        wrong=multibranch_ledger(arcs,0,1,Q(1,2))
        self.assertEqual((correct.start_count,correct.end_count),(1,1))
        self.assertEqual((wrong.start_count,wrong.end_count),(0,0))
        self.assertTrue(correct.balances and wrong.balances)

    def test_reflected_tent_switches_square_orientation(self):
        a=polygon_example_two_roots()
        reflected=tuple(MovingArc(tuple(pt(v.x,-v.y) for v in arc.vertices),
                                  pt(arc.velocity.x,-arc.velocity.y)) for arc in a)
        good=multibranch_ledger(reflected,0,1,Q(1,1024),orientation=-1)
        opposite=multibranch_ledger(reflected,0,1,Q(1,1024),orientation=1)
        self.assertEqual((good.start_count,good.end_count,len(good.internal_events)),(2,0,1))
        self.assertEqual((opposite.start_count,opposite.end_count),(0,0))

    def test_optional_v220_tagged_notches_are_not_orbit_count(self):
        try:
            from accumulating_notches_v220 import notch,side_squared
        except ImportError:
            self.skipTest('v220 only in full repository checkout')
        from cutoff_incidence_v221 import shrinking_notch_tagged_count
        for n in range(1,9):
            self.assertEqual(side_squared(notch(n)),Q(1,16**(n+1)))
            self.assertEqual(shrinking_notch_tagged_count(Q(1,16**(n+1))),n)

    def test_report(self):
        r=report()
        self.assertEqual(r['two_root_model']['start'],2)
        self.assertEqual(r['two_root_model']['end'],0)
        self.assertEqual(r['two_root_model']['internal_pairs'],1)
        self.assertEqual(r['dihedral_orbit_labels'],8)
        self.assertTrue(r['two_root_model']['balances'])

    def test_optional_v218_polygon_orbit_crosscheck(self):
        try:
            from configuration_symmetry_v218 import isolated_distinct_edge_squares, pt as p218
        except ImportError:
            self.skipTest('v218 only in full repository checkout')
        P=tuple(p218(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        orbits=isolated_distinct_edge_squares(P)
        self.assertEqual(len(orbits),1)
        self.assertEqual(len(next(iter(orbits.values()))),8)

if __name__=='__main__':unittest.main()
