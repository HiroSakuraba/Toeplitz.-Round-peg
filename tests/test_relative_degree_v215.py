"""Exact-rational regression tests for v215's boundary-certified PL degree."""
import random
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from relative_degree_v215 import (
    Point, Arc, arc, point, arcs_from_periodic_windows, Segment, _cell, _solve, _det, _rank,
    _feasible_box, _residual, source_lift_degree,
    BoundarySquare, CertificateInconclusive,
)


def seed_arcs():
    a=point(Q(1,2),0);e=point(Q(1,8),Q(1,8));r=point(-e.y,e.x)
    p=(a,point(a.x+e.x,a.y+e.y),
       point(a.x+e.x+r.x,a.y+e.y+r.y),point(a.x+r.x,a.y+r.y))
    slopes=(Q(-2),Q(3),Q(7,2),Q(-1,3));w=Q(1,32)
    arcs=tuple(arc((v.x-w,v.y-s*w),(v.x+w,v.y+s*w))
               for v,s in zip(p,slopes))
    return arcs,p,slopes,w


def change_arc(arcs,index,vertices):
    v=list(arcs);v[index]=arc(*vertices);return tuple(v)


def shift_y(arcs,index,offset):
    v=list(arcs)
    v[index]=arc(*((p.x,p.y+offset) for p in arcs[index].vertices))
    return tuple(v)


class DegreeV215Tests(unittest.TestCase):
    def test_periodic_graph_window_adapter(self):
        # Duck-typed PeriodicGraph fixture with exact graph interpolation.
        a,_,_,_=seed_arcs()
        class Graph:
            def __init__(self, vertices):
                self.vertices=vertices
            def height(self,x):
                for v,w in zip(self.vertices,self.vertices[1:]):
                    if v.x<=x<=w.x:
                        return v.y+(w.y-v.y)*(x-v.x)/(w.x-v.x)
                raise ValueError("outside graph")
        graphs=tuple(Graph((point(Q(0),c.vertices[0].y),)
                                +c.vertices
                                +(point(Q(1),c.vertices[-1].y),)) for c in a)
        windows=tuple((c.vertices[0].x,c.vertices[-1].x) for c in a)
        self.assertEqual(arcs_from_periodic_windows(graphs,windows),a)
        self.assertEqual(source_lift_degree(arcs_from_periodic_windows(graphs,windows)).degree,-1)

    def test_seed_transverse_square(self):
        arcs,p,s,w=seed_arcs()
        c=source_lift_degree(arcs)
        self.assertEqual(c.degree,-1)
        self.assertEqual(len(c.preimages_of_regular_target),1)
        self.assertEqual(c.opposite_vertical_gap,Q(5,64))
        self.assertGreater(c.boundary_margin,0)
        self.assertGreater(c.vertical_perturbation_bound,0)
        self.assertGreater(c.checked_boundary_cells,0)
        self.assertEqual(_residual(p,1),(0,0,0,0))

    def test_singular_but_empty_cells_still_have_margin(self):
        a=tuple(arc((Q(7,16),h),(Q(9,16),h)) for h in (1,0,-1,0))
        cert=source_lift_degree(a)
        self.assertEqual(cert.degree,0)
        self.assertGreater(cert.boundary_margin,0)
        self.assertGreater(cert.vertical_perturbation_bound,0)

    def test_square_at_internal_knot_is_counted(self):
        arcs,p,slopes,w=seed_arcs()
        center=p[2]; slope=slopes[2]; eta=Q(1,100)
        a=change_arc(arcs,2,((center.x-w,center.y-(slope-eta)*w),
                             (center.x,center.y),
                             (center.x+w,center.y+(slope+eta)*w)))
        cert=source_lift_degree(a)
        self.assertEqual(cert.degree,-1)
        self.assertEqual(len(cert.preimages_of_regular_target),1)
        self.assertGreater(cert.boundary_margin,0)
        self.assertEqual(_residual(p,1),(0,0,0,0))

    def test_two_squares_can_cancel(self):
        arcs,p,s,w=seed_arcs()
        l,c,r=p[2].x-w,p[2].x,p[2].x+w
        delta=Q(1,100)
        target=lambda x:Q(1,4)-2*(x-Q(1,2))
        modified=change_arc(arcs,2,((l,target(l)-delta),
                                     (c,target(c)+delta),
                                     (r,target(r)-delta)))
        cert=source_lift_degree(modified)
        self.assertEqual(cert.degree,0)
        self.assertEqual(len(cert.preimages_of_regular_target),2)
        self.assertEqual({1 if z.determinant>0 else -1
                          for z in cert.preimages_of_regular_target},{-1,1})

    def test_boundary_root_detected_even_at_knot(self):
        a=tuple(arc((Q(0),h),(Q(1,2),h),(Q(1),h))
                for h in (1,0,-1,0))
        with self.assertRaises(BoundarySquare):
            source_lift_degree(a,shifts=(0,-1,0,1))

    def test_lift_specificity(self):
        # The projected windows are all identical; the correct square
        # uses deck offsets not represented in the source's physical lift.
        a=tuple(arc((Q(7,16),h),(Q(9,16),h)) for h in (1,0,-1,0))
        empty=source_lift_degree(a,shifts=(0,0,0,0))
        self.assertEqual(empty.degree,0)
        # A full square on the shifted windows may have a whole family,
        # including boundary zeros; this is a correctly rejected degree.
        with self.assertRaises(BoundarySquare):
            source_lift_degree(a,shifts=(0,-1,0,1))

    def test_negative_orientation_is_separate(self):
        a,_,_,_=seed_arcs()
        self.assertEqual(source_lift_degree(a,orientation=-1).degree,0)
        self.assertEqual(source_lift_degree(a,orientation=1).degree,-1)

    def test_exact_linear_algebra(self):
        A=((Q(1),Q(2),Q(0),Q(0)),(Q(2),Q(4),Q(0),Q(0)),
           (Q(0),Q(0),Q(1),Q(0)),(Q(0),Q(0),Q(0),Q(1)))
        self.assertEqual(_rank(A),3)
        self.assertEqual(_det(A),0)
        self.assertTrue(_feasible_box(A,(Q(0),)*4))
        self.assertFalse(_feasible_box(A,(Q(0),Q(1),Q(0),Q(0))))

    def test_strict_robustness_under_random_vertical_shifts(self):
        a,_,_,_=seed_arcs()
        cert=source_lift_degree(a)
        epsilon=cert.vertical_perturbation_bound/Q(3)
        rng=random.Random(215)
        for _ in range(20):
            cur=a
            for j in range(4):
                cur=shift_y(cur,j,Q(rng.randrange(-3,4),3)*epsilon)
            other=source_lift_degree(cur)
            self.assertEqual(other.degree,cert.degree)

    def test_degree_is_invariant_under_collinear_refinement(self):
        a,p,slopes,w=seed_arcs()
        c=p[2];s=slopes[2]
        refined=change_arc(a,2,((c.x-w,c.y-s*w),
                                (c.x,c.y),
                                (c.x+w,c.y+s*w)))
        one=source_lift_degree(a)
        two=source_lift_degree(refined)
        self.assertEqual(one.degree,two.degree)
        self.assertEqual(one.degree,-1)
        # The extra internal face contains a true zero at target 0,
        # so v214's sign-at-0 approach would reject this subdivision.
        self.assertEqual(len(two.preimages_of_regular_target),1)

    def test_exact_boundary_wall_crossing(self):
        # At t=55/624 the planted transverse square reaches the
        # endpoint of source arc 1. Its degree is -1 before and 0 after.
        a,_,_,_=seed_arcs()
        wall=Q(55,624)
        before=source_lift_degree(shift_y(a,1,wall-Q(1,10000)))
        self.assertEqual(before.degree,-1)
        with self.assertRaises(BoundarySquare):
            source_lift_degree(shift_y(a,1,wall))
        after=source_lift_degree(shift_y(a,1,wall+Q(1,10000)))
        self.assertEqual(after.degree,0)

    def test_degree_invariant_when_internal_knot_moves(self):
        a,p,slopes,w=seed_arcs()
        c=p[2];s=slopes[2];eta=Q(1,100)
        # A shifted piecewise kink crosses the actual square's parameter.
        for offset in (-Q(1,10000),Q(0),Q(1,10000)):
            px=c.x+offset
            y=c.y+s*offset
            mod=change_arc(a,2,((c.x-w,c.y-s*w),
                                (px,y),
                                (c.x+w,c.y+s*w)))
            cert=source_lift_degree(mod)
            self.assertEqual(cert.degree,-1)
            self.assertGreater(cert.boundary_margin,0)

    def test_invalid_input(self):
        a,_,_,_=seed_arcs()
        with self.assertRaises(ValueError):source_lift_degree(a,shifts=(1,0,0,0))
        with self.assertRaises(ValueError):source_lift_degree(a,orientation=0)
        with self.assertRaises(ValueError):source_lift_degree(a[:3])
        with self.assertRaises(ValueError):arc((Q(1,2),0),(Q(1,2),1))

    def test_degeneracy_separation_guard(self):
        a,_,_,_=seed_arcs()
        bad=shift_y(a,2,-Q(1,4))
        with self.assertRaises(ValueError):source_lift_degree(bad)

    def test_resource_caps_are_explicit(self):
        a,_,_,_=seed_arcs()
        with self.assertRaises(CertificateInconclusive):
            source_lift_degree(a,max_subdivision_nodes=1)
        with self.assertRaises(CertificateInconclusive):
            source_lift_degree(a,max_target_attempts=1)


if __name__=='__main__':
    unittest.main()
