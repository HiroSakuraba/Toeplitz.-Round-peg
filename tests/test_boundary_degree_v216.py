"""Rational v216 regression suite, including global source-window selection limits."""
import random
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from boundary_degree_v216 import (
    pt,arc,Arc,add,sub,mul,Point,boundary_degree,
    BoundaryRoot,DegreeInconclusive,select_certified_four_arcs,
    vertex_neighborhood,segment_distance_squared,arc_distance_squared,
    feasible_box,_solve,_det,_rank,_affine_pieces,validate_simple_polygon,
)


def planted_arcs(slopes=(Q(-2),Q(3),Q(7,2),Q(-1,3)),hw=Q(1,32),dy=(Q(0),)*4,
                 edge=None,shifts=(0,0,0,0)):
    origin=pt(Q(1,2),0)
    e=edge or pt(Q(1,8),Q(1,8))
    j=pt(-e.y,e.x)
    p=(origin,add(origin,e),add(add(origin,e),j),add(origin,j))
    result=[]
    for k,v,s,y in zip(shifts,p,slopes,dy):
        x=v.x-k
        result.append(arc((x-hw,v.y-s*hw+y),(x+hw,v.y+s*hw+y)))
    return tuple(result)


def shift_arc(a,d):
    return Arc(tuple(pt(p.x,p.y+d) for p in a.vertices))


class BoundaryDegreeTests(unittest.TestCase):
    def test_affine_linear_algebra_and_rank_deficient_box(self):
        self.assertEqual(_det(((1,2),(3,4))),Q(-2))
        self.assertEqual(_rank(((1,2),(2,4))),1)
        self.assertTrue(feasible_box(((1,1),(2,2)),(Q(1),Q(2)),2))
        self.assertFalse(feasible_box(((1,1),(2,2)),(Q(1),Q(3)),2))
        self.assertEqual(_solve(((1,2),(3,4)),(Q(5),Q(11))),(Q(1),Q(2)))

    def test_nonmonotone_arcs_supported(self):
        a=planted_arcs()
        # Move a point backwards in x outside the chosen root, while keeping
        # the source endpoint geometry and a fixed original source arc.
        a=list(a)
        p,q=a[0].vertices
        extra=pt(p.x-Q(1,128),p.y-Q(1,1000))
        a[0]=Arc((p,extra,q))
        c=boundary_degree(a)
        self.assertEqual(c.degree,-1)
        self.assertGreater(c.boundary_margin,0)

    def test_base_degree_without_interior_root_enumeration(self):
        c=boundary_degree(planted_arcs())
        self.assertEqual(c.degree,-1)
        self.assertGreaterEqual(c.intersections,1)
        self.assertGreater(c.boundary_margin,0)
        self.assertGreater(c.opposite_distance_sq,0)
        self.assertGreater(c.perturbation_bound_linf,0)

    def test_wrong_orientation_is_zero(self):
        self.assertEqual(boundary_degree(planted_arcs(),orientation=-1).degree,0)

    def test_shifts_do_not_masquerade_as_physical_lifts(self):
        a=planted_arcs(edge=pt(-1,-1),shifts=(0,-1,0,1))
        self.assertNotEqual(boundary_degree(a,shifts=(0,-1,0,1)).degree,0)
        self.assertEqual(boundary_degree(a,shifts=(0,0,0,0)).degree,0)

    def test_internal_exact_collinear_knot_remains_same(self):
        a=list(planted_arcs())
        p,q=a[0].vertices
        a[0]=Arc((p,pt((p.x+q.x)/2,(p.y+q.y)/2),q))
        self.assertEqual(boundary_degree(a).degree,-1)

    def test_internal_genuine_kink_through_square(self):
        a=list(planted_arcs())
        p,q=a[2].vertices
        mid=pt((p.x+q.x)/2,(p.y+q.y)/2)
        # Make a genuine kink preserving the root and both source endpoints.
        a[2]=Arc((p,pt((p.x+mid.x)/2,(p.y+mid.y)/2+Q(1,256)),mid,q))
        self.assertEqual(boundary_degree(a).degree,-1)

    def test_negative_degree_cancellation_two_squares(self):
        a=list(planted_arcs())
        lo,center,hi=Q(15,32),Q(1,2),Q(17,32)
        delta=Q(1,100)
        target=lambda x:Q(1,4)-2*(x-Q(1,2))
        a[2]=arc((lo,target(lo)-delta),(center,target(center)+delta),
                 (hi,target(hi)-delta))
        self.assertEqual(boundary_degree(a).degree,0)
        lowered=list(a);lowered[2]=shift_arc(lowered[2],-2*delta)
        self.assertEqual(boundary_degree(lowered).degree,0)

    def test_boundary_escape_reports_explicit_failure(self):
        a=planted_arcs()
        t=Q(55,624)
        before=list(a);before[1]=shift_arc(before[1],t-Q(1,10000))
        boundary=list(a);boundary[1]=shift_arc(boundary[1],t)
        after=list(a);after[1]=shift_arc(after[1],t+Q(1,10000))
        self.assertEqual(boundary_degree(before).degree,-1)
        with self.assertRaises(BoundaryRoot):boundary_degree(boundary)
        self.assertEqual(boundary_degree(after).degree,0)

    def test_perturbation_bound_rational_family(self):
        rng=random.Random(216)
        a=planted_arcs()
        c=boundary_degree(a)
        eps=c.perturbation_bound_linf
        for _ in range(10):
            vals=[Q(rng.randint(-3,3),4)*eps for i in range(4)]
            shifted=tuple(shift_arc(ar,t) for ar,t in zip(a,vals))
            self.assertEqual(boundary_degree(shifted).degree,c.degree)

    def test_segment_distances_exact_and_zero_intersection(self):
        self.assertEqual(segment_distance_squared(pt(0,0),pt(1,0),pt(2,-1),pt(2,1)),Q(1))
        self.assertEqual(segment_distance_squared(pt(0,0),pt(2,2),pt(0,2),pt(2,0)),Q(0))
        self.assertEqual(segment_distance_squared(pt(0,0),pt(1,0),pt(0,1),pt(1,1)),Q(1))

    def test_degree_ignores_opposite_height_projection_requirement(self):
        # Choose a rotated seed where opposite arc y ranges overlap while
        # their Euclidean distance is strictly positive.
        a=planted_arcs(edge=pt(Q(1,8),-Q(1,8)))
        # With this seed, opposite y intervals need not be disjoint,
        # whereas Euclidean separation is checked without projection.
        c=boundary_degree(a,orientation=1)
        self.assertNotEqual(c.degree,0)

    def test_rational_result_is_deterministic(self):
        a=planted_arcs()
        c=boundary_degree(a)
        self.assertEqual(c,boundary_degree(a))
        self.assertEqual(c.ray,(Q(1),Q(c.ray[1]),Q(c.ray[2]),Q(c.ray[3])))

    def test_vertex_arc_selection_is_a_priori(self):
        poly=(pt(0,0),pt(2,0),pt(2,2),pt(0,2))
        arcs=tuple(vertex_neighborhood(poly,i) for i in range(4))
        self.assertEqual(len(arcs),4)
        self.assertEqual(arcs[0].vertices[1],pt(0,0))
        self.assertGreater(arc_distance_squared(arcs[0],arcs[2]),0)

    def test_global_search_returns_only_conditional_results(self):
        poly=(pt(0,0),pt(2,0),pt(2,2),pt(0,2))
        result,attempts,err=select_certified_four_arcs(poly,max_candidates=2)
        self.assertEqual(attempts,2)
        # It can fail even on curves known to inscribe squares: negative
        # output is not a theorem of square nonexistence.
        self.assertTrue(result is None or result[2].degree!=0)

    def test_boundary_count_against_independent_interior_linear_solve(self):
        # This test computes the only-cell affine root and sign directly;
        # boundary_degree() itself never enumerates interior square roots.
        base=planted_arcs()
        rng=random.Random(216216)
        seen={}
        for _ in range(72):
            deltas=tuple(Q(rng.randint(-35,35),512) for _ in range(4))
            arcs=tuple(shift_arc(a,t) for a,t in zip(base,deltas))
            st=tuple((a.vertices[0],a.vertices[1]) for a in arcs)
            f0,cols=_affine_pieces(st,(0,0,0,0),1)
            A=[[cols[j][i] for j in range(4)] for i in range(4)]
            det=_det(A)
            self.assertNotEqual(det,0)
            solution=_solve(A,[-f for f in f0])
            if any(v in (0,1) for v in solution):
                with self.assertRaises(BoundaryRoot):boundary_degree(arcs)
                continue
            expected=(1 if det>0 else -1) if all(0<v<1 for v in solution) else 0
            if arc_distance_squared(arcs[0],arcs[2])==0:
                continue  # v216 requires positive opposite-arc gap
            c=boundary_degree(arcs)
            self.assertEqual(c.degree,expected)
            seen[expected]=seen.get(expected,0)+1
        self.assertIn(-1,seen)
        self.assertIn(0,seen)

    def test_parameter_reversal_changes_degree_sign(self):
        a=list(planted_arcs())
        self.assertEqual(boundary_degree(a).degree,-1)
        a[0]=Arc(tuple(reversed(a[0].vertices)))
        self.assertEqual(boundary_degree(a).degree,1)

    def test_rational_rotation_and_translation_of_marked_arcs(self):
        a=planted_arcs()
        def turn(p):
            return pt(-p.y+Q(9,7),p.x-Q(2,5))
        b=tuple(Arc(tuple(turn(p) for p in ar.vertices)) for ar in a)
        self.assertEqual(boundary_degree(a).degree,boundary_degree(b).degree)

    def test_certified_convex_quadrilateral_neighborhood(self):
        poly=tuple(pt(*x) for x in ((0,0),(3,1),(2,3),(-1,2)))
        self.assertGreater(validate_simple_polygon(poly),0)
        r,attempts,_=select_certified_four_arcs(poly)
        self.assertIsNotNone(r)
        family,indices,orientation,cert=r
        self.assertEqual((family,indices,orientation),('vertex',(0,1,2,3),1))
        self.assertEqual(cert.degree,-1)
        self.assertEqual(cert.boundary_margin,Q(1,8))
        self.assertEqual(cert.perturbation_bound_linf,Q(1,32))
        self.assertEqual(attempts,1)

    def test_certified_quadrilateral_robust_to_two_dimensional_vertex_motion(self):
        poly=tuple(pt(*x) for x in ((0,0),(3,1),(2,3),(-1,2)))
        base=boundary_degree(tuple(vertex_neighborhood(poly,i) for i in range(4)))
        rng=random.Random(2160)
        self.assertEqual(base.perturbation_bound_linf,Q(1,32))
        for _ in range(15):
            # Corresponding points on both half-edges move by at most
            # 1/64 in each planar coordinate under vertex perturbations.
            perturbed=tuple(pt(v.x+Q(rng.randint(-1,1),64),
                               v.y+Q(rng.randint(-1,1),64)) for v in poly)
            validate_simple_polygon(perturbed)
            cert=boundary_degree(tuple(vertex_neighborhood(perturbed,i)
                                       for i in range(4)))
            self.assertEqual(cert.degree,base.degree)

    def test_edge_selection_succeeds_when_vertex_count_cancels(self):
        poly=tuple(pt(*x) for x in ((0,0),(2,0),(3,3),(0,2)))
        r,n,errors=select_certified_four_arcs(poly)
        self.assertIsNotNone(r)
        self.assertEqual((r[0],r[1],r[2],r[3].degree),
                         ('edge',(0,1,2,3),1,-1))
        self.assertEqual(r[3].perturbation_bound_linf,Q(1,64))

    def test_square_polygon_scan_is_inconclusive_not_false(self):
        poly=tuple(pt(*x) for x in ((0,0),(1,0),(1,1),(0,1)))
        r,n,errors=select_certified_four_arcs(poly)
        self.assertIsNone(r)
        self.assertEqual(n,4)
        self.assertTrue(any(e[-1]=='BoundaryRoot' for e in errors))
        # The polygon itself is a square! No certificate != no square.

    def test_exact_simple_polygon_validation_rejects_crossings(self):
        with self.assertRaises(ValueError):validate_simple_polygon(
            tuple(pt(*p) for p in ((0,0),(2,2),(0,2),(2,0))))
        with self.assertRaises(ValueError):validate_simple_polygon(
            tuple(pt(*p) for p in ((0,0),(1,0),(0,0),(0,1))))
        with self.assertRaises(ValueError):select_certified_four_arcs(
            tuple(pt(*p) for p in ((0,0),(2,2),(0,2),(2,0))))

    def test_periodic_shifts_apply_to_size_guard(self):
        # The zero-side guard must use the *lifted* arcs; quotient arcs
        # coinciding cannot be mistaken for physical coincidence.
        a=list(planted_arcs())
        self.assertGreater(boundary_degree(a).opposite_distance_sq,0)
        # Force opposite input arcs to coincide geometrically, but move the
        # physical lift of arc 2 far enough to separate them.
        a[2]=a[0]
        with self.assertRaises(ValueError):boundary_degree(a)
        c=boundary_degree(a,shifts=(0,0,5,0))
        self.assertGreater(c.opposite_distance_sq,0)

    def test_validation(self):
        a=planted_arcs()
        with self.assertRaises(ValueError):boundary_degree(a[:3])
        with self.assertRaises(ValueError):boundary_degree(a,shifts=(1,0,0,0))
        with self.assertRaises(ValueError):boundary_degree(a,orientation=0)
        with self.assertRaises(ValueError):Arc((pt(0,0),pt(0,0)))

if __name__=='__main__':unittest.main()