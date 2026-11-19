"""Exact algebraic-time semialgebraic rank/feasibility regression tests for v228."""
import sys
import unittest
import random
from fractions import Fraction as Q
from pathlib import Path
from itertools import product
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from singular_event_feasibility_v228 import (
    AlgebraicQuotient,AlgebraicFeasibilityInconclusive,algebraic_box_vertices,
    rank_and_basis_rows_at,singular_polygon_chart_at_event,
    audit_singular_event_slices,report,
)
from local_orbit_germs_v227 import AlgebraicTime,sign_at_time
from polygon_orbit_census_v223 import (
    pt,convex_vertex_alignment_polygon,census_polygon,_box_vertices,_rank,
)
from persistent_strata_v225 import (
    affine_system_polynomials,singular_chart_snapshot,
)
from moving_polygon_events_v224 import interpolated_polygon


R=AlgebraicTime((-Q(1,2),1),Q(1,2),Q(1,2))
I=AlgebraicTime((-1,2,1),Q(2,5),Q(1,2))   # sqrt(2)-1


def system(b_first):
    # 2x2 block [[t+1,1],[2,t+1]] drops rank at sqrt(2)-1.
    A=(((1,1),(1,), (0,), (0,)),
       ((2,), (1,1),(0,), (0,)),
       ((0,), (0,), (1,), (0,)),
       ((0,), (0,), (0,), (1,)))
    b=tuple(b_first)+((Q(1,2),),(Q(1,2),))
    return A,b


class ExactAlgebraicBoxTests(unittest.TestCase):
    def test_irrational_rank_loss_is_exact(self):
        A,b=system(((0,),(0,)))
        rank,rows=rank_and_basis_rows_at(A,I)
        self.assertEqual(rank,3)
        self.assertEqual(len(rows),3)
        self.assertEqual(rank_and_basis_rows_at(A,R)[0],4)

    def test_irrational_singular_isolated_vertex(self):
        A,b=system(((0,),(0,)))
        rank,vertices=algebraic_box_vertices(A,b,I)
        self.assertEqual(rank,3)
        self.assertEqual(len(vertices),1)
        for j in (0,1):self.assertEqual(vertices[0][j].sign(),0)
        for j in (2,3):
            self.assertTrue(vertices[0][j].equal(AlgebraicQuotient((Q(1,2),),(1,),I)))

    def test_irrational_singular_positive_dimension(self):
        A,b=system(((Q(1,2),Q(1,2)), (1,)))
        rank,vertices=algebraic_box_vertices(A,b,I)
        self.assertEqual(rank,3)
        self.assertEqual(len(vertices),2)
        self.assertFalse(all(x.equal(y) for x,y in zip(vertices[0],vertices[1])))

    def test_irrational_singular_inconsistent(self):
        A,b=system(((1,),(0,)))
        rank,vertices=algebraic_box_vertices(A,b,I)
        self.assertEqual(rank,3)
        self.assertEqual(vertices,())

    def test_algebraic_quotient_cross_multiplication(self):
        x=AlgebraicQuotient((1,1),(2,),I)
        y=AlgebraicQuotient((2,2),(4,),I)
        z=AlgebraicQuotient((1,),(1,),I)
        self.assertTrue(x.equal(y))
        self.assertFalse(x.equal(z))
        self.assertEqual(x.sign(),1)
        self.assertEqual(AlgebraicQuotient((-1,),(2,),I).sign(),-1)
        with self.assertRaises(ValueError):x.rational_value()

    def test_zero_denominator_rejected(self):
        with self.assertRaises(ZeroDivisionError):
            AlgebraicQuotient((1,),(-1,2,1),I)

    def test_budget_protection(self):
        A,b=system(((0,),(0,)))
        with self.assertRaises(AlgebraicFeasibilityInconclusive):
            algebraic_box_vertices(A,b,I,max_candidates=0)

    def test_invalid_shape_and_event(self):
        with self.assertRaises(ValueError):algebraic_box_vertices(((1,),),((1,),),I)
        with self.assertRaises(TypeError):algebraic_box_vertices(*system(((0,),(0,))),Q(1,2))

    def test_rational_full_rank_exact_values(self):
        A=tuple(tuple(((1,) if i==j else (0,)) for j in range(4)) for i in range(4))
        b=((Q(1,4),),(Q(1,2),),(Q(3,4),),(Q(1,3),))
        rank,vertices=algebraic_box_vertices(A,b,R)
        self.assertEqual(rank,4)
        self.assertEqual(len(vertices),1)
        self.assertEqual(tuple(x.rational_value() for x in vertices[0]),
                         (Q(1,4),Q(1,2),Q(3,4),Q(1,3)))

    def test_independent_random_rational_box_vertices(self):
        # Independent v223 static rational polytope implementation.
        rng=random.Random(228)
        for _ in range(36):
            A0=tuple(tuple(Q(rng.randint(-2,2)) for j in range(4)) for i in range(4))
            b0=tuple(Q(rng.randint(-2,2),2) for i in range(4))
            A=tuple(tuple((x,) for x in row) for row in A0)
            b=tuple((x,) for x in b0)
            rank,vertices=algebraic_box_vertices(A,b,R)
            expected=set(tuple(Q(y) for y in x) for x in _box_vertices(A0,b0))
            actual=set(tuple(z.rational_value() for z in x) for x in vertices)
            self.assertEqual(rank,_rank(A0))
            self.assertEqual(actual,expected)

    def test_rational_full_rank_infeasible(self):
        A=tuple(tuple(((1,) if i==j else (0,)) for j in range(4)) for i in range(4))
        b=((2,),(0,),(0,),(0,))
        rank,vertices=algebraic_box_vertices(A,b,R)
        self.assertEqual(rank,4)
        self.assertEqual(vertices,())

    def test_rational_rank_zero_all_box_corners(self):
        A=tuple(((0,),)*4 for _ in range(4))
        b=((0,),)*4
        rank,vertices=algebraic_box_vertices(A,b,R)
        self.assertEqual(rank,0)
        self.assertEqual(len(vertices),16)
        self.assertEqual(set(tuple(x.rational_value() for x in v) for v in vertices),
                         set(product((Q(0),Q(1)),repeat=4)))

    def test_rank_zero_inconsistent_detected(self):
        A=tuple(((0,),)*4 for _ in range(4))
        b=((1,),)+((0,),)*3
        rank,vertices=algebraic_box_vertices(A,b,R)
        self.assertEqual((rank,len(vertices)),(0,0))

    def test_rational_box_boundary_is_inclusive(self):
        A=tuple(tuple(((1,) if i==j else (0,)) for j in range(4)) for i in range(4))
        b=((0,),(1,),(Q(1,2),),(Q(1,2),))
        _,vertices=algebraic_box_vertices(A,b,R)
        self.assertEqual(len(vertices),1)
        self.assertEqual(tuple(v.rational_value() for v in vertices[0]),
                         (Q(0),Q(1),Q(1,2),Q(1,2)))


class PolygonEventSliceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.R=tuple(pt(*p) for p in ((0,0),(2,0),(2,1),(0,1)))
        cls.T=tuple(pt(*p) for p in ((0,0),(2,0),(2,Q(3,2)),(0,Q(1,2))))

    def test_rectangle_continuum(self):
        s=singular_polygon_chart_at_event(self.R,self.R,(0,0,2,2),-1,R,Q(1,100))
        self.assertEqual(s.rank,3)
        self.assertEqual(s.status,'positive-dimensional-square-family')
        self.assertEqual(len(s.parameter_vertices),2)
        self.assertTrue(s.has_positive_square)
        self.assertIsNotNone(s.witness)
        self.assertEqual(len(s.witness),4)

    def test_rectangle_cutoff_excludes_continuum(self):
        s=singular_polygon_chart_at_event(self.R,self.R,(0,0,2,2),-1,R,1)
        self.assertEqual(s.status,'no-positive-size-square')
        self.assertFalse(s.has_positive_square)

    def test_rectangle_infeasible_chart(self):
        s=singular_polygon_chart_at_event(self.R,self.R,(0,2,0,2),1,R,Q(1,100))
        self.assertEqual(s.status,'infeasible')
        self.assertEqual(len(s.parameter_vertices),0)

    def test_independent_v225_rectangle_snapshot(self):
        for ids,sigma in (((0,0,2,2),-1),((0,2,0,2),1),((0,0,2,2),1)):
            x=singular_polygon_chart_at_event(self.R,self.R,ids,sigma,R,Q(1,100))
            y=singular_chart_snapshot(self.R,ids,sigma,Q(1,100))
            self.assertEqual(x.rank,y.rank)
            self.assertEqual(len(x.parameter_vertices),len(y.vertices))
            if x.status=='positive-dimensional-square-family':
                self.assertEqual(y.status,'positive-dimensional-family')

    def test_all_trapezoid_singular_charts_independent_snapshot(self):
        # Complete check across the rank-deficient edge charts of one
        # rational nonrectangular quadrilateral.
        checked=0
        for ids in product(range(4),repeat=4):
            if any(ids.count(i)>=3 for i in set(ids)):continue
            for sigma in (-1,1):
                x=singular_polygon_chart_at_event(self.T,self.T,ids,sigma,R,Q(1,100))
                if x.rank==4:continue
                y=singular_chart_snapshot(self.T,ids,sigma,Q(1,100))
                checked+=1
                self.assertEqual(len(x.parameter_vertices),len(y.vertices))
                self.assertEqual(x.rank,y.rank)
        self.assertEqual(checked,28)

    def test_independent_v225_trapezoid_singular_snapshot(self):
        for ids,sigma in (((0,1,0,1),-1),((0,1,0,3),1),((1,0,3,0),1)):
            x=singular_polygon_chart_at_event(self.T,self.T,ids,sigma,R,Q(1,100))
            y=singular_chart_snapshot(self.T,ids,sigma,Q(1,100))
            self.assertEqual(x.rank,y.rank)
            self.assertEqual(len(x.parameter_vertices),len(y.vertices))

    def test_convex_event_all_singular_charts_have_no_square_above_cutoff(self):
        P0=convex_vertex_alignment_polygon(0)
        P1=convex_vertex_alignment_polygon(Q(9,8))
        a=AlgebraicTime((-Q(8,9),1),Q(8,9),Q(8,9))
        r=audit_singular_event_slices(P0,P1,a,Q(1,100))
        self.assertEqual(r.checked_singular_charts,28)
        self.assertEqual(dict(r.statuses),{'infeasible':20,'no-positive-size-square':8})
        self.assertEqual(len(r.witnesses),0)
        self.assertEqual(r.unresolved_chart_count,0)
        # The two physical event squares found in v227 are thus represented
        # by nonsingular charts; no hidden positive-size singular square.
        self.assertEqual(census_polygon(interpolated_polygon(P0,P1,Q(8,9)),Q(1,100)).count,2)

    def test_irrational_scaling_event_singular_chart_audit(self):
        P0=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        P1=tuple(pt(2*p.x,2*p.y) for p in P0)
        r=audit_singular_event_slices(P0,P1,I,Q(144,25))
        self.assertEqual(r.checked_singular_charts,4)
        self.assertEqual(dict(r.statuses),{'no-positive-size-square':4})
        self.assertEqual(r.unresolved_chart_count,0)

    def test_cutoff_sign_exactly_at_threshold(self):
        # If cutoff equals the largest side-squared on the rectangle chart,
        # no strictly larger square survives the filter.
        s=singular_polygon_chart_at_event(self.R,self.R,(0,0,2,2),-1,R,Q(1))
        self.assertFalse(s.has_positive_square)

    def test_invalid_polygon_chart_args(self):
        with self.assertRaises(ValueError):
            singular_polygon_chart_at_event(self.R,self.R,(0,2,0,5),1,R,0)
        with self.assertRaises(ValueError):
            singular_polygon_chart_at_event(self.R,self.R,(0,0,2,2),7,R,0)
        with self.assertRaises(ValueError):
            singular_polygon_chart_at_event(self.R,self.R,(0,0,2,2),-1,R,-1)

    def test_audit_budget_fails_closed(self):
        with self.assertRaises(AlgebraicFeasibilityInconclusive):
            audit_singular_event_slices(self.R,self.R,R,Q(1,100),max_charts=0)

    def test_report_exact_consistency(self):
        r=report()
        self.assertEqual(r['rectangle_status'],'positive-dimensional-square-family')
        self.assertEqual(r['rectangle_rank'],3)
        self.assertEqual(r['feasible_vertices'],2)


if __name__=='__main__':unittest.main()
