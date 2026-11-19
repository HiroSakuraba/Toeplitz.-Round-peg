"""Exact algebraic event-cover tests for general linear polygon vertex motion."""
import random,sys,unittest
from pathlib import Path
from fractions import Fraction as Q
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from polygon_orbit_census_v223 import (
    pt,sub,add,mul,cross,census_polygon,convex_vertex_alignment_polygon,
    validate_jordan_polygon,_affine_system,_solve,_evaluate,_side_sq,
)
from moving_polygon_events_v224 import (
    ZERO,ONE,T,polynomial,iszero,plus,minus,neg,times,divide,
    derivative,evaluate,poly_gcd,squarefree,sturm_sequence,
    count_roots_open,isolate_roots,determinant_poly,
    chart_polynomials,isolate_chart_events,polygon_event_cover,
    rational_vertex_square,interpolated_polygon,certify_jordan_isotopy,
    three_point_orientation_polynomial,
)

class PolynomialTests(unittest.TestCase):
    def test_fraction_polynomial_arithmetic(self):
        f=polynomial((Q(2),-Q(3),Q(1)))  # (t-1)(t-2)
        self.assertEqual(evaluate(f,1),0)
        self.assertEqual(evaluate(f,2),0)
        self.assertEqual(derivative(f),polynomial((-3,2)))
        self.assertEqual(times(polynomial((-1,1)),polynomial((-2,1))),f)
        q,r=divide(f,polynomial((-1,1)))
        self.assertEqual(q,polynomial((-2,1)))
        self.assertTrue(iszero(r))
        self.assertEqual(plus(f,neg(f)),ZERO)
        self.assertEqual(minus(ONE,ONE),ZERO)
        self.assertEqual(poly_gcd(f,derivative(f)),ONE)

    def test_exact_rational_sturm_roots(self):
        f=times(polynomial((-Q(1,2),1)),polynomial((-Q(1,3),1)))
        r=isolate_roots(f,0,1,width=Q(1,2**28))
        self.assertEqual(len(r),2)
        for q in (Q(1,3),Q(1,2)):
            self.assertTrue(any(e.lower<=q<=e.upper for e in r))
        self.assertTrue(any(e.rational and e.lower==Q(1,2) for e in r))
        self.assertEqual(count_roots_open(sturm_sequence(f),0,1),2)

    def test_irrational_and_negative_roots_sturm(self):
        r=isolate_roots(polynomial((-2,0,1)),-2,2,width=Q(1,2**32))
        self.assertEqual(len(r),2)
        self.assertTrue(r[0].upper<0 and r[1].lower>0)
        for e in r:
            self.assertFalse(e.rational)
            self.assertLess(e.upper-e.lower,Q(1,2**31))
            self.assertTrue(evaluate(polynomial((-2,0,1)),e.lower)*
                            evaluate(polynomial((-2,0,1)),e.upper)<0)

    def test_repeated_roots_squarefree(self):
        f=times(times(polynomial((-Q(1,3),1)),polynomial((-Q(1,3),1))),
                polynomial((Q(1,2),1)))
        self.assertEqual(len(squarefree(f)),3)
        rs=isolate_roots(f,-1,1,width=Q(1,2**20))
        self.assertEqual(len(rs),2)
        self.assertTrue(any(r.lower==r.upper==-Q(1,2) for r in rs))
        self.assertTrue(any(r.lower<=Q(1,3)<=r.upper for r in rs))

    def test_even_multiplicity_dyadic_root(self):
        f=times(polynomial((-Q(1,2),1)),polynomial((-Q(1,2),1)))
        self.assertEqual(tuple((r.lower,r.upper) for r in isolate_roots(f,0,1)),
                         ((Q(1,2),Q(1,2)),))

    def test_endpoint_roots_do_not_contaminate_open_count(self):
        f=polynomial((0,-1,1))
        st=sturm_sequence(f)
        self.assertEqual(count_roots_open(st,0,1),0)
        self.assertEqual(tuple((r.lower,r.upper) for r in isolate_roots(f,0,1)),
                         ((Q(0),Q(0)),(Q(1),Q(1))))

    def test_no_roots_identically_zero_and_budget(self):
        self.assertEqual(isolate_roots(polynomial((1,0,1)),-3,3),())
        self.assertEqual(isolate_roots(ONE,0,1),())
        with self.assertRaises(ValueError):isolate_roots(ZERO,0,1)
        with self.assertRaises(ValueError):isolate_roots(T,1,0)
        with self.assertRaises(RuntimeError):isolate_roots(polynomial((-2,0,1)),0,2,max_nodes=1)

    def test_polynomial_determinant(self):
        M=((ONE,T),(T,ONE))
        self.assertEqual(determinant_poly(M),polynomial((1,0,-1)))
        self.assertEqual(determinant_poly(((T,ONE),(T,ONE))),ZERO)


class MovingPolygonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=convex_vertex_alignment_polygon(Q(0))
        cls.R=convex_vertex_alignment_polygon(Q(9,8))

    def test_cramer_polynomials_match_independent_v223_solver(self):
        rng=random.Random(2240)
        done=0
        for _ in range(80):
            ids=tuple(rng.randrange(4) for __ in range(4))
            sigma=rng.choice((-1,1))
            c=chart_polynomials(self.P,self.R,ids,sigma,Q(1,40))
            t=Q(rng.randint(0,12),12)
            D=evaluate(c.determinant,t)
            if not D:continue
            S=interpolated_polygon(self.P,self.R,t)
            edges=tuple((S[i],S[(i+1)%4]) for i in ids)
            A,b=_affine_system(edges,sigma)
            u=_solve(A,b)
            self.assertIsNotNone(u)
            self.assertEqual(u,c.algebraic_parameters(t))
            q=_side_sq(_evaluate(edges,u))
            self.assertEqual(evaluate(c.cutoff,t),
                             (q-c.cutoff_squared)*D*D)
            done+=1
        self.assertGreater(done,40)

    def test_degree_bounds_determinant_and_cutoff(self):
        c=chart_polynomials(self.P,self.R,(0,0,1,2),1,Q(1,100))
        self.assertLessEqual(len(c.determinant)-1,4)
        self.assertLessEqual(max(len(v)-1 for v in c.numerators),4)
        self.assertLessEqual(len(c.cutoff)-1,10)
        self.assertFalse(c.identically_singular)

    def test_known_convex_vertex_event_is_exactly_certified(self):
        c=chart_polynomials(self.P,self.R,(0,0,1,2),1,Q(1,100))
        t=Q(8,9)
        self.assertEqual(evaluate(c.upper_differences[1],t),0)
        w=rational_vertex_square(self.P,self.R,c,t)
        self.assertEqual(w.squared_side,1)
        self.assertEqual(w.endpoint_dimensions,(1,))
        self.assertEqual(w.parameters,(Q(1,2),Q(1),Q(2,3),Q(1,2)))
        self.assertEqual(w.points,(pt(1,0),pt(2,0),pt(2,1),pt(1,1)))
        cover=isolate_chart_events(c,0,1,width=Q(1,2**24))
        self.assertTrue(any(e.kind=='edge-1-end' and
                            e.bracket.lower<=t<=e.bracket.upper for e in cover.events))

    def test_convex_polygon_remains_jordan_across_event(self):
        cert=certify_jordan_isotopy(self.P,self.R)
        self.assertEqual(cert.triples_checked,4)
        self.assertEqual(census_polygon(interpolated_polygon(self.P,self.R,Q(0))).count,1)
        self.assertEqual(census_polygon(interpolated_polygon(self.P,self.R,Q(8,9))).count,2)
        self.assertEqual(census_polygon(interpolated_polygon(self.P,self.R,Q(11,12))).count,3)

    def test_vertex_event_requires_repeated_edge_charts(self):
        selected=chart_polynomials(self.P,self.R,(0,0,1,2),1,Q(1,100))
        self.assertEqual(len(set(selected.edges)),3)
        w=rational_vertex_square(self.P,self.R,selected,Q(8,9))
        self.assertEqual(w.squared_side,1)
        # Distinct-open-edge-only diagnostics miss this endpoint event.
        covers,unresolved=polygon_event_cover(self.P,self.R,Q(1,100),
            start=Q(3,4),end=1,width=Q(1,2**12),include_repeated_edges=False)
        self.assertEqual(len(covers),48)
        self.assertFalse(unresolved)
        self.assertFalse(any(e.bracket.lower<=Q(8,9)<=e.bracket.upper
                             for c in covers for e in c.events))

    def test_global_cover_reports_persistent_chart_face_identities(self):
        covers,unresolved=polygon_event_cover(self.P,self.R,Q(1,100),
            start=Q(3,4),end=1,width=Q(1,2**12))
        self.assertGreater(len(covers),200)
        self.assertTrue(unresolved)
        self.assertTrue(any(e.bracket.lower<=Q(8,9)<=e.bracket.upper
                             for c in covers for e in c.events))
        # We don't claim a complete global cobordism in presence of
        # persistently-on-edge solutions or singular-in-time charts.
        self.assertGreater(sum(len(c.events) for c in covers),0)

    def test_wrong_chart_event_rejected_as_not_a_square(self):
        c=chart_polynomials(self.P,self.R,(0,1,2,3),1,Q(1,100))
        with self.assertRaises(ValueError):
            rational_vertex_square(self.P,self.R,c,Q(8,9))

    def test_persistently_singular_chart_reported_unresolved(self):
        base=tuple(pt(*q) for q in ((0,0),(2,0),(2,1),(0,1)))
        c=chart_polynomials(base,base,(0,0,0,0),1,Q(1,16))
        self.assertTrue(c.identically_singular)
        result=isolate_chart_events(c)
        self.assertIn('singular',result.unresolved)
        with self.assertRaises(ValueError):result.certified_chambers(0,1)

    def test_similarity_size_cutoff_with_irrational_time(self):
        # First compute a square from the exact v223 static polygon census.
        p=tuple(pt(*v) for v in ((0,0),(2,0),(3,3),(0,2)))
        s=census_polygon(p).squares[0]
        self.assertEqual(s.side_squared,Q(72,25))
        ids=tuple(int(v) for v in s.curve_positions)
        end=tuple(mul(v,2) for v in p)
        k=2*s.side_squared
        chart=chart_polynomials(p,end,ids,1,k)
        # Its square has side²=(1+t)²*72/25, so cutoff at sqrt(2)-1.
        r=isolate_roots(chart.cutoff,0,1,width=Q(1,2**28))
        # Additional algebraic factor at t=-1 cannot enter [0,1].
        self.assertEqual(len(r),1)
        self.assertFalse(r[0].rational)
        self.assertLess(r[0].lower,Q(1,2))
        self.assertGreater(r[0].upper,Q(1,3))
        self.assertEqual(chart.eligible(Q(0)),False)
        self.assertEqual(chart.eligible(Q(1)),True)
        self.assertEqual(census_polygon(p,k).count,0)
        self.assertEqual(census_polygon(end,k).count,1)

    def test_jordan_isotopy_rejects_self_intersection_deformation(self):
        a=tuple(pt(*x) for x in ((0,0),(1,0),(1,1),(0,1)))
        b=tuple(pt(*x) for x in ((0,0),(0,1),(1,1),(1,0)))
        validate_jordan_polygon(a);validate_jordan_polygon(b)
        with self.assertRaisesRegex(ValueError,'collinear'):
            certify_jordan_isotopy(a,b)

    def test_global_cover_rejects_non_jordan_intermediate_path(self):
        a=tuple(pt(*x) for x in ((0,0),(1,0),(1,1),(0,1)))
        b=tuple(pt(*x) for x in ((0,0),(0,1),(1,1),(1,0)))
        with self.assertRaisesRegex(ValueError,'collinear'):
            polygon_event_cover(a,b,Q(1,100),max_charts=20)

    def test_jordan_isotopy_rejects_persistent_collinearity(self):
        a=tuple(pt(*x) for x in ((0,0),(1,0),(2,0),(2,1),(0,1)))
        validate_jordan_polygon(a)
        with self.assertRaisesRegex(ValueError,'permanently collinear'):
            certify_jordan_isotopy(a,a)

    def test_reference_orientation_polynomial(self):
        p=self.P;q=self.R
        for i,j,k in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)):
            f=three_point_orientation_polynomial(p,q,i,j,k)
            for t in (Q(0),Q(1,3),Q(7,9),Q(1)):
                P=interpolated_polygon(p,q,t)
                self.assertEqual(evaluate(f,t),
                    cross(sub(P[j],P[i]),sub(P[k],P[i])))

    def test_certified_chambers_exclude_candidate_roots(self):
        c=chart_polynomials(self.P,self.R,(0,0,1,2),1,Q(1,100))
        cov=isolate_chart_events(c,0,1,width=Q(1,2**15))
        self.assertFalse(cov.unresolved)
        chambers=cov.certified_chambers(0,1)
        self.assertTrue(chambers)
        for a,b in chambers:
            mid=(a+b)/2
            self.assertTrue(all(evaluate(f,mid)!=0 for _,f in c.candidate_polynomials()))

    def test_inputs_and_budgets(self):
        c=chart_polynomials(self.P,self.R,(0,0,1,2),1,Q(1,100))
        with self.assertRaises(ValueError):chart_polynomials(self.P,self.R,(0,1),1,0)
        with self.assertRaises(ValueError):chart_polynomials(self.P,self.R,(0,0,1,2),0,0)
        with self.assertRaises(ValueError):chart_polynomials(self.P,self.R,(0,0,1,2),1,-1)
        with self.assertRaises(RuntimeError):polygon_event_cover(self.P,self.R,Q(1,100),max_charts=1)
        with self.assertRaises(ValueError):certify_jordan_isotopy(self.P,self.R,1,0)
        with self.assertRaises(ValueError):rational_vertex_square(self.P,self.R,chart_polynomials(self.P,self.R,(0,1,2,3),1,Q(1,100)),Q(8,9))

if __name__=='__main__':unittest.main()
