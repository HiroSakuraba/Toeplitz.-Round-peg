import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import _evaluate, _square_residual, census_polygon, pt
from balanced_sweep_v300 import cumulative
from tilted_sweep_v301 import (
    CORNER_RANGES, SAMPLE_TIMES, WINNER, audit_annulus, audit_charts, audit_ranges,
    audit_tube, corner_x, cramer, end_packet_cdf, mass_bounds, red_fraction, seed,
    strict_sign, winning_tags,
)


def edges_for(vertices, chart):
    return [(vertices[i], vertices[(i+1) % 4]) for i in chart]


class TiltedSweepTests(unittest.TestCase):
    def test_winner_formula_as_polynomial_identity(self):
        # Multiplying each residual by D gives degree <= 3. Four distinct
        # zeros suffice. The cyclic-order determinant has degree <= 4.
        for t in SAMPLE_TIMES:
            d = t/100
            points = _evaluate(edges_for(seed(d), WINNER), winning_tags(d))
            self.assertEqual(_square_residual(points, 1), (0,0,0,0))
            det = cramer(d, (3,1,2,3))[0]
            self.assertEqual(det, 9*(4*d*d+12*d+17)/256)
            self.assertEqual(corner_x(d), tuple(p.x for p in (points[3], *points[:3])))

    def test_all_distinct_chart_escapes(self):
        for d in (Q(1,1000000), Q(1,250), Q(1,100)):
            den, *nums = cramer(d, (0,1,2,3))
            self.assertEqual(den, -9*d*d/32)
            self.assertEqual(tuple(n/den for n in nums),
                             ((2*d+3)/(2*d), 3/(2*d), (4*d+3)/(4*d), 3/(4*d)))
        self.assertEqual(cramer(Q(0), (0,1,2,3))[0], 0)

    def test_interval_census_and_independent_full_label_solver(self):
        records = audit_charts()
        self.assertEqual(len(records), 13)
        self.assertEqual([r["edges"] for r in records if "interior" in r["certificate"]], [list(WINNER)])
        for d in (Q(1,100), Q(1,257), Q(1,1000000)):
            expected = _evaluate(edges_for(seed(d), WINNER), winning_tags(d))
            for sign in (1,-1):
                cs = census_polygon(seed(sign*d))
                self.assertEqual(cs.count, 1)
                self.assertEqual(cs.squares[0].uses_repeated_edge, True)
                points = expected if sign == 1 else [pt(1-p.x,p.y-d) for p in expected]
                self.assertEqual(set(cs.squares[0].vertices), set(points))

    def test_distinct_one_sided_square_limits(self):
        points = set(_evaluate(edges_for(seed(0), WINNER), winning_tags(Q(0))))
        self.assertIn(pt(0,0), points)
        self.assertNotIn(pt(1,0), points)
        self.assertNotEqual(points, {pt(1-p.x,p.y) for p in points})
        self.assertEqual(corner_x(Q(0)), (0,Q(9,17),Q(15,17),Q(6,17)))
        for z in (Q(-4,17),Q(-1,17),Q(0),Q(4,17)):
            s, t, w = Q(4,17)+z, Q(4,17)-z/2, Q(9,17)
            P = (pt(s,-3*s/4),pt(s+w,-3*(1-s-w)/4),
                 pt(t+w,3*(1-t-w)/2),pt(t,3*t/2))
            self.assertEqual(_square_residual(P, 1), (0,0,0,0))

    def test_punctured_bernstein_endpoint_behavior(self):
        self.assertEqual(strict_sign([t*t for t in SAMPLE_TIMES]), 1)
        self.assertEqual(strict_sign([-t for t in SAMPLE_TIMES]), -1)
        self.assertEqual(strict_sign([1-t for t in SAMPLE_TIMES]), 0)
        self.assertEqual(strict_sign([t-Q(1,2) for t in SAMPLE_TIMES]), 0)

    def test_original_sweep_misses_and_replacement_balances(self):
        audit_ranges()
        wide, tall, end = cumulative(0), cumulative(1), end_packet_cdf()
        self.assertGreater(mass_bounds(wide)[0], Q(1,2))
        self.assertGreater(mass_bounds(tall)[0], Q(1,2))
        self.assertLess(mass_bounds(end)[1], Q(1,2))
        for d in (Q(1,1000000),Q(1,250),Q(1,150),Q(1,100)):
            xs = corner_x(d)
            for x, (lo,hi) in zip(xs,CORNER_RANGES):
                self.assertLessEqual(lo,x)
                self.assertLessEqual(x,hi)
            q0, q1 = red_fraction(end,xs), red_fraction(tall,xs)
            u = (Q(1,2)-q0)/(q1-q0)
            self.assertTrue(0 < u < 1)
            self.assertEqual((1-u)*q0+u*q1,Q(1,2))
        xs = corner_x(Q(1,100))
        self.assertEqual(red_fraction(wide,xs),Q(58807,114136))
        self.assertEqual(red_fraction(tall,xs),Q(242015,342408))

    def test_degree_tube_at_both_shape_faces(self):
        from polygon_orbit_census_v223 import _solve
        for d in (Q(1,1000000),Q(1,100)):
            edges = edges_for(seed(d),(3,1,2,3))
            a,b,c,e = winning_tags(d)
            center = (e,a,b,c)
            for k in (1-d/20000,1+d/20000):
                def f(tags):
                    p,q,r,s = _evaluate(edges,tags)
                    dx,dy = q.x-p.x,q.y-p.y
                    return (r.x-q.x+k*dy,r.y-q.y-k*dx,s.x-p.x+k*dy,s.y-p.y-k*dx)
                origin = f([0]*4)
                cols = [f([int(i==j) for i in range(4)]) for j in range(4)]
                matrix = [[c[i]-origin[i] for c in cols] for i in range(4)]
                tags = _solve(matrix,[-v for v in origin])
                self.assertEqual(f(tags),(0,0,0,0))
                self.assertLess(max(abs(x-y) for x,y in zip(tags,center)),d/8)
                self.assertTrue(all(0 < v < 1 for v in tags))
        self.assertEqual(audit_tube()["preconditioned_error_per_drift"], "11/250")

    def test_rough_displacement_and_moving_packet(self):
        self.assertEqual(audit_annulus()["square_side_squared_lower"], "841/5000")
        # Distinct unaligned arcs, including a packet wrapping around the seam.
        red = [(Q(1,10),Q(7,20)),(Q(11,20),Q(4,5))]
        q = sum(b-a for a,b in red)
        def mass(start):
            end = start+Q(1,12)
            overlap = sum(max(0,min(b+k,end)-max(a+k,start)) for a,b in red for k in (0,1))
            return q/2 + 6*overlap
        self.assertGreater(mass(Q(1,8)),Q(1,2))
        self.assertLess(mass(Q(23,24)),Q(1,2))
        self.assertEqual(mass(Q(1,10)-Q(1,24)),Q(1,2))


if __name__ == "__main__":
    unittest.main()
