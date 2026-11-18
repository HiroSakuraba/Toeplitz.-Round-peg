"""Exact tests for source/lift-aware finite-dimensional relative degree (v214)."""
import random
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from periodic_bridge_v212 import Point, PeriodicGraph, rotate90
from periodic_square_cert_v213 import (
    segments_of, _affine_square_system, _solve_square, find_labeled_square)
from relative_degree_v214 import (
    relative_degree, _determinant, _cell_face_residual_bound,
    InadmissibleDegreeDomain,
)


def make_example(vertical_offsets=(F(0),)*4, period_shifts=(0,0,0,0),
                 edge=Point(F(1,8), F(1,8))):
    """Four graph windows supporting one transverse, noncollapsed square.

    `period_shifts` select lifted x positions for the four square vertices,
    while the source windows live in the quotient x mod 1.
    """
    origin = Point(F(1,2), F(0))
    r = rotate90(edge)
    planted = (origin, origin+edge, origin+edge+r, origin+r)
    slopes = (F(-2),F(3),F(7,2),F(-1,3))
    hw = F(1,32)
    graphs = []
    windows = []
    for i,(p,m,dy) in enumerate(zip(planted,slopes,vertical_offsets)):
        x = p.x - F(period_shifts[i])
        lo,hi = x-hw, x+hw
        if not (0 < lo < hi < 1):
            raise ValueError('choose window center strictly inside one period')
        base = p.y - m*hw + dy
        graphs.append(PeriodicGraph((Point(F(0),base),
                                     Point(lo,p.y-m*hw+dy),
                                     Point(hi,p.y+m*hw+dy),
                                     Point(F(1),base))))
        windows.append((lo,hi))
    return tuple(graphs),tuple(windows)


def with_vertical_shift(graphs, index, t):
    out=list(graphs)
    g=out[index]
    out[index]=PeriodicGraph(tuple(Point(p.x,p.y+t) for p in g.vertices))
    return tuple(out)


class RelativeDegreeTests(unittest.TestCase):
    def test_seeded_independent_exhaustive_solver_crosscheck(self):
        # Crosscheck isolated-degree enumerator with v213's independent
        # complete feasibility search on 48 perturbations. A degree-zero
        # configuration may still have two opposite-sign squares.
        original,ws=make_example()
        rng=random.Random(214214)
        for _ in range(48):
            offsets=tuple(F(rng.randint(-4,4),160) for i in range(4))
            gs=tuple(with_vertical_shift(original,i,v)[i]
                     for i,v in enumerate(offsets))
            cert=relative_degree(gs,ws,orientation=1)
            feasible=find_labeled_square(gs,ws,orientations=(1,),
                                         relative_lifts=(0,0,0,0))
            self.assertEqual(bool(cert.roots), not feasible.no_square)
            if cert.degree:
                self.assertIsNotNone(feasible.witness)

    def test_determinants_are_exact(self):
        self.assertEqual(_determinant(((F(1),F(2)),(F(3),F(4)))),F(-2))
        self.assertEqual(_determinant(((F(1),F(2)),(F(2),F(4)))),F(0))
        self.assertEqual(_determinant(((F(0),F(1)),(F(1),F(0)))),F(-1))

    def test_unique_interior_root_nonzero_degree(self):
        gs,ws=make_example()
        c=relative_degree(gs,ws)
        self.assertIn(c.degree,(-1,1))
        self.assertEqual(len(c.roots),1)
        self.assertEqual(c.roots[0].parameters,(F(1,2),)*4)
        self.assertEqual(c.roots[0].side_squared,F(1,32))
        self.assertGreater(c.opposite_height_gap,0)
        self.assertGreater(c.boundary_residual_lower_bound,0)
        self.assertGreater(c.vertical_perturbation_bound,0)
        self.assertEqual(c.singular_empty_cells,0)

    def test_two_orientations_separate(self):
        gs,ws=make_example()
        a=relative_degree(gs,ws,orientation=1)
        b=relative_degree(gs,ws,orientation=-1)
        self.assertIn(a.degree,(-1,1))
        self.assertEqual(b.degree,0)

    def test_small_vertical_perturbations_preserve_degree(self):
        gs,ws=make_example()
        c=relative_degree(gs,ws)
        eps=c.vertical_perturbation_bound/F(3)
        for shifts in ((eps,eps,eps,eps),
                       (eps,-eps,eps,-eps),
                       (-eps,eps,F(0),F(0))):
            perturbed=tuple(with_vertical_shift(gs,i,t)[i] for i,t in enumerate(shifts))
            d=relative_degree(perturbed,ws)
            self.assertEqual(d.degree,c.degree)
            self.assertEqual(len(d.roots),1)
            self.assertEqual(d.roots[0].shifts,(0,0,0,0))

    def test_deck_coherence_is_part_of_the_domain(self):
        # Construct a square of side sqrt(2) whose projected arc parameter
        # is the same for all four vertices, but whose physical lifts differ.
        gs,ws=make_example(edge=Point(F(-1),F(-1)),
                           period_shifts=(0,-1,0,1))
        c=relative_degree(gs,ws,relative_lifts=(0,-1,0,1))
        self.assertIn(c.degree,(-1,1))
        self.assertEqual(len(c.roots),1)
        self.assertEqual(c.roots[0].shifts,(0,-1,0,1))
        d=relative_degree(gs,ws,relative_lifts=(0,0,0,0))
        self.assertEqual(d.degree,0)
        self.assertEqual(len(d.roots),0)


    def test_singular_no_root_yields_degree_but_no_robustness_bound(self):
        gs=tuple(PeriodicGraph((Point(F(0),F(h)),Point(F(1),F(h))))
                 for h in (1,0,-1,0))
        ws=((F(7,16),F(9,16)),)*4
        cert=relative_degree(gs,ws)
        self.assertEqual(cert.degree,0)
        self.assertEqual(len(cert.roots),0)
        self.assertGreater(cert.singular_empty_cells,0)
        self.assertIsNone(cert.vertical_perturbation_bound)

    def test_internal_knot_root_is_rejected_not_double_counted(self):
        gs,ws=make_example()
        gs=list(gs)
        l,c,r=F(15,32),F(1,2),F(17,32)
        d=F(1,100)
        line=lambda x:F(1,4)-2*(x-F(1,2))
        yl=line(l)-d
        yc=line(c)
        yr=line(r)-d
        gs[2]=PeriodicGraph((Point(F(0),yl),Point(l,yl),
                             Point(c,yc),Point(r,yr),Point(F(1),yl)))
        with self.assertRaises(InadmissibleDegreeDomain):
            relative_degree(gs,ws)

    def test_two_squares_cancel_exact_signed_degree(self):
        # The first three source arcs remain unchanged; graph 2 has a
        # two-segment tent crossing the required p2 locus twice. Both
        # roots are interior and transverse, with opposite Jacobian signs.
        # Thus degree=0 does NOT mean there are no inscribed squares.
        gs,ws=make_example()
        gs=list(gs)
        l,c,r=F(15,32),F(1,2),F(17,32)
        delta=F(1,100)
        target=lambda x:F(1,4)-2*(x-F(1,2))
        yl=target(l)-delta
        yc=target(c)+delta
        yr=target(r)-delta
        gs[2]=PeriodicGraph((Point(F(0),yl),Point(l,yl),
                             Point(c,yc),Point(r,yr),Point(F(1),yl)))
        cert=relative_degree(tuple(gs),ws)
        self.assertEqual(cert.degree,0)
        self.assertEqual(len(cert.roots),2)
        self.assertEqual(tuple(r.segment_indices[2] for r in cert.roots),(0,1))
        self.assertEqual({1 if r.jacobian_det>0 else -1 for r in cert.roots}, {1,-1})
        self.assertGreater(cert.boundary_residual_lower_bound,0)
        self.assertGreater(cert.vertical_perturbation_bound,0)
        # Lowering the tent by 2 delta yields no squares on the four
        # marked source windows: both pairs leave via an internal PL knot.
        lowered=with_vertical_shift(tuple(gs),2,-2*delta)
        empty=relative_degree(lowered,ws)
        self.assertEqual(empty.degree,0)
        self.assertEqual(len(empty.roots),0)

    def test_no_source_region_square_zero_degree(self):
        gs,ws=make_example()
        separated=with_vertical_shift(gs,1,F(2))
        c=relative_degree(separated,ws)
        self.assertEqual(c.degree,0)
        self.assertEqual(c.roots,())

    def test_boundary_zero_requires_explicit_failure(self):
        gs,ws=make_example()
        # Move one graph through the first point at which the unique square
        # exits an original segment chart. The algorithm must not silently
        # assign an invariant when a zero reaches the source-domain boundary.
        segs=tuple(segments_of(g,w)[0] for g,w in zip(gs,ws))
        A,b=_affine_square_system(segs,(0,0,0,0),1)
        u0=_solve_square(A,b)
        moved=with_vertical_shift(gs,1,F(1))
        s1=tuple(segments_of(g,w)[0] for g,w in zip(moved,ws))
        A1,b1=_affine_square_system(s1,(0,0,0,0),1)
        u1=_solve_square(A1,b1)
        assert u0 is not None and u1 is not None
        crossings=[]
        for x,y in zip(u0,u1):
            delta=y-x
            for endpoint in (F(0),F(1)):
                if delta:
                    t=(endpoint-x)/delta
                    if t>0: crossings.append(t)
        t=min(crossings)
        boundary=with_vertical_shift(gs,1,t)
        with self.assertRaises(InadmissibleDegreeDomain):
            relative_degree(boundary,ws)
        after=relative_degree(with_vertical_shift(gs,1,t+F(1,50)),ws)
        self.assertEqual(after.degree,0)

    def test_invalid_domain_and_degenerate_guard(self):
        gs,ws=make_example()
        with self.assertRaises(ValueError):
            relative_degree(gs,ws,relative_lifts=(1,0,0,0))
        with self.assertRaises(ValueError):
            relative_degree(gs,ws,orientation=0)
        with self.assertRaises(ValueError):
            relative_degree(gs,((F(1,2),F(1,2)),)+ws[1:])
        g=list(gs)
        # Force 0 and 2 to overlap in their y ranges, independently of
        # whether a genuine nondegenerate square happens to exist.
        g[2]=with_vertical_shift(gs,2,-F(1,4))[2]
        with self.assertRaises(InadmissibleDegreeDomain):
            relative_degree(tuple(g),ws)

    def test_row_inverse_lower_bound(self):
        gs,ws=make_example()
        segs=tuple(segments_of(g,w)[0] for g,w in zip(gs,ws))
        A,b=_affine_square_system(segs,(0,0,0,0),1)
        m=_cell_face_residual_bound(A,b)
        self.assertGreater(m,F(0))
        self.assertLessEqual(m,F(1))


if __name__=='__main__':
    unittest.main()
