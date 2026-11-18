"""Independent exact checks for the finite PL periodic square decision theorem."""
import random
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from periodic_bridge_v212 import (  # noqa: E402
    Point, PeriodicGraph, SOURCE_WINDOWS, four_strands,
    is_labeled_square, rotate90, squared_length,
)
from periodic_square_cert_v213 import (  # noqa: E402
    Segment, _extreme_solutions, _intervals_allow_square,
    find_labeled_square, segments_of, solve_four_segments,
    source_window_exclusion, window_height_range, window_horizontal_gap,
    certify_transverse_square,
)


def graph_at_height(height):
    return PeriodicGraph((Point(F(0),F(height)),Point(F(1),F(height))))


def make_segments_around_square(origin, edge, orientation):
    re = rotate90(edge)
    re = Point(orientation*re.x,orientation*re.y)
    p0=origin
    p1=origin+edge
    p2=p1+re
    p3=p0+re
    slopes=(F(-2), F(3), F(7,2), F(-1,3))
    halfwidth=F(1,12)
    pieces=[]
    for p,m in zip((p0,p1,p2,p3),slopes):
        pieces.append(Segment(Point(p.x-halfwidth,p.y-m*halfwidth),
                              Point(p.x+halfwidth,p.y+m*halfwidth)))
    return tuple(pieces)


class FiniteSquareCertificateTests(unittest.TestCase):
    def test_linear_slice_rank_zero_only_degenerate(self):
        A=[[F(0)]*4 for _ in range(4)]
        b=[F(0)]*4
        vertices=list(_extreme_solutions(A,b))
        self.assertEqual(len(vertices),16)
        self.assertEqual(set(vertices),set(__import__('itertools').product((F(0),F(1)),repeat=4)))
        self.assertEqual(list(_extreme_solutions(A,[F(1),0,0,0])),[])

    def test_random_constructions_both_orientations(self):
        rng=random.Random(213)
        for sign in (1,-1):
            for _ in range(24):
                a,b=(F(rng.randint(-9,9),7) for _ in range(2))
                if a==0 and b==0: a=F(1,7)
                x,y=(F(rng.randint(-20,20),9) for _ in range(2))
                arcs=make_segments_around_square(Point(x,y),Point(a,b),sign)
                w=solve_four_segments(arcs,orientation=sign)
                self.assertIsNotNone(w)
                self.assertGreater(w.side_squared,F(0))
                # The four segments can support a whole family of squares;
                # the returned vertex need not equal the planted midpoint.
                planted=tuple(segment.point(F(1,2)) for segment in arcs)
                self.assertEqual(squared_length(planted[1]-planted[0]),a*a+b*b)
                self.assertTrue(all(F(0)<=u<=F(1) for u in w.parameters))
                if sign==1:
                    self.assertTrue(is_labeled_square(w.points))
                else:
                    # Original v212 checker selects only the positive orientation.
                    self.assertFalse(is_labeled_square(w.points))

    def test_rank_deficient_horizontal_family(self):
        rows=(graph_at_height(1),graph_at_height(0),
              graph_at_height(-1),graph_at_height(0))
        w=find_labeled_square(rows)
        self.assertIsNotNone(w.witness)
        self.assertEqual(w.witness.side_squared,F(2))
        self.assertTrue(is_labeled_square(w.witness.points))
        self.assertGreater(w.rejected_cases,0)

    def test_incompatible_height_levels_exhaustively_rejected(self):
        rows=tuple(graph_at_height(h) for h in (0,0,0,1))
        r=find_labeled_square(rows)
        self.assertTrue(r.no_square)
        self.assertGreater(r.rejected_cases,0)

    def test_zero_height_case_and_singleton_windows(self):
        flat=(graph_at_height(0),)*4
        self.assertTrue(find_labeled_square(flat).no_square)
        rows=(graph_at_height(1),graph_at_height(0),graph_at_height(-1),graph_at_height(0))
        at_point=((F(1,2),F(1,2)),)*4
        # Horizontal-lifted square has all four period parameters equal.
        self.assertEqual(len(segments_of(rows[0],at_point[0])),1)
        w=find_labeled_square(rows,at_point)
        self.assertIsNotNone(w.witness)
        self.assertEqual(w.witness.side_squared,F(2))
        self.assertRaises(ValueError,find_labeled_square,rows,None,(1,0))

    def test_own_windows_are_square_free_exactly(self):
        rows=four_strands()
        result=find_labeled_square(rows,SOURCE_WINDOWS)
        self.assertTrue(result.no_square)
        self.assertGreater(result.rejected_cases,0)
        self.assertEqual(result.solved_cases,0)
        self.assertEqual(result.lift_bound,4)
        cert=source_window_exclusion(rows,SOURCE_WINDOWS)
        self.assertTrue(cert.proves_absence)
        self.assertEqual(cert.horizontal_separation,F(7,16))
        self.assertEqual(cert.maximum_opposite_height_difference,F(1,16))
        self.assertEqual(cert.margin,F(3,8))

    def test_common_windows_do_contain_known_square(self):
        rows=four_strands()
        common=((F(7,8),F(7,8)),)*4
        result=find_labeled_square(rows,common)
        self.assertFalse(result.no_square)
        self.assertEqual(result.witness.side_squared,F(2))
        self.assertFalse(source_window_exclusion(rows,common).proves_absence)

    def test_rational_window_margin_for_amplitude_family(self):
        for amp in (F(0),F(1,128),F(1,16),F(1,8)):
            graphs=four_strands(amp)
            cert=source_window_exclusion(graphs,SOURCE_WINDOWS)
            self.assertTrue(cert.proves_absence)
            self.assertEqual(cert.margin,F(7,16)-amp)
            self.assertTrue(find_labeled_square(graphs,SOURCE_WINDOWS).no_square)

    def test_window_range_extrema_attained(self):
        g=four_strands()[1]
        self.assertEqual(window_height_range(g,SOURCE_WINDOWS[1]),(-F(1,16),F(0)))
        self.assertEqual(window_horizontal_gap(SOURCE_WINDOWS[0],SOURCE_WINDOWS[2]),F(7,16))
        self.assertEqual(window_horizontal_gap((F(0),F(1,4)),(F(7,8),F(1))),F(0))
        self.assertEqual(window_horizontal_gap((F(1,8),F(3,8)),(F(5,8),F(7,8))),F(1,4))

    def test_interval_screen_never_rejects_known_square(self):
        points=(Point(F(7,8),F(1)),Point(-F(1,8),F(0)),
                Point(F(7,8),-F(1)),Point(F(15,8),F(0)))
        arcs=make_segments_around_square(points[0],points[1]-points[0],1)
        self.assertTrue(_intervals_allow_square(arcs,(0,0,0,0),1))
        self.assertFalse(_intervals_allow_square(arcs,(0,0,30,0),1))

    def test_projected_source_localization_is_not_lift_coherence(self):
        # Every projected vertex lies in the same marked parameter window.
        # The periodic square uses deck lifts 0,-1,0,+1; no planar square
        # exists when all original arcs occupy their common fundamental cell.
        graphs=tuple(graph_at_height(h) for h in (1,0,-1,0))
        windows=((F(7,16),F(9,16)),)*4
        periodic=find_labeled_square(graphs,windows)
        self.assertFalse(periodic.no_square)
        self.assertEqual(periodic.witness.side_squared,F(2))
        self.assertNotEqual(periodic.witness.shifts,(0,0,0,0))
        planar=find_labeled_square(graphs,windows,relative_lifts=(0,0,0,0))
        self.assertTrue(planar.no_square)
        specified=find_labeled_square(graphs,windows,relative_lifts=(0,-1,0,1))
        self.assertFalse(specified.no_square)
        self.assertEqual(specified.witness.shifts,(0,-1,0,1))
        self.assertRaises(ValueError,find_labeled_square,graphs,windows,(1,-1),(1,0,0,0))

    def test_14_planted_periodic_squares_with_fixed_lifts(self):
        # Constructs genuinely piecewise-linear nonconstant height graphs
        # passing through a square whose x coordinates span multiple lifts.
        rng=random.Random(77)
        for _ in range(14):
            p0=Point(F(rng.choice([1,3,5,7]),8), F(rng.randint(-3,3),8))
            edge=Point(F(rng.randint(-3,3),4), F(rng.randint(-3,3),4))
            if squared_length(edge)==0:
                edge=Point(F(1,4),F(1,4))
            re=rotate90(edge)
            planted=(p0,p0+edge,p0+edge+re,p0+re)
            graphs=[];shifts=[]
            for p in planted:
                q=p.x%1;shifts.append(int(p.x-q))
                vs=(Point(F(0),p.y),Point(F(1),p.y)) if q==0 else (
                    Point(F(0),p.y+F(1,10)),Point(q,p.y),
                    Point(F(1),p.y+F(1,10)))
                graphs.append(PeriodicGraph(vs))
            r=find_labeled_square(graphs,relative_lifts=tuple(shifts))
            self.assertIsNotNone(r.witness)
            self.assertTrue(is_labeled_square(r.witness.points))

    def test_large_integer_lifts_are_finitely_bounded(self):
        graphs=tuple(graph_at_height(h) for h in (3,0,-3,0))
        r=find_labeled_square(graphs)
        self.assertFalse(r.no_square)
        self.assertEqual(r.lift_bound,7)
        self.assertEqual(r.witness.side_squared,F(18))
        self.assertTrue(any(abs(k)>1 for k in r.witness.shifts))

    def test_transverse_persistence_bound_and_exact_perturbations(self):
        arcs=make_segments_around_square(Point(F(7,8),F(1)),
                                         Point(-F(1),-F(1)),1)
        cert=certify_transverse_square(arcs)
        self.assertEqual(cert.parameter_margin,F(1,2))
        self.assertEqual(cert.inverse_inf_norm,F(819,55))
        self.assertEqual(cert.opposite_height_gap,F(37,24))
        self.assertEqual(cert.epsilon_bound,F(55,3276))
        rng=random.Random(214)
        eps=cert.epsilon_bound/F(2)
        for _ in range(24):
            perturbed=[]
            for s in arcs:
                dy0=eps*F(rng.randint(-50,50),50)
                dy1=eps*F(rng.randint(-50,50),50)
                perturbed.append(Segment(Point(s.start.x,s.start.y+dy0),
                                         Point(s.end.x,s.end.y+dy1)))
            w=solve_four_segments(tuple(perturbed))
            self.assertIsNotNone(w)
            self.assertTrue(all(0<u<1 for u in w.parameters))
            self.assertGreater(w.side_squared,F(0))
        # Rank-deficient families do not satisfy the transverse criterion.
        flat=tuple(Segment(Point(F(0),F(h)),Point(F(1),F(h)))
                   for h in (1,0,-1,0))
        with self.assertRaises(ValueError):
            certify_transverse_square(flat)

    def test_pruning_does_not_discard_48_planted_squares(self):
        rng=random.Random(2131)
        for sign in (1,-1):
            for _ in range(24):
                edge=Point(F(rng.randint(-7,7),9),F(rng.randint(-7,7),9))
                if squared_length(edge)==0:
                    edge=Point(F(1,9),F(0))
                origin=Point(F(rng.randint(-7,7),11),F(rng.randint(-7,7),11))
                arcs=make_segments_around_square(origin,edge,sign)
                self.assertTrue(_intervals_allow_square(arcs,(0,0,0,0),sign))

    def test_input_validation(self):
        rows=four_strands()
        with self.assertRaises(ValueError):
            find_labeled_square(rows[:3])
        with self.assertRaises(ValueError):
            find_labeled_square(rows,((-F(1),F(0)),)*4)
        with self.assertRaises(ValueError):
            find_labeled_square(rows,orientations=())
        with self.assertRaises(ValueError):
            window_horizontal_gap((F(0),F(2)),(F(0),F(1)))


if __name__=='__main__':
    unittest.main()
