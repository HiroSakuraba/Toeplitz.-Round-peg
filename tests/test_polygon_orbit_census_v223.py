"""v223 exact, complete finite square-orbit census on ONE embedded Jordan polygon."""
import sys
import unittest
import random
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from polygon_orbit_census_v223 import (
    pt,add,sub,mul,dot,Point,PositiveDimensionalSquares,
    validate_jordan_polygon,census_polygon,canonical_orbit,_square_residual,
    _affine_system,_rank,_solve,_box_vertices, _can_share_edge,
    similarity_polygon,similarity_cutoff_ledger,convex_vertex_alignment_polygon,
    vertex_alignment_snapshots,report,
)

GENERIC=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
TRAPEZOID=tuple(pt(*xy) for xy in ((0,0),(2,0),(2,Q(3,2)),(0,Q(1,2))))
PENTAGON=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,1),(2,2),(0,1)))
UNIT=tuple(pt(*xy) for xy in ((0,0),(1,0),(1,1),(0,1)))
RECTANGLE=tuple(pt(*xy) for xy in ((0,0),(2,0),(2,1),(0,1)))


class PolygonSquareCensusTests(unittest.TestCase):
    def test_generic_polygon_complete_one_geometric_orbit(self):
        c=census_polygon(GENERIC)
        self.assertEqual(c.count,1)
        self.assertEqual(c.parity,1)
        self.assertEqual(c.graceful_count,1)
        self.assertEqual(c.nongraceful_count,0)
        self.assertEqual(c.squares[0].side_squared,Q(72,25))
        self.assertEqual(c.squares[0].labelings_checked,8)
        self.assertFalse(c.squares[0].has_polygon_vertex)
        self.assertFalse(c.squares[0].uses_repeated_edge)
        self.assertEqual(c.squares[0].cyclic_order_sector,(0,1,2,3))
        self.assertEqual(len(set(c.squares[0].vertices)),4)
        self.assertTrue(c.singular_charts>0)
        self.assertEqual(c.redundant_labelings,7)

    def test_two_isolated_orbits_on_trapezoid_repeated_edges_and_endpoints(self):
        c=census_polygon(TRAPEZOID)
        self.assertEqual(c.count,2)
        self.assertEqual(c.parity,0)
        self.assertEqual({x.side_squared for x in c.squares},{Q(45,49),Q(1)})
        self.assertTrue(all(x.uses_repeated_edge for x in c.squares))
        self.assertEqual([x.labelings_checked for x in c.squares],[8,16])
        self.assertEqual(sum(x.has_polygon_vertex for x in c.squares),1)
        self.assertTrue(all(x.graceful for x in c.squares))

    def test_pentagon_vertex_root_duplicated_in_edge_charts_not_orbits(self):
        c=census_polygon(PENTAGON)
        self.assertEqual(c.count,1)
        self.assertEqual(c.squares[0].side_squared,Q(20,9))
        self.assertEqual(c.squares[0].labelings_checked,16)
        self.assertTrue(c.squares[0].has_polygon_vertex)
        self.assertTrue(c.squares[0].uses_repeated_edge)

    def test_unit_square_has_continuum_and_cannot_have_finite_parity(self):
        with self.assertRaises(PositiveDimensionalSquares) as cm:
            census_polygon(UNIT)
        self.assertGreater(cm.exception.max_side_sq,0)
        self.assertGreater(cm.exception.vertex_count,1)
        self.assertEqual(_square_residual(cm.exception.sample_square,
                                          cm.exception.orientation),(0,0,0,0))
        self.assertEqual(cm.exception.max_side_sq,
                         dot(sub(cm.exception.sample_square[1],cm.exception.sample_square[0]),
                             sub(cm.exception.sample_square[1],cm.exception.sample_square[0])))

    def test_rectangle_has_continuum_on_two_repeated_edges(self):
        with self.assertRaises(PositiveDimensionalSquares) as cm:
            census_polygon(RECTANGLE)
        self.assertLess(len(set(cm.exception.edge_tuple)),4)

    def test_large_cutoff_excludes_all_even_on_continuum_polygons(self):
        self.assertEqual(census_polygon(UNIT,Q(10)).count,0)
        self.assertEqual(census_polygon(RECTANGLE,Q(100)).count,0)
        self.assertEqual(census_polygon(GENERIC,Q(72,25)).count,0)
        self.assertEqual(census_polygon(GENERIC,Q(72,25)+Q(1,100)).count,0)

    def test_strict_cutoff_keeps_only_large_trapezoid_orbit(self):
        c=census_polygon(TRAPEZOID,Q(19,20))
        self.assertEqual(c.count,1)
        self.assertEqual(c.squares[0].side_squared,Q(1))

    def test_cyclic_polygon_starting_vertex_is_not_a_new_geometric_square(self):
        for k in range(4):
            c=census_polygon(GENERIC[k:]+GENERIC[:k])
            self.assertEqual(c.count,1)
            self.assertEqual(c.squares[0].side_squared,Q(72,25))
            self.assertEqual(c.squares[0].vertices,census_polygon(GENERIC).squares[0].vertices)

    def test_polygon_parameter_reversal_changes_sector_not_existence(self):
        c=census_polygon(tuple(reversed(GENERIC)))
        self.assertEqual(c.count,1)
        self.assertEqual(c.squares[0].cyclic_order_sector,(0,3,2,1))
        self.assertTrue(c.squares[0].graceful)

    def test_rigid_and_rational_similarity_covariance(self):
        base=census_polygon(GENERIC).squares
        def transform(p):return pt(-3*p.y+Q(7,5),3*p.x-Q(13,11))
        c=census_polygon(tuple(transform(p) for p in GENERIC))
        self.assertEqual(c.count,1)
        self.assertEqual(c.squares[0].side_squared,9*base[0].side_squared)
        self.assertEqual(c.squares[0].labelings_checked,8)
        reflected=census_polygon(tuple(pt(-p.x,p.y) for p in GENERIC))
        self.assertEqual(reflected.count,1)
        self.assertEqual(reflected.squares[0].side_squared,base[0].side_squared)
        self.assertEqual(reflected.graceful_count,1)

    def test_exact_affine_square_system_solves_plant(self):
        edges=tuple((GENERIC[i],GENERIC[(i+1)%4]) for i in range(4))
        A,b=_affine_system(edges,1)
        self.assertEqual(_rank(A),4)
        sol=_solve(A,b)
        self.assertIsNotNone(sol)
        self.assertTrue(all(0<u<1 for u in sol))
        self.assertTrue(all(sum(t*u for t,u in zip(row,sol))==v for row,v in zip(A,b)))

    def test_rank_deficient_feasible_polytope_vertices_exact(self):
        A=[[Q(1),0,0,0],[Q(2),0,0,0],[0,0,0,0],[0,0,0,0]]
        b=[Q(1,2),Q(1),Q(0),Q(0)]
        V=set(_box_vertices(A,b))
        self.assertEqual(len(V),8)
        self.assertEqual({v[0] for v in V},{Q(1,2)})
        self.assertEqual(_rank(A),1)
        self.assertFalse(tuple(_box_vertices(A,[Q(1,2),Q(0),Q(0),Q(0)])))

    def test_known_non_square_order_fails_canonical_normalization(self):
        vertices=(pt(0,0),pt(1,1),pt(1,0),pt(0,1))
        with self.assertRaises(ValueError):canonical_orbit(vertices,tuple(map(Q,range(4))))

    def test_repeated_edge_classification_stable_under_vertex_assignment(self):
        self.assertTrue(_can_share_edge((Q(0),Q(1),Q(2),Q(3)),4))
        self.assertFalse(_can_share_edge((Q(1,4),Q(5,4),Q(9,4),Q(13,4)),4))
        self.assertTrue(_can_share_edge((Q(1,4),Q(3,4),Q(9,4),Q(13,4)),4))

    def test_invalid_polygons_are_rejected_before_square_search(self):
        bad=[((0,0),(2,2),(0,2),(2,0)), # bowtie
             ((0,0),(1,0),(1,1),(0,0)), # repeat
             ((0,0),(2,0),(1,0),(1,1),(0,1)), # adjacent backtracking
             ((0,0),(1,0),(2,0))] # zero area
        for pairs in bad:
            with self.subTest(pairs=pairs):
                with self.assertRaises(ValueError):census_polygon(tuple(pt(*x) for x in pairs))

    def test_similarity_family_is_jordan_and_has_fixed_unfiltered_orbits(self):
        for t in (Q(0),Q(1,7),Q(1,2),Q(1)):
            P=similarity_polygon(GENERIC,t,Q(1,2),Q(1),
                                 (Q(3,11),-Q(1,9)),(Q(2,3),Q(4,5)))
            validate_jordan_polygon(P)
            c=census_polygon(P)
            self.assertEqual(c.count,1)
            self.assertEqual(c.squares[0].side_squared,Q(72,25)*(Q(1,2)+t)**2)

    def test_similarity_one_irrational_cutoff_event_certified_by_rational_bracket(self):
        ledger=similarity_cutoff_ledger(GENERIC,0,1,Q(1),Q(1,2),Q(1))
        self.assertEqual((ledger.initial_orbits,ledger.final_orbits),(0,1))
        self.assertEqual(len(ledger.events),1)
        self.assertTrue(ledger.balances)
        event=ledger.events[0]
        self.assertIsNone(event.exact_time)
        self.assertLess(event.upper-event.lower,Q(1,2**39))
        lo,hi=event.lower,event.upper
        q=event.base_side_squared
        self.assertLess(q*(Q(1,2)+lo)**2,1)
        self.assertGreater(q*(Q(1,2)+hi)**2,1)

    def test_similarity_rational_event_time(self):
        led=similarity_cutoff_ledger(TRAPEZOID,0,1,Q(1,4),Q(1,4),Q(1))
        self.assertEqual((led.initial_orbits,led.final_orbits),(0,2))
        self.assertEqual(len(led.events),2)
        self.assertTrue(led.balances)
        self.assertIn(Q(1,4),tuple(e.exact_time for e in led.events))
        self.assertEqual(sum(e.exact_time is None for e in led.events),1)

    def test_similarity_time_reversal(self):
        a=similarity_cutoff_ledger(GENERIC,0,1,Q(1),Q(1,2),Q(1))
        b=similarity_cutoff_ledger(GENERIC,0,1,Q(1),Q(3,2),-Q(1))
        self.assertEqual((a.initial_orbits,a.final_orbits),(b.final_orbits,b.initial_orbits))
        self.assertEqual(len(a.events),len(b.events))
        self.assertEqual(a.parity_change,b.parity_change)

    def test_similarity_constant_translation_has_no_events(self):
        a=similarity_cutoff_ledger(GENERIC,0,1,Q(1),Q(1),Q(0),
                                   (Q(10),Q(12)),(-Q(7),Q(1)))
        self.assertEqual((a.initial_orbits,a.final_orbits),(1,1))
        self.assertEqual(a.events,())

    def test_similarity_generic_time_slices_enforced(self):
        with self.assertRaisesRegex(ValueError,'time slice'):
            similarity_cutoff_ledger(TRAPEZOID,0,1,Q(1,4),Q(1,2),1)
        with self.assertRaises(ValueError):
            similarity_cutoff_ledger(GENERIC,0,1,0)
        with self.assertRaises(ValueError):
            similarity_cutoff_ledger(GENERIC,0,1,Q(1),Q(0),1)

    def test_orbit_keys_match_independent_v218_distinct_open_edge_auditor(self):
        try:
            from configuration_symmetry_v218 import (pt as vpt,isolated_distinct_edge_squares)
        except ImportError:
            self.skipTest('v218 available only in full checkout')
        P=tuple(vpt(p.x,p.y) for p in GENERIC)
        independent=isolated_distinct_edge_squares(P)
        c=census_polygon(GENERIC)
        self.assertEqual(len(independent),1)
        self.assertEqual(c.count,len(independent))
        self.assertEqual(next(iter(independent.values())).__len__(),8)

    def test_polygon_validator_matches_independent_v216(self):
        try:
            from boundary_degree_v216 import pt as vpt,validate_simple_polygon
        except ImportError:
            self.skipTest('v216 available only in full checkout')
        for P in (GENERIC,TRAPEZOID,PENTAGON):
            self.assertNotEqual(validate_simple_polygon(tuple(vpt(p.x,p.y) for p in P)),0)
            self.assertEqual(census_polygon(P).count,census_polygon(tuple(P)).count)

    def test_convex_vertex_alignment_challenges_nongeneric_orbit_parity(self):
        samples=vertex_alignment_snapshots()
        self.assertEqual(tuple((t,c.count) for t,c in samples),
                         ((Q(0),1),(Q(63,64),1),(Q(1),2),(Q(33,32),3)))
        self.assertEqual(sum(s.has_polygon_vertex for s in samples[2][1].squares),1)
        self.assertEqual(sum(s.has_polygon_vertex for s in samples[3][1].squares),0)
        self.assertTrue(all(s.graceful for _,c in samples for s in c.squares))
        for t in (Q(0),Q(17,25),Q(1),Q(5,4),Q(79,60)):
            P=convex_vertex_alignment_polygon(t)
            turns=[]
            for i in range(4):
                a=sub(P[(i+1)%4],P[i])
                b=sub(P[(i+2)%4],P[(i+1)%4])
                turns.append(a.x*b.y-a.y*b.x)
            self.assertTrue(all(x>0 for x in turns))
        with self.assertRaises(ValueError):convex_vertex_alignment_polygon(Q(4,3))

    def test_report_uses_complete_polygon_census_not_distinct_edge_only(self):
        r=report()
        self.assertEqual(r['geometric_square_orbits'],1)
        self.assertEqual(r['labeling_duplications'],7)
        self.assertEqual(r['convex_vertex_alignment_snapshots'],
                         (('0',1),('63/64',1),('1',2),('33/32',3)))

if __name__=='__main__':
    unittest.main()
