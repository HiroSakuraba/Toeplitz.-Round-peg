"""v264: generic exact convex event ledger, refinement stability and fail-closed tests."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt,add,mul,census_polygon
from big_peg_pair_death_v254 import (
    family_vertices,ROOT_LEFT,ROOT_RIGHT,QPOLY,
)
from big_peg_multichart_v249 import ppoint_at,peval
from convex_square_event_atlas_v264 import (
    InconclusiveAtlas,_strict_sign,_valid_endpoints,
    compress_collinear_refinement,certify_strict_convex_motion,
    certify_uniform_side_floor,classify_chart,
    audit_moving_convex_polygon,
)


def moving_pentagon():
    P=family_vertices()
    return (tuple(ppoint_at(v,ROOT_LEFT) for v in P),
            tuple(ppoint_at(v,ROOT_RIGHT) for v in P))


def translated(P,offset):
    return tuple(add(v,offset) for v in P)


def scaled(P,s):
    return tuple(mul(v,s) for v in P)


def subdivided(P,fraction=Q(1,2)):
    n=len(P)
    out=[]
    for i,p in enumerate(P):
        nxt=P[(i+1)%n]
        out.append(p)
        out.append(add(p,mul(pt(nxt.x-p.x,nxt.y-p.y),fraction)))
    return tuple(out)


class ReusableConvexAtlas(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.before,cls.after=moving_pentagon()
        cls.pentagon=audit_moving_convex_polygon(cls.before,cls.after)
        cls.static_translated=audit_moving_convex_polygon(
            cls.before,translated(cls.before,pt(Q(1,7),-Q(2,11))))
        cls.static_scaled=audit_moving_convex_polygon(
            cls.after,scaled(cls.after,Q(3,2)))
        cls.refined=audit_moving_convex_polygon(
            subdivided(cls.before),subdivided(cls.after))

    def test_nonpentagon_static_convex_quadrilateral(self):
        P=tuple(pt(x,y) for x,y in (
            (0,0),(4,0),(3,2),(-1,3)))
        B=translated(P,pt(Q(1,10),-Q(1,7)))
        r=audit_moving_convex_polygon(P,B)
        self.assertEqual(r["initial_vertex_count"],4)
        self.assertEqual(r["compressed_vertex_count"],4)
        self.assertEqual(r["all_weakly_increasing_support_charts"],35)
        self.assertEqual(r["possible_event_clusters"],[])
        count=census_polygon(P).count
        self.assertEqual(r["complete_all_real_time_geometric_square_orbits"],
                         count)
        self.assertGreater(count,0)
        self.assertEqual(r["complete_all_real_time_quotient_parity"],count%2)

    def test_linear_translation_is_complete_for_all_real_times(self):
        r=self.static_translated
        self.assertEqual(r["possible_event_clusters"],[])
        self.assertEqual(r["complete_all_real_time_geometric_square_orbits"],3)
        self.assertEqual(r["endpoint_geometric_square_orbits"],[3,3])
        self.assertEqual(r["complete_all_real_time_quotient_parity"],1)
        self.assertIn("complete all-real-time",r["status"])

    def test_positive_homothety_is_complete_for_all_real_times(self):
        r=self.static_scaled
        self.assertEqual(r["possible_event_clusters"],[])
        self.assertEqual(r["complete_all_real_time_geometric_square_orbits"],1)
        self.assertEqual(r["endpoint_geometric_square_orbits"],[1,1])
        self.assertEqual(r["complete_all_real_time_quotient_parity"],1)

    def test_genuine_v254_algebraic_event_is_localized(self):
        r=self.pentagon
        self.assertEqual(r["initial_vertex_count"],5)
        self.assertEqual(r["all_weakly_increasing_support_charts"],70)
        self.assertEqual(r["unresolved_charts"],{})
        self.assertEqual(r["endpoint_geometric_square_orbits"],[3,1])
        self.assertIsNone(r["complete_all_real_time_geometric_square_orbits"])
        self.assertGreater(len(r["possible_event_clusters"]),0)
        # Check that some exact rational isolator strictly encloses
        # the existing algebraically proved irrational collision.
        found=False
        for e in r["possible_event_clusters"]:
            lo=ROOT_LEFT+(ROOT_RIGHT-ROOT_LEFT)*Q(e["lower"])
            hi=ROOT_LEFT+(ROOT_RIGHT-ROOT_LEFT)*Q(e["upper"])
            if peval(QPOLY,lo)*peval(QPOLY,hi)<0:
                found=True
        self.assertTrue(found)

    def test_v263_irrational_event_slice_never_guessed(self):
        r=self.pentagon
        self.assertTrue(any(
            e["irrational_or_merged_event_slice_unresolved"] for
            e in r["possible_event_clusters"]))
        self.assertTrue(any(
            e["exact_event_geometric_square_orbits"] is None for
            e in r["possible_event_clusters"]))
        self.assertIn("event-free chamber",r["status"])

    def test_v254_chambers_include_three_before_and_one_after(self):
        r=self.pentagon
        chambers=r["rational_event_free_chambers"]
        self.assertGreaterEqual(len(chambers),2)
        self.assertEqual(chambers[0]["geometric_square_orbits"],3)
        self.assertEqual(chambers[-1]["geometric_square_orbits"],1)
        self.assertEqual(chambers[0]["quotient_parity"],1)
        self.assertEqual(chambers[-1]["quotient_parity"],1)

    def test_every_pentagon_chart_covered_by_exclusion_or_event(self):
        r=self.pentagon
        self.assertEqual(sum(r["chart_classification"].values()),70)
        self.assertEqual(len(r["full_chart_ledger"]),70)
        self.assertEqual(r["unresolved_charts"],{})

    def test_exact_refinement_collapse_to_original_five_edges(self):
        A,B,removed=compress_collinear_refinement(
            subdivided(self.before),subdivided(self.after))
        self.assertEqual(A,self.before)
        self.assertEqual(B,self.after)
        self.assertEqual(len(removed),5)
        self.assertEqual(self.refined["compressed_vertex_count"],5)
        self.assertEqual(self.refined["initial_vertex_count"],10)
        self.assertEqual(self.refined["all_weakly_increasing_support_charts"],70)

    def test_refinement_does_not_change_gap_or_event_chambers(self):
        a,b=self.pentagon,self.refined
        self.assertEqual(a["uniform_positive_square_size_certificate"],
                         b["uniform_positive_square_size_certificate"])
        self.assertEqual(a["endpoint_geometric_square_orbits"],
                         b["endpoint_geometric_square_orbits"])
        self.assertEqual(a["possible_event_clusters"],
                         b["possible_event_clusters"])
        self.assertEqual(a["rational_event_free_chambers"],
                         b["rational_event_free_chambers"])

    def test_refinement_with_nonmidpoint_fraction_still_exact(self):
        A,B,removed=compress_collinear_refinement(
            subdivided(self.before,Q(1,3)),
            subdivided(self.after,Q(1,3)))
        self.assertEqual((A,B),(self.before,self.after))
        self.assertEqual(len(removed),5)

    def test_tiny_bent_vertex_is_not_compressed(self):
        P=subdivided(self.before)
        bent=list(P)
        target=bent[1]
        bent[1]=pt(target.x,target.y+Q(1,100000))
        A,B,removed=compress_collinear_refinement(tuple(bent),tuple(bent))
        self.assertIn(bent[1],A)
        self.assertNotIn(1,removed)
        self.assertGreaterEqual(len(A),6)

    def test_exact_uniform_nonadjacent_edge_gap_is_positive(self):
        for r in (self.pentagon,self.static_translated,self.static_scaled):
            gap=r["uniform_positive_square_size_certificate"]
            self.assertGreater(Q(gap[
                "edge_pair_distance_squared_lower_bound"]),0)
            self.assertEqual(
                Q(gap["every_positive_square_side_squared_lower_bound"]),
                Q(gap["edge_pair_distance_squared_lower_bound"])/2)
            self.assertEqual(
                Q(gap["safe_side_squared_cutoff"]),
                Q(gap["edge_pair_distance_squared_lower_bound"])/4)

    def test_all_chamber_squares_exceed_certified_side_floor(self):
        for r in (self.pentagon,self.static_translated,
                  self.static_scaled,self.refined):
            floor=Q(r["uniform_positive_square_size_certificate"][
                "every_positive_square_side_squared_lower_bound"])
            for ch in r["rational_event_free_chambers"]:
                self.assertTrue(all(Q(x)>=floor
                    for x in ch["sample_square_side_squared"]))

    def test_regular_signed_local_index_available_on_generic_chambers(self):
        r=self.pentagon
        chambers=r["rational_event_free_chambers"]
        self.assertTrue(chambers[0][
            "all_sample_squares_regular_in_open_positive_sector"])
        self.assertEqual(chambers[0]["regular_oriented_index_sum"],1)
        self.assertEqual(chambers[-1]["regular_oriented_index_sum"],1)
        self.assertEqual(self.static_translated[
            "rational_event_free_chambers"][0]["regular_oriented_index_sum"],1)

    def test_strict_polynomial_signs_fail_on_roots(self):
        self.assertEqual(_strict_sign((Q(1),)),1)
        self.assertEqual(_strict_sign((Q(-1),)), -1)
        self.assertEqual(_strict_sign((Q(0),)),0)
        self.assertEqual(_strict_sign((Q(1),Q(0),Q(1))),1)
        self.assertIsNone(_strict_sign((Q(1),Q(-2),Q(1))))
        self.assertIsNone(_strict_sign((Q(-1),Q(2))))
        self.assertIsNone(_strict_sign((Q(0),Q(1))))

    def test_general_noncanonical_convex_input_is_accepted(self):
        P=tuple(pt(x,y) for x,y in ((0,0),(4,0),(3,2),(-1,3)))
        A,B,removed=compress_collinear_refinement(P,P)
        self.assertEqual(A,P)
        self.assertEqual(removed,())
        self.assertEqual(certify_strict_convex_motion(A,B),4)
        gap=certify_uniform_side_floor(A,B)
        self.assertGreater(gap["safe_side_squared_cutoff"],0)

    def test_nonconvex_or_self_intersecting_family_fails_closed(self):
        bad=tuple(pt(x,y) for x,y in (
            (0,0),(2,0),(1,Q(1,2)),(0,2)))
        with self.assertRaises((InconclusiveAtlas,ValueError)):
            audit_moving_convex_polygon(bad,bad)

    def test_unequal_counts_and_wrong_fraction_inputs_rejected(self):
        with self.assertRaises(ValueError):
            audit_moving_convex_polygon(self.before,self.after[:-1])
        with self.assertRaises(ValueError):
            audit_moving_convex_polygon(self.before,
                tuple(pt(v.x,v.y) for v in self.after[:-1]))

    def test_invalid_pair_separation_subdivision_depth_rejected(self):
        with self.assertRaises(ValueError):
            certify_uniform_side_floor(self.before,self.after,
                                       max_depth=-1)
        with self.assertRaises(ValueError):
            certify_uniform_side_floor(self.before,self.after,
                                       max_depth=1.5)

    def test_uncertified_size_cutoff_rejected_in_chart_classifier(self):
        with self.assertRaises(ValueError):
            classify_chart(self.before,self.after,
                (0,1,2,4),Q(0))
        with self.assertRaises(ValueError):
            classify_chart(self.before,self.after,
                (4,3,2,1),Q(1,100))

    def test_report_does_not_claim_global_jordan_proof(self):
        self.assertIn("No nonconvex/all-Jordan theorem",
                      self.pentagon["scope"])


if __name__=="__main__":
    unittest.main()
