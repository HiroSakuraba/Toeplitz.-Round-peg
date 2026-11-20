"""v265 exact irrational vertex-gluing / all-real-time singular census tests."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt, add, mul
from big_peg_pair_death_v254 import (
    family_vertices,ROOT_LEFT,ROOT_RIGHT,QPOLY,
    COLLIDING_LEFT,COLLIDING_RIGHT,
)
from big_peg_multichart_v249 import ppoint_at,peval
from moving_polygon_events_v224 import divide,iszero
from square_relative_degree_v260 import (
    BOUNDARY_RESIDUAL_FLOOR as V260_FLOOR,
    _normalized_event_polynomial,
)
from convex_square_event_atlas_v264 import InconclusiveAtlas, audit_moving_convex_polygon
from square_irrational_vertex_gluing_v265 import (
    UnresolvedVertexEvent,_positive_min,_strict_sign,
    certify_vertex_pair,complete_v264_vertex_event,
)


def original_motion():
    P=family_vertices()
    return tuple(ppoint_at(v,ROOT_LEFT) for v in P),tuple(
        ppoint_at(v,ROOT_RIGHT) for v in P)


def add_motion(A,B,v_start,v_end):
    return (tuple(add(p,v_start) for p in A),
            tuple(add(p,v_end) for p in B))


def scale_motion(A,B,factor):
    return (tuple(mul(p,factor) for p in A),
            tuple(mul(p,factor) for p in B))


class ExactIrrationalVertexGluing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.A,cls.B=original_motion()
        cls.pair=certify_vertex_pair(
            cls.A,cls.B,COLLIDING_LEFT,COLLIDING_RIGHT,1,
            Q(51780781,400000000))
        cls.all=complete_v264_vertex_event(
            cls.A,cls.B,COLLIDING_LEFT,COLLIDING_RIGHT,1)

    def test_common_contact_polynomial_exactly_same_as_v254(self):
        cert=self.pair
        q=cert.exact_event_polynomial
        w=_normalized_event_polynomial()
        _,r=divide(q,w)
        self.assertTrue(iszero(r))
        _,r2=divide(w,q)
        self.assertTrue(iszero(r2))
        self.assertLess(peval(q,Q(0))*peval(q,Q(1)),0)
        self.assertEqual(cert.exact_divisibilities_checked,5)

    def test_exactly_one_strictly_interior_irrational_event(self):
        lo,hi=self.pair.event_bracket
        self.assertLess(Q(0),lo)
        self.assertLess(lo,hi)
        self.assertLess(hi,Q(1))
        self.assertLess(hi-lo,Q(1,2**24)+Q(1,2**24))
        actual_lo=ROOT_LEFT+(ROOT_RIGHT-ROOT_LEFT)*lo
        actual_hi=ROOT_LEFT+(ROOT_RIGHT-ROOT_LEFT)*hi
        self.assertLess(peval(QPOLY,actual_lo)*peval(QPOLY,actual_hi),0)

    def test_two_charts_share_exact_same_vertex_and_other_corners(self):
        c=self.pair
        self.assertEqual(c.left_edges,COLLIDING_LEFT)
        self.assertEqual(c.right_edges,COLLIDING_RIGHT)
        self.assertEqual(c.corner_index,1)
        self.assertEqual(c.shared_vertex_index,2)
        self.assertEqual(c.exact_divisibilities_checked,5)

    def test_local_indices_cancel(self):
        c=self.pair
        self.assertEqual((c.left_orientation,c.right_orientation),(-1,1))
        self.assertEqual(c.left_orientation+c.right_orientation,0)
        self.assertEqual(self.all["glued_pair_relative_degree"],0)
        self.assertEqual(self.all["pair_signed_regular_indices"],[-1,1])

    def test_uniform_extreme_4d_outer_faces_zero_free(self):
        self.assertGreater(self.pair.source_root_face_gap,0)
        self.assertGreater(self.pair.external_defect_floor,0)
        self.assertEqual(self.pair.exact_strict_face_signs_checked,16)
        self.assertGreater(Q(self.all[
            "uniform_outer_boundary_square_defect_floor"]),0)
        self.assertGreater(V260_FLOOR,0)

    def test_whole_pentagon_completed_at_irrational_event(self):
        r=self.all
        self.assertEqual(r["status"],
            "PASS: complete three-slice whole-polygon positive-square census from exact generic irrational vertex gluing")
        self.assertEqual(
            [r["event_before_geometric_squares"],
             r["event_exact_geometric_squares"],
             r["event_after_geometric_squares"]],[3,2,1])
        self.assertTrue(r["complete_whole_polygon_all_time_ledger"])
        self.assertEqual(r["relative_unlabeled_parity_before"],1)
        self.assertEqual(r["relative_unlabeled_parity_after"],1)
        self.assertFalse(r["unrestricted_toeplitz_proved"])

    def test_complete_nonsingular_chamber_counts(self):
        self.assertEqual(self.all["all_real_time_pre_event_count"],3)
        self.assertEqual(self.all["all_real_time_post_event_count"],1)
        self.assertGreater(Q(self.all[
            "all_real_time_uniform_positive_square_side_floor"]),0)

    def test_affine_translation_in_time_preserves_event_certificate(self):
        A,B=add_motion(self.A,self.B,pt(Q(1,7),Q(-1,9)),
                       pt(Q(7,10),Q(1,3)))
        c=certify_vertex_pair(A,B,COLLIDING_LEFT,COLLIDING_RIGHT,1,
                              Q(51780781,400000000))
        self.assertEqual(c.exact_event_polynomial,
                         self.pair.exact_event_polynomial)
        self.assertEqual(c.event_bracket,self.pair.event_bracket)
        self.assertEqual((c.left_orientation,c.right_orientation),(-1,1))
        self.assertGreater(c.external_defect_floor,0)

    def test_independent_positively_scaled_motion(self):
        A,B=scale_motion(self.A,self.B,Q(3,2))
        c=certify_vertex_pair(A,B,COLLIDING_LEFT,COLLIDING_RIGHT,1)
        self.assertEqual(c.exact_event_polynomial,
                         self.pair.exact_event_polynomial)
        self.assertEqual((c.left_orientation,c.right_orientation),(-1,1))
        self.assertGreater(c.external_defect_floor,0)

    def test_wrong_vertex_and_edge_alignment_rejected(self):
        with self.assertRaises(ValueError):
            certify_vertex_pair(self.A,self.B,
                COLLIDING_LEFT,COLLIDING_RIGHT,2)
        with self.assertRaises(ValueError):
            certify_vertex_pair(self.A,self.B,
                COLLIDING_LEFT,COLLIDING_LEFT,1)
        with self.assertRaises(ValueError):
            certify_vertex_pair(self.A,self.B,
                COLLIDING_RIGHT,COLLIDING_LEFT,1)

    def test_same_geometry_reverse_motion_is_a_birth_not_death(self):
        with self.assertRaises(UnresolvedVertexEvent):
            certify_vertex_pair(self.B,self.A,
                COLLIDING_LEFT,COLLIDING_RIGHT,1)

    def test_uncertified_event_on_static_convex_quadrilateral(self):
        P=tuple(pt(x,y) for x,y in (
            (0,0),(4,0),(3,2),(-1,3)))
        with self.assertRaises((UnresolvedVertexEvent,ValueError)):
            certify_vertex_pair(P,P,(0,1,2,3),(0,2,2,3),1)

    def test_nonconvex_input_fails_closed(self):
        p=tuple(pt(x,y) for x,y in (
            (0,0),(2,0),(1,Q(1,2)),(0,2)))
        with self.assertRaises((InconclusiveAtlas,ValueError)):
            certify_vertex_pair(p,p,(0,1,2,3),(0,2,2,3),1)

    def test_bad_cutoff_type_rejected(self):
        for x in (Q(-1),Q(-1,2),0.05,None):
            with self.subTest(x=x):
                with self.assertRaises(ValueError):
                    certify_vertex_pair(
                        self.A,self.B,COLLIDING_LEFT,
                        COLLIDING_RIGHT,1,x)

    def test_broken_event_polynomial_cannot_be_called_positive(self):
        self.assertEqual(_strict_sign((Q(1),)),1)
        self.assertEqual(_strict_sign((Q(-1),)),-1)
        self.assertIsNone(_strict_sign((Q(-1,2),Q(1))))
        with self.assertRaises(UnresolvedVertexEvent):
            _positive_min((Q(-1,2),Q(1)))

    def test_v264_third_party_vertex_locked_root_must_be_rejected(self):
        atlas=audit_moving_convex_polygon(self.A,self.B)
        ledger=dict(atlas)
        ledger["chart_classification"]=dict(atlas["chart_classification"])
        ledger["chart_classification"]["stable_boundary"]=1
        with patch(
            "square_irrational_vertex_gluing_v265.audit_moving_convex_polygon",
            return_value=ledger,
        ):
            with self.assertRaises(UnresolvedVertexEvent):
                complete_v264_vertex_event(
                    self.A,self.B,COLLIDING_LEFT,COLLIDING_RIGHT,1)

    def test_remaining_global_root_matches_stable_interior_chart(self):
        self.assertEqual(self.all["stable_nonparticipating_square_orbits"],1)
        self.assertEqual(self.all["regular_oriented_degree_in_both_chambers"],1)

    def test_perturbed_polygon_does_not_silently_inherit_exact_event(self):
        altered=list(self.B)
        v=altered[4]
        altered[4]=pt(v.x+Q(1,1000),v.y)
        with self.assertRaises(UnresolvedVertexEvent):
            complete_v264_vertex_event(
                self.A,tuple(altered),
                COLLIDING_LEFT,COLLIDING_RIGHT,1)


if __name__=="__main__":
    unittest.main()
