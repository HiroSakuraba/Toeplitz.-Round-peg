"""v256-v257: exact nonconvex central-square and near-symmetric boundary witnesses."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import (
    pt,add,sub,dot,cross,_on_seg,_square_residual,
    validate_jordan_polygon
)
from central_symmetry_square_v256 import (
    J,minus,point_location,distance_squared_point_segment,
    disk_certificate,_segment_intersections,quarter_turn_pairs,
    is_centrally_symmetric,certify_centered_square,
    certify_near_square,NOTCHED_12,audit,
)


class QuarterTurnJordanCertificates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=NOTCHED_12
        cls.Q=tuple((add(p,pt(Q(1,20),Q(1,40)))
                    if i==3 else p)
                    for i,p in enumerate(NOTCHED_12))

    def test_polygon_is_jordan_with_reflex_corners(self):
        self.assertEqual(validate_jordan_polygon(self.P),self.P)
        turns=[cross(sub(self.P[i],self.P[(i-1)%12]),
                     sub(self.P[(i+1)%12],self.P[i]))
               for i in range(12)]
        self.assertTrue(any(x<0 for x in turns))
        self.assertTrue(any(x>0 for x in turns))

    def test_nonstar_shape_about_origin(self):
        self.assertEqual(point_location(self.P,pt(0,0)),1)
        self.assertEqual(point_location(self.P,pt(4,4)),0)
        self.assertEqual(point_location(self.P,pt(Q(5,2),Q(5,2))),-1)

    def test_exact_centered_disk_inside_nonconvex_polygon(self):
        d2=disk_certificate(self.P,Q(4))
        self.assertGreater(d2,Q(4))
        self.assertEqual(point_location(self.P,pt(0,0)),1)

    def test_source_is_centrally_symmetric_without_ninety_degree_symmetry(self):
        self.assertTrue(is_centrally_symmetric(self.P))
        self.assertEqual(len(self.P),12)
        for i in range(6):
            self.assertEqual(self.P[i+6],minus(self.P[i]))
        self.assertNotEqual(tuple(self.P),tuple(J(p) for p in self.P))

    def test_turning_origin_square_corners_are_exact(self):
        c=certify_centered_square(self.P,Q(4))
        self.assertTrue(c.polygon_centrally_symmetric)
        self.assertGreaterEqual(c.side_squared,Q(8))
        self.assertEqual(_square_residual(c.vertices,1),(0,0,0,0))
        self.assertEqual(len(set(c.vertices)),4)
        self.assertEqual(c.vertices[2],minus(c.vertices[0]))
        self.assertEqual(c.vertices[1],J(c.vertices[0]))

    def test_four_square_vertices_all_on_actual_boundary(self):
        from central_symmetry_square_v256 import _on_boundary
        cert=certify_centered_square(self.P,Q(4))
        self.assertTrue(all(_on_boundary(self.P,p) for p in cert.vertices))

    def test_every_quarter_turn_intersection_is_valid(self):
        pairs=quarter_turn_pairs(self.P)
        self.assertGreater(len(pairs),0)
        for w in pairs:
            self.assertEqual(w.rotated_source,minus(J(w.point)))
            self.assertGreaterEqual(w.first_parameter,0)
            self.assertLessEqual(w.first_parameter,1)
            self.assertGreaterEqual(w.other_parameter,0)
            self.assertLessEqual(w.other_parameter,1)
            self.assertTrue(_on_seg(
                w.point,self.P[w.edge_index],
                self.P[(w.edge_index+1)%12]))
            self.assertTrue(_on_seg(
                w.rotated_source,self.P[w.other_edge_index],
                self.P[(w.other_edge_index+1)%12]))

    def test_segment_crossings_and_overlaps_are_exact(self):
        a,b,c,d=pt(0,0),pt(2,2),pt(0,2),pt(2,0)
        x=_segment_intersections(a,b,c,d)
        self.assertEqual(x,((pt(1,1),Q(1,2),Q(1,2)),))
        overlaps=_segment_intersections(pt(0,0),pt(4,0),pt(2,0),pt(5,0))
        self.assertEqual({q[0] for q in overlaps},{pt(2,0),pt(4,0)})
        self.assertEqual(_segment_intersections(pt(0,0),pt(1,0),pt(0,1),pt(1,1)),())

    def test_disk_rejects_missing_origin_or_too_large_radius(self):
        with self.assertRaises(ValueError):
            disk_certificate(self.P,Q(25))
        shifted=tuple(add(p,pt(10,0)) for p in self.P)
        with self.assertRaises(ValueError):
            disk_certificate(shifted,Q(1))

    def test_central_square_rejects_asymmetric_shape(self):
        with self.assertRaises(ValueError):
            certify_centered_square(self.Q,Q(4))

    def test_boundary_rotation_rejects_origin_outside(self):
        moved=tuple(add(p,pt(10,0)) for p in self.P)
        with self.assertRaises(ValueError):
            quarter_turn_pairs(moved)


class NearSymmetryTheorem(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=NOTCHED_12
        cls.Q=tuple(add(p,pt(Q(1,20),Q(1,40))) if i==3 else p
                    for i,p in enumerate(cls.P))
        cls.result=certify_near_square(cls.Q,Q(4))

    def test_perturbed_polygon_is_valid_and_non_symmetric(self):
        self.assertEqual(validate_jordan_polygon(self.Q),self.Q)
        self.assertFalse(is_centrally_symmetric(self.Q))
        self.assertGreater(disk_certificate(self.Q,Q(4)),Q(4))

    def test_exact_vertex_symmetry_defect(self):
        self.assertEqual(self.result.paired_vertex_defect_squared,Q(1,320))
        self.assertLess(self.result.paired_vertex_defect_squared,Q(4)/4)

    def test_near_square_has_four_distinct_actual_boundary_points(self):
        from central_symmetry_square_v256 import _on_boundary
        pts=self.result.exact_boundary_points
        self.assertEqual(len(pts),4)
        self.assertEqual(len(set(pts)),4)
        self.assertTrue(all(_on_boundary(self.Q,p) for p in pts))

    def test_near_square_is_an_ideal_square_plus_certified_corner_errors(self):
        exact=self.result.ideal_centered_square
        approximate=self.result.exact_boundary_points
        self.assertEqual(_square_residual(exact,1),(0,0,0,0))
        self.assertEqual(exact[0],approximate[0])
        self.assertEqual(exact[3],approximate[3])
        for a,b in zip(exact,approximate):
            d=sub(a,b)
            self.assertLessEqual(dot(d,d),Q(1,320))

    def test_ideal_square_size_is_large_relative_to_centered_disk(self):
        self.assertGreaterEqual(self.result.ideal_side_squared,Q(8))
        self.assertGreater(self.result.quarter_turn_intersection_count,0)
        self.assertLessEqual(self.result.max_corner_error_squared,
                             self.result.paired_vertex_defect_squared)

    def test_no_false_claim_of_exact_square_for_asymmetric_shape(self):
        self.assertFalse(self.result.exact_square_claimed)
        self.assertNotEqual(self.result.paired_vertex_defect_squared,0)

    def test_method_rejects_defect_too_large_relative_to_declared_disk(self):
        # Same perfectly good polygon but a very small claimed disk
        # cannot prove separation of four approximate corners.
        with self.assertRaises(ValueError):
            certify_near_square(self.Q,Q(1,100))

    def test_near_square_rejects_odd_vertex_count(self):
        # Delete a non-essential notch vertex, keep a Jordan polygon.
        P=self.P[:3]+self.P[4:]
        self.assertEqual(len(P),11)
        validate_jordan_polygon(P)
        with self.assertRaises(ValueError):
            certify_near_square(P,Q(1))

    def test_symmetric_limit_has_zero_corner_error(self):
        res=certify_near_square(self.P,Q(4))
        self.assertEqual(res.paired_vertex_defect_squared,Q(0))
        self.assertEqual(res.max_corner_error_squared,Q(0))
        self.assertEqual(res.exact_boundary_points,res.ideal_centered_square)
        self.assertFalse(res.exact_square_claimed) # generic API does not claim it

    def test_report_scope_and_numerical_values(self):
        data=audit()
        self.assertTrue(data["central_example_nonconvex_and_not_star_shaped_about_centre"])
        self.assertEqual(data["near_example_max_antipodal_vertex_mismatch_squared"],"1/320")
        self.assertFalse(data["near_example_exact_square_claimed"])
        self.assertEqual(data["unrestricted_Toeplitz"],"OPEN")


if __name__=="__main__":
    unittest.main()
