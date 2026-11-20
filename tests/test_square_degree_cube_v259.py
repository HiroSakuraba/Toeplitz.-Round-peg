"""v259: exact 4D Poincare-Miranda square degree across an actual kink."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt,add,sub,dot,cross,validate_jordan_polygon,census_polygon
from big_peg_persistence_v248 import square_corner
from big_peg_corner_v250 import family_vertices,prove_jordan_family
from big_peg_multichart_v249 import HALF,ppoint_at
from square_degree_cube_v259 import (
    RADIUS,EPSILON,ASYMMETRIC_EDGE_SHIFT,
    PiecewiseArc,square_defect,column_for_arc,
    _invert4,_precondition,_average_preconditioner,
    _face_vertices,_perturbation_response,
    certify_cube_degree,build_v250_kink_arcs,
    _linear_square_witness,audit,
)


class PiecewiseBoundaryDegree(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base,cls.arcs=build_v250_kink_arcs()
        cls.perturbed,cls.arcs2=build_v250_kink_arcs(shifted=True)
        cls.c=certify_cube_degree(cls.arcs,RADIUS,EPSILON)
        cls.c2=certify_cube_degree(cls.arcs2,RADIUS,EPSILON)

    def test_v250_reference_is_jordan_for_all_time(self):
        proof=prove_jordan_family(family_vertices())
        self.assertTrue(proof["strictly_convex_and_jordan_for_all_times"])
        self.assertEqual(validate_jordan_polygon(self.base),self.base)

    def test_reference_cube_crosses_actual_noncollinear_vertex(self):
        a0=self.arcs[0]
        self.assertEqual(a0.center,self.base[1])
        self.assertNotEqual(a0.negative_slope,a0.positive_slope)
        self.assertGreater(cross(a0.negative_slope,a0.positive_slope),0)
        self.assertEqual(a0.point(Q(0)),self.base[1])
        self.assertNotEqual(a0.point(-RADIUS),a0.point(RADIUS))

    def test_reference_center_four_points_form_an_exact_square(self):
        centers=[a.center for a in self.arcs]
        self.assertEqual(square_defect(centers),(Q(0),)*4)
        self.assertEqual(len(set(centers)),4)

    def test_actual_asymmetric_center_points_are_only_approximate(self):
        self.assertTrue(any(square_defect([a.center for a in self.arcs2])))
        self.assertEqual(
            sub(self.arcs2[1].center,self.arcs[1].center),
            ASYMMETRIC_EDGE_SHIFT)
        self.assertEqual(validate_jordan_polygon(self.perturbed),self.perturbed)

    def test_exact_reference_face_margin(self):
        self.assertEqual(self.c.minimum_face_margin,Q(17387,2696700))
        self.assertEqual(len(self.c.exact_face_margins),8)
        self.assertEqual(self.c.checked_face_vertices,216)
        self.assertTrue(all(x>0 for x in self.c.exact_face_margins))

    def test_asymmetric_face_margin_exact(self):
        self.assertEqual(self.c2.minimum_face_margin,
                         Q(11853176941,1785754740000))
        self.assertGreater(self.c2.minimum_face_margin,0)
        self.assertEqual(self.c2.checked_face_vertices,216)

    def test_operator_bound_for_arbitrary_coordinatewise_continuous_error(self):
        self.assertEqual(self.c.perturbation_response,Q(89890,3311))
        self.assertEqual(self.c.maximum_transformed_perturbation,
                         Q(8989,33110000))
        self.assertLess(self.c.maximum_transformed_perturbation,
                        self.c.minimum_face_margin)
        self.assertLess(self.c2.maximum_transformed_perturbation,
                        self.c2.minimum_face_margin)

    def test_degree_one_uses_no_single_four_edge_chart(self):
        self.assertEqual(self.c.local_degree,1)
        self.assertEqual(self.c2.local_degree,1)
        self.assertNotEqual(
            self.arcs[0].negative_slope,
            self.arcs[0].positive_slope)

    def test_preconditioner_inverse_is_genuine_four_dimensional(self):
        matrix=_average_preconditioner(self.arcs)
        inverse=_invert4(matrix)
        for row in range(4):
            for col in range(4):
                actual=sum(inverse[row][k]*matrix[k][col] for k in range(4))
                self.assertEqual(actual,Q(int(row==col)))

    def test_square_residual_vanishes_for_any_genuine_square(self):
        p=pt(Q(11,10),Q(7,10))
        z=(p.x,p.y,Q(2),Q(1))
        self.assertEqual(square_defect([square_corner(z,k) for k in range(4)]),
                         (Q(0),)*4)

    def test_face_grid_contains_the_kink_and_both_faces(self):
        samples=list(_face_vertices(1,-1,RADIUS))
        self.assertEqual(len(samples),27)
        self.assertEqual({s[0] for s in samples},{-RADIUS,Q(0),RADIUS})
        self.assertTrue(all(s[1]==-RADIUS for s in samples))
        for axis in range(4):
            for sign in (-1,1):
                self.assertEqual(len(list(_face_vertices(axis,sign,RADIUS))),27)

    def test_rational_side_floor_stronger_than_nine(self):
        self.assertGreater(self.c.certified_square_side_squared_lower,Q(9))
        self.assertGreater(self.c2.certified_square_side_squared_lower,Q(9))

    def test_exact_witness_from_piecewise_square_constraints(self):
        coordinates,pts=_linear_square_witness(self.arcs2,self.c2)
        self.assertEqual(square_defect(pts),(Q(0),)*4)
        self.assertTrue(all(-RADIUS<t<RADIUS for t in coordinates))
        self.assertGreater(coordinates[0],0)
        self.assertEqual(len(set(pts)),4)

    def test_independent_all_square_census_contains_actual_witness(self):
        _,pts=_linear_square_witness(self.arcs2,self.c2)
        census=census_polygon(self.perturbed)
        self.assertTrue(any(set(pts)==set(s.vertices) for s in census.squares))
        self.assertGreater(census.count,0)

    def test_base_independent_piecewise_witness_is_exact(self):
        parameters,pts=_linear_square_witness(self.arcs,self.c)
        self.assertEqual(square_defect(pts),(Q(0),)*4)
        self.assertTrue(all(-RADIUS<t<RADIUS for t in parameters))

    def test_wrong_small_radius_fails_perturbation_budget(self):
        with self.assertRaises(ValueError):
            certify_cube_degree(self.arcs,Q(1,10000),EPSILON)

    def test_excessive_arbitrary_continuous_error_fails_closed(self):
        with self.assertRaises(ValueError):
            certify_cube_degree(self.arcs,RADIUS,Q(1,1000))

    def test_zero_perturbation_admitted(self):
        z=certify_cube_degree(self.arcs,RADIUS,Q(0))
        self.assertEqual(z.local_degree,1)
        self.assertEqual(z.maximum_transformed_perturbation,0)

    def test_singular_preconditioner_cannot_be_inverted(self):
        singular=((Q(1),Q(0),Q(0),Q(0)),)*4
        with self.assertRaises(ValueError):
            certify_cube_degree(self.arcs,RADIUS,Q(0),preconditioner=singular)

    def test_rank_deficient_piecewise_chart_is_rejected_if_pm_faces_fail(self):
        # Four parallel displacement directions give a singular average
        # Jacobian. We do NOT infer any existence theorem just because
        # the four central points happen to form a square.
        diamond=[pt(1,0),pt(0,1),pt(-1,0),pt(0,-1)]
        arcs=tuple(PiecewiseArc(p,pt(1,0),pt(1,0)) for p in diamond)
        with self.assertRaises(ValueError):
            certify_cube_degree(arcs,RADIUS,Q(0))

    def test_shifted_near_square_fails_if_far_outside_cube(self):
        arcs=list(self.arcs)
        original=arcs[1]
        arcs[1]=PiecewiseArc(add(original.center,pt(1,1)),
                             original.negative_slope,original.positive_slope)
        with self.assertRaisesRegex(ValueError,"face sign"):
            certify_cube_degree(tuple(arcs),RADIUS,Q(0))

    def test_degenerate_reference_side_fails_positive_side_floor(self):
        # Artificial input: if face signs fail there is no false
        # positive side guarantee even when center points collapse.
        arcs=tuple(PiecewiseArc(pt(0,0),a.negative_slope,a.positive_slope)
                   for a in self.arcs)
        with self.assertRaises(ValueError):
            certify_cube_degree(arcs,RADIUS,EPSILON)

    def test_invalid_float_radius_and_error_rejected(self):
        with self.assertRaises(ValueError):
            certify_cube_degree(self.arcs,0.01,EPSILON)
        with self.assertRaises(ValueError):
            certify_cube_degree(self.arcs,RADIUS,0.00001)
        with self.assertRaises(ValueError):
            certify_cube_degree(self.arcs,Q(0),EPSILON)
        with self.assertRaises(ValueError):
            certify_cube_degree(self.arcs,RADIUS,Q(-1,1000))

    def test_invalid_arc_parameter_type_rejected(self):
        with self.assertRaises(ValueError):
            self.arcs[0].point(0.1)

    def test_not_a_global_jordan_degree_or_toeplitz_result(self):
        report=audit()
        self.assertEqual(report["local_degree"],1)
        self.assertIn("not a global square-orbit invariant",report["scope"])
        self.assertIn("No universal positive Big Peg",report["scope"])
        self.assertIn("Toeplitz proof",report["scope"])


if __name__=="__main__":
    unittest.main()
