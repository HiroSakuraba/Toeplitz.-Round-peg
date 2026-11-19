"""v248 exact fixed-point persistence and adversarial negative controls."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import pt, add, sub, dot, validate_jordan_polygon
from big_peg_exact_v243 import read_polygon
from big_peg_persistence_v248 import (
    OCTAGON, OCTAGON_ASSIGNMENT, ORIGINAL_ASSIGNMENT,
    V247_ASSIGNMENT, CORNER_COEFFS, square_corner,
    invert_and_solve, certify_persistence, report, audit,
)


class ExactSquarePersistence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.octagon=certify_persistence(OCTAGON,OCTAGON_ASSIGNMENT,Q(1,512))
        cls.original,_,_=read_polygon(
            ROOT/"research"/"v244-pochuev-original-ten-gon.json")
        cls.v247,_,_=read_polygon(
            ROOT/"research"/"v247-refined-ten-gon.json")

    def test_robust_jordan_octagon(self):
        self.assertEqual(len(OCTAGON),8)
        self.assertEqual(validate_jordan_polygon(OCTAGON),OCTAGON)
        self.assertTrue(all(dot(v,v)<9 for v in OCTAGON))

    def test_exact_original_square_root(self):
        c=self.octagon
        self.assertEqual(c.root,(
            Q(21561,15136),Q(808131,529760),
            Q(-44719,15136),Q(-20205,105952)))
        self.assertEqual(len(set(c.assignment)),4)
        self.assertTrue(all(0<t<1 for t in c.parameters))

    def test_edge_parameters_replayed(self):
        c=self.octagon
        for k,i in enumerate(c.assignment):
            A=OCTAGON[i]
            E=sub(OCTAGON[(i+1)%8],A)
            p=square_corner(c.root,k)
            self.assertEqual(dot(E,sub(p,A))/dot(E,E),c.parameters[k])

    def test_uniform_rough_curve_side_bound_exact(self):
        c=self.octagon
        self.assertEqual(c.cube_radius,Q(367709,3390464))
        self.assertEqual(c.minimum_square_side_squared,
                         Q(46593827703305,5747623067648))
        self.assertGreater(c.minimum_square_side_squared,Q(27,10)**2)
        self.assertTrue(c.nondegenerate)

    def test_uniform_graph_image_cube_self_map(self):
        c=self.octagon
        self.assertLess(c.epsilon*c.normal_response_bound,c.cube_radius)
        self.assertTrue(all(
            c.cube_radius*l<min(t,1-t) for t,l in zip(
                c.parameters,c.parameter_sensitivities)))

    def test_the_restricted_big_peg_bound(self):
        # Every Jordan curve inside B_3 containing the four normal
        # graphs has r_in<=3 and an inscribed square side > 2.7.
        self.assertGreater(self.octagon.minimum_square_side_squared,Q(9,10)**2*9)

    def test_more_displacement_can_break_certificate(self):
        with self.assertRaises(ValueError):
            certify_persistence(OCTAGON,OCTAGON_ASSIGNMENT,Q(1,256))

    def test_invalid_or_float_epsilon_rejected(self):
        for bad in (Q(0),Q(-1,100),0.001,float(1/512),"1/512"):
            with self.subTest(epsilon=bad):
                with self.assertRaises(ValueError):
                    certify_persistence(OCTAGON,OCTAGON_ASSIGNMENT,bad)

    def test_bad_assignments_rejected(self):
        for indices in ((0,2,4),(0,0,4,6),(-1,2,4,6),(0,2,4,8)):
            with self.subTest(indices=indices):
                with self.assertRaises(ValueError):
                    certify_persistence(OCTAGON,indices,Q(1,512))

    def test_singular_exact_system_rejected(self):
        symmetric=tuple(pt(*xy) for xy in (
            (2,1),(1,2),(-1,2),(-2,1),
            (-2,-1),(-1,-2),(1,-2),(2,-1)))
        with self.assertRaises(ValueError):
            certify_persistence(symmetric,OCTAGON_ASSIGNMENT,Q(1,10000))

    def test_original_pochuev_stability(self):
        c=certify_persistence(self.original,ORIGINAL_ASSIGNMENT,Q(1,10**6))
        self.assertTrue(c.nondegenerate)
        self.assertEqual(c.assignment,ORIGINAL_ASSIGNMENT)
        self.assertLess(min(min(t,1-t) for t in c.parameters),Q(1,10**4))

    def test_v247_very_small_nonsmooth_stability(self):
        c=certify_persistence(self.v247,V247_ASSIGNMENT,Q(1,10**9))
        self.assertTrue(c.nondegenerate)
        self.assertLess(min(min(t,1-t) for t in c.parameters),Q(1,10**7))

    def test_v247_larger_perturbation_not_certified(self):
        with self.assertRaises(ValueError):
            certify_persistence(self.v247,V247_ASSIGNMENT,Q(1,10**8))

    def test_constant_normal_graphs_have_explicit_squares(self):
        # Choose four *constant* graph normal offsets with exact
        # h_i = epsilon * e_i.x. Since |e_x|<=||e|| this is admissible.
        # Solve new linear equations and check all four corners.
        c=self.octagon
        E=[];A=[];M=[];b=[];h=[]
        for k,i in enumerate(c.assignment):
            a=OCTAGON[i];e=sub(OCTAGON[(i+1)%8],a)
            n=pt(-e.y,e.x)
            cx,cy=CORNER_COEFFS[k]
            M.append(tuple(n.x*cx[j]+n.y*cy[j] for j in range(4)))
            b.append(dot(n,a))
            h.append(c.epsilon*e.x*(1 if k%2 else -1))
            A.append(a);E.append(e)
        inverse,root=invert_and_solve(M,b)
        z=tuple(root[j]+sum(inverse[j][i]*h[i] for i in range(4))
                for j in range(4))
        self.assertLessEqual(max(abs(z[j]-root[j]) for j in range(4)),
                             c.cube_radius)
        for k,(a,e,offset) in enumerate(zip(A,E,h)):
            p=square_corner(z,k)
            t=dot(e,sub(p,a))/dot(e,e)
            n=pt(-e.y,e.x)
            graph=add(add(a,pt(t*e.x,t*e.y)),
                      pt(offset*n.x/dot(e,e),offset*n.y/dot(e,e)))
            self.assertEqual(p,graph)
            self.assertTrue(0<t<1)
        self.assertGreaterEqual(z[2]**2+z[3]**2,
                                c.minimum_square_side_squared)

    def test_report_is_exact_and_scoped(self):
        data=audit()
        self.assertIn("LOCAL theorem only",data["scope"])
        self.assertIn("EVERY",data["robust_octagon"]["guarantee"])
        self.assertIn("9/10",data["restricted_rough_curve_class_big_peg_bound"])


if __name__=="__main__":
    unittest.main()
