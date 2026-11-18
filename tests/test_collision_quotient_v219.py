"""Exact v219 quotient / collision-link tests. No universal existence claims."""
import random
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from collision_quotient_v219 import (
    B, apply, norm_squared, conf, center, from_center, cyclic_shift,
    equivariant_retract, orbit, real_rotation_basis,
    equivariant_section_formal, section_equivariance_identity,
    leading_collision_residual, leading_link_lower_bound_squared,
    pt, is_simple_polygon, notch_polygon, notch_square,
    is_oriented_square, notch_square_side_squared, audit,
)


class QuotientCollisionTests(unittest.TestCase):
    def test_reject_invalid_gap_simplex(self):
        for gaps in ((0,0,0,0),(Q(1,2),Q(1,2),Q(1,2),-Q(1,2)),
                     (1,0,0,0,0)):
            with self.assertRaises(ValueError):conf(0,gaps)

    def test_closed_simplex_boundary_allowed(self):
        for k in range(4):
            gaps=tuple(Q(int(i==k)) for i in range(4))
            self.assertEqual(conf(0,gaps).gaps,gaps)

    def test_equivariant_center_exact_shift(self):
        rng=random.Random(219)
        for _ in range(100):
            v=[rng.randrange(0,21) for _ in range(4)]
            if sum(v)==0:v[0]=1
            g=tuple(Q(x,sum(v)) for x in v)
            x=conf(Q(rng.randrange(-100,100),17),g)
            self.assertEqual((center(cyclic_shift(x))-center(x))%1,Q(1,4))

    def test_center_inverse(self):
        x=conf(Q(-13,29),(Q(0),Q(1,2),Q(0),Q(1,2)))
        self.assertEqual(from_center(center(x),x.gaps),x)

    def test_no_stabilizers_even_at_collision_vertices(self):
        points=[conf(Q(1,7),tuple(Q(int(i==k)) for i in range(4)))
                for k in range(4)]
        points.append(conf(Q(5,13),(Q(1,4),)*4))
        points.append(conf(Q(2,9),(Q(0),Q(1,3),Q(1,3),Q(1,3))))
        for x in points:
            states=orbit(x)
            self.assertEqual(len(set(states)),4)
            self.assertEqual(cyclic_shift(states[-1]),x)
            self.assertEqual(tuple((center(y)-center(x))%1 for y in states),
                             (Q(0),Q(1,4),Q(1,2),Q(3,4)))

    def test_retract_is_deformation_retraction(self):
        x=conf(Q(1,7),(Q(1,8),Q(1,3),Q(1,6),Q(3,8)))
        self.assertEqual(equivariant_retract(x,0),x)
        y=equivariant_retract(x,1)
        self.assertEqual(y.gaps,(Q(1,4),)*4)
        self.assertEqual(center(y),center(x))
        self.assertEqual(equivariant_retract(y,Q(1,2)),y)

    def test_retract_equivariant_all_stages(self):
        rng=random.Random(17)
        for _ in range(50):
            raw=[rng.randrange(0,11) for _ in range(4)]
            if sum(raw)==0:raw[0]=1
            x=conf(Q(rng.randrange(100),97),tuple(Q(k,sum(raw)) for k in raw))
            s=Q(rng.randrange(0,20),19)
            self.assertEqual(equivariant_retract(cyclic_shift(x),s),
                             cyclic_shift(equivariant_retract(x,s)))

    def test_retract_reject_invalid_times(self):
        x=conf(0,(Q(1,4),)*4)
        for s in (Q(-1,5),Q(6,5)):
            with self.assertRaises(ValueError):equivariant_retract(x,s)

    def test_B_quarter_turn_both_orientations(self):
        for s in (-1,1):
            u,v=real_rotation_basis(s)
            self.assertEqual(apply(B(s),u),v)
            self.assertEqual(apply(B(s),v),tuple(-x for x in u))
            for k in range(4):
                e=tuple(Q(int(i==k)) for i in range(4))
                z=e
                for _ in range(4):z=apply(B(s),z)
                self.assertEqual(z,e)

    def test_section_exact_quarter_turn(self):
        for s in (-1,1):
            for cos,sin in ((Q(1),Q(0)),(Q(0),Q(1)),(Q(3,5),Q(4,5)),
                            (Q(-5,13),Q(12,13)),(Q(7,25),Q(-24,25))):
                self.assertTrue(section_equivariance_identity(s,cos,sin))
                z=equivariant_section_formal(s,cos,sin)
                self.assertEqual(norm_squared(z),Q(2))
                self.assertNotEqual(z,(Q(0),)*4)

    def test_section_symbolic_norm_relation_all_rationals(self):
        for sigma in (-1,1):
            for a in range(-4,5):
                for b in range(-4,5):
                    p=equivariant_section_formal(sigma,Q(a,7),Q(b,11))
                    self.assertEqual(norm_squared(p),
                                     2*(Q(a,7)**2+Q(b,11)**2))

    def test_collision_leading_term_orthogonal_lower_bound(self):
        rng=random.Random(778)
        for _ in range(120):
            a,b,c=[rng.randrange(0,15) for _ in range(3)]
            if a+b+c==0:a=1
            n=a+b+c
            aa,bb,cc=Q(a,n),Q(b,n),Q(c,n)
            vx,vy=Q(rng.randrange(-5,6),3),Q(rng.randrange(-7,8),4)
            if vx==vy==0:vx=Q(1)
            for s in (-1,1):
                val=leading_collision_residual(vx,vy,aa,bb,cc,s)
                self.assertEqual(val[2]**2+val[3]**2,
                                 (vx**2+vy**2)*(1+aa**2))
                self.assertGreaterEqual(val[2]**2+val[3]**2,
                                        leading_link_lower_bound_squared(vx,vy))
                self.assertGreater(norm_squared(val),0)

    def test_collision_first_order_is_exact_for_straight_line(self):
        v=(Q(5,3),Q(-7,4))
        for a,b,c in ((Q(0),Q(0),Q(1)),(Q(1),Q(0),Q(0)),
                      (Q(1,3),Q(1,2),Q(1,6))):
            eps=Q(1,131)
            p=((Q(0),Q(0)), (eps*a*v[0],eps*a*v[1]),
               (eps*(a+b)*v[0],eps*(a+b)*v[1]),(eps*v[0],eps*v[1]))
            for sigma in (-1,1):
                e=(p[1][0]-p[0][0],p[1][1]-p[0][1])
                f=(p[2][0]-p[1][0]+sigma*e[1],
                   p[2][1]-p[1][1]-sigma*e[0],
                   p[3][0]-p[0][0]+sigma*e[1],
                   p[3][1]-p[0][1]-sigma*e[0])
                expect=leading_collision_residual(*v,a,b,c,sigma)
                self.assertEqual(f,tuple(eps*x for x in expect))

    def test_collision_invalid_weights_and_orientation(self):
        with self.assertRaises(ValueError):leading_collision_residual(1,0,1,1,1)
        with self.assertRaises(ValueError):leading_collision_residual(1,0,1,0,0,2)

    def test_shrinking_notches_are_simple(self):
        for n in (2,3,5,10,17,128,10000):
            self.assertTrue(is_simple_polygon(notch_polygon(n)))

    def test_notch_gives_exact_nondegenerate_squares(self):
        for n in (2,3,5,10,17,128,10000):
            sq=notch_square(n)
            self.assertTrue(is_oriented_square(sq,1))
            self.assertFalse(is_oriented_square(sq,-1))
            self.assertEqual(notch_square_side_squared(n),Q(1,n*n))

    def test_notch_converges_to_rectangle_as_sets(self):
        for n in (2,5,10,100,10000):
            p=notch_polygon(n);e=Q(1,n)
            self.assertEqual(p[0],pt(0,0))
            self.assertEqual(p[-1],pt(0,1))
            self.assertEqual(p[2],pt(1,-e))
            self.assertEqual(p[3],pt(1+e,-e))
            # All extra points lie in [1,1+e] x [-e,0]; all other
            # curve segments lie on the original rectangle boundary.
            self.assertTrue(all(1<=z.x<=1+e and -e<=z.y<=0
                                for z in p[1:5]))
            self.assertTrue(all(0<=z.x<=2 and 0<=z.y<=1
                                for z in (p[0],p[5],p[6],p[7])))

    def test_self_crossing_polygon_is_not_simple(self):
        self.assertFalse(is_simple_polygon((pt(0,0),pt(2,2),pt(0,2),pt(2,0))))
        self.assertFalse(is_simple_polygon((pt(0,0),pt(1,0),pt(0,0),pt(0,1))))

    def test_audit(self):
        r=audit()
        self.assertTrue(r['equivariant_retraction'])
        self.assertTrue(r['section_equivariance_both_signs'])
        self.assertFalse(r['global_general_Jordan_square_peg_proved'])
        self.assertEqual(r['notch_side_squared_n64'],'1/4096')

    def test_independent_v218_residual_representation(self):
        try:
            from configuration_symmetry_v218 import (
                residual_cyclic_matrix, CYCLIC_CHART_JACOBIAN, determinant)
        except ImportError:
            self.skipTest('v218 implementation absent outside full checkout')
        for sigma in (-1,1):
            self.assertEqual(B(sigma),residual_cyclic_matrix(sigma))
            self.assertEqual(determinant(B(sigma)),1)
        self.assertEqual(determinant(CYCLIC_CHART_JACOBIAN),-1)

    def test_independent_v216_polygon_validator(self):
        try:
            from boundary_degree_v216 import pt as p216,validate_simple_polygon
        except ImportError:
            self.skipTest('v216 implementation absent outside full checkout')
        for n in (2,3,10,100):
            p=notch_polygon(n)
            validate_simple_polygon(tuple(p216(v.x,v.y) for v in p))


if __name__=='__main__':unittest.main()
