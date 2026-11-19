"""v249: exact time-cover certificate through a moving polygon vertex."""
import sys
from fractions import Fraction as Q
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from big_peg_multichart_v249 import (
    Z,O,HALF,ALPHA,S2_LIMIT,OCTAGON,
    pol,padd,psub,pmul,peval,bcoeff,certify_sign,
    family_vertices,ppoint_at,chart_polynomials,root_at,params_at,
    certify_chart,prove_jordan_family,certify_cover,audit,
)


class PolynomialExactCertificates(unittest.TestCase):
    def test_power_to_bernstein_identity(self):
        p=pol(Q(1,3),Q(-3,4),Q(1,2))
        b=bcoeff(p,Q(1,4),Q(3,4))
        self.assertEqual(peval(p,Q(1,4)),b[0])
        self.assertEqual(peval(p,Q(3,4)),b[-1])

    def test_subdivision_proves_nonnegative_quadratic(self):
        p=pol(Q(1,4),Q(-1),Q(1))  # (t - 1/2)^2
        self.assertGreater(certify_sign(p,Q(0),Q(1)),1)
        with self.assertRaises(ValueError):
            certify_sign(p,Q(0),Q(1),strict=True,max_depth=8)

    def test_negative_polynomial_fails_closed(self):
        with self.assertRaises(ValueError):
            certify_sign(pol(1,-2),Q(0),Q(1),max_depth=5)

    def test_invalid_polynomial_intervals(self):
        for a,b in ((0,0),(Q(1,2),Q(1,2)),(Q(-1),0),(0,2),(Q(3,4),Q(1,4))):
            with self.subTest(a=a,b=b):
                with self.assertRaises(ValueError):
                    certify_sign(O,Q(a),Q(b))


class ChangingJordanFamily(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.vertices=family_vertices()
        cls.report=audit()
        cls.atlas=certify_cover(cls.vertices,
                                ((Q(0),HALF,(1,3,5,7)),
                                 (HALF,Q(1),(0,3,5,7))))

    def test_family_has_nine_vertices_and_redundant_split(self):
        self.assertEqual(len(self.vertices),9)
        self.assertEqual(ppoint_at(self.vertices[0],Q(0)),OCTAGON[0])
        self.assertEqual(ppoint_at(self.vertices[2],Q(1)),OCTAGON[1])
        self.assertNotEqual(ppoint_at(self.vertices[4],Q(0)),
                            ppoint_at(self.vertices[4],Q(1)))

    def test_every_time_jordan_certificate(self):
        proof=prove_jordan_family(self.vertices)
        self.assertTrue(proof["all_polygons_jordan"])
        self.assertTrue(proof["all_polygons_inside_radius_3"])
        self.assertEqual(proof["minimal_nonadjacent_coordinate_separation"],"9/10")
        self.assertGreater(proof["turn_and_edge_sign_certificates"],0)

    def test_unverified_vertex_motion_cannot_pass_jordan_cover(self):
        # The global separation proof requires a verified L-infinity
        # displacement envelope for EVERY supplied moving vertex.
        tampered=list(self.vertices)
        vertex=tampered[4]
        tampered[4]=(padd(vertex[0],pol(Q(1,10))),vertex[1])
        with self.assertRaises(ValueError):
            prove_jordan_family(tuple(tampered))

    def test_exact_jordan_sample_regressions(self):
        for t in (Q(0),Q(1,10),Q(1,2),Q(2,3),Q(1)):
            with self.subTest(t=t):
                polygon=tuple(ppoint_at(v,t) for v in self.vertices)
                self.assertEqual(validate_jordan_polygon(polygon),polygon)

    def test_both_subcharts_have_strict_square_size_lower_bound(self):
        self.assertEqual(len(self.atlas),2)
        for c in self.atlas:
            self.assertEqual(len(c.assignment),4)
            self.assertGreater(c.sign_checks,0)
            for t,pars,root in ((c.left,c.left_parameters,c.left_root),
                                (c.right,c.right_parameters,c.right_root)):
                self.assertTrue(all(0<=v<=1 for v in pars))
                self.assertGreater(root[2]**2+root[3]**2,S2_LIMIT)

    def test_joint_root_matches_exactly(self):
        low,high=self.atlas
        self.assertEqual(low.right,high.left)
        self.assertEqual(low.right_root,high.left_root)
        self.assertEqual(low.right_parameters[0],Q(0))
        self.assertEqual(high.left_parameters[0],Q(1))
        self.assertEqual(pt(low.right_root[0],low.right_root[1]),
                         ppoint_at(self.vertices[1],HALF))

    def test_switch_is_essential_to_these_edge_charts(self):
        # At t=0 the first-edge version is outside its segment.
        D,N=chart_polynomials(self.vertices,(0,3,5,7))
        z=root_at(D,N,Q(0))
        self.assertGreater(params_at(self.vertices,(0,3,5,7),z,Q(0))[0],1)
        # At t=1 the second-edge version is outside its segment.
        D,N=chart_polynomials(self.vertices,(1,3,5,7))
        z=root_at(D,N,Q(1))
        self.assertLess(params_at(self.vertices,(1,3,5,7),z,Q(1))[0],0)

    def test_wrong_chart_full_half_interval_fails_closed(self):
        for wrong in ((Q(0),HALF,(0,3,5,7)),
                      (HALF,Q(1),(1,3,5,7))):
            with self.subTest(wrong=wrong):
                with self.assertRaises(ValueError):
                    certify_chart(self.vertices,wrong[2],wrong[0],wrong[1])

    def test_symbolic_reparametrization_of_same_square_branch(self):
        low,high=self.atlas
        for a,b in zip(low.numerators,high.numerators):
            self.assertEqual(pmul(a,high.denominator),
                             pmul(b,low.denominator))

    def test_cover_with_gap_is_rejected(self):
        with self.assertRaises(ValueError):
            certify_cover(self.vertices,(
                (Q(0),Q(2,5),(1,3,5,7)),
                (HALF,Q(1),(0,3,5,7))))

    def test_cover_missing_endpoint_is_rejected(self):
        with self.assertRaises(ValueError):
            certify_cover(self.vertices,((Q(1,10),HALF,(1,3,5,7)),
                                          (HALF,Q(1),(0,3,5,7))))

    def test_large_square_cutoff_is_not_certified(self):
        with self.assertRaises(ValueError):
            certify_cover(self.vertices,((Q(0),HALF,(1,3,5,7)),
                                          (HALF,Q(1),(0,3,5,7))),
                          side2_limit=Q(10))

    def test_negative_cutoff_not_accepted(self):
        with self.assertRaises(ValueError):
            certify_chart(self.vertices,(1,3,5,7),Q(0),HALF,
                          side2_limit=Q(0))

    def test_rational_event_location_and_monotone_split(self):
        self.assertEqual(ALPHA,Q(39906,73777))
        low,high=self.atlas
        self.assertEqual(low.right_root,high.left_root)
        q0=ALPHA-Q(1,4)
        q1=ALPHA+Q(1,4)
        self.assertTrue(0<q0<ALPHA<q1<1)

    def test_restricted_big_peg_bound_proved(self):
        self.assertEqual(S2_LIMIT,Q(3)**2*Q(9,10)**2) # (2.7)^2 = 9*(0.9)^2
        self.assertEqual(S2_LIMIT,Q(27,10)**2)
        self.assertIn("9/10",self.report["restricted_family_big_peg_ratio_strictly_above"])

    def test_scope_disallows_general_continuation_claim(self):
        self.assertTrue(self.report["handoff_corner_on_moving_vertex"])
        self.assertTrue(self.report["geometric_square_branch_symbolically_identical_across_charts"])
        self.assertIn("not general singular-orbit continuation",self.report["scope"])


if __name__=="__main__":
    unittest.main()
