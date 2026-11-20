"""v254: exact pair death, algebraic one-sided contacts and complete snapshots."""
import sys
from fractions import Fraction as F
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

from polygon_orbit_census_v223 import (
    pt,sub,cross,validate_jordan_polygon,census_polygon,
)
from big_peg_multichart_v249 import (
    pol,padd,psub,pmul,pneg,peval,ppoint_at,
    chart_polynomials,corner_numerators,
)
from big_peg_cutoff_ledger_v252 import _edge_parameter_inequalities
from moving_polygon_events_v224 import divide,iszero
from singular_continuation_v229 import root_sign,side_signs
from big_peg_pair_death_v254 import (
    START,END,ROOT_LEFT,ROOT_RIGHT,CUT,EVENT,QPOLY,
    P0,PERSISTENT,COLLIDING_LEFT,COLLIDING_RIGHT,CRITICAL,
    family_vertices,certify_jordan_family,_strict_sign,
    _contact_polynomials,certify_branch,certify_collision,
    _sample_snapshots,audit,
)


class ExactPositiveSizePairDeath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.P=family_vertices()
        cls.events=certify_collision(cls.P)
        cls.report=audit()

    def test_all_real_time_jordan_and_nonadjacent_edge_gap(self):
        r=certify_jordan_family(self.P)
        self.assertTrue(r["convex_Jordan_for_all_real_times"])
        self.assertGreater(r["global_nonadjacent_edge_rect_gap"],F(0))
        self.assertGreater(r["strict_turn_sign_checks"],0)

    def test_exact_cyclic_pentagon_sample_endpoints(self):
        self.assertEqual(len(self.P),5)
        for t in (START,ROOT_LEFT,ROOT_RIGHT,END):
            polygon=tuple(ppoint_at(v,t) for v in self.P)
            self.assertEqual(validate_jordan_polygon(polygon),polygon)

    def test_critical_event_quadratic_and_unique_root(self):
        self.assertEqual(QPOLY,pol(291143052416293,-496738335230000,72468080000000))
        self.assertEqual(side_signs(QPOLY,EVENT),(1,0,-1))
        self.assertGreater(peval(QPOLY,START),0)
        self.assertLess(peval(QPOLY,END),0)
        self.assertGreater(peval(QPOLY,ROOT_LEFT),0)
        self.assertLess(peval(QPOLY,ROOT_RIGHT),0)
        self.assertEqual(root_sign(QPOLY,EVENT),0)

    def test_event_root_is_algebraic_not_rounded_to_rational(self):
        self.assertFalse(EVENT.rational)
        self.assertEqual(EVENT.lo,ROOT_LEFT)
        self.assertEqual(EVENT.hi,ROOT_RIGHT)
        self.assertLess(ROOT_LEFT,ROOT_RIGHT)

    def test_both_collision_chart_determinants_nonzero_over_whole_interval(self):
        for ids in (COLLIDING_LEFT,COLLIDING_RIGHT,PERSISTENT):
            with self.subTest(ids=ids):
                D,N,_=_contact_polynomials(self.P,ids)
                sign=_strict_sign(D)
                self.assertEqual(sign,certify_branch(self.P,ids).determinant_sign)
                self.assertNotEqual(root_sign(D,EVENT),0)

    def test_both_colliding_squares_have_seven_strict_contacts(self):
        for ids in (COLLIDING_LEFT,COLLIDING_RIGHT):
            c=certify_branch(self.P,ids)
            self.assertEqual(c.other_strict_contact_checks,7)
            self.assertEqual(c.critical_contact_index,CRITICAL[ids])
            self.assertEqual(c.event_factor_sign,1)

    def test_critical_inequalities_share_exact_algebraic_factor(self):
        for ids in (COLLIDING_LEFT,COLLIDING_RIGHT):
            D,N,ineq=_contact_polynomials(self.P,ids)
            idx=CRITICAL[ids]
            q,r=divide(ineq[idx],QPOLY)
            self.assertTrue(iszero(r))
            self.assertEqual(_strict_sign(q),1)
            self.assertGreater(peval(ineq[idx],ROOT_LEFT),0)
            self.assertLess(peval(ineq[idx],ROOT_RIGHT),0)
            self.assertEqual(root_sign(ineq[idx],EVENT),0)

    def test_root_coordinates_agree_at_event_in_all_four_components(self):
        da,na,_=_contact_polynomials(self.P,COLLIDING_LEFT)
        db,nb,_=_contact_polynomials(self.P,COLLIDING_RIGHT)
        for a,b in zip(na,nb):
            q,r=divide(psub(pmul(a,db),pmul(b,da)),QPOLY)
            self.assertTrue(iszero(r))
            self.assertEqual(root_sign(psub(pmul(a,db),pmul(b,da)),EVENT),0)

    def test_two_square_branches_are_distinct_strictly_before_event(self):
        da,na,_=_contact_polynomials(self.P,COLLIDING_LEFT)
        db,nb,_=_contact_polynomials(self.P,COLLIDING_RIGHT)
        q,r=divide(psub(pmul(na[0],db),pmul(nb[0],da)),QPOLY)
        self.assertTrue(iszero(r))
        self.assertNotEqual(_strict_sign(q),0)

    def test_collision_square_second_corner_is_exact_pentagon_vertex(self):
        d,n,_=_contact_polynomials(self.P,COLLIDING_LEFT)
        x,y=corner_numerators(n,1)
        for num,coord in ((x,self.P[2][0]),(y,self.P[2][1])):
            q,r=divide(psub(num,pmul(d,coord)),QPOLY)
            self.assertTrue(iszero(r))

    def test_shared_square_side_squared_strictly_positive_at_event(self):
        self.assertEqual(certify_branch(self.P,COLLIDING_LEFT).minimum_side_squared_exceeds,CUT)
        self.assertEqual(certify_branch(self.P,COLLIDING_RIGHT).minimum_side_squared_exceeds,CUT)
        self.assertTrue(self.events["determinants_nonzero_at_event"])
        self.assertTrue(self.events["shared_second_corner_at_pentagon_vertex_2"])

    def test_persistent_third_square_stays_strictly_on_four_edges(self):
        c=certify_branch(self.P,PERSISTENT)
        self.assertEqual(c.other_strict_contact_checks,8)
        self.assertIsNone(c.critical_contact_index)
        self.assertIsNone(c.event_factor)
        self.assertGreater(c.minimum_side_squared_exceeds,0)

    def test_persistent_square_distinct_at_event(self):
        a,na,_=_contact_polynomials(self.P,COLLIDING_LEFT)
        c,nc,_=_contact_polynomials(self.P,PERSISTENT)
        diff=psub(pmul(na[0],c),pmul(nc[0],a))
        self.assertNotEqual(root_sign(diff,EVENT),0)

    def test_independent_full_positive_orbit_counts_both_sides(self):
        snapshot=self.report["total_positive_geometric_square_counts_at_rational_times"]
        self.assertEqual([v["all_positive_geometric_square_orbits"]
                          for v in snapshot],[3,3,1,1])
        self.assertEqual([v["time"] for v in snapshot],
                         list(map(str,(START,ROOT_LEFT,ROOT_RIGHT,END))))

    def test_unfiltered_census_checks_repeated_edges_and_orientations(self):
        for t,count in ((START,3),(END,1)):
            polygon=tuple(ppoint_at(p,t) for p in self.P)
            result=census_polygon(polygon)
            self.assertEqual(result.count,count)
            self.assertGreater(result.checked_edge_charts,100)
            self.assertTrue(all(s.side_squared>1 for s in result.squares))

    def test_negative_control_rejects_wrong_geometry(self):
        changed=list(self.P)
        changed[1]=(padd(changed[1][0],pol(F(1,100000))),changed[1][1])
        for checker in (certify_jordan_family,certify_collision):
            with self.assertRaises(ValueError):
                checker(tuple(changed))

    def test_negative_control_rejects_invalid_square_charts(self):
        for invalid in ((0,0,1,2),(0,1,2,3),(-1,1,3,4)):
            with self.assertRaises(ValueError):
                certify_branch(self.P,invalid)

    def test_event_is_not_rank_drop_or_size_cutoff(self):
        self.assertEqual(self.report["singularity_type"],
                         "two nonsingular Cramer charts glued at a polygon vertex; NO rank drop of their incidence matrices")
        self.assertIn("GEOMETRIC square pair death",self.report["status"])
        self.assertEqual(self.report["local_mod2_balance"],
                         "two geometric orbits disappear as a pair (2 mod 2 = 0)")

    def test_no_universal_or_full_time_total_orbit_claim(self):
        self.assertIn("NOT an all-time complete orbit census",self.report["scope"])
        self.assertIn("Toeplitz",self.report["scope"])


if __name__=="__main__":
    unittest.main()
