"""v218 exact finite checks of cyclic symmetry and collision strata."""
import random
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from configuration_symmetry_v218 import (
    pt, square_residual, residual_cyclic_matrix, matmul, determinant,
    CYCLIC_CHART_JACOBIAN, positive_gap_chart, cyclic_chart_shift,
    orbit_local_signs, order_sectors, graceful_sector, closed_gap_chart,
    total_collision_vertex, chart_points, is_nondegenerate_square,
    isolated_distinct_edge_squares, run_audit,
)


class ConfigurationSymmetryTests(unittest.TestCase):
    def test_four_dimensional_domain_orientation_reversal(self):
        self.assertEqual(determinant(CYCLIC_CHART_JACOBIAN),-1)

    def test_target_orientation_preserved_both_orientations(self):
        for sigma in (1,-1):
            self.assertEqual(determinant(residual_cyclic_matrix(sigma)),1)

    def test_cyclic_residual_identity_random_exact(self):
        rng=random.Random(218)
        for sigma in (-1,1):
            B=residual_cyclic_matrix(sigma)
            for _ in range(100):
                p=tuple(pt(Q(rng.randint(-5,5),3),Q(rng.randint(-7,7),4))
                        for i in range(4))
                self.assertEqual(square_residual(p[1:]+p[:1],sigma),
                                 matmul(B,square_residual(p,sigma)))

    def test_four_turns_act_as_identity_on_target(self):
        for sigma in (-1,1):
            B=residual_cyclic_matrix(sigma)
            for e in range(4):
                v=tuple(Q(j==e) for j in range(4))
                w=v
                for _ in range(4):w=matmul(B,w)
                self.assertEqual(w,v)

    def test_four_turns_act_as_identity_on_gap_chart(self):
        c=positive_gap_chart(Q(7,8),Q(1,8),Q(3,16),Q(1,4))
        d=c
        for _ in range(4):d=cyclic_chart_shift(d)
        self.assertEqual(d,c)

    def test_cyclic_chart_preserves_positive_gaps(self):
        c=positive_gap_chart(Q(5,7),Q(1,9),Q(2,9),Q(1,3))
        for _ in range(16):
            self.assertGreater(1-sum(c[1:]),0)
            c=cyclic_chart_shift(c)
            self.assertTrue(all(t>0 for t in c[1:]))

    def test_four_signs_cancel_both_starting_indices(self):
        for sign in (-1,1):
            self.assertEqual(orbit_local_signs(sign),
                             (sign,-sign,sign,-sign))
            self.assertEqual(sum(orbit_local_signs(sign)),0)

    def test_six_cyclic_order_sectors_two_graceful(self):
        sectors=order_sectors()
        self.assertEqual(len(sectors),6)
        self.assertEqual(len(set(sectors)),6)
        self.assertEqual([s for s in sectors if graceful_sector(s)],
                         [(0,1,2,3),(0,3,2,1)])

    def test_opposite_cyclic_order_sector_not_silently_discarded(self):
        self.assertFalse(graceful_sector((0,2,1,3)))
        self.assertFalse(graceful_sector((0,1,3,2)))

    def test_closed_simplex_has_four_total_collisions(self):
        for i in range(4):
            g=tuple(Q(int(j==i)) for j in range(4))
            self.assertTrue(total_collision_vertex(g))
        self.assertFalse(total_collision_vertex(tuple([Q(1,4)]*4)))
        self.assertFalse(total_collision_vertex((Q(0),Q(1,3),Q(1,3),Q(1,3))))

    def test_boundary_stratum_only_all_collision_is_residual_zero(self):
        # Simple square polygon, collision and non-total-collision faces.
        P=tuple(pt(*q) for q in ((0,0),(2,0),(2,2),(0,2)))
        for j in range(4):
            gaps=tuple(Q(int(k==j)) for k in range(4))
            pts=chart_points(P,Q(1,13),gaps)
            self.assertEqual(len(set(pts)),1)
            for sigma in (-1,1):
                self.assertEqual(square_residual(pts,sigma),(0,0,0,0))
        for gaps in ((0,Q(1,3),Q(1,3),Q(1,3)),
                     (Q(1,2),0,Q(1,4),Q(1,4)),
                     (0,0,Q(1,2),Q(1,2)),
                     (0,Q(1,2),0,Q(1,2))):
            pts=chart_points(P,Q(1,13),closed_gap_chart(*gaps))
            self.assertGreater(len(set(pts)),1)
            for sigma in (-1,1):
                self.assertNotEqual(square_residual(pts,sigma),(0,0,0,0))

    def test_cyclic_rotation_preserves_square_not_size_collapse(self):
        p=(pt(0,0),pt(2,0),pt(2,2),pt(0,2))
        for j in range(4):
            q=p[j:]+p[:j]
            self.assertTrue(is_nondegenerate_square(q,1))
            self.assertFalse(is_nondegenerate_square(q,-1))

    def test_invalid_charts_rejected(self):
        with self.assertRaises(ValueError):positive_gap_chart(0,0,Q(1,4),Q(1,4))
        with self.assertRaises(ValueError):positive_gap_chart(0,Q(1,2),Q(1,2),Q(1,2))
        with self.assertRaises(ValueError):closed_gap_chart(0,0,0,0)
        with self.assertRaises(ValueError):total_collision_vertex((Q(1),Q(1),0,0))

    def test_one_generic_open_edge_square_eight_labelings(self):
        P=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        orbits=isolated_distinct_edge_squares(P)
        self.assertEqual(len(orbits),1)
        records=next(iter(orbits.values()))
        self.assertEqual(len(records),8)
        self.assertEqual({sig for _,sig,_ in records},{-1,1})
        for sig in (-1,1):
            labeled=[(inds,sgn) for inds,s,sng in records if s==sig
                     for sgn in (sng,)]
            self.assertEqual(len(labeled),4)
            self.assertEqual(sum(sgn for _,sgn in labeled),0)

    def test_square_orbit_cyclic_relations_in_finite_polygon(self):
        P=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        records=next(iter(isolated_distinct_edge_squares(P).values()))
        for inds,sigma,sign in records:
            rotated=inds[1:]+inds[:1]
            found=[sgn for id2,sig2,sgn in records
                   if id2==rotated and sig2==sigma]
            self.assertEqual(found,[-sign])

    def test_full_orbit_label_count_is_even_in_finite_example(self):
        P=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
        orbits=isolated_distinct_edge_squares(P)
        self.assertEqual(sum(len(v) for v in orbits.values())%2,0)
        self.assertEqual(len(orbits)%2,1)  # Quotient parity differs from labeled parity.

    def test_main_audit_report(self):
        report=run_audit()
        self.assertEqual(report['domain_shift_det'],'-1')
        self.assertEqual(report['target_shift_det_plus'],'1')
        self.assertEqual(report['orbit_sign_sum'],0)
        self.assertEqual(report['graceful_sectors'],2)
        self.assertEqual(report['distinct_open_edge_square_orbits'],1)
        self.assertEqual(report['distinct_edge_labelings'],(8,))

    def test_independent_v216_cyclic_local_degrees(self):
        try:
            from boundary_degree_v216 import (
                pt as p216, vertex_neighborhood, boundary_degree,
            )
        except ImportError:
            self.skipTest('v216 module absent outside full checkout')
        P=tuple(p216(*xy) for xy in ((0,0),(3,1),(2,3),(-1,2)))
        arcs=[vertex_neighborhood(P,i) for i in range(4)]
        degrees=[]
        for k in range(4):
            cert=boundary_degree(arcs[k:]+arcs[:k],max_rays=24)
            degrees.append(cert.degree)
        self.assertEqual(degrees,[-1,1,-1,1])


if __name__=='__main__':
    unittest.main()