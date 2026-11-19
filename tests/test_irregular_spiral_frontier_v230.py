"""v230 exact threshold and conditional positive-scale checks."""
from fractions import Fraction as F
from pathlib import Path
import random
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from irregular_spiral_frontier_v230 import (
    profile,finite_q_variation,conditional_box_a_bounds,protected_by_box_a,
    proof_inequality_checks,report,
)


class SpiralThresholds(unittest.TestCase):
    def test_rectifiable_p_below_one(self):
        self.assertTrue(profile(F(9,10)).rectifiable)
        self.assertFalse(profile(1).rectifiable)
    def test_liouville_endpoint(self):
        self.assertTrue(profile(F(1999,1000)).liouville_regular)
        self.assertFalse(profile(2).liouville_regular)
    def test_non_rectifiable_but_liouville(self):
        a=profile(F(3,2))
        self.assertFalse(a.rectifiable)
        self.assertTrue(a.liouville_regular)
        self.assertTrue(a.finite_variation_below_two)
    def test_irregular_examples(self):
        for p in (2,F(5,2),3,10,100):
            x=profile(p)
            self.assertFalse(x.rectifiable)
            self.assertFalse(x.liouville_regular)
            self.assertFalse(x.finite_variation_below_two)
            self.assertTrue(x.box_a_applies)
    def test_variation_threshold_critical(self):
        self.assertTrue(finite_q_variation(F(1,2),1))
        self.assertFalse(finite_q_variation(1,1))
        self.assertTrue(finite_q_variation(1,F(5,4)))
        self.assertFalse(finite_q_variation(F(3,2),F(3,2)))
        self.assertTrue(finite_q_variation(F(3,2),F(7,4)))
        self.assertFalse(finite_q_variation(2,F(1999,1000)))
        self.assertTrue(finite_q_variation(2,F(2001,1000)))
    def test_variation_seeded_exponents(self):
        rng=random.Random(230)
        for _ in range(60):
            p=F(rng.randint(1,160),40)
            q=F(rng.randint(40,160),40)
            self.assertEqual(finite_q_variation(p,q), p<1 or q>p)
    def test_invalid(self):
        for p in (0,-1):
            with self.assertRaises(ValueError):profile(p)
            with self.assertRaises(ValueError):finite_q_variation(p,2)
        with self.assertRaises(ValueError):finite_q_variation(2,F(99,100))


class ProtectedWitness(unittest.TestCase):
    def test_archimedean_rational_inequalities(self):
        self.assertTrue(proof_inequality_checks())
    def test_box_a_exact_bounds(self):
        for p in (2,F(5,2),3,7,1000):
            x=conditional_box_a_bounds(p)
            self.assertEqual(x.min_vertex_radius_strictly_greater_than,F(1,5))
            self.assertEqual(x.side_length_strictly_greater_than,F(1,4))
            self.assertEqual(x.max_vertex_radius_strictly_less_than,F(1))
            self.assertEqual(x.cutoff_radius_that_preserves_all_vertices,F(1,5))
            self.assertTrue(x.needs_external_face_certificate)
    def test_only_box_a(self):
        with self.assertRaises(ValueError):conditional_box_a_bounds(F(1999,1000))
    def test_truncation_certifiable(self):
        for p in (2,3,5,20):
            for e in (F(1,100),F(1,6),F(1,5)):
                self.assertTrue(protected_by_box_a(p,e))
            self.assertFalse(protected_by_box_a(p,F(1,4)))
    def test_truncation_cannot_claim_absence(self):
        # Negative return means no *uniform v230* certificate, not absence.
        self.assertFalse(protected_by_box_a(F(3,2),F(1,10)))
    def test_invalid_truncation(self):
        self.assertFalse(protected_by_box_a(2,0))
        self.assertFalse(protected_by_box_a(2,-1))
    def test_report_labels_proof_dependency(self):
        d=report()
        self.assertEqual(d['conditional_uniform_square_side_floor'],'>1/4')
        self.assertEqual(d['conditional_outer_vertex_radius_ceiling'],'<1')
        self.assertIn('not computed by this script',d['source_of_square_existence'])
        self.assertEqual(d['profiles']['2']['liouville_regular'],False)


if __name__=='__main__':unittest.main()
