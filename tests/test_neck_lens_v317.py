from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from neck_lens_v317 import (Graph, add_folds, affine_refinement, components, local_test, scan_tubes, select,
                            staircase_fixture, validate_grid_loop)


class LocalLensTests(unittest.TestCase):
    def test_genuine_tips(self):
        f = Graph([(0,0),(1,-1),(2,0)])
        g = Graph([(0,0),(1,1),(2,0)])
        self.assertIsNone(local_test(f,g,0,2)["minimum_margin"])

    def test_clipped_diamond_fails_neck_test(self):
        e = Q(1,1024)
        f = Graph([(e,-e),(1,-1),(2-e,-e)])
        g = Graph([(e,e),(1,1),(2-e,e)])
        self.assertIsNone(local_test(f,g,e,2-e))

    def test_reports_missing_geometry(self):
        with self.assertRaises(ValueError):
            Graph([(0,0),(1,2)])
        f = Graph([(0,0),(2,0)])
        g = Graph([(0,1),(1,0),(2,1)])
        with self.assertRaisesRegex(ValueError,"positive open lens"):
            local_test(f,g,0,2)
        self.assertIsNone(local_test(f,Graph([(0,1),(2,2)]),0,2))

    def test_level_contacts_and_plateaus(self):
        profile = [(Q(i),Q(z)) for i,z in enumerate((0,2,1,1,3,1,2,0))]
        self.assertEqual(list(components(profile,Q(1))),
                         [(Q(1,2),2),(3,5),(5,Q(13,2))])

    def test_unmarked_selection(self):
        curve = staircase_fixture(64)
        for loop in (curve, curve[11:]+curve[:11], list(reversed(curve))):
            cert = select(loop)
            self.assertEqual(cert["interval"],["0","8"])
            self.assertEqual(cert["height"],"771/256")

    def test_exact_transverse_box_and_singular_control(self):
        e = Q(1,256)
        f = Graph([(0,0),(4,-1),(8,1)])
        g = Graph([(0,e),(4,2+e),(8,1+e)])
        cert = affine_refinement(f,g,Q(0),Q(8))
        self.assertEqual(2*Q(cert["side_squared_floor"]),Q(7470893,802816))
        self.assertIsNone(affine_refinement(Graph([(0,0),(8,0)]),
                                           Graph([(0,3),(8,3)]),Q(0),Q(8)))

    def test_unmarked_folded_traversals_and_missing_support(self):
        e = Q(1,256)
        f = Graph([(0,0),(4,-1),(8,1)])
        g = Graph([(0,e),(4,2+e),(8,1+e)])
        box = affine_refinement(f,g,Q(0),Q(8))
        target = add_folds(staircase_fixture(64),64)
        self.assertEqual(len(scan_tubes(target,box)),4)
        self.assertIsNone(scan_tubes([(100,100),(102,100),(102,101),(100,101)],box))

    def test_grid_loop_simplicity(self):
        validate_grid_loop([(0,0),(1,0),(2,0),(2,1),(0,1)])
        with self.assertRaises(ValueError):
            select([(0,0),(3,0),(3,3),(1,3),(1,-1),(0,-1)])
        with self.assertRaises(ValueError):
            validate_grid_loop([(0,0),(3,0),(2,0),(2,2),(0,2)])
