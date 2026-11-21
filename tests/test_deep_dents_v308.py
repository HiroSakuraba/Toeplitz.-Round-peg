from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from deep_dents_v308 import certify, covers_segment, example, pt
from polygon_orbit_census_v223 import _on_seg


class DeepDentTests(unittest.TestCase):
    def test_opposing_deep_dents(self):
        poly, lanes, low, high, cap = example("opposite_dents")
        for loop in (poly,poly[::-1],poly[3:]+poly[:3]):
            record = certify(loop,lanes,low,high,cap)
            self.assertEqual(record["occupied_cells"],[1,1,1])
            self.assertEqual(record["forced_side_squared_at_least"],"225/32")
            self.assertEqual(record["forbidden_half_diagonal_squared_lower_open"],"1/16")
        self.assertEqual(poly[2].y-low,Q(255,64))
        self.assertEqual(high-poly[8].y,Q(255,64))

    def test_even_class_can_have_large_square(self):
        poly,*args = example("two_runs")
        r = certify(poly,*args)
        self.assertEqual(r["runs"],2)
        self.assertIsNone(r["forced_side_squared_at_least"])
        for x in (-6,-2):
            for y in (-2,2):
                self.assertTrue(any(_on_seg(pt(x,y),a,b) for a,b in zip(poly,poly[1:]+poly[:1])))

    def test_three_runs(self):
        record = certify(*example("three_runs"))
        self.assertEqual(record["occupied_cells"],[1,0,1,0,1])
        self.assertEqual(record["runs"],3)
        self.assertEqual(record["parity"],1)

    def test_whole_segment_coverage(self):
        boxes = [(0,1,0,1),(2,3,0,1)]
        self.assertFalse(covers_segment(pt(0,0),pt(3,1),boxes))
        self.assertTrue(covers_segment(pt(0,0),pt(3,1),boxes+[(1,2,0,1)]))

    def test_reject_strips_too_wide_and_bad_edges(self):
        poly,lanes,low,high,cap = example("opposite_dents")
        with self.assertRaisesRegex(ValueError,"excluded"):
            certify(poly,lanes,low,high,1)
        with self.assertRaisesRegex(ValueError,"uncovered"):
            certify((pt(-6,-2),pt(6,-2),pt(6,2)),lanes,low,high,cap)

    def test_perturbed_strip_allowance(self):
        poly,lanes,low,high,cap = example("opposite_dents")
        eta = Q(1,8)
        # Rigid translation attains the registered perturbation budget.
        moved = tuple(pt(p.x+eta,p.y) for p in poly)
        enlarged = [(a-eta,b+eta) for a,b in lanes]
        record = certify(moved,enlarged,low-eta,high+eta,cap+2*eta)
        self.assertEqual(record["occupied_cells"],[1,1,1])
        self.assertEqual(record["forced_side_squared_at_least"],"49/8")

    def test_interior_cross_control(self):
        square = (pt(1,0),pt(0,1),pt(-1,0),pt(0,-1))
        self.assertTrue(all(p.x==0 or p.y==0 for p in square))
        self.assertTrue(all(p.x**2+p.y**2==1 for p in square))
        self.assertTrue(all((a.x-b.x)**2+(a.y-b.y)**2==2
                            for a,b in zip(square,square[1:]+square[:1])))
