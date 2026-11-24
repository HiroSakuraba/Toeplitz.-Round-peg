from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from cutoff_interactions_v363 import CENTERS, EPS, SIDE, WIDTH, cutoff_chart_cover, moving_square, triangle
from free_contact_v360 import _inside_clear
from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from three_patch_sweep_v364 import (
    completion_boxes, completion_side, crossing_times, fixture_boxes, fixture_cutoff, target_roof,
)


class ThreePatchSweepTests(unittest.TestCase):
    def test_every_mask_and_all_cutoff_charts(self):
        ids = (0,0,1,2)
        expected = {ids[k:]+ids[:k] for k in range(4)}
        self.assertEqual(set(cutoff_chart_cover(SIDE, upper_cutoff=fixture_cutoff(0))),expected)
        for n in (0,1,4,9):
            s,e,w = fixture_cutoff(n), EPS/2**n, WIDTH/2**n
            self.assertEqual(set(cutoff_chart_cover(s,e,w)),expected)
            for mask in range(8):
                triangle(mask,w,e)
                self.assertEqual(moving_square(mask,error=e)[1]>s*s,mask==7)
            corners, sq = moving_square(7,Q(13,14),error=e)
            self.assertEqual(sq,s*s)
            self.assertEqual(corners[3],CENTERS[3])

    def test_target_sides_and_actual_limiting_contacts(self):
        target = target_roof()
        self.assertIsNone(completion_side(target,fixture_boxes(0,0),fixture_cutoff(0)))
        for n in (1,3,6):
            s = fixture_cutoff(n)
            self.assertIs(completion_side(target,fixture_boxes(n,0),s),True)
            self.assertIs(completion_side(target,fixture_boxes(n,1),s),False)
            self.assertIsNone(completion_side(target,fixture_boxes(n,Q(1,2)),s))
        x = Q(23,28)
        self.assertEqual(x,CENTERS[3][0]+Q(1,4))
        poly = validate_jordan_polygon([pt(*p) for p in target])
        edges = tuple(zip(poly,poly[1:]+poly[:1]))
        for corner in ((x,0),(x+SIDE,0),(x+SIDE,SIDE),(x,SIDE)):
            self.assertIsNone(_inside_clear(corner,edges,0))

    def test_fold_pair_birth_and_death(self):
        for lam,count in ((0,1),(Q(2,7),2),(Q(1,2),3),(Q(11,14),2),(1,1)):
            self.assertEqual(len(crossing_times(lam)),count)
        self.assertGreater(crossing_times(0)[0],Q(2,3))
        self.assertLess(crossing_times(1)[0],Q(1,3))

    def test_false_positive_with_three_near_corners(self):
        poly = triangle(8)
        edges = tuple(zip(poly,poly[1:]+poly[:1]))
        self.assertTrue(_inside_clear(CENTERS[3],edges,EPS/2))
        self.assertIsNone(_inside_clear(CENTERS[3],edges,EPS))
        self.assertEqual(len(cutoff_chart_cover(SIDE)),4)
        self.assertGreater(moving_square(8)[1],SIDE**2)
        for x,y in CENTERS[:3]:
            self.assertIsNone(_inside_clear((x,y),edges,0))

    def test_ambiguous_assignments_and_closed_contacts(self):
        boxes = [(0,0,0,0),(1,1,0,0),(1,1,1,1)]
        self.assertEqual(completion_boxes(boxes,1),((0,0,1,1),))
        self.assertEqual(completion_boxes(boxes,2),())
        self.assertIsNone(completion_side(((0,0),(2,0),(0,1)),boxes,1))
        big = [(-2,2,-2,2)]*3
        self.assertEqual(len(completion_boxes(big,1)),3)
