from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from cutoff_interactions_v363 import (
    CUTOFF, EPS, cutoff_chart_cover, moving_square, reconstruction_masks, triangle,
)
from polygon_orbit_census_v223 import _on_seg, census_polygon, pt


class CutoffInteractionTests(unittest.TestCase):
    def test_parity_formula_on_every_basis_monomial(self):
        for n in range(9):
            for order in range(5):
                masks = list(reconstruction_masks(n,order))
                for k in range(min(n,order)+1):
                    for subset in combinations(range(n),k):
                        term = sum(1 << i for i in subset)
                        self.assertEqual(sum(m & term == term for m in masks)%2,1)
        self.assertEqual(list(reconstruction_masks(8,1)), [0,1,2,4,8,16,32,64,128])

    def test_global_cutoff_cover_and_seed(self):
        ids = (0,0,1,2)
        self.assertEqual(set(cutoff_chart_cover()),{ids[k:]+ids[:k] for k in range(4)})
        ledger = census_polygon(triangle(0))
        self.assertEqual(sorted(s.side_squared for s in ledger.squares),
                         [Q(72,25),Q(144,49),Q(360,121)])
        self.assertEqual(sum(s.side_squared>CUTOFF**2 for s in ledger.squares),1)

    def test_fourth_order_is_attained(self):
        for mask in range(16):
            poly = triangle(mask)
            points, size = moving_square(mask)
            self.assertEqual(size>CUTOFF**2,mask==15)
            for p in points:
                self.assertTrue(any(_on_seg(pt(*p),a,b) for a,b in zip(poly,poly[1:]+poly[:1])))
        self.assertEqual(moving_square(15,Q(19,20))[1],CUTOFF**2)
        self.assertLess(moving_square(15,Q(19,20)-EPS)[1],CUTOFF**2)
        self.assertGreater(moving_square(15,Q(19,20)+EPS)[1],CUTOFF**2)

    def test_distant_events_on_unequal_bays(self):
        from cutoff_interactions_v363 import distant_bays
        events, clearance = distant_bays()
        self.assertEqual([site for _,site in events],[3,2,1,0])
        self.assertEqual(len({t for t,_ in events}),4)
        self.assertTrue(all(0<t<1 for t,_ in events))
        self.assertEqual(clearance,Q(2,5))

    def test_box_spectra_and_cutoff_selection(self):
        from cutoff_interactions_v363 import box_size_intervals, clear_squared_cutoff
        corners = [(0,0,0,0),(1,1,0,0),(1,1,1,1),(0,0,1,1)]
        self.assertEqual(box_size_intervals(corners),[(1,1)])
        self.assertEqual(box_size_intervals(corners[:2]),[(Q(1,2),Q(1,2)),(1,1)])
        self.assertEqual(box_size_intervals([(0,1,0,1)]*2),[(0,2)])
        self.assertEqual(clear_squared_cutoff((Q(1,2),Q(3,2)),[corners]),Q(3,4))
        self.assertIsNone(clear_squared_cutoff((Q(1,2),1),[[(0,1,0,1)]*2]))

    def test_pair_distance_can_overestimate_three_way_contact(self):
        from cutoff_interactions_v363 import box_size_intervals, clear_squared_cutoff
        e = Q(1,100)
        boxes = [(x-e,x+e,y-e,y+e) for x,y in ((0,0),(1,0),(Q(1,2),Q(7,8)))]
        for pair in combinations(boxes,2):
            self.assertTrue(any(a<=1<=b for a,b in box_size_intervals(pair)))
        self.assertEqual(box_size_intervals(boxes),[])
        self.assertEqual(clear_squared_cutoff((Q(99,100),Q(101,100)),[boxes]),1)
