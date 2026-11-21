from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from center_walls_v309 import certify_runs, example, perturbed_fixture, pt
from polygon_orbit_census_v223 import _on_seg


class CenterWallTests(unittest.TestCase):
    def test_original_and_quarter_unit_tube(self):
        data = example("two_runs")
        for eta,inset,expected_x,expected_side in [(0,1,3,"225/32"),(Q(1,4),0,2,"169/32")]:
            r = certify_runs(*data,displacement=eta,inset=inset)
            w, = r["walls"]
            self.assertEqual(Q(w["right_wall"]),expected_x)
            self.assertEqual(Q(w["left_wall"]),-expected_x)
            self.assertLess(Q(w["wall_radius_squared_upper"]),Q(r["gap_upper"])**2)
            self.assertEqual(r["side_squared_floor"],expected_side)
            self.assertEqual(r["large_squares_at_least"],2)

    def test_census_polygon_registration(self):
        p,*_ = example("two_runs")
        self.assertTrue(all((p0.x-p1.x)**2+(p0.y-p1.y)**2 < Q(1,1024)**2
                            for p0,p1 in zip(p,perturbed_fixture())))
        self.assertEqual(Q(1,1024)+Q(255,1024),Q(1,4))

    def test_midline_is_not_a_separator(self):
        p,*_ = example("two_runs")
        for x in (-2,2):
            for y in (-2,2):
                self.assertTrue(any(_on_seg(pt(x,y),a,b) for a,b in zip(p,p[1:]+p[:1])))

    def test_opposite_bridges(self):
        _,lanes,low,high,cap = example("three_runs")
        upper,lower = high-cap/2,low+cap/2
        vertices = [(-10,-2),(-6,-2),(-6,upper),(-2,upper),(-2,-2),(10,-2),
                    (10,2),(6,2),(6,lower),(2,lower),(2,2),(-10,2)]
        p = tuple(pt(x,y) for x,y in vertices)
        for loop in (p,p[::-1],p[2:]+p[:2]):
            r = certify_runs(loop,lanes,low,high,cap,inset=1)
            self.assertEqual([w["bridge"] for w in r["walls"]],["top","bottom"])
            self.assertEqual(r["large_squares_at_least"],3)

    def test_reject_two_caps_even_with_same_cell_marks(self):
        p,lanes,low,high,cap = example("two_runs")
        p = (p[0],pt(-1,-2),pt(-1,-2+cap/2),pt(-2,-2+cap/2),*p[2:])
        with self.assertRaisesRegex(ValueError,"both caps"):
            certify_runs(p,lanes,low,high,cap)

    def test_allowance_and_inset_fail_closed(self):
        for kw in ({"displacement":Q(1,2)},{"inset":2},{"displacement":-1}):
            with self.assertRaises(ValueError):
                certify_runs(*example("two_runs"),**kw)
