from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from half_turn_degree_v313 import (
    cut_crossings,ear_census,mul,norm_bounds,polygon,profile,pt,
    tube_certificate,winding_witness,
)
from boundary_or_winding_v313 import contact_interval, window_alternative
from half_turn_degree_v313 import grid_witnesses
from polygon_orbit_census_v223 import _evaluate,_side_sq,_square_residual


class HalfTurnTests(unittest.TestCase):
    def test_dyadic_both_alternatives(self):
        boundary,local=grid_witnesses()
        self.assertEqual(boundary["side_squared"],"4")
        self.assertEqual(local["winding"],1)
        self.assertGreater(Q(local["separation_squared"]),Q(4,25))

    def test_collinear_intervals_and_constants(self):
        a,b,c,d=pt(0,0),pt(4,0),pt(1,0),pt(2,0)
        self.assertEqual(contact_interval(a,b,d,c),(Q(1,4),Q(1,2)))
        self.assertEqual(contact_interval(c,c,a,b),(0,1))
        self.assertEqual(contact_interval(b,b,a,c),None)
        self.assertEqual(contact_interval(a,b,b,pt(4,3)),(1,1))

    def test_degenerate_rays_and_separation(self):
        p=[pt(0,0),pt(4,0),pt(4,1),pt(0,1)]
        # The full moving completion runs along a horizontal boundary edge.
        result=window_alternative(p,[pt(1,0),pt(2,0)],[pt(1,1),pt(2,1)])
        self.assertEqual(result["outcome"],"boundary square")
        self.assertEqual(result["side_squared"],"1")

    def test_vertex_hit_can_still_have_defined_winding(self):
        from rough_arc_winding_v312 import whole_curve_winding
        p=[pt(0,0),pt(4,0),pt(1,3)]
        j=[pt(Q(1,2),0),pt(1,0)]
        k=[pt(2,0),pt(4,0)]
        with self.assertRaisesRegex(ValueError,"segment endpoint"):
            whole_curve_winding(p,j,k)
        self.assertEqual(window_alternative(p,j,k),{"outcome":"winding","winding":-1})

    def test_regular_profile_census_and_crossing(self):
        rows = ear_census()
        self.assertEqual(len(rows),2)
        edges=list(zip(profile(),profile()[1:]))
        for size,ids,tags in rows:
            points=_evaluate([edges[i] for i in ids],tags)
            self.assertEqual(_side_sq(points),size)
            self.assertEqual(_square_residual(points,1),(0,0,0,0))
        events=cut_crossings()
        event=Q(events["event_scale_squared"])
        self.assertEqual(event*rows[-1][0],Q(4,25))
        self.assertTrue(Q(81,100)<event<1)
        self.assertLess(rows[0][0],Q(4,25))

    def test_half_turn_and_broken_pair(self):
        for top,bottom,symmetric in [(1,1,True),(Q(9,10),Q(9,10),True),(1,Q(9,10),False)]:
            points=set(polygon(top,bottom))
            self.assertEqual(points=={mul(p,-1) for p in points},symmetric)

    def test_certified_tubes(self):
        small=tube_certificate(polygon())
        large=tube_certificate(polygon(),Q(7,5),Q(1,32),upper_side=Q(8,5))
        self.assertEqual(large["side_interval"],["7/5","8/5"])
        self.assertTrue(all(min(r["excluded_leaves"].values())>0 for r in (small,large)))
        with self.assertRaisesRegex(ValueError,"unresolved"):
            tube_certificate(polygon(),Q(2),Q(0),max_depth=8)

    def test_invalid_intervals_and_norms(self):
        with self.assertRaises(ValueError):
            tube_certificate(polygon(),upper_side=Q(1,10))
        self.assertEqual(norm_bounds((-1,2),(3,4)),(9,20))
        self.assertEqual(norm_bounds((-3,-1),(-2,1)),(1,13))

    def test_nonzero_window_and_explicit_grid_contact(self):
        self.assertEqual(winding_witness()["winding"],1)
        from polygon_orbit_census_v223 import _on_seg
        boundary=polygon()
        edges=list(zip(boundary,boundary[1:]+boundary[:1]))
        points=[pt(x,y) for x,y in [(Q(-7,2),-1),(Q(-3,2),-1),(Q(-3,2),1),(Q(-7,2),1)]]
        self.assertTrue(all(any(_on_seg(p,a,b) for a,b in edges) for p in points))
        self.assertEqual(_square_residual(points,1),(0,0,0,0))
