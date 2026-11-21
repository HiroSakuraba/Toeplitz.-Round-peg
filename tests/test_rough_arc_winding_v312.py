from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from rough_arc_winding_v312 import (
    boundary_bit, complete, corridor, edges, folded_arcs, full_side_window_census,
    loop, marked_census, model_arcs, notch_strip, pt, roughening_diamonds,
    triangle_example, whole_curve_winding,
)


class RoughArcTests(unittest.TestCase):
    def test_winding_and_census_on_four_folds(self):
        arcs = folded_arcs()
        report = boundary_bit(arcs)
        self.assertEqual(report["winding"],1)
        self.assertEqual(sum(report["positive_ray_hits_by_face"]),3)
        sizes = marked_census(arcs)
        self.assertEqual(len(sizes),31)
        self.assertGreaterEqual(min(sizes), Q(7,8)**2)
        p = loop(arcs)
        self.assertEqual(len(p),30)
        self.assertEqual(whole_curve_winding(p,arcs[1],arcs[2])["winding"],-1)
        self.assertEqual(full_side_window_census(p,arcs[1],arcs[2]),sizes)

    def test_nonuniform_face_and_cancellation(self):
        arcs = model_arcs("winding_only")
        ps = [complete(arcs[1][0],b)[0].x for b in arcs[2]]
        self.assertLess(min(ps),0)
        self.assertGreater(max(ps),0)
        for name,expected in [("winding_only",1),("two",2),("empty",0)]:
            arcs = model_arcs(name)
            self.assertEqual(len(marked_census(arcs)),expected)
            self.assertEqual(boundary_bit(arcs)["bit"],expected%2)

    def test_reversed_traversal_changes_signed_winding(self):
        i,j,k,l = folded_arcs()
        self.assertEqual(boundary_bit((i,j[::-1],k,l))["winding"],-1)
        self.assertEqual(boundary_bit((i[::-1],j,k,l[::-1]))["winding"],1)

    def test_sharp_width_transition(self):
        e = Q(1,16)
        for h,expected in ((Q(1,8),"no uniform guarantee"),
                           (3*e,"closed-arc existence"),(Q(1,4),"odd degree")):
            self.assertEqual(corridor(h,e)["status"],expected)
            i = [pt(e,-Q(5,8)),pt(e,Q(5,8))]
            j = [pt(1-e,-h),pt(1-e,h)]
            k = [pt(1-h,1+e),pt(1+h,1+e)]
            l = [pt(-Q(5,8),1+e),pt(Q(5,8),1+e)]
            if h<3*e:
                self.assertEqual(marked_census((i,j,k,l)),[])
            elif h==3*e:
                q,b = j[-1],pt(1-e,1+e)
                self.assertEqual(complete(q,b),(pt(e,3*e),pt(e,1+e)))
            else:
                self.assertEqual(len(marked_census((i,j,k,l))),1)
                k,l = ([pt(p.x,1-e) for p in arc] for arc in (k,l))
                self.assertEqual(marked_census((i,j,k,l)),[(1-2*e)**2])

    def test_end_guards_and_domain_checks(self):
        with self.assertRaises(ValueError):
            corridor(a=Q(1,2))
        with self.assertRaises(ValueError):
            notch_strip(Q(1,64))
        i,j,k,l = model_arcs("winding_only")
        with self.assertRaisesRegex(ValueError,"artificial tail"):
            boundary_bit(([pt(0,-Q(1,16)),pt(0,Q(1,16))],j,k,l))
        with self.assertRaisesRegex(ValueError,"polygon"):
            whole_curve_winding([pt(0,0),pt(4,0),pt(1,3)],j,k)

    def test_triangle_and_infinite_series_domains(self):
        result = triangle_example()
        self.assertEqual((result["winding"],result["side_floor"]),(-1,"883/560"))
        self.assertEqual(roughening_diamonds(folded_arcs())["marked_edges"],20)
        self.assertEqual(notch_strip()["side_floor"],"39/40")

    def test_circle_diagonal_and_exact_squares(self):
        points = [pt((1-t*t)/(1+t*t),2*t/(1+t*t)) for t in map(Q,range(-9,10))]
        for q in points:
            for b in points:
                p,d = complete(q,b)
                self.assertEqual(p.x*p.x+p.y*p.y,d.x*d.x+d.y*d.y)
            p,d = complete(q,pt(-q.y,q.x))
            self.assertEqual((p.x*p.x+p.y*p.y,d.x*d.x+d.y*d.y),(1,1))
