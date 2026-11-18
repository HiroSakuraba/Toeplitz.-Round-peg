"""v220 exact-rational tests for one fixed rectifiable Jordan curve."""
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from accumulating_notches_v220 import (
    point,notch,flat_control,polygon,assert_simple_polygon,side_squared,
    squared_diameter,square_residual,signed_area,rectilinear_perimeter,
    predicted_finite_area,predicted_finite_length,limiting_area,limiting_length,
    collision_chart_data,parameter_span,normalized_residual_at_collision,
    omitted_tail_vertical,report,
)


class InfiniteRectifiableJordanTests(unittest.TestCase):
    def test_all_finite_polygons_simple_exact(self):
        for n in range(0,11):
            p=polygon(n)
            self.assertEqual(len(p),4*n+4)
            self.assertTrue(assert_simple_polygon(p))

    def test_every_notch_is_exact_positive_square(self):
        for n in list(range(1,18)) + [32,64]:
            p=notch(n)
            s=Q(1,4**(n+1))
            self.assertEqual(side_squared(p),s*s)
            self.assertEqual(squared_diameter(p),2*s*s)
            self.assertEqual(square_residual(p,+1),(0,0,0,0))
            self.assertNotEqual(square_residual(p,-1),(0,0,0,0))
            self.assertEqual(len(set(p)),4)

    def test_notch_intervals_separated_from_each_other(self):
        for n in range(1,65):
            self.assertLess(notch(n+1)[3].x, notch(n)[0].x)
            self.assertLess(notch(n+1)[3].x,flat_control(n)[0].x)

    def test_notches_accumulate_at_origin(self):
        for n in range(2,24):
            self.assertLess(max(abs(p.x)+abs(p.y) for p in notch(n+1)),
                            max(abs(p.x)+abs(p.y) for p in notch(n)))
        self.assertLess(max(abs(p.x)+abs(p.y) for p in notch(24)),Q(1,100000))

    def test_rectifiable_exact_perimeter_formula(self):
        for N in range(0,13):
            self.assertEqual(rectilinear_perimeter(polygon(N)),
                             predicted_finite_length(N))
        self.assertEqual(limiting_length(),Q(37,6))
        for N in range(2,13):
            self.assertLess(predicted_finite_length(N),limiting_length())
            self.assertEqual(limiting_length()-predicted_finite_length(N),
                             Q(1,6*4**N))

    def test_area_formula(self):
        for N in range(0,13):
            self.assertEqual(signed_area(polygon(N)),predicted_finite_area(N))
        self.assertEqual(limiting_area(),Q(481,240))
        self.assertEqual(limiting_area()-predicted_finite_area(7),
                         Q(1,240*16**7))

    def test_omitted_notch_Hausdorff_upper_bound(self):
        for N in range(0,12):
            self.assertEqual(omitted_tail_vertical(N),notch(N+1)[3].x-notch(N+1)[0].x)
        with self.assertRaises(ValueError):omitted_tail_vertical(-1)

    def test_equal_parameter_gap_ratios_on_two_distinct_sequences(self):
        for n in range(1,25):
            sq=collision_chart_data(n,'square')
            fl=collision_chart_data(n,'flat')
            self.assertEqual(sq[1],(Q(1,3),)*3)
            self.assertEqual(fl[1],(Q(1,3),)*3)
            self.assertEqual(sq[2],fl[2])
            self.assertEqual(sq[2],parameter_span(n))
            self.assertEqual(sq[0]-fl[0],parameter_span(n))
            self.assertGreater(fl[0],0)
            self.assertLess(sq[0],Q(1,2**(n-1)))

    def test_normalized_residual_two_incompatible_limits(self):
        L=limiting_length()
        expected=(L/3,-L/3,L,-L/3)
        for n in list(range(1,15))+[40,100]:
            self.assertEqual(normalized_residual_at_collision(n,'square'),(0,0,0,0))
            self.assertEqual(normalized_residual_at_collision(n,'flat'),expected)
            self.assertEqual(square_residual(flat_control(n),+1),
                             (Q(1,4**(n+1)),-Q(1,4**(n+1)),
                              3*Q(1,4**(n+1)),-Q(1,4**(n+1))))
        self.assertNotEqual(expected,(0,0,0,0))

    def test_no_zero_free_parameter_ratio_collision_link(self):
        # At the SAME blown-up ratio point (t=0; a=b=c=1/3), one
        # sequence has normalized square residual zero; the other has
        # nonzero constant residual. A continuous target-valued link
        # retaining F/epsilon cannot exist on that bare parameter chart.
        self.assertEqual(collision_chart_data(40,'square')[1],
                         collision_chart_data(40,'flat')[1])
        self.assertEqual(normalized_residual_at_collision(40,'square'),(0,0,0,0))
        self.assertNotEqual(normalized_residual_at_collision(40,'flat'),(0,0,0,0))

    def test_finite_polygon_crosscheck_against_v216(self):
        try:
            from boundary_degree_v216 import pt as vpt,validate_simple_polygon
        except ImportError:
            self.skipTest('v216 independent polygon checker available in full checkout')
        for n in range(0,7):
            p=tuple(vpt(v.x,v.y) for v in polygon(n))
            self.assertGreater(validate_simple_polygon(p),0)

    def test_input_rejections(self):
        for bad in (-4,0,1.25,True):
            with self.assertRaises(ValueError):notch(bad)
        for bad in (-2,1.5,True):
            with self.assertRaises(ValueError):polygon(bad)
        with self.assertRaises(ValueError):collision_chart_data(1,'other')
        with self.assertRaises(ValueError):square_residual(notch(1),0)

    def test_report_consistency(self):
        x=report(8)
        self.assertEqual(x['smallest_square_residual'],('0','0','0','0'))
        self.assertEqual(x['polygon_vertices'],36)
        self.assertEqual(x['limiting_Jordan_curve_perimeter'],'37/6')
        self.assertEqual(x['normalized_gap_direction'],('1/3','1/3','1/3'))

if __name__=='__main__':
    unittest.main()