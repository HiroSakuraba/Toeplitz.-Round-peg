from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cutoff_payment_v315 import SIDE, bracket, displacement, event_data, pt, repair, replay


class CutoffPaymentTests(unittest.TestCase):
    def test_unique_crossing_and_regular_support(self):
        critical, _, size = event_data()
        self.assertEqual(critical, F(2312, 2665))
        self.assertEqual(critical * size, SIDE**2)

    def test_nested_brackets_and_repaired_squares(self):
        intervals = []
        for depth in (4, 8, 12):
            data = replay(depth)
            a, b = F(data["lower_scale"]), F(data["upper_scale"])
            self.assertEqual(b-a, F(1, 10 * 2**depth))
            self.assertGreater(F(data["repaired_square_side_squared"]), SIDE**2)
            self.assertEqual(F(data["remote_crossing_distance_lower_bound"]), SIDE)
            intervals.append((a,b))
        self.assertTrue(all(a <= c < d <= b for (a,b),(c,d) in zip(intervals,intervals[1:])))

    def test_plateau_taper_and_support(self):
        p, q = pt(0,0), pt(1,-1)
        for x, expected in ((0,1), (2,1), (3,F(1,2)), (4,0), (5,0)):
            self.assertEqual(displacement(pt(x,0),p,q,F(1)),pt(expected,-expected))

    def test_repair_at_polygon_vertices(self):
        curve = (pt(0,0),pt(4,0),pt(4,4),pt(0,4))
        shift = F(1,32)
        square = tuple(pt(p.x + shift, p.y) for p in curve)
        self.assertGreaterEqual(len(repair(curve, curve, square, shift)),4)
        with self.assertRaisesRegex(ValueError,"overlapping"):
            repair(curve, curve, square, F(1))
        with self.assertRaisesRegex(ValueError,"displacement"):
            repair(curve, curve, square, shift/2)
        with self.assertRaisesRegex(ValueError,"square"):
            repair(curve, curve, square[:3]+(pt(0,5),), shift)

    def test_invalid_brackets(self):
        for q, depth in ((F(1), 3), (F(2312,2665), -1)):
            with self.assertRaises(ValueError):
                bracket(q, depth)
