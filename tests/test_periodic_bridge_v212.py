"""Exact-rational regression tests for the v212 conditional lemma and obstruction."""
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from periodic_bridge_v212 import (  # noqa: E402
    Point, PeriodicGraph, SOURCE_WINDOWS, alternating_area, audit,
    four_strands, is_labeled_square, minimum_mod1_separation,
    no_square_inside_all_windows, rotate90, squared_length,
    square_at_unmodified_period, strand,
)


class PeriodicBridgeTests(unittest.TestCase):
    def test_periodic_graph_necessities(self):
        s = four_strands()
        self.assertEqual(len(s), 4)
        for g in s:
            self.assertEqual(g.height(F(0)), g.height(F(1)))
            for x in (F(0), F(1, 5), F(7, 8), F(19, 7), -F(4, 3)):
                self.assertEqual(g.height(x), g.height(x + 1))

    def test_exact_area_balance(self):
        s = four_strands()
        self.assertEqual(tuple(g.area_y_dx() for g in s),
                         (F(1), F(0), F(-1), F(0)))
        self.assertEqual(alternating_area(s), 0)
        # Equal + and - tent areas cancel for arbitrary rational amplitude.
        for a in (F(0), F(1, 128), F(1, 16), F(1, 8)):
            self.assertEqual(alternating_area(four_strands(a)), 0)

    def test_nonintersection_by_vertical_gap(self):
        for x in (F(i, 127) for i in range(128)):
            g0, _, g2, _ = four_strands()
            self.assertGreaterEqual(g0.height(x), F(15, 16))
            self.assertLessEqual(g2.height(x), -F(15, 16))
        # The bounds are analytic from |tent| <= 1 (sampling isn't proof).

    def test_exact_square_on_periodic_completion(self):
        s = four_strands()
        w = square_at_unmodified_period(s)
        self.assertTrue(is_labeled_square(w))
        self.assertEqual(squared_length(w[1] - w[0]), F(2))
        self.assertEqual(tuple(p.y for p in w), (F(1), F(0), F(-1), F(0)))
        self.assertEqual(rotate90(w[1] - w[0]), w[2] - w[1])
        self.assertEqual(tuple((p.x % 1) for p in w), (F(7, 8),) * 4)
        self.assertTrue(all(not(lo <= F(7, 8) <= hi) for lo, hi in SOURCE_WINDOWS))

    def test_no_square_on_all_original_windows_exact(self):
        self.assertEqual(minimum_mod1_separation(0, 2), F(7, 16))
        self.assertTrue(no_square_inside_all_windows())
        self.assertTrue(no_square_inside_all_windows(F(1, 8)))
        self.assertFalse(no_square_inside_all_windows(F(1, 4)))

    def test_perturbed_unbalanced_control(self):
        base = list(four_strands())
        shifted = PeriodicGraph(tuple(Point(p.x, p.y + F(1, 8))
                                      for p in base[0].vertices))
        base[0] = shifted
        self.assertEqual(alternating_area(tuple(base)), F(1, 8))
        self.assertFalse(is_labeled_square(square_at_unmodified_period(tuple(base))))

    def test_exact_audit(self):
        r = audit()
        self.assertEqual(r["alternating_area"], "0")
        self.assertTrue(r["is_exact_square"])
        self.assertTrue(r["no_all_window_square"])
        self.assertEqual(r["minimum_window_0_2_x_separation_mod1"], "7/16")


if __name__ == "__main__":
    unittest.main()
