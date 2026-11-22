from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from angular_windows_v328_v330 import (
    INF, arctan_enclosure, carrier_exclusion, circle_values, ellipse_window, fixed_target_bounds,
    gadget, original_family_bounds, pi_enclosure, plateau_values, primitive, retuned_turns,
    spiral_fixture, square_with_gadgets, standard_mesh, tagged_return_integral, window_enclosure,
)


class Enclosures(unittest.TestCase):
    def test_pi_and_arctan_brackets(self):
        lo, hi = pi_enclosure(30)
        self.assertTrue(Q(314159265358979323846, 10**20) < lo < hi < Q(314159265358979323847, 10**20))
        self.assertLess(hi-lo, Q(1, 10**40))
        a, b = arctan_enclosure(Q(1), 400)          # contains pi/4
        self.assertTrue(a <= hi/4 and lo/4 <= b)
        with self.assertRaises(ValueError):
            arctan_enclosure(Q(2), 5)


class WindowAverages(unittest.TestCase):
    def test_circle_model_is_contained(self):
        mesh = standard_mesh(32)
        for m in (0, 1, 4):
            lo, hi = window_enclosure(mesh, circle_values(mesh, Q(4)), m, Q(4))
            self.assertTrue(lo <= 4 <= hi)

    def test_refinement_tightens(self):
        coarse, fine = standard_mesh(8), standard_mesh(64)
        w1 = window_enclosure(coarse, circle_values(coarse, Q(1)), 2, Q(1))
        w2 = window_enclosure(fine, circle_values(fine, Q(1)), 2, Q(1))
        self.assertLess(w2[1]-w2[0], w1[1]-w1[0])

    def test_plateau_goes_to_zero_and_three_points_cannot_resolve(self):
        mesh = standard_mesh(64)
        hi = [window_enclosure(mesh, plateau_values(mesh, Q(1), Q(1, 2)), m, Q(1))[1] for m in (1, 4, 16)]
        self.assertTrue(hi[0] > hi[1] > hi[2] and hi[2] < Q(1, 10**4))
        lo, up = window_enclosure([Q(0), Q(1), INF], circle_values([Q(0), Q(1), INF], Q(9)), 4, Q(9))
        self.assertEqual(lo, 0)
        self.assertGreater(up, 18)

    def test_mesh_and_data_validation(self):
        with self.assertRaises(ValueError):
            window_enclosure([Q(0), Q(2), INF], [(0, 0)]*3, 1, 1)
        mesh = standard_mesh(4)
        bad = [(Q(5), Q(5))] + circle_values(mesh, Q(1))[1:]
        with self.assertRaises(ValueError):
            window_enclosure(mesh, bad, 1, Q(1))


class Ellipses(unittest.TestCase):
    def test_mean_slope_at_m_zero_and_circle(self):
        self.assertEqual(ellipse_window(0, Q(3), Q(1))["D_single"], 3)    # A/pi = ab
        w = ellipse_window(7, Q(2), Q(2))
        self.assertEqual(w["D_single"], 4)
        self.assertEqual(w["D_low"], (4, 4))

    def test_bounds_bracket_and_converge_to_square_radius(self):
        for a, b in ((2, 1), (8, 1), (32, 1), (6, 3)):
            q2 = Q(2*a*a*b*b, a*a+b*b)
            gaps = []
            for m in (1, 16, 256):
                w = ellipse_window(m, Q(a), Q(b))
                self.assertTrue(w["D_low"][1] <= w["D_single"] <= w["D_high"][0])
                self.assertTrue(w["D_low"][1] < q2 < w["D_high"][0])
                gaps.append(w["D_high"][1]-w["D_low"][0])
            self.assertTrue(gaps[0] > gaps[1] > gaps[2])

    def test_scale_invariance(self):
        w1, w3 = ellipse_window(9, Q(8), Q(1)), ellipse_window(9, Q(24), Q(3))
        self.assertEqual(w3["D_single"], 9*w1["D_single"])


class Spirals(unittest.TestCase):
    def test_general_turn_count_return_integral(self):
        for k, N in ((2, 4), (3, 9), (3, 50), (5, 7)):
            P = spiral_fixture(k, N)
            e = Q(1, k)
            start, end = P.index((Q(0), Q(0))), P.index((-e, Q(0)))
            back, stop = P.index((-e, e)), P.index((Q(0), e))
            # Return paths of v324 measured with the origin at the attachment point.
            joint = primitive(P[start:end+1])[-1] + primitive(P[back:stop+1])[-1]
            self.assertEqual(joint, tagged_return_integral(k, N))

    def test_retuned_limit_is_half_the_area(self):
        self.assertEqual(retuned_turns(8), 647)
        self.assertLess(abs(tagged_return_integral(64, retuned_turns(64))-8), Q(1, 1000))

    def test_original_passes_and_retuned_fails(self):
        self.assertIsNotNone(carrier_exclusion(spiral_fixture(8, 64))["window_half_width"])
        self.assertIsNone(carrier_exclusion(spiral_fixture(8, retuned_turns(8)))["window_half_width"])
        self.assertIsNone(carrier_exclusion(spiral_fixture(4, 16))["window_half_width"])

    def test_analytic_bounds_cover_exact_values(self):
        for k in (3, 6, 9):
            c, b = carrier_exclusion(spiral_fixture(k, k*k)), original_family_bounds(k)
            self.assertLessEqual(c["oscillation"], b["osc_upper"])
            self.assertLessEqual(c["radius_sq"], b["radius_sq_upper"])
            self.assertGreaterEqual(c["area"], 16)
        self.assertGreater(original_family_bounds(7)["margin"], 0)
        self.assertLess(original_family_bounds(6)["margin"], 0)

    def test_fixed_target_truncation(self):
        gadgets = [(4**(j+3), Q(1, 2**(j+3)), Q(1, 2**j)) for j in (2, 1)]
        c = carrier_exclusion(square_with_gadgets(gadgets))
        b = fixed_target_bounds()
        self.assertLessEqual(c["oscillation"], b["osc_upper"])
        self.assertEqual(c["radius_sq"], 8)
        self.assertGreater(b["margin"], 0)

    def test_gadget_validation(self):
        with self.assertRaises(ValueError):
            gadget(0)


if __name__ == "__main__":
    unittest.main()
