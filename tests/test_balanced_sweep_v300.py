import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from balanced_sweep_v300 import (
    LinearMap, atlas_band_red, atlas_cost, balanced_slice, balancing_weight,
    cumulative, physical_squares, profile, profile_value, zero_intervals,
)


class BalancedSweepTests(unittest.TestCase):
    def test_mass_normalization_and_cap(self):
        for u in (0, Q(1, 7), Q(1, 2), 1):
            mass = cumulative(u)
            self.assertEqual((mass(0), mass(1)), (0, 1))
            for a, b in zip(mass.knots, mass.knots[1:]):
                slope = (mass(b) - mass(a)) / (b - a)
                self.assertLessEqual(slope, Q(9, 2))
                self.assertGreaterEqual(slope, Q(1, 4))
                x = (a + 2*b) / 3
                self.assertEqual(mass.inverse()(mass(x)), x)
        self.assertEqual(cumulative(1)(Q(9, 16)) - cumulative(1)(Q(7, 16)), Q(9, 16))

    def test_atlas_source_reconstruction(self):
        edges = [((0, 0), (4, 0)), ((4, 0), (1, 1)),
                 ((1, 1), (0, 3)), ((0, 3), (0, 0))]
        points = [(9, 0), (17, 9), (8, 17), (0, 8)]
        tags = []
        for i, ((a, b), p) in enumerate(zip(edges, points)):
            axis = 0 if b[0] != a[0] else 1
            loc = (Q(p[axis], 11) - a[axis]) / (b[axis] - a[axis])
            self.assertEqual(tuple(A + loc*(B-A) for A, B in zip(a, b)), tuple(Q(v, 11) for v in p))
            tags.append((i + loc) / 4)
        self.assertEqual(atlas_cost(tags), (Q(409, 528), Q(264, 119)))
        self.assertEqual(atlas_cost([Q(0), Q(1, 16), Q(1, 2), Q(11, 12)]), (Q(23, 48), Q(24, 23)))

    def test_reference_band_against_four_edge_parameters(self):
        for x in (Q(2, 5), Q(1, 2), Q(3, 5), Q(4, 7)):
            H = 1 + 3*x + x*x
            local = [(1 - 3*x/4 - x*x/4)/H, (4*x + x*x)/H,
                     (1-x*x)/H, (10*x/3 + 5*x*x/3)/H]
            q, cap = atlas_cost([(i + v)/4 for i, v in enumerate(local)])
            self.assertEqual(q, atlas_band_red(x))
            self.assertLessEqual(cap, Q(316, 103))

    def test_square_coordinates_and_balancing_are_independent(self):
        for depth in range(5):
            h = profile(depth)
            roots = physical_squares(h)
            self.assertEqual(len(roots), 1)
            (s, end), = roots
            t = s + 2*h(s)
            self.assertEqual(s, end)
            self.assertEqual(h(s), h(t))
            self.assertGreaterEqual(t-s, Q(1, 9))
            u, v = balancing_weight(s, t)
            self.assertEqual(u, v)
            mass = cumulative(u)
            self.assertEqual(mass(t)-mass(s), Q(1, 2))
            balanced, _ = balanced_slice(h, u)
            self.assertTrue(any(lo <= mass(s) <= hi for lo, hi in balanced))
            corners = [(s, -2*h(s)/3), (t, -2*h(t)/3), (t, 4*h(t)/3), (s, 4*h(s)/3)]
            squared = [(b[0]-a[0])**2 + (b[1]-a[1])**2 for a,b in zip(corners, corners[1:]+corners[:1])]
            self.assertEqual(squared, [(t-s)**2]*4)
        self.assertNotEqual(physical_squares(profile(0)), physical_squares(profile(1)))

    def test_all_finite_gate_zeros_have_the_required_aspect(self):
        for depth in range(5):
            h = profile(depth)
            self.assertGreaterEqual(balanced_slice(h, 0)[1][0], Q(5, 32))
            self.assertLessEqual(balanced_slice(h, 1)[1][1], -Q(1, 4))

    def test_refinement_retains_nodes_and_adds_variation(self):
        for n in range(1, 5):
            coarse = profile(n-1)
            for x in coarse.knots:
                self.assertEqual(coarse(x), profile_value(n, x))
            for i in range(4**n):
                a, b = Q(i, 4**n), Q(i+1, 4**n)
                middle = (a+b)/2
                excess = profile_value(n+1, middle) - (profile_value(n+1, a)+profile_value(n+1, b))/2
                expected = Q(1 if middle < Q(1, 2) else 2, 2**(n+6))
                self.assertEqual(excess, expected)

    def test_multiple_roots_at_one_weight(self):
        left = [(0, 0), (2, 2), (4, 4), (5, 11), (6, 8), (7, 11),
                (8, 4), (9, 12), (14, 16), (16, 18)]
        points = left + [(32-x,y) for x,y in reversed(left[:-1])]
        h = LinearMap([Q(x, 32) for x,y in points], [Q(y, 32) for x,y in points])
        roots, bounds = balanced_slice(h, Q(1, 4))
        self.assertEqual(len(roots), 5)
        self.assertEqual(bounds, (Q(7, 400), Q(49, 300)))
        self.assertGreaterEqual(balanced_slice(h, 0)[1][0], Q(5, 32))
        self.assertLessEqual(balanced_slice(h, 1)[1][1], -Q(1, 4))

    def test_zero_plateaus_and_knot_ownership(self):
        h = LinearMap([0,Q(1,8),Q(7,8),1], [0,Q(1,8),Q(1,8),0])
        self.assertEqual(physical_squares(h), [(Q(1,8),Q(5,8))])
        roots, bounds = balanced_slice(h, 1)
        self.assertEqual(roots, [(Q(1,16), Q(7,16))])
        self.assertEqual(bounds, (-Q(5,36), Q(5,36)))
        self.assertEqual(zero_intervals([0,1,2,3], lambda x: [1,0,0,-1][x]), [(1,2)])

    def test_independent_midline_drift_removes_vertical_rectangles(self):
        h = profile(1)
        (s, _), = physical_squares(h)
        t = s + 2*h(s)
        delta = Q(1, 10**6)
        def mid(x):
            return h(x)/3 + delta*x
        self.assertEqual(mid(t) - mid(s), delta * (t-s))
        self.assertNotEqual(mid(s), mid(t))

    def test_elliptic_component_weight_crossing(self):
        # On the circle, rotate a square through a quarter turn. The two
        # source densities alternate between 3/2 and 1/2 on quarter arcs.
        def signed_red(start, phase):
            knots = [Q(j,4) + phase for j in range(-1,9)]
            mass = Q(0)
            for a in (start, start + Q(1,2)):
                b = a + Q(1,4)
                for j, (l, r) in enumerate(zip(knots, knots[1:]), -1):
                    density = Q(3,2) if j % 2 == 0 else Q(1,2)
                    mass += density * max(0, min(b,r) - max(a,l))
            return mass - Q(1,2)
        for j in range(17):
            theta = Q(1,8) + Q(j,128)
            a, b = signed_red(theta, 0), signed_red(theta, Q(1,8))
            self.assertEqual(signed_red(theta+Q(1,4), 0), -a)
            self.assertEqual(signed_red(theta+Q(1,4), Q(1,8)), -b)
            self.assertEqual(-a/(b-a), Q(j,16))

    def test_invalid_distribution_and_interpolation(self):
        with self.assertRaises(ValueError):
            cumulative(Q(3,2))
        with self.assertRaises(ValueError):
            LinearMap([0,0], [1,2])
        with self.assertRaises(ValueError):
            LinearMap([0,1],[0,0]).inverse()
        with self.assertRaises(ValueError):
            profile(0)(-1)


if __name__ == "__main__":
    unittest.main()
