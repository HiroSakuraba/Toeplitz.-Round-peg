from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from direction_search_v303 import arc_reserve, broken_symmetry, direction_certificates, example, moved
from polygon_orbit_census_v223 import pt


class DirectionSearch(unittest.TestCase):
    def test_similarity_and_boundary_labels(self):
        points = tuple(map(moved, example()))
        expected = [{"start": pt(7, -2), "end": pt(10, 2),
                     "length_squared": Q(25), "normalized_tube": Q(1, 1000)}]
        for loop in (points, points[3:]+points[:3], points[::-1]):
            self.assertEqual(direction_certificates(loop), expected)

    def test_crest_interior_must_be_checked(self):
        arc = (pt(0, 0), pt(Q(1, 4), Q(3, 8)), pt(Q(3, 4), Q(3, 8)), pt(1, 0))
        self.assertIsNone(arc_reserve(arc, True))

    def test_zero_reserve_is_accepted(self):
        points = (pt(0, 0), pt(Q(1, 2), Q(-3, 8)+Q(11, 1000)), pt(1, 0),
                  pt(Q(1, 2), Q(3, 4)))
        self.assertEqual(direction_certificates(points)[0]["normalized_tube"], 0)

    def test_smooth_obstruction_has_centered_rectangles_on_both_sides(self):
        # The Pythagorean triple gives rational unit radii directions near 45 degrees.
        for x, y in ((Q(120, 169), Q(119, 169)), (Q(119, 169), Q(120, 169))):
            def radius(c, s):
                return (1 + (c**4-6*c*c*s*s+s**4)/100)/2
            self.assertEqual(x*x+y*y, 1)
            r = radius(x, y)
            self.assertGreater(r, Q(29, 100))
            self.assertLess(r, Q(69, 100))
            for u, v in ((-x, y), (-x, -y), (x, -y), (-y, x)):
                self.assertEqual(radius(u, v), r)
            self.assertLess(Q(99, 100), y/x)
            self.assertLess(y/x, Q(100, 99))

    def test_unique_centered_square_after_symmetry_breaking(self):
        c2, s2 = Q(3, 5), Q(4, 5)  # Double angle of the diagonal direction (2,1).
        second_harmonic = (4*c2-3*s2)/5
        fourth_harmonic = c2*c2-s2*s2
        radius = (1 + fourth_harmonic/100 + second_harmonic/1000)/2
        self.assertEqual(second_harmonic, 0)
        self.assertEqual(radius, Q(2493, 5000))
        self.assertNotEqual((4*s2+3*c2)/5, 0)
        self.assertGreater(Q(broken_symmetry()["derivative_margin"]), 0)

    def test_no_certificate_and_invalid_input(self):
        square = (pt(0, 0), pt(1, 0), pt(1, 1), pt(0, 1))
        self.assertEqual(direction_certificates(square), [])
        with self.assertRaises(ValueError):
            direction_certificates((square[0], square[2], square[1], square[3]))
