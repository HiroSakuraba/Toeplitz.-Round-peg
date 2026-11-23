"""v331 integer Floer grading. Needs numpy; set REQUIRE_NUMERIC_DEPS=1 to make a missing numpy a
failure instead of a skip."""
from fractions import Fraction as Q
from pathlib import Path
import os
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
try:
    from floer_grading_v331 import Grader, ellipse_grading, grading_filter
    from floer_parity_v331 import calibration, square_generators
    from floer_actions_v331 import PolygonParam, square_actions
    from residue_census_v331 import certify, float_search
    from angular_windows_v328_v330 import area, pi_enclosure, radius_squared
    HAVE_NUMPY = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_NUMPY = False
REQUIRE = os.environ.get("REQUIRE_NUMERIC_DEPS") == "1"
needs_numpy = unittest.skipUnless(HAVE_NUMPY or REQUIRE, "numpy not installed")

PENTAGON = [(0, 0), (7, 1), (8, 6), (3, 9), (-2, 5)]
TWELVE = [(0, 0), (5, 1), (9, 0), (11, 4), (8, 5), (10, 9), (6, 8), (4, 11), (2, 7), (-2, 9), (-1, 4), (-3, 2)]


def squares(P):
    quads, *_ = float_search(P, qmax2=1e9)
    found, _ = certify(P, quads, qmax2=Q(10**9))
    return [v[0] for v in found.values()]


@needs_numpy
class Calibration(unittest.TestCase):
    def test_ellipse_branches_are_two_and_one(self):
        for a, b in ((2.0, 1.0), (5.0, 1.0), (1.3, 1.0)):
            for theta in (0.3, 0.8, 1.2, 1.5):
                self.assertEqual((ellipse_grading(a, b, theta, "major"), ellipse_grading(a, b, theta, "minor")), (2, 1))

    def test_pentagon_carrier_has_grading_two(self):
        P = [(Q(x), Q(y)) for x, y in PENTAGON]
        A, top = area(P), pi_enclosure(30)[1]*radius_squared(P)/2
        for loop in (P, P[::-1]):
            G = Grader(loop)
            (sq,) = squares(loop)
            acts = square_actions(PolygonParam(loop), sq)
            grades = [G.grading(*g) for g in G.square_generators(sq)]
            self.assertEqual(sorted(grades), [1, 2])
            carrier = grades[[A/2 <= a <= top for a in acts].index(True)]
            self.assertEqual(carrier, 2)


@needs_numpy
class Consistency(unittest.TestCase):
    def test_parity_and_duality_on_the_twelve_gon(self):
        even = calibration()
        P = [(Q(x), Q(y)) for x, y in TWELVE]
        for loop in (P, P[::-1]):
            G = Grader(loop)
            all_grades = []
            for sq in squares(loop):
                grades = [G.grading(*g) for g in G.square_generators(sq)]
                self.assertEqual(sum(grades), 3)                      # duality: k <-> 3 - k
                for (_, sign), g in zip(square_generators(loop, sq), grades):
                    self.assertEqual(g % 2 == 0, sign == even)        # parity reproduces the v331 signs
                all_grades += grades
            self.assertEqual(sorted(all_grades), [0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3])

    def test_smoothed_grading_agrees_on_generic_contacts(self):
        P = [(Q(x), Q(y)) for x, y in TWELVE]
        G = Grader(P)
        for sq in squares(P):
            for gen in G.square_generators(sq):
                self.assertEqual(G.possible_gradings(*gen), {G.grading(*gen)})


@needs_numpy
class Spiral(unittest.TestCase):
    def test_no_grading_two_candidate_at_twenty_four_turns(self):
        row = grading_filter(24, workers=2)
        self.assertEqual(row["generic_grading_two"], 0)
        self.assertEqual(row["vertex_cases_admitting_two"], 0)
        self.assertEqual(row["vertex_contact_cases"], 43)
        self.assertEqual(row["duality_grading_sums"], {"3": 1001})
        self.assertEqual(max(map(int, row["generic_window_gradings"])), -35)


if __name__ == "__main__":
    unittest.main()
