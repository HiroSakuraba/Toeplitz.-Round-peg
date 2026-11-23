from fractions import Fraction as Q
from importlib.util import find_spec
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import _side_sq, pt
from terrain_component_degree_v337 import GLUED, TERRAIN, chart_cells, shear_witnesses


def orbit(square):
    return tuple(sorted(square))


def endpoints(*squares):
    return {orbit([pt(*p) for p in square]) for square in squares}


class ComponentDegree(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.terrain = dict(chart_cells(TERRAIN))
        cls.glued = dict(chart_cells(GLUED))

    def test_complete_terrain_components(self):
        big = endpoints(((2, 0), (2, 2), (0, 2), (0, 0)),
                        ((1, 0), (2, 1), (1, 2), (0, 1)))
        small = endpoints(((2, 2), (2, 3), (1, 3), (1, 2)),
                          ((2, 6), (2, 7), (1, 7), (1, 6)))
        self.assertEqual(len(self.terrain), 14)
        segments = {ids: {orbit(p) for p in corners}
                    for ids, corners in self.terrain.items() if len(corners) > 1}
        self.assertEqual(segments, {(0, 1, 4, 5): big, (1, 1, 3, 3): small})
        for corners in self.terrain.values():
            self.assertLessEqual({orbit(p) for p in corners}, big | small)

    def test_glued_component_and_size_gap(self):
        expected = endpoints(((3, 0), (7, 3), (4, 7), (0, 4)),
                             ((4, 0), (7, 4), (3, 7), (0, 3)))
        # Check the ordered segment, as matching unordered ends alone
        # would allow a spurious cyclic relabeling along the segment.
        ids = (0, 5, 6, 9)
        a, b = self.glued[ids]
        self.assertEqual({orbit(a), orbit(b)}, expected)
        self.assertEqual(tuple(pt((p.x+q.x)/2, (p.y+q.y)/2) for p, q in zip(a, b)),
                         tuple(pt(*p) for p in ((Q(7,2),0),(7,Q(7,2)),
                                                 (Q(7,2),7),(0,Q(7,2)))))
        for chart, corners in self.glued.items():
            if chart in (ids, (0, 4, 6, 9), (0, 5, 7, 9)):
                self.assertLessEqual({orbit(p) for p in corners}, expected)
            else:
                self.assertLess(max(map(_side_sq, corners)), 16)

    @unittest.skipUnless(find_spec("sympy"), "symbolic shear replay requires SymPy")
    def test_regular_shears(self):
        import sympy as sp

        t = sp.Symbol('t', positive=True)
        for polygon, cells, cutoff, direction, count, expected in (
            (TERRAIN, self.terrain, Q(3,2), 1, 9, {(0,1,4,5): (-56*t, 2)}),
            (GLUED, self.glued, 24, -1, 3, {(0,5,6,9): (1536*t, Q(49,2))}),
        ):
            candidates = [ids for ids, corners in cells.items()
                          if max(map(_side_sq, corners)) >= cutoff]
            self.assertEqual(len(candidates), count)
            self.assertEqual(shear_witnesses(polygon, candidates, direction), expected)

    def test_complete_census_after_finite_shear(self):
        for polygon, shear, cutoff, expected in (
            (TERRAIN, Q(1,100), Q(3,2), 2), (GLUED, Q(-1,100), 24, Q(49,2))
        ):
            deformed = [(x + shear*y, y) for x, y in polygon]
            large = [corners for _, corners in chart_cells(deformed)
                     if max(map(_side_sq, corners)) > cutoff]
            self.assertEqual(len(large), 1)
            square, = large[0]
            self.assertEqual(_side_sq(square), expected)

    def test_strict_crossing_terrain(self):
        roof = ((0,0), (3,0), (3,5), (1,5), (1,2), (0,2))
        square, = dict(chart_cells(roof))[0,1,3,5]
        self.assertEqual(square, tuple(pt(*p) for p in ((2,0),(3,2),(1,3),(0,1))))
        self.assertEqual(_side_sq(square), 5)


if __name__ == '__main__':
    unittest.main()
