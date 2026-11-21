from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from polygon_orbit_census_v223 import _affine_system, _evaluate, _solve, add, dot, mul, sub
from rough_crossing_v302 import ORIGIN, PACKET, SIDE, T_MIN, anchor, annulus, orientation_gate, packet_start
from tilted_sweep_v301 import seed


class RoughCrossing(unittest.TestCase):
    def test_margins(self):
        self.assertGreater(annulus()["side_squared"], SIDE**2)
        self.assertEqual(orientation_gate()["excess"], Q(131, 8000))
        self.assertEqual(anchor()["sorted_charts"], 13)

    def test_actual_endpoint_rectangles_on_both_shears_and_zero(self):
        # Independent contact search, with no call to winning_tags or fixed winning chart.
        for d, k in product((Q(-1, 100), Q(0), Q(1, 100)), (T_MIN, 1/T_MIN)):
            polygon = seed(d)
            edges = tuple(zip(polygon, polygon[1:]+polygon[:1]))
            for chart in product(edges, repeat=4):
                matrix, rhs = _affine_system(chart, k)
                tags = _solve(matrix, rhs)
                if tags is None or not all(0 <= u <= 1 for u in tags):
                    continue
                a, b, c, e = _evaluate(chart, tags)
                v = sub(b, a)
                if v.x <= abs(v.y):
                    continue
                center = mul(add(a, c), Q(1, 2))
                s, h = sub(center, ORIGIN), sub(center, a)
                if dot(s, s) >= dot(h, h):
                    continue
                jv = sub(c, b)
                self.assertEqual(sub(e, a), jv)
                self.assertEqual(dot(v, jv), 0)
                self.assertEqual(dot(jv, jv), k*k*dot(v, v))
                self.assertGreater(min(dot(v, v), dot(jv, jv)), SIDE**2)
                break
            else:
                self.fail(f"no actual endpoint witness at d={d}, k={k}")

    def test_packet_mass_on_each_side_of_the_seam(self):
        lengths = (Q(1, 5), Q(3, 20), Q(2, 5), Q(1, 4))
        q = lengths[0]+lengths[2]
        for origin in (Q(0), Q(49, 50), Q(81, 100), Q(1, 3)):
            s1 = (origin+lengths[0]) % 1
            phase = packet_start(s1, q)
            red = [(origin, origin+lengths[0]),
                   (origin+sum(lengths[:2]), origin+sum(lengths[:3]))]
            overlap = sum(max(0, min(b, phase+shift+PACKET)-max(a, phase+shift))
                          for a, b in red for shift in (-1, 0, 1, 2))
            self.assertEqual(q/2+6*overlap, Q(1, 2))
            self.assertEqual(packet_start(s1+1, q), phase)
