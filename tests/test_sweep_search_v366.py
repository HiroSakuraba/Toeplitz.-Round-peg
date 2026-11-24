from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from sweep_search_v366 import (  # noqa: E402
    RegisteredPolygon, certificate, search, symmetric_circle,
)
from three_patch_sweep_v364 import target_roof  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DECAGON = [tuple(map(Q, v)) for v in json.loads((ROOT / 'research/v247-refined-ten-gon.json').read_text())['vertices']]
L = Q(12, 7)


def is_square(sq):
    P, Qp, R, D = sq
    u = (Qp[0] - P[0], Qp[1] - P[1])
    return (R[0] - Qp[0], R[1] - Qp[1]) in ((-u[1], u[0]), (u[1], -u[0])) and \
        (D[0] - P[0], D[1] - P[1]) == (R[0] - Qp[0], R[1] - Qp[1])


class RoofRecovery(unittest.TestCase):
    poly = RegisteredPolygon(target_roof())

    def test_v364_labeling_is_found_without_marks(self):
        c = certificate(self.poly, (Q(1, 32), Q(1, 32)), (Q(3, 32), Q(1, 32)), (Q(14, 32), Q(1, 32)), 1)
        self.assertEqual(c['degree'], 1)
        x = Q(23, 28)
        self.assertEqual(c['squares'], [((x, 0), (x + L, 0), (x + L, L), (x, L))])

    def test_blind_search_finds_the_same_square(self):
        res = search(self.poly, bins=16, spans=(1,), gap=1)
        best = res[0]
        self.assertEqual((best['arcs'], best['orient'], best['side_sq_floor']), ((1, 7, 11), 1, L * L))
        x = Q(23, 28)
        self.assertEqual(set(best['squares'][0]), {(x, 0), (x + L, 0), (x + L, L), (x, L)})
        for cert in res:
            for sq in cert['squares']:
                self.assertTrue(is_square(sq))
                self.assertTrue(all(self.poly.side_of(p) is None for p in sq))


class Decagon(unittest.TestCase):
    poly = RegisteredPolygon(DECAGON)

    def test_best_certificate_and_its_square(self):
        res = search(self.poly, bins=16, spans=(1, 2), gap=1)
        self.assertEqual(len(res), 20)
        best = res[0]
        self.assertEqual((best['arcs'], best['span'], best['orient']), ((7, 3, 10), 1, -1))
        self.assertGreater(best['side_sq_floor'], Q(196, 100))       # side floor > 1.4
        (sq,) = best['squares']
        self.assertTrue(is_square(sq))
        self.assertTrue(all(self.poly.side_of(p) is None for p in sq))

    def test_persists_under_a_small_registered_redraw(self):
        arcs = [(Q(x, 16), Q(1, 16)) for x in (7, 3, 10)]
        moved = [(x + Q((-1) ** i, 4096), y + Q(i % 3 - 1, 4096)) for i, (x, y) in enumerate(DECAGON)]
        c = certificate(RegisteredPolygon(moved), *arcs, -1)
        self.assertIsNotNone(c)
        self.assertEqual(c['degree'], 1)


class CircleControl(unittest.TestCase):
    def test_quarter_turn_symmetric_polygons_have_no_certificate(self):
        for m in (2, 3, 4):
            self.assertEqual(search(RegisteredPolygon(symmetric_circle(m)), bins=16, spans=(1, 2)), [])

    def test_one_bump_on_the_octagon_gives_clearance_proportional_to_it(self):
        clear = []
        for bump in (Q(1, 64), Q(1, 256), Q(1, 1024)):
            res = search(RegisteredPolygon(symmetric_circle(2, bump=bump)), bins=16, spans=(1, 2))
            self.assertEqual(len(res), 8)
            clear.append(max(c['clearance_sq'] for c in res))
        # quadrupling the bump multiplies the clearance by four, up to five percent
        for a, b in zip(clear, clear[1:]):
            self.assertLess(abs(a / b - 16), Q(16 * 11, 100))


if __name__ == '__main__':
    unittest.main()
