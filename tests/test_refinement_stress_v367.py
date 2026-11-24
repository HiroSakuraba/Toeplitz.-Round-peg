from fractions import Fraction as Q
from pathlib import Path
import json
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from refinement_stress_v367 import bump_across, run, weakest_endpoint  # noqa: E402
from sweep_search_v366 import RegisteredPolygon, certificate, search  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DECAGON = [tuple(map(Q, v)) for v in json.loads((ROOT / 'research/v247-refined-ten-gon.json').read_text())['vertices']]


class Adversary(unittest.TestCase):
    def test_bump_flips_the_weakest_fourth_corner_and_keeps_registration(self):
        poly = RegisteredPolygon(DECAGON)
        best = search(poly, bins=16, spans=(1, 2))[0]
        c2, D = weakest_endpoint(poly, best, 16)
        new = bump_across(poly, D)
        self.assertIsNotNone(new)
        self.assertNotEqual(new.side_of(D), poly.side_of(D))
        self.assertEqual(len(new.V), len(poly.V) + 3)
        # original vertices keep their times; inserted ones lie inside their edge's interval
        self.assertTrue(set(poly.times) <= set(new.times))
        # the attacked certificate dies, its neighbour across the attacked face becomes odd, and the
        # certificate on the union of the two Q-bins keeps its degree and its clearance
        arc = lambda b, n=1: (Q(b, 16), Q(n, 16))
        self.assertEqual(best['arcs'], (7, 3, 10))
        self.assertEqual(certificate(new, arc(7), arc(3), arc(10), -1)['degree'], 0)
        self.assertEqual(certificate(new, arc(7), arc(2), arc(10), -1)['degree'], 1)
        u0 = certificate(poly, arc(7), arc(2, 2), arc(10), -1)
        u1 = certificate(new, arc(7), arc(2, 2), arc(10), -1)
        self.assertEqual((u0['degree'], u1['degree']), (1, 1))
        self.assertEqual(u0['clearance_sq'], u1['clearance_sq'])

    def test_two_stages_of_each_policy(self):
        for policy in ('best', 'robust'):
            recs, poly = run(RegisteredPolygon(DECAGON), 2, policy=policy)
            self.assertEqual([r['certificates'] for r in recs][0], 20)
            self.assertTrue(all(r['certificates'] > 0 for r in recs))
            self.assertLess(recs[1]['target_clearance_sq'], recs[0]['target_clearance_sq'])

    def test_recorded_greedy_push(self):
        data = json.loads((ROOT / 'results' / 'v367' / 'decagon_greedy_push.json').read_text())
        self.assertEqual(data['path'], [[7, 3, 10], [7, 2, 10], [6, 2, 10], [6, 2, 11], [6, 1, 11], [6, 1, 12]])
        self.assertEqual(data['union_before']['degree'], 1)
        self.assertEqual(data['union_after']['degree'], 1)
        self.assertEqual(data['union_before']['clearance_sq'], data['union_after']['clearance_sq'])
        final = RegisteredPolygon([tuple(map(Q, v)) for v in data['final_vertices']], data['final_times'])
        union = [tuple(map(Q, a)) for a in data['union']]
        again = certificate(final, *union, -1)
        self.assertEqual((again['degree'], str(again['clearance_sq'])), (1, data['union_after']['clearance_sq']))
        last = [(Q(b, 16), Q(1, 16)) for b in data['path'][-1]]
        self.assertEqual(certificate(final, *last, -1)['degree'], 1)

    def test_recorded_runs(self):
        for name in ('decagon_attack_best', 'decagon_attack_robust'):
            f = ROOT / 'results' / 'v367' / f'{name}.json'
            data = json.loads(f.read_text())
            stages = data['stages']
            self.assertTrue(all(s['certificates'] > 0 for s in stages))
            floors = [Q(s['best_floor_sq']) for s in stages]
            self.assertGreaterEqual(min(floors), Q(196, 100))   # the best side floor stays above 1.4
            margins = [Q(s['best_margin_sq']) for s in stages]
            self.assertLess(margins[-1], margins[0])


if __name__ == '__main__':
    unittest.main()
