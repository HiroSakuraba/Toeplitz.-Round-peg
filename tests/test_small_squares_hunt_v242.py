"""v242 small-squares hunt: finder validation, recorded-result integrity, controls.

Needs numpy and scipy. Set REQUIRE_NUMERIC_DEPS=1 to make missing packages a
failure instead of a skip (the dedicated CI job does this).
"""
import json
import math
import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
try:
    import numpy as np
    import small_squares_hunt_v242 as h
except ImportError:  # pragma: no cover - environment dependent
    h = None

REQUIRE = os.environ.get('REQUIRE_NUMERIC_DEPS') == '1'
SQRT2 = math.sqrt(2)


class RecordedData(unittest.TestCase):
    """Standard library only: the committed JSON is internally consistent."""

    def setUp(self):
        self.data = json.loads((ROOT / 'research' / 'v242-small-squares-hunt.json').read_text())

    def test_no_verified_record_below_sqrt2(self):
        vals = [r['verified_J'] for r in self.data['runs'] if r['verified_J'] is not None]
        self.assertGreaterEqual(len(vals), 12)
        self.assertGreaterEqual(min(vals), SQRT2 - 1e-5)

    def test_discarded_runs_are_the_overflow_ones(self):
        for r in self.data['runs']:
            if r['verified_J'] is None:
                self.assertTrue(r['role'].startswith('discarded'))
                self.assertGreater(r['log_radius_spread'], 100)

    def test_search_values_below_sqrt2_were_artefacts(self):
        for r in self.data['runs']:
            if r['verified_J'] is not None and r['search_J'] < SQRT2 - 1e-4:
                self.assertTrue(r['role'].startswith('resolution trap'))
                self.assertGreater(r['verified_J'], 1.6)

    def test_reference_shapes_and_families_stay_above_sqrt2(self):
        self.assertTrue(all(v >= SQRT2 - 1e-6 for v in self.data['reference_shapes'].values()))
        self.assertTrue(all(f['J'] >= SQRT2 - 1e-6 for f in self.data['families']))
        self.assertEqual(self.data['random_nonstar_sample']['below_sqrt2'], 0)


@unittest.skipIf(h is None and not REQUIRE, 'numpy/scipy not installed')
class Finder(unittest.TestCase):
    def setUp(self):
        if h is None:
            self.fail('REQUIRE_NUMERIC_DEPS=1 but numpy/scipy are missing')

    def test_ellipse_exact(self):
        e = h.ellipse(2.0, 1.0)
        sq, side = h.find_squares(e, 256)
        self.assertEqual(len(sq), 1)
        self.assertAlmostEqual(side, 4 / math.sqrt(5), places=4)
        self.assertAlmostEqual(h.inradius(e), 1.0, places=4)

    def test_regular_polygons_exact(self):
        self.assertAlmostEqual(h.verified_ratio(h.regular_polygon(3), high=False), 6 * (2 - math.sqrt(3)), places=3)
        self.assertAlmostEqual(h.verified_ratio(h.regular_polygon(6), high=False), 2 * math.sqrt(3) - 2, places=3)

    def test_two_parametrisations_agree_and_counts_are_odd(self):
        rng = np.random.default_rng(242)
        for _ in range(3):
            K = 6
            k = np.r_[np.arange(1, K + 1), np.arange(1, K + 1)]
            c = h.Radial.from_vector(rng.normal(0, 0.3, 2 * K) / k)
            sq1, s1 = h.find_squares(c, 256, 'diag')
            sq2, s2 = h.find_squares(c, 256, 'adj')
            self.assertAlmostEqual(s1, s2, places=8)
            self.assertEqual(len(sq1) % 2, 1)

    def test_polygon_and_radial_finders_agree(self):
        rng = np.random.default_rng(7)
        c = h.Radial.from_vector(rng.normal(0, 0.15, 8) / np.r_[1, 2, 3, 4, 1, 2, 3, 4])
        P = h.Polygon(c.gamma(np.linspace(0, 2 * np.pi, 1500, endpoint=False)))
        self.assertAlmostEqual(h.find_squares(c, 256)[1], h.find_squares(P, 256)[1], places=3)

    def test_small_perturbations_of_a_circle_stay_above_sqrt2(self):
        rng = np.random.default_rng(3)
        for _ in range(3):
            x = rng.normal(0, 1, 10)
            x *= 0.03 / np.abs(x).max()
            self.assertGreater(h.verified_ratio(h.Radial.from_vector(x), high=False), SQRT2 - 1e-6)


@unittest.skipIf(h is None and not REQUIRE, 'numpy/scipy not installed')
class Controls(unittest.TestCase):
    def setUp(self):
        if h is None:
            self.fail('REQUIRE_NUMERIC_DEPS=1 but numpy/scipy are missing')
        self.data = h.load_data()

    def test_floor_shape_touches_its_incircle_at_a_square(self):
        r = next(r for r in self.data['runs'] if r['role'] == 'floor shape')
        c = h.Radial.from_vector(r['x'])
        self.assertGreater(c.log_radius_spread(), 4)          # far from round
        sq, side = h.find_squares(c, 768)
        rad, cen = h.incircle(c, N=12000, G=250)
        self.assertAlmostEqual(side / rad, SQRT2, places=4)
        V = max(sq, key=lambda V: abs(V[1] - V[0]))
        self.assertTrue(np.allclose(np.abs(V - V.mean()), rad, rtol=1e-3))

    def test_resolution_trap_is_caught_by_recheck(self):
        r = next(r for r in self.data['runs'] if r['role'] == 'resolution trap')
        c = h.Radial.from_vector(r['x'])
        ir = h.inradius(c, N=12000, G=250)
        self.assertLess(h.find_squares(c, 192)[1] / ir, 1.2)
        self.assertGreater(h.find_squares(c, 768)[1] / ir, 1.6)


if __name__ == '__main__':
    unittest.main()
