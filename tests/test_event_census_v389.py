from pathlib import Path
import json
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from event_census_v389 import continue_branch, half_turn_partner, pair_events, seeds_at  # noqa: E402
from smooth_cores_v386 import sector_two_arm  # noqa: E402
from twisted_classification_v387 import twisted_corners  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


class Mechanics(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(5)
        self.found = []
        while not self.found:
            self.found = seeds_at(0.0, 0.02, 40, rng, (0.03, 0.15))

    def test_partner_is_the_half_turn_image(self):
        lam, w, ks = self.found[0]
        wp, ksp = half_turn_partner(w, ks)
        self.assertAlmostEqual(wp[3], w[3], places=10)
        Z, Zp = twisted_corners(w), twisted_corners(wp)
        self.assertAlmostEqual(abs(sum(Z)), abs(sum(Zp)), places=10)
        self.assertEqual(sector_two_arm(list(ks), Z), sector_two_arm(list(ksp), Zp))

    def test_branch_moves_little(self):
        lam, w, ks = self.found[0]
        br = continue_branch(0.0, w, ks, steps=20)
        sides = [x[1][3] for x in br]
        self.assertEqual(br[0][0], 0.0)
        self.assertLess((max(sides) - min(sides)) / sides[0], 0.3)

    def test_pair_events(self):
        legs = [dict(sector="graceful", sides=[1.0, 1.1]), dict(sector="graceful", sides=[1.0, 0.95])]
        self.assertEqual(pair_events([dict(legs=legs)], 1.05), (1, 0))
        self.assertEqual(pair_events([dict(legs=legs)], 0.97), (1, 0))
        self.assertEqual(pair_events([dict(legs=legs)], 1.2), (0, 0))


class Results(unittest.TestCase):
    def test_results_file(self):
        d = json.loads((ROOT / "results" / "v389" / "event_census.json").read_text())
        for s0, v in d.items():
            c = v["combined"]
            self.assertEqual(c["partners_same_sector"], c["pairs"])
            self.assertLess(c["move_median"], 0.01)
            self.assertEqual(v["distinct_pairs_observed"], c["pairs"])
            self.assertEqual(v["complete_pairs_observed"] + v["incomplete_pairs_observed"], c["pairs"])
            self.assertNotIn("population_lower_bound", v)          # the product of sample sizes is not a bound
            if v["overlap"] == 0:
                self.assertIsNone(v["lincoln_petersen_heuristic"])
            self.assertGreater(v["population_scale_heuristic"], 1000)
            self.assertGreater(v["estimated_events_per_cutoff_heuristic"], 10)
        self.assertGreater(d["0.02"]["combined"]["odd_fraction"], 0.3)


if __name__ == "__main__":
    unittest.main()
