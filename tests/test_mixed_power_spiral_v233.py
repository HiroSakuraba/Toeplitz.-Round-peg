"""Regression controls for exact mixed-power square witnesses."""

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from mixed_power_spiral_v233 import FLOORS, LAWS, WITNESSES, certify, discover, report
from spiral_pruned_v231 import _side_floor_holds


class MixedPowerSpiralTests(unittest.TestCase):
    def test_three_exact_witnesses_replay_without_discovery(self):
        with patch("mixed_power_spiral_v233.discover", side_effect=AssertionError):
            results = report()
        self.assertEqual(set(results), set(LAWS))
        for name, law in LAWS.items():
            with self.subTest(name=name):
                proof = certify(law, *WITNESSES[name], *FLOORS[name])
                self.assertTrue(proof.certified, proof.status)
                self.assertLess(4*proof.coordinate_error, proof.boundary_margin)

    def test_impossible_size_floors_rejected(self):
        for name, law in LAWS.items():
            orientation, endpoints = WITNESSES[name]
            self.assertEqual(certify(law, orientation, endpoints, Q(2)).status,
                             "insufficient-side-bound")
            self.assertEqual(certify(law, orientation, endpoints,
                                     radius_floor=Q(1)).status, "insufficient-radius-bound")

    def test_wrong_orientation_and_different_law(self):
        first, second, _ = LAWS
        orientation, endpoints = WITNESSES[first]
        self.assertFalse(certify(LAWS[first], -orientation, endpoints).certified)
        self.assertFalse(certify(LAWS[second], orientation, endpoints).certified)

    def test_invalid_or_overlapping_intervals_rejected(self):
        name = next(iter(LAWS))
        orientation, endpoints = WITNESSES[name]
        with self.assertRaises(ValueError):
            certify(LAWS[name], orientation, endpoints[:2]*4)
        with self.assertRaises(ValueError):
            certify(LAWS[name], orientation, (0, *endpoints[1:]))

    def test_side_bound_accounts_for_movement_along_chords(self):
        starts = ((Q(0), Q(0)), (Q(1), Q(0)))
        ends = ((Q(1), Q(0)), (Q(2), Q(0)))
        root = (Q(1, 2), Q(1, 2))
        self.assertFalse(_side_floor_holds(starts, ends, root, Q(0), Q(1, 2)))

    def test_sample_and_pair_budget_exhaustion_are_inconclusive(self):
        law = next(iter(LAWS.values()))
        sample_limited = discover(law, radii=(.55,), max_segments=1)
        self.assertEqual(sample_limited["status"], "inconclusive")
        self.assertEqual(sample_limited["ledger"][0]["status"], "sample-budget")
        pair_limited = discover(law, radii=(.55,), max_pairs=1)
        self.assertEqual(pair_limited["status"], "inconclusive")
        self.assertEqual(pair_limited["ledger"][0]["pairs"], 1)
        self.assertTrue(pair_limited["ledger"][0]["pair_budget_exhausted"])
