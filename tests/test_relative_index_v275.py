"""Tests of the relative sign system and the certified graph class."""

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from relative_index_v275 import ARTIFACT, local_protection, quotient_index, relative_obstructions


class RelativeIndexTests(unittest.TestCase):
    def test_interval_distinguishes_boundary_signs(self):
        edges = [(0, 1, 0), (1, 2, 0)]
        self.assertIsNone(relative_obstructions(range(3), edges, {0: 0, 2: 0})[0]["witness"])
        self.assertEqual(relative_obstructions(range(3), edges, {0: 0, 2: 1})[0]["witness"],
                         {"kind": "relative_boundary", "path": [0, 1, 2]})

    def test_odd_cycle_and_frame_change(self):
        edges = [(0, 1, 0), (1, 2, 0), (2, 0, 1)]
        changed = [(a, b, bit ^ (a == 1) ^ (b == 1)) for a, b, bit in edges]
        for transitions in (edges, changed):
            witness = relative_obstructions(range(3), transitions, {})[0]["witness"]
            self.assertEqual(witness["kind"], "odd_monodromy")
            self.assertEqual(witness["path"][0], witness["path"][-1])
        self.assertIsNone(relative_obstructions(range(3), [(a, b, 0) for a, b, _ in edges],
                                                {})[0]["witness"])

    def test_atlas_class_has_a_protected_square(self):
        atlas = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        info = quotient_index(atlas)
        witnesses = [c["witness"] for c in info["components"] if c["witness"]]
        self.assertEqual(len(witnesses), 2)
        self.assertTrue(all(w["square_nodes"] for w in witnesses))
        self.assertTrue(all(Fraction(w["squared_radius_floor"]) > Fraction(9, 20) for w in witnesses))
        box = local_protection(atlas, info)
        self.assertEqual(box["node"], 337)
        self.assertEqual(box["edge_ids"], [0, 1, 2, 3])

    def test_invalid_relative_boundary_rejected(self):
        with self.assertRaisesRegex(ValueError, "relative boundary"):
            relative_obstructions((0, 1), [(0, 1, 1)], {2: 1})
