"""Exact regressions for the refinement controls."""

from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from relative_witness_refinement_v282_v283 import (
    ARTIFACT, diagonal_dot, notch_rectangle, report, spiral_audit, truncation_audit,
)


class RelativeWitnessRefinementTests(unittest.TestCase):
    def test_exact_record(self):
        self.assertEqual(report(), json.loads(ARTIFACT.read_text(encoding="utf-8")))

    def test_notch_flat_sign_and_square_apex(self):
        eps = Q(1, 32)
        for h in range(4):
            bottom = h*eps/4
            for k in range(h, 5):
                top = k*eps/4
                corners = notch_rectangle(32, bottom, top)
                self.assertEqual(diagonal_dot(corners), (top-bottom)**2-eps**2)
                self.assertEqual(diagonal_dot(corners) == 0, h == 0 and k == 4)
        with self.assertRaises(ValueError):
            notch_rectangle(32, 0, 2*eps)

    def test_spiral_frame_and_square(self):
        result = spiral_audit()
        self.assertEqual(result["flat_dot_products"], ["4", "-1/64"])
        self.assertEqual(result["squared_radius_lower_on_entire_witness"], "1/512")
        self.assertEqual(Q(result["one_shared_square_squared_side"]), Q(101252, 58081))

    def test_fixed_contacts_survive_deeper_truncation(self):
        for counterwind in (False, True):
            record = truncation_audit(5, counterwind)
            self.assertGreater(record["changed_arc_vertices"], 2000)
