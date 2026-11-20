"""Exact replay and independent angle check for the moving concave atlas."""

from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The v280 moving atlas uses SageMath")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from sage.all import QQ
from moving_rectangle_atlas_v280 import (
    ARTIFACT, END, certify_atlas, expected_bands, points_at, representative_charts,
)
from adversarial_rectangles_v278 import concave_polygon
from rectangle_atlas_v274 import independent_slice


class MovingConcaveAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.audit = certify_atlas()

    def test_projection_and_complete_replay(self):
        self.assertEqual(self.audit, json.loads(ARTIFACT.read_text(encoding="utf-8")))
        self.assertEqual(self.audit["ordered_contact_charts"], 168)
        self.assertEqual(self.audit["projection_event_times"], ["-1/100", "0", "1/100"])

    def test_relative_indices_at_the_fold(self):
        self.assertEqual(self.audit["geometric_squares_negative_zero_positive"], [1, 2, 3])
        self.assertEqual(self.audit["framed_relative_intervals_negative_zero_positive"],
                         [[1], [1, 1], [1, 1, 1]])
        self.assertEqual(self.audit["moving_components"][1]["regular_time_total_parity"], 0)

    def test_independent_non_square_angle_slice(self):
        time = END/2
        angle = QQ(1)/3
        polygon = tuple(tuple(QQ(v) for v in p) for p in concave_polygon(Fraction(str(time))))
        by_ids = {c.ids: c for c in representative_charts()}
        bands = expected_bands(time)
        actual = {tuple(sorted(points_at(by_ids[ids], time, x)))
                  for ids, ranges in bands.items()
                  for a, b, _, _ in ranges for x in (angle, 1-angle) if a <= x <= b}
        self.assertEqual(len(actual), 6)
        self.assertEqual(sorted(actual), sorted(independent_slice(polygon, angle)))
