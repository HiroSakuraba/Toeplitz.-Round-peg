"""Check the phase envelopes and failure controls under SageMath."""

import importlib.util
import json
from math import ceil, floor
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The v235 interval tests run in the SageMath workflow")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from sage.all import QQ
from rough_spiral_family_v235 import (
    MANIFEST, _scaled_tent, phase_enclosure, replay, tent_chord_error, verify_rough_box,
)


class RoughSpiralTests(unittest.TestCase):
    def test_tent_chord_error_includes_interior_half_integer_peak(self):
        a, b = QQ(1)/3, QQ(2)/3
        self.assertEqual(tent_chord_error(a, b, 1), (a, a, QQ(1)/6))

    def test_shortcut_matches_exhaustive_rational_knot_check(self):
        for a, b in ((QQ(1)/3, QQ(7)/10), (QQ(417)/1024, QQ(429)/1024),
                     (QQ(25)/64, QQ(27)/64)):
            for frequency in (1, 2, 32, 1024):
                with self.subTest(a=a, b=b, frequency=frequency):
                    start, end, error = tent_chord_error(a, b, frequency)
                    points = (QQ(j)/(2*frequency) for j in range(
                        ceil(2*frequency*a), floor(2*frequency*b)+1))
                    exhaustive = max((abs(_scaled_tent(t, frequency) -
                                          (start+(end-start)*(t-a)/(b-a)))
                                      for t in points), default=QQ(0))
                    self.assertEqual(exhaustive, error)

    def test_enclosure_rejects_range_crossing_cutoff(self):
        self.assertIsNone(phase_enclosure(QQ(1)/4, QQ(1)/2))

    def test_coverage_and_orientation_controls(self):
        rows = json.loads(MANIFEST.read_text(encoding="utf-8"))
        a, b, orientation, *box = rows[0]
        self.assertEqual(verify_rough_box(QQ(a), QQ(b), orientation, *box), "certified")
        self.assertNotEqual(verify_rough_box(QQ(a), QQ(b), -orientation, *box), "certified")
        with self.assertRaisesRegex(ValueError, "does not start at 0"):
            replay(rows[1:], verify_rough_box)
