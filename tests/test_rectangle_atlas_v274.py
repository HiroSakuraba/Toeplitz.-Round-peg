"""End-to-end atlas coverage, geometry and independent slice checks."""

from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The v274 algebraic atlas uses SageMath")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from sage.all import QQ
from rectangle_atlas_v274 import (
    ARTIFACT, AtlasInconclusive, artifact, atlas_slice, build_atlas, independent_slice,
)


class RectangleAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = build_atlas()
        cls.result = artifact(cls.data)

    def test_complete_chart_count_and_replay(self):
        self.assertEqual(self.result, json.loads(ARTIFACT.read_text(encoding="utf-8")))
        self.assertEqual(self.result["ordered_charts"], 168)
        self.assertEqual(len(self.result["squares"]), 2)

    def test_ratio_census_at_rational_and_irrational_angles(self):
        algebraic = min(x for (x, _), (_, kind) in self.data["nodes"].items()
                        if kind == "genuine" and x.minpoly().degree() > 1)
        for x in (QQ(1)/3, QQ(2)/3, algebraic):
            with self.subTest(x=x):
                expected = independent_slice(self.data["points"], x)
                self.assertTrue(expected)
                self.assertEqual(atlas_slice(self.data, x), expected)

    def test_explicit_gluing_includes_source_vertex_and_flat_end(self):
        incidence = Counter(i for b in self.result["branches"] for i in b["ends"])
        self.assertTrue(any(n["kind"] == "two_pair_flat" and incidence[i] > 0
                            for i, n in enumerate(self.result["nodes"])))
        self.assertTrue(any(n["kind"] == "square" and "2" in n["parameters"]
                            and incidence[i] > 0 for i, n in enumerate(self.result["nodes"])))

    def test_parallel_edge_continuum_fails_closed(self):
        with self.assertRaisesRegex(AtlasInconclusive, "parallel"):
            build_atlas(((0, 0), (3, 0), (3, 2), (0, 2)))
