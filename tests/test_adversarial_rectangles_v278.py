"""Exact controls for v278's falsified candidate hypotheses."""

import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from adversarial_rectangles_v278 import MANIFEST, concave_polygon, moving_notch, report
from polygon_orbit_census_v223 import pt, validate_jordan_polygon


class AdversarialLedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = report()
        cls.cases = {case["id"]: case["evidence"] for case in cls.report["specimens"]}

    def test_exact_artifact_replays_and_keeps_general_claim_open(self):
        canonical = json.dumps(self.report, indent=2, sort_keys=True) + "\n"
        self.assertEqual(MANIFEST.read_text(encoding="utf-8"), canonical)
        self.assertEqual(self.report["unrestricted_toeplitz"], "OPEN")
        self.assertTrue(all(case["status"] == "DISPROVED_CANDIDATE"
                            for case in self.report["specimens"]))

    def test_even_raw_singular_count_and_independent_orbit_sizes(self):
        case = self.cases["concave-vertex-pair"]
        samples = case["independent_whole_polygon_slices"]
        self.assertEqual([s["orbit_count"] for s in samples], [1, 2, 3])
        geom = samples[1]["geometric_orbits"]
        self.assertEqual([o["size_squared"] for o in geom], ["1", "145/121"])
        self.assertEqual([o["support_labels"] for o in geom], [32, 8])
        self.assertTrue(all(o["diagonal_pairing"] == [[0, 2], [1, 3]] for o in geom))
        self.assertTrue(all(o["source_order_sector"] == [0, 1, 2, 3] for o in geom))
        self.assertTrue(geom[0]["vertex_contact"])
        self.assertFalse(geom[1]["vertex_contact"])

    def test_two_selected_branches_have_distinct_contact_patterns(self):
        case = self.cases["concave-vertex-pair"]
        self.assertEqual(case["selected_branch_counts_negative_zero_positive"], [0, 1, 2])
        self.assertEqual([b["edge_choices"] for b in case["exact_selected_branches"]],
                         [(3, 0, 2, 3), (0, 0, 2, 3)])
        self.assertEqual(case["geometric_double_vertex_event_square"],
                         [["0", "0"], ["1", "0"], ["1", "1"], ["0", "1"]])

    def test_reject_outside_window_and_nonjordan_polygon(self):
        with self.assertRaisesRegex(ValueError, "outside"):
            concave_polygon(Q(1, 10))
        crossing = tuple(pt(x, y) for x, y in ((0, 0), (3, 0), (1, 2), (2, -1)))
        with self.assertRaisesRegex(ValueError, "nonadjacent"):
            validate_jordan_polygon(crossing)

    def test_scalar_closure_and_finite_continuum_are_rejected(self):
        closure = self.cases["false-closure"]
        self.assertEqual(closure["both_squared_radii"], "25")
        self.assertNotEqual(closure["start"], closure["end"])
        self.assertIsNone(self.cases["circle-continuum"]["square_zero_count"])
        self.assertIsNone(self.cases["bump-cell"]["entire_square_orbit_count"])
        self.assertFalse(self.cases["bump-cell"]["whole_rectangle_atlas_verified_here"])

    def test_collapse_and_independent_large_witness_coexist(self):
        fine = moving_notch(64)
        self.assertEqual(fine["uniform_coordinate_error_upper"], "1/32")
        self.assertEqual(fine["selected_square_side_squared"], "1/4096")
        self.assertEqual(fine["other_square_side_squared"], "1")
        collision = self.cases["fixed-collision-link"]
        self.assertEqual(collision["gap_direction"], ["1/3"]*3)
        self.assertNotEqual(collision["square_limit"], collision["flat_control_limit"])


if __name__ == "__main__":
    unittest.main()
