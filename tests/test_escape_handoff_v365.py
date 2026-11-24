from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from cutoff_interactions_v363 import CENTERS, SIDE, WIDTH, cutoff_chart_cover, moving_square  # noqa: E402
from escape_handoff_v365 import (  # noqa: E402
    ROOT, SEAM, crossings, cube_difference, four_patch_hit, levels, seam_events, stage_counts,
    union_bit_difference, union_separation,
)
from three_patch_sweep_v364 import completion_side, fixture_boxes, fixture_cutoff, target_roof  # noqa: E402


class EscapeAndReentry(unittest.TestCase):
    def test_single_crossing_moves_between_stages(self):
        self.assertEqual(seam_events(), {'t=1/3': Q(3, 7), 't=2/3': Q(1, 4)})
        table = {Q(0): Q(11, 15), Q(1, 4): SEAM, Q(3, 7): Q(1, 3), Q(1): Q(13, 42)}
        for lam, t in table.items():
            self.assertEqual(crossings(lam), (t,))
        for k in range(57):
            lam = Q(k, 56)
            self.assertEqual(len(crossings(lam)), 1)
            self.assertTrue(all(0 <= z <= 1 for z in levels(lam)))
        self.assertEqual(stage_counts(0), (0, 1))
        self.assertEqual(stage_counts(Q(1, 4)), (1, 1))
        self.assertEqual(stage_counts(1), (1, 0))

    def test_cutoff_side_at_the_root_and_sides_of_the_fourth_corner(self):
        for n in (0, 2, 5):
            e = Q(1, 2**18) / 2**n
            corners, sq = moving_square(7, ROOT, error=e)
            self.assertEqual(sq, fixture_cutoff(n) ** 2)
            self.assertEqual(corners[3], CENTERS[3])
        for n in (2, 3, 6):
            s = fixture_cutoff(n)
            self.assertIs(completion_side(target_roof(), fixture_boxes(n, 0), s), True)
            self.assertIs(completion_side(target_roof(), fixture_boxes(n, Q(1, 4)), s), True)
            self.assertIs(completion_side(target_roof(), fixture_boxes(n, Q(3, 4)), s), False)
            self.assertIs(completion_side(target_roof(), fixture_boxes(n, 1), s), False)
        self.assertIs(completion_side(target_roof(), fixture_boxes(3, Q(3, 7)), fixture_cutoff(3)), True)

    def test_stage_counts(self):
        # inside-portion counts: start, end, seam, as described in the note
        inside = lambda lam: lam < Q(1, 2)
        s1_start = sum(1 for t in crossings(0) if t <= SEAM)
        s2_start = sum(1 for t in crossings(0) if t >= SEAM)
        s1_end_inside = sum(1 for t in crossings(1) if t <= SEAM and inside(Q(1)))
        seam = [lam for lam in (Q(1, 4),) if SEAM in crossings(lam)]
        self.assertEqual((s1_start, s2_start, s1_end_inside), (0, 1, 0))
        self.assertEqual(seam, [Q(1, 4)])
        self.assertTrue(inside(seam[0]))
        # the chain: odd start, even end, seam interior
        self.assertEqual((s1_start + s2_start + s1_end_inside) % 2, 1)


class ChangedStretch(unittest.TestCase):
    def test_boolean_identity_for_every_parity_function(self):
        for code in range(1 << 16):
            p = lambda S: (code >> sum(1 << b for b in S)) & 1
            self.assertEqual(union_bit_difference(p, 0, 1, 2, 3),
                             cube_difference(p, (0, 1, 2)) ^ cube_difference(p, (0, 1, 3))
                             ^ cube_difference(p, (0, 1, 2, 3)))

    def test_triangle_swap_fails_separation(self):
        boxes = [(x - 3 * WIDTH, x + 3 * WIDTH, y - 3 * WIDTH, y + 3 * WIDTH) for x, y in CENTERS]
        s = fixture_cutoff(0)
        self.assertTrue(four_patch_hit(boxes, s))
        self.assertIsNone(union_separation(*boxes, s))
        p = lambda S: 0 if len(S) == 4 else 1   # v363: odd on proper subsets, even on the full set
        self.assertEqual((cube_difference(p, (0, 1, 2)), cube_difference(p, (0, 1, 3)),
                          union_bit_difference(p, 0, 1, 2, 3)), (0, 0, 1))

    def test_separated_patches_give_a_positive_bound(self):
        # the base two corners of the triangle's square and two candidate third stretches:
        # the true upper-right contact and a far patch no side-s square through the base can reach
        b = lambda c, r: (c[0] - r, c[0] + r, c[1] - r, c[1] + r)
        k1, k2, k3 = (b(CENTERS[i], 3 * WIDTH) for i in range(3))
        far = b((Q(16, 7), -SIDE), 3 * WIDTH)
        d2 = union_separation(k1, k2, k3, far, fixture_cutoff(0))
        self.assertIsNotNone(d2)
        self.assertGreater(d2, Q(1, 2))
        self.assertFalse(four_patch_hit([k1, k2, k3, far], fixture_cutoff(0)))

    def test_chart_cover_is_schedule_independent(self):
        ids = (0, 0, 1, 2)
        self.assertEqual(set(cutoff_chart_cover(SIDE, upper_cutoff=fixture_cutoff(0))),
                         {ids[k:] + ids[:k] for k in range(4)})


if __name__ == '__main__':
    unittest.main()
