from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from convex_shell_v307 import certify, covered, fixture, pt, replay


class ConvexShellTests(unittest.TestCase):
    def test_macroscopic_root(self):
        r = replay()
        self.assertEqual(Q(r["cutoff_polynomial"]),Q(-1,128))
        self.assertEqual(Q(r["large_root_discriminant"]),Q(321,1024))

    def test_notch_also_has_a_small_square(self):
        p = fixture(True)
        self.assertEqual(certify(p[::-1],pt(0,0),1,5,Q(1,32),Q(3,4)),replay())
        self.assertEqual(p[3].y-p[4].y,p[4].x-p[5].x)
        self.assertEqual(p[3].y-p[4].y,Q(1,32))

    def test_uncovered_segment_interior(self):
        hull = [pt(0,0),pt(1,0),pt(1,1),pt(0,1)]
        self.assertFalse(covered(pt(0,Q(1,2)),pt(1,Q(1,2)),hull,Q(1,32)))
        self.assertTrue(covered(pt(0,Q(1,8)),pt(Q(1,8),0),hull,Q(1,16)))

    def test_missing_enclosure(self):
        p = tuple(pt(x,y) for x,y in [(-4,-1),(4,-1),(4,1),(-4,1),(-4,Q(63,64)),
                                      (Q(255,64),Q(63,64)),(Q(255,64),Q(-63,64)),(-4,Q(-63,64))])
        hull = [pt(-4,-1), pt(4,-1), pt(4,1), pt(-4,1)]
        for a,b in zip(p,p[1:]+p[:1]):
            self.assertTrue(covered(a,b,hull,Q(1,32)))
        with self.assertRaisesRegex(ValueError,"enclose"):
            certify(p,pt(0,0),1,5,Q(1,32),Q(3,4))

    def test_equality_is_inconclusive(self):
        with self.assertRaisesRegex(ValueError,"strict"):
            certify(fixture(),pt(0,0),1,5,Q(3,92),Q(3,4))

    def test_small_circle_in_unit_shell(self):
        d = Q(1,32)
        c, radius = 1-d/2, d/4
        self.assertTrue(1-d <= c-radius < c+radius <= 1)
        self.assertLess(radius,Q(3,4))
