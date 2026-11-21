from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import sys
import unittest
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from accumulating_notches_v220 import notch
from separated_arc_parity_v311 import detecting_windows, four_cut, notch_certificate, notched_census, pt


class ArcTests(unittest.TestCase):
    def test_all_single_cycle_placements(self):
        # Include windows with zero, one, two, three or four contacts.
        for tags in combinations(range(1,16,2),4):
            four_cut([tags], (0,6,8,14), 16)
            four_cut([tags], (12,2,4,10), 16)

    def test_cyclic_relabeling_and_many_local_squares(self):
        cuts = (0, 2, 4, 6)
        square = (1,5,7,9)
        counts = four_cut([square],cuts,20)
        self.assertEqual(counts["bit"],1)
        extra = [(7,8,9,10), (9,10,11,12)]
        self.assertEqual(four_cut([square,*extra],cuts,20)["bit"],1)
        shifted = lambda row: tuple((v+7)%20 for v in row)
        self.assertEqual(four_cut([shifted(square)],shifted(cuts),20),counts)

    def test_even_shared_source_edge(self):
        cuts = (0,2,4,6)
        # Two cycles use the same consecutive contacts 1 and 5.
        self.assertEqual(four_cut([(1,5,7,9),(1,5,8,10)],cuts,12)["bit"],0)

    def test_strip_exposes_an_edge_even_with_shared_contacts(self):
        points = [pt(0,0),pt(1,0),pt(2,0),pt(2,1),pt(1,1),pt(0,1)]
        for order in ([0,1,2,3,4,5], [0,1,2,3,5,4]):
            squares = [SimpleNamespace(vertices=[points[i] for i in ids],
                                       curve_positions=[Q(order[i]) for i in ids])
                       for ids in ([0,1,4,5],[1,2,3,4])]
            result = detecting_windows(squares,6)
            self.assertEqual(result["bit"],1)
            self.assertGreaterEqual(Q(result["contact_distance_squared"]),1)

    def test_large_side_canceled_by_smaller_diagonal(self):
        a,b,c,d = pt(0,0),pt(2,0),pt(2,2),pt(0,2)
        e,f = pt(1,1),pt(1,-1)
        large = SimpleNamespace(vertices=(a,b,c,d),curve_positions=(0,1,4,5))
        # A and B occur consecutively along the source, as do E and F.
        small = SimpleNamespace(vertices=(a,f,b,e),curve_positions=(0,3,1,2))
        result = detecting_windows([large,small],6)
        self.assertEqual(result["contact_distance_squared"],"4")
        self.assertNotEqual(result["contact_pair"],["0","1"])

    def test_cut_at_corner_and_interleaving_fail(self):
        with self.assertRaises(ValueError):
            four_cut([(0,1,2,3)],(0,2,4,6),8)
        with self.assertRaises(ValueError):
            four_cut([], (0,4,2,6),8)

    def test_uniform_certificate_and_arbitrarily_deep_notches(self):
        self.assertEqual(notch_certificate()["side_floor"], "255/256")
        with self.assertRaises(ValueError):
            notch_certificate(Q(1,256))
        for n in (1,2,16,100):
            pts = notch(n)
            self.assertTrue(all(-Q(1,16)<=p.y<=0 and 0<p.x<=Q(9,16) for p in pts))

    def test_individual_entries_can_change(self):
        first, second = (notched_census(1,k) for k in (1,4))
        self.assertNotEqual(first["matrix"],second["matrix"])
        self.assertEqual((first["bit"],second["bit"]), (1,1))
