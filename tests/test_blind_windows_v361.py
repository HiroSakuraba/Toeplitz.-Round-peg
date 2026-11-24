from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from blind_windows_v361 import Frame, find_windows, scan_frames, _clear_intervals
from free_contact_v360 import free_contact_certificate
from test_free_contact_v360 import fixture


def replay(vertices, hit):
    # Insert the automatically selected endpoints into the original polygon
    # and run v360's independently implemented Euclidean probe checker.
    n = len(vertices)
    ends = [(w.start%n,w.end%n) for w in hit.windows]
    cuts = sorted(set(range(n)).union(t for pair in ends for t in pair))
    u,v = hit.frame.direction
    f = hit.frame
    polygon = []
    for t in cuts:
        i = int(t)
        a,b = vertices[i],vertices[(i+1)%n]
        x,y = [Q(lo)+(t-i)*(Q(hi)-lo) for lo,hi in zip(a,b)]
        polygon.append(((u*x+v*y-f.left)/f.width,(u*y-v*x-f.bottom)/f.width))
    arcs = [(cuts.index(a),cuts.index(b)) for a,b in ends]
    return free_contact_certificate(polygon,arcs,1,hit.relative_error)


class BlindWindowTests(unittest.TestCase):
    def test_all_contact_exchanges_without_marks(self):
        for height in (12,11,Q(21,2)):
            for folded in (False,True):
                p,_ = fixture(height,folded)
                hit = find_windows(p,Q(1,2816),Q(777,100))
                self.assertEqual(hit.side_squared_floor,Q(17816841,247808))
                self.assertGreater(hit.side_squared_floor,Q(847,100)**2)
                self.assertIsNotNone(replay(p,hit))

    def test_fixed_frame_rediscovers_folded_arcs_after_registration_changes(self):
        p,_ = fixture(11,True)
        p = [(Q(x),Q(y)) for x,y in p]
        split = []
        for a,b in zip(p,p[1:]+p[:1]):
            split.extend((a,tuple((x+y)/2 for x,y in zip(a,b))))
        for polygon in (p,list(reversed(p)),split[12:]+split[:12]):
            hit = scan_frames(polygon,[Frame((1,0),0,0,11)],Q(1,2816))
            self.assertGreater(hit.side_squared_floor,Q(777,100)**2)
            self.assertIsNotNone(replay(polygon,hit))

    def test_similarity_and_miss(self):
        p,_ = fixture(12,True)
        # Rotate by a rational unit complex number, scale by 2, translate.
        moved = [(Q(6,5)*x-Q(8,5)*y+2,Q(8,5)*x+Q(6,5)*y-3) for x,y in p]
        hit = find_windows(moved,Q(1,2816),Q(169,10))
        self.assertGreater(hit.side_squared_floor,Q(169,10)**2)
        triangle = [(0,0),(4,0),(2,3)]
        self.assertIsNone(find_windows(triangle,0))
        # The equidistant upper contacts at height 12/7 are 12/7 apart.
        self.assertEqual(4-2*Q(2,3)*Q(12,7),Q(12,7))

    def test_tube_slice_includes_collinear_and_point_contacts(self):
        edges = [((Q(1,4),Q(1)),(Q(1,2),Q(1))),
                 ((Q(3,4),Q(2)),(Q(3,4),Q(1)))]
        self.assertEqual(list(_clear_intervals(edges,Q(0))),
                         [(0,Q(1,4)),(Q(1,2),Q(3,4)),(Q(3,4),1)])
