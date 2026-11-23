from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"code"))
from polygon_component_v341 import PolygonComponents
from polygon_orbit_census_v223 import pt
from separated_arcs_v359 import certified_pairs

OFFSET = Q(1, 128)
NOTCH = ((0,0),(10,0),(10,9),(9,9),(9,8),(5,8),(5,9),(2,9),(2,4),(1,4),(1,9),(0,9))
STAIRS = ((0,0),(15,0),(15,12),(11,12),(11,4),(8,4),(8,12),(5,12),(5,9),(2,9),(2,12),(0,12))
FIFTH = ((6,0),(11,6),(5,11),(0,5))


def knots(poly):
    return [Q(i,len(poly)) for i in range(len(poly)+1)]


def fold_contacts(tiles):
    # Re-register each replacement inside the same original edge interval.
    poly, times = [], []
    contacts = dict(zip((0,3,7,11), FIFTH))
    for i, a in enumerate(STAIRS):
        b = STAIRS[(i+1) % len(STAIRS)]
        poly.append(a)
        times.append(Q(i,12))
        if i not in contacts:
            continue
        axis = int(a[0] == b[0])
        length = abs(b[axis]-a[axis])
        direction = (b[axis]-a[axis])//length
        mid = direction*(contacts[i][axis]-a[axis])
        pattern = ((0,0),(3,0),(3,1),(2,1),(2,2),(4,2),(4,0))
        for j in range(tiles):
            for k,(along,across) in enumerate(pattern):
                s = mid-Q(1,2)+(j+Q(along,5))/tiles
                point = list(a)
                point[axis] += direction*s
                point[1-axis] += Q(across,128*tiles)
                poly.append(tuple(point))
                registered = mid-Q(1,2)+(j+Q(k,7))/tiles
                times.append((i+registered/length)/12)
        point = list(a)
        point[axis] += direction*(mid+Q(1,2))
        poly.append(tuple(point))
        times.append((i+(mid+Q(1,2))/length)/12)
    return poly, times+[Q(1)]


class ArcTests(unittest.TestCase):
    def test_search_tie_and_decagon_under_local_edits(self):
        dec = json.loads((ROOT/'research/v247-refined-ten-gon.json').read_text())
        for poly, seed, span, pair in (
            (NOTCH, ((2,7),(3,0),(10,1),(9,8)), 2, (5,9)),
            (dec['vertices'], None, 1, (3,7)),
        ):
            poly = [tuple(map(Q,p)) for p in poly]
            engine = PolygonComponents(poly)
            seed = (engine.cell(tuple(dec['claimed_exact_invariants']['maximizing_edge_assignment'])).squares[0]
                    if seed is None else tuple(pt(*p) for p in seed))
            cells = engine.trace(seed)
            times = knots(poly)
            if span == 2:
                self.assertFalse(certified_pairs(cells,times,16,OFFSET,1))
            self.assertIn(pair, certified_pairs(cells,times,16,OFFSET,span))
            untouched = next(i for i in range(len(poly)) if all(i not in ids for ids in cells))
            a, b = poly[untouched], poly[(untouched+1) % len(poly)]
            dx, dy = b[0]-a[0], b[1]-a[1]
            for scale in (128, 512):
                kink = ((a[0]+b[0])/2-dy/scale, (a[1]+b[1])/2+dx/scale)
                index = untouched+1
                refined = poly[:index]+[kink]+poly[index:]
                registered = times[:index]+[(times[index-1]+times[index])/2]+times[index:]
                self.assertIn(pair, certified_pairs(PolygonComponents(refined).trace(seed),
                                                    registered,16,OFFSET,span))

    def test_four_folded_contacts(self):
        seed = tuple(pt(*p) for p in FIFTH)
        for poly, times in [(STAIRS, knots(STAIRS)), fold_contacts(2), fold_contacts(4)]:
            self.assertIn((4,15), certified_pairs(PolygonComponents(poly).trace(seed),
                                                  times,16,OFFSET))
        lost = tuple(pt(*p) for p in ((4,0),(9,4),(5,9),(0,5)))
        self.assertFalse(certified_pairs(PolygonComponents(STAIRS).trace(lost),
                                        knots(STAIRS),16,OFFSET,2))

    def test_pair_stays_fixed_when_witness_changes(self):
        seed = tuple(pt(*p) for p in ((4,0),(9,4),(5,9),(0,5)))
        for g in (Q(89,10), Q(9), Q(91,10)):
            poly = list(STAIRS)
            poly[8], poly[9] = (5,g), (2,g)
            chosen = seed if g < 9 else tuple(pt(*p) for p in FIFTH)
            self.assertIn((0,9), certified_pairs(PolygonComponents(poly).trace(chosen),
                                                knots(poly),16,OFFSET,2))
