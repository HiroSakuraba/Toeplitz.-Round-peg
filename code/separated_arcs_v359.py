"""Fixed coarse-arc certificates on point/segment cells of the v341 tracer."""
from fractions import Fraction as Q
from itertools import combinations
from math import floor

from polygon_component_v341 import UnsupportedComponent, degree


def compulsory_bins(cells, knots, count=16, offset=Q(0), span=1):
    """Arcs of span adjacent bins whose OPEN interiors are visited by every square in the component.

    knots give persistent traversal times on the polygon, one per vertex plus 1.
    A bin j begins at offset+j/count modulo one.
    """
    knots = tuple(map(Q, knots))
    if (not 1 <= span < count or count < 4 or not cells or not knots or knots[0] != 0 or knots[-1] != 1
            or any(a >= b for a, b in zip(knots, knots[1:]))):
        raise ValueError("need traced cells, increasing unit-circle knots and at least four bins")
    common = set(range(count))
    for edges, cell in cells.items():
        if len(cell.parameters) not in (1, 2):
            raise UnsupportedComponent("coarse label test requires point or segment cells")
        tracks = []
        events = {Q(0), Q(1)}
        for edge, a, b in zip(edges, cell.parameters[0], cell.parameters[-1]):
            width = knots[edge+1]-knots[edge]
            start = count*(knots[edge]+a*width-Q(offset))
            change = count*(b-a)*width
            tracks.append((start, change))
            if change:
                low, high = sorted((start, start+change))
                for boundary in range(floor(low)+1, floor(high)+1):
                    events.add((boundary-start)/change)
        ordered = sorted(events)
        samples = ordered + [(a+b)/2 for a, b in zip(ordered, ordered[1:])]
        for t in samples:
            positions = [start+t*change for start, change in tracks]
            common.intersection_update(j for j in tuple(common)
                                       if any(0 < (x-j) % count < span for x in positions))
    return common


def certified_pairs(cells, knots, count=16, offset=Q(0), span=1):
    if sum(sign for _, sign in degree(cells)[1]) % 2 == 0:
        return set()
    bins = compulsory_bins(cells, knots, count, offset, span)
    return {(a,b) for a,b in combinations(sorted(bins), 2) if span < b-a < count-span}
