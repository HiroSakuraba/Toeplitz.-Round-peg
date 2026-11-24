"""Check lifted paths of approximate squares on registered rational polygons."""
from bisect import bisect_right
from collections import deque
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement

from polygon_orbit_census_v223 import add, dot, mul, pt, sub, validate_jordan_polygon


@dataclass(frozen=True)
class PathReplay:
    relabeling: int
    max_residual_squared: Q
    separated_parameters: tuple
    separated_pair: tuple
    parameter_distance: Q
    chord_squared: Q


def _gaps(row):
    return tuple(b-a for a,b in zip(row,row[1:]+(row[0]+1,)))


def _balanced_interval(gaps0, gaps1):
    lo,hi = Q(0),Q(1)
    for a,b in zip(gaps0,gaps1):
        if a == b:
            if a > Q(1,2):
                return None
        else:
            cross = (Q(1,2)-a)/(b-a)
            if b > a:
                hi = min(hi,cross)
            else:
                lo = max(lo,cross)
    return (lo,hi) if lo <= hi else None


def _residual_squared(points):
    p,q,r,s = points
    side = sub(q,p)
    turn = pt(-side.y,side.x)
    er = sub(sub(r,q),turn)
    es = sub(sub(s,p),turn)
    return max(dot(er,er),dot(es,es))


def replay_path(vertices, knots, rows):
    """Rows are four unwrapped, cyclically ordered times. Interpolate linearly.

    The path must end at a nontrivial cyclic choice of the starting label,
    up to a common integer turn. Knots include both 0 and 1.
    """
    poly = validate_jordan_polygon([pt(*p) for p in vertices])
    knots = tuple(map(Q,knots))
    rows = [tuple(map(Q,row)) for row in rows]
    if (len(knots) != len(poly)+1 or knots[0] != 0 or knots[-1] != 1
            or any(a >= b for a,b in zip(knots,knots[1:]))):
        raise ValueError("invalid registration")
    if len(rows) < 2 or any(len(row) != 4 or min(_gaps(row)) < 0 for row in rows):
        raise ValueError("need lifted cyclically ordered rows")
    def rotate(k):
        return rows[0][k:]+tuple(x+1 for x in rows[0][:k])
    def matches(k):
        shifts = {b-a for a,b in zip(rotate(k),rows[-1])}
        return len(shifts) == 1 and next(iter(shifts)).denominator == 1
    k = next((k for k in (1,2,3) if matches(k)), None)
    if k is None:
        raise ValueError("endpoint does not cyclically relabel the starting row")

    def position(t):
        t %= 1
        i = bisect_right(knots,t)-1
        f = (t-knots[i])/(knots[i+1]-knots[i])
        return add(poly[i],mul(sub(poly[(i+1)%len(poly)],poly[i]),f))

    maximum = Q(0)
    witness = None
    for start,end in zip(rows,rows[1:]):
        breaks = {Q(0),Q(1)}
        for a,b in zip(start,end):
            if a != b:
                low,high = sorted((a,b))
                for lap in range(low//1,high//1+1):
                    for knot in knots[:-1]:
                        t = lap+knot
                        if low < t < high:
                            breaks.add((t-a)/(b-a))
        balanced = _balanced_interval(_gaps(start),_gaps(end))
        witness_time = None
        if balanced is not None and witness is None:
            witness_time = sum(balanced)/2
            breaks.add(witness_time)
        for t in sorted(breaks):
            row = tuple(a+t*(b-a) for a,b in zip(start,end))
            points = [position(s) for s in row]
            maximum = max(maximum,_residual_squared(points))
            if t == witness_time:
                def separation(pair):
                    i,j = pair
                    gap = row[j]-row[i]
                    return min(gap,1-gap)
                pair = max(combinations(range(4),2),key=separation)
                chord = sub(points[pair[0]],points[pair[1]])
                witness = row,pair,separation(pair),dot(chord,chord)
    if witness is None or witness[2] < Q(1,3):
        raise ArithmeticError("balanced configuration missing")
    return PathReplay(k,maximum,*witness)


def find_vertex_path(vertices, max_residual_squared):
    """Find a relabeling path on the vertex grid with uniform edge times.

    Moves advance or retreat one corner, or all four together, by one edge.
    Return lifted time rows for independent replay. A grid miss is inconclusive.
    """
    poly = validate_jordan_polygon([pt(*p) for p in vertices])
    limit = Q(max_residual_squared)
    if limit < 0:
        raise ValueError("negative squared residual")
    n = len(poly)

    def canonical(row):
        shift = (row[0]//n)*n
        return tuple(x-shift for x in row)

    allowed = set()
    for i in range(n):
        for rest in combinations_with_replacement(range(i,i+n+1),3):
            row = (i,)+rest
            if _residual_squared([poly[z%n] for z in row]) <= limit:
                allowed.add(row)

    moves = [(1,1,1,1)] + [tuple(int(i == j) for i in range(4)) for j in range(4)]
    moves += [tuple(-x for x in move) for move in moves]

    def neighbors(row):
        for move in moves:
            new = canonical(tuple(x+y for x,y in zip(row,move)))
            if new in allowed:
                yield new

    def explore(root, goal=None):
        queue = deque([root])
        parents = {root: None}
        while queue:
            row = queue.popleft()
            if row == goal:
                break
            for new in neighbors(row):
                if new not in parents:
                    parents[new] = row
                    queue.append(new)
        return parents

    remaining = allowed.copy()
    while remaining:
        component = explore(min(remaining))
        remaining.difference_update(component)
        for row in component:
            for k in (1,2,3):
                end = canonical(row[k:]+tuple(x+n for x in row[:k]))
                if end not in component:
                    continue
                parents = explore(row,end)
                chain = []
                cursor = end
                while cursor is not None:
                    chain.append(cursor)
                    cursor = parents[cursor]
                path = []
                for node in reversed(chain):
                    start = node[0]
                    if path:
                        # n >= 3: the one-step lift nearest the previous start is unique.
                        previous = path[-1][0]
                        start += ((previous-start+n//2)//n)*n
                    path.append(tuple(x+start-node[0] for x in node))
                return [tuple(Q(x,n) for x in node) for node in path]
    return None
