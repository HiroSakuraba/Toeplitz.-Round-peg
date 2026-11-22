"""Exact local two-graph certificates on an ordered square-grid boundary."""
from bisect import bisect_left, bisect_right, insort
from fractions import Fraction as Q
from itertools import product

from polygon_orbit_census_v223 import _solve
import json

RIFFORD = Q(9, 500)


class Graph:
    def __init__(self, points):
        if len(points) < 2:
            raise ValueError("graph needs two abscissas")
        self.x, self.y = zip(*[(Q(x), Q(y)) for x, y in points])
        for a, b, c, d in zip(self.x, self.x[1:], self.y, self.y[1:]):
            if b <= a or abs(d-c) > b-a:
                raise ValueError("expected an increasing 1-Lipschitz graph")

    def at(self, x):
        if not self.x[0] <= x <= self.x[-1]:
            raise ValueError("outside graph")
        i = min(bisect_right(self.x, x)-1, len(self.x)-2)
        if self.x[i] == x:
            return self.y[i]
        return self.y[i]+(x-self.x[i])*(self.y[i+1]-self.y[i])/(self.x[i+1]-self.x[i])


def gap_profile(f, g, left, right):
    knots = sorted({left, right} | {x for graph in (f, g) for x in graph.x if left < x < right})
    return [(x, g.at(x)-f.at(x)) for x in knots]


def local_test(f, g, left, right):
    left, right = Q(left), Q(right)
    if left >= right:
        raise ValueError("empty interval")
    profile = gap_profile(f, g, left, right)
    gaps = [w for _, w in profile]
    if min(gaps) < 0 or any(w == 0 for w in gaps[1:-1]) or max(gaps) == 0:
        raise ValueError("graphs must bound a positive open lens")
    height = max(gaps)
    lag = height/80
    if lag > right-left:
        return None
    margins = []
    for x, y, w in ((left, left+lag, gaps[0]), (right, right-lag, gaps[-1])):
        if w:
            margins.extend((f.at(y)-g.at(x)+lag, f.at(x)+lag-g.at(y)))
    if margins and min(margins) <= 0:
        return None
    return {"interval": [str(left), str(right)], "height": str(height), "lag": str(lag),
            "endpoint_gaps": [str(gaps[0]), str(gaps[-1])],
            "minimum_margin": str(min(margins)) if margins else None,
            "side_floor": str(RIFFORD*height)}


def components(profile, level):
    """Components of the strict superlevel set with their closures, including plateaus."""
    start = profile[0][0] if profile[0][1] > level else None
    for (a, v), (b, w) in zip(profile, profile[1:]):
        if v <= level < w:
            start = a + (b-a)*(level-v)/(w-v)
        elif v > level >= w:
            end = a + (b-a)*(level-v)/(w-v)
            yield start, end
            start = None
    if start is not None:
        yield start, profile[-1][0]


def validate_grid_loop(curve):
    """Count all closed perpendicular intersections and exclude collinear overlaps."""
    if len(curve) < 4 or len(set(curve)) != len(curve):
        raise ValueError("expected distinct cyclic polygon vertices")
    groups = {}
    events = []
    directions = []
    for (x,y), (a,b) in zip(curve, curve[1:]+curve[:1]):
        if (x == a) == (y == b):
            raise ValueError("expected nonzero orthogonal edges")
        horizontal = y == b
        directions.append(horizontal)
        if horizontal:
            lo, hi = sorted((x,a))
            events.extend(((lo,0,y,y),(hi,2,y,y)))
        else:
            lo, hi = sorted((y,b))
            events.append((x,1,lo,hi))
        groups.setdefault((horizontal, y if horizontal else x), []).append((lo,hi))
    for group in groups.values():
        ordered = sorted(group)
        if any(b > c for (_,b),(c,_) in zip(ordered,ordered[1:])):
            raise ValueError("overlapping grid edges")
    active = []
    total = 0
    for _, kind, lo, hi in sorted(events):
        if kind == 0:
            insort(active, lo)
        elif kind == 1:
            total += bisect_right(active, hi)-bisect_left(active, lo)
        else:
            active.pop(bisect_left(active, lo))
    if total != sum(d != directions[i-1] for i,d in enumerate(directions)):
        raise ValueError("self-contact in grid polygon")


def monotone_runs(curve, frame):
    if frame not in (0, 1) or len(curve) < 4:
        raise ValueError("need a square-grid loop and a diagonal frame")
    points = []
    for (x, y), (a, b) in zip(curve, curve[1:]+curve[:1]):
        if (x == a) == (y == b):
            raise ValueError("expected nonzero horizontal or vertical grid edges")
        u, v = Q(x+y)/2, Q(y-x)/2
        points.append((u, v) if frame == 0 else (v, -u))
    direction = [p[0] < q[0] for p, q in zip(points, points[1:]+points[:1])]
    cuts = [i for i, sign in enumerate(direction) if sign != direction[i-1]]
    runs = []
    for i, j in zip(cuts, cuts[1:]+cuts[:1]):
        chain = points[i:j+1] if i < j else points[i:]+points[:j+1]
        runs.append(Graph(chain if direction[i] else chain[::-1]))
    return runs


def dyadic_levels(profile):
    positive = [w for _, w in profile if w > 0]
    if not positive:
        return
    top, bottom = max(positive), min(positive)/2
    level = Q(2)**(top.numerator.bit_length()-top.denominator.bit_length())
    while level >= bottom:
        if level < top:
            yield level
        level /= 2


def select(curve, frame=0):
    """Input is a simple polygon in cyclic order. Return the strongest hit in this frame."""
    validate_grid_loop(curve)
    runs = monotone_runs(curve, frame)
    best = None
    for i, first in enumerate(runs):
        for j in range(i+1, len(runs)):
            second = runs[j]
            left, right = max(first.x[0], second.x[0]), min(first.x[-1], second.x[-1])
            if left >= right:
                continue
            mid = (left+right)/2
            f, g = (first, second) if first.at(mid) < second.at(mid) else (second, first)
            profile = gap_profile(f, g, left, right)
            if min(w for _, w in profile) < 0 or any(w <= 0 for _, w in profile[1:-1]):
                raise ValueError("crossing or touching graph interiors")
            # First try full overlap, then look for narrow sections inside larger lobes.
            largest_possible_gap = max(w for _, w in profile)/40
            trials = {(left, right)}
            for level in dyadic_levels(profile):
                trials.update(components(profile, level))
            for a, b in sorted(trials):
                # The endpoint cone tests imply 2d > gap at either artificial seam.
                # No subinterval can reach a larger height than the full overlap.
                if max(g.at(a)-f.at(a), g.at(b)-f.at(b)) >= largest_possible_gap:
                    continue
                result = local_test(f, g, a, b)
                if result and (best is None or Q(result["height"]) > Q(best["height"])):
                    best = dict(result, frame=frame, run_pair=[i, j], run_count=len(runs))
    return best



def square_points(x, y, a, b):
    return ((x,y), (x+a,y+b), (x+a-b,y+a+b), (x-b,y+a))


def affine_refinement(f, g, left, right):
    """Try both secants on each graph and bound a four-dimensional degree box."""
    mid = (left+right)/2
    halves = [(left,mid), (mid,right)]
    padding = (right-left)/512
    supports = []
    for graph in (f, g):
        lines = []
        for l, r in halves:
            slope = (graph.at(r)-graph.at(l))/(r-l)
            intercept = graph.at(l)-slope*l
            errors = [abs(y-slope*x-intercept) for x,y in zip(graph.x,graph.y) if l < x < r]
            lines.append((slope,intercept,l,r,max(errors, default=Q(0))))
        supports.append(lines)

    best = None
    basis = [square_points(*(int(i==j) for i in range(4))) for j in range(4)]
    for lines in product(supports[0], supports[0], supports[1], supports[1]):
        if any(error >= padding for _,_,_,_,error in lines):
            continue
        matrix = [[basis[j][i][1]-line[0]*basis[j][i][0] for j in range(4)]
                  for i,line in enumerate(lines)]
        rhs = [line[1] for line in lines]
        root = _solve(matrix, rhs)
        if root is None:
            continue
        points = square_points(*root)
        if root[2] <= 0 or any(not l < x < r for (x,_),(_,_,l,r,_) in zip(points,lines)):
            continue
        vertices = [_solve(matrix,[c+sign*padding for c,sign in zip(rhs,signs)])
                    for signs in product((-1,1),repeat=4)]
        if any(not l < x < r for z in vertices
               for (x,_),(_,_,l,r,_) in zip(square_points(*z),lines)):
            continue
        ranges = [(min(z[j] for z in vertices),max(z[j] for z in vertices)) for j in (2,3)]
        if ranges[0][0] <= 0:
            continue
        floor = sum(max(Q(0), lo, -hi)**2 for lo,hi in ranges)
        if best is None or floor > Q(best["side_squared_floor"]):
            best = {"padding":str(padding), "vertical_error":str(max(line[-1] for line in lines)),
                    "contact_intervals":[list(map(str,(line[2],line[3]))) for line in lines],
                    "tube_lines":[[str(m),str(c)] for m,c,_,_,_ in lines],
                    "abscissa_ranges":[[str(min(square_points(*z)[j][0] for z in vertices)),
                                        str(max(square_points(*z)[j][0] for z in vertices))]
                                       for j in range(4)],
                    "horizontal_run_and_rise":[list(map(str,pair)) for pair in ranges],
                    "side_squared_floor":str(floor)}
    return best


def staircase_fixture(n):
    """World grid loop in (x,y). The two marked subarcs are not returned."""
    if type(n) is not int or n < 1:
        raise ValueError("need positive integer mesh factor")
    e = Q(1,256)
    def roof(upper):
        points = [(Q(0), e if upper else Q(0))]
        for m in ((Q(1,2), -Q(1,4)) if upper else (-Q(1,4), Q(1,2))):
            for _ in range(4*n):
                x, y = points[-1]
                change = 1 if upper else -1
                first = (1+change*m)/(2*n)
                points.extend(((x+first, y+change*first), (x+Q(1,n), y+m/n)))
        return points

    def detour(x, y, sign, scale):
        half = [(0,0), (2,-2), (1,-3), (2,-4), (6,0)]
        lower = [(x+sign*scale*a, y+scale*b) for a,b in half]
        upper = [(x+sign*scale*a, y+e-scale*b) for a,b in reversed(half)]
        return lower + [(x+sign*(6*scale+e/2), y+e/2)] + upper

    low, high = roof(False), roof(True)
    right = detour(*low[-1], 1, Q(1))
    left = detour(*low[0], -1, Q(1,2))
    uv = low + right[1:] + list(reversed(high))[1:] + list(reversed(left))[1:-1]
    return [(u-v, u+v) for u,v in uv]



def add_folds(curve, n):
    """Replace the middle third of each short grid edge with a small rectangular detour."""
    grid = Q(1,4*n)
    depth = grid/24
    result = []
    for (x,y),(a,b) in zip(curve,curve[1:]+curve[:1]):
        ex, ey = (a>x)-(a<x), (b>y)-(b<y)
        cells = (abs(a-x)+abs(b-y))/grid
        if cells.denominator != 1:
            raise ValueError("edge is off the fixture grid")
        remaining = int(cells)
        while remaining:
            take = min(4,remaining)
            size = take*grid
            for along, offset in ((0,0),(size/3,0),(size/3,depth),
                                  (2*size/3,depth),(2*size/3,0)):
                result.append((x+ex*along-ey*offset, y+ey*along+ex*offset))
            x, y = x+ex*size, y+ey*size
            remaining -= take
    return result


def scan_tubes(curve, box):
    """Check a previously certified box on an unmarked simple target in frame zero."""
    validate_grid_loop(curve)
    uv = [(Q(x+y)/2, Q(y-x)/2) for x,y in curve]
    p = Q(box["padding"])
    if p <= 0 or len(box["tube_lines"]) != 4 or len(box["abscissa_ranges"]) != 4:
        raise ValueError("expected four supports from a certified positive-width box")
    hits = []
    for coeffs, bounds in zip(box["tube_lines"],box["abscissa_ranges"],strict=True):
        m,c = map(Q,coeffs)
        left,right = map(Q,bounds)
        outside = next((i for i,(u,v) in enumerate(uv) if abs(v-m*u-c) >= p), None)
        start = 0 if outside is None else (outside+1)%len(uv)
        low = high = None
        for k in range(len(uv)):
            i = (start+k)%len(uv)
            u,v = uv[i]
            if abs(v-m*u-c) >= p:
                low = high = None
                continue
            if low is None or u < uv[low][0]:
                low = i
            if high is None or u > uv[high][0]:
                high = i
            if uv[low][0] < left and uv[high][0] > right:
                hits.append([low,high])
                break
        else:
            return None
    return hits


def report():
    rows = []
    for n in (64, 128):
        curve = staircase_fixture(n)
        result = select(curve)
        if result is None or Q(result["height"]) < 3+Q(1,256):
            raise ArithmeticError("failed to recover a central lobe")
        runs = monotone_runs(curve, result["frame"])
        left, right = map(Q, result["interval"])
        f, g = sorted((runs[i] for i in result["run_pair"]), key=lambda graph: graph.at((left+right)/2))
        refined = affine_refinement(f, g, left, right)
        if refined is None or 2*Q(refined["side_squared_floor"]) <= 9:
            raise ArithmeticError("lost the selected local degree box")
        rows.append(dict(n=n, vertices=len(curve), refinement=refined,
                         refined_physical_floor_squared=str(2*Q(refined["side_squared_floor"])),
                         physical_side_squared_floor=str(2*Q(result["side_floor"])**2), **result))
    n = 64
    folded = add_folds(staircase_fixture(n),n)
    box = rows[-1]["refinement"]  # Independently selected on the other refinement.
    hits = scan_tubes(folded,box)
    if hits is None:
        raise ArithmeticError("lost odd degree on the folded target")
    spans = [max(run.x[-1]-run.x[0] for run in monotone_runs(folded,i)) for i in (0,1)]
    return {"staircases": rows, "folded_target":{"n":n,"vertices":len(folded),
            "max_run_spans":list(map(str,spans)),"traversal_end_indices":hits,
            "physical_side_squared_floor":str(2*Q(box["side_squared_floor"]))}}


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
