"""Exact unmarked three-window certificates on rational Jordan polygons."""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product

from free_contact_v360 import _inside_clear
from polygon_orbit_census_v223 import pt, validate_jordan_polygon


@dataclass(frozen=True)
class Frame:
    # Coordinates before normalization: X = p.u, Y = p.Ju.
    direction: tuple
    left: Q
    bottom: Q
    width: Q


@dataclass(frozen=True)
class Traversal:
    # Increasing original polygon parameter, extended cyclically across index zero.
    start: Q
    end: Q
    low: Q
    high: Q


@dataclass(frozen=True)
class Certificate:
    frame: Frame
    relative_error: Q
    windows: tuple  # bottom, right, left
    outside: Q      # normalized probe abscissas, at ordinate 1
    inside: Q
    side_squared_floor: Q


def _clip(lo, hi, start, end, lower, upper):
    if start == end:
        return (lo, hi) if lower <= start <= upper else None
    a, b = sorted(((lower-start)/(end-start), (upper-start)/(end-start)))
    lo, hi = max(lo, a), min(hi, b)
    return (lo, hi) if lo <= hi else None


def _traversals(edges, normal, level, error, pad):
    """Clip cyclically and take the two projection extrema on each component."""
    runs = []
    for i, (p, q) in enumerate(edges):
        interval = _clip(Q(0), Q(1), p[normal], q[normal], level-error, level+error)
        if interval is None:
            continue
        s, t = p[1-normal], q[1-normal]
        interval = _clip(*interval, s, t, pad, 1-pad)
        if interval is None:
            continue
        samples = [(i+u, s+u*(t-s)) for u in interval]
        if runs and runs[-1][-1][0] == samples[0][0]:
            runs[-1].append(samples[1])
        else:
            runs.append(samples)
    n = len(edges)
    if len(runs) >= 2 and runs[0][0][0] == 0 and runs[-1][-1][0] == n:
        tail = runs.pop()
        runs[0] = tail + [(i+n, x) for i, x in runs[0][1:]]
    spans = []
    for run in runs:
        a, low = min(run, key=lambda item: item[1])
        b, high = max(run, key=lambda item: item[1])
        if low < high:
            spans.append(Traversal(min(a,b), max(a,b), low, high))
    return spans


def _clear_intervals(edges, radius):
    # The forbidden line intervals are exact horizontal sections of the
    # closed max-norm tube. This also handles horizontal tangencies and overlaps.
    blocked = []
    for (x,y), (xx,yy) in edges:
        interval = _clip(Q(0), Q(1), y, yy, 1-radius, 1+radius)
        if interval is not None:
            a,b = sorted(x+t*(xx-x) for t in interval)
            blocked.append((a-radius, b+radius))
    cursor = Q(0)
    for a,b in sorted(blocked):
        if cursor < min(a,1):
            yield cursor, min(a,1)
        cursor = max(cursor,b)
        if cursor >= 1:
            return
    if cursor < 1:
        yield cursor, Q(1)


def _project(poly, direction):
    u,v = direction
    return [(u*p.x+v*p.y, u*p.y-v*p.x) for p in poly]


def _frame_test(projected, frame, error, pad):
    L = frame.width
    vertices = [((x-frame.left)/L, (y-frame.bottom)/L) for x,y in projected]
    edges = list(zip(vertices, vertices[1:]+vertices[:1]))
    areas = [_traversals(edges, axis, level, error, pad)
             for axis,level in ((1,0),(0,1),(0,0))]
    if not all(areas):
        return None
    points = [pt(*p) for p in vertices]
    point_edges = list(zip(points, points[1:]+points[:1]))
    clear = []
    for a,b in _clear_intervals(edges, 8*error):
        sign = _inside_clear(((a+b)/2,Q(1)), point_edges, 0)
        if sign is not None:
            clear.append((a,b,sign))
    if len({sign for _,_,sign in clear}) < 2:
        return None
    for bottom, right, left in product(*areas):
        lo = max(left.low, 1-bottom.high+4*error, 1-right.high+8*error)
        hi = min(left.high, 1-bottom.low-4*error, 1-right.low-8*error)
        probes = {}
        for a,b,sign in clear:
            start, end = max(a,lo), min(b,hi)
            if start < end:
                probes[sign] = (start+end)/2
            if len(probes) == 2:
                u,v = frame.direction
                floor = (L*(1-2*error))**2/(2*(u*u+v*v))
                return Certificate(frame, error, (bottom,right,left),
                                   probes[False], probes[True], floor)
    return None


def _settings(relative_error, padding):
    e = Q(relative_error)
    pad = max(2*e, Q(1,1024)) if padding is None else Q(padding)
    if not 0 <= e < pad < Q(1,2):
        raise ValueError("require 0 <= relative_error < padding < 1/2")
    return e,pad


def scan_frames(vertices, frames, relative_error, padding=None):
    """Try a supplied iterable of frames, selecting all arcs and probes blindly."""
    poly = validate_jordan_polygon([pt(*v) for v in vertices])
    e,pad = _settings(relative_error,padding)
    for candidate in frames:
        frame = Frame(tuple(map(Q,candidate.direction)), Q(candidate.left),
                      Q(candidate.bottom), Q(candidate.width))
        if frame.width <= 0 or frame.direction == (0,0):
            raise ValueError("frame requires nonzero direction and positive width")
        hit = _frame_test(_project(poly,frame.direction), frame, e,pad)
        if hit is not None:
            return hit
    return None


def find_windows(vertices, relative_error, minimum_side=0, padding=None):
    """Try frames supported by one polygon edge and two perpendicular parallels."""
    poly = validate_jordan_polygon([pt(*v) for v in vertices])
    e,pad = _settings(relative_error,padding)
    cutoff = Q(minimum_side)
    if cutoff < 0:
        raise ValueError("negative minimum side")
    directions = set()
    for p,q in zip(poly, poly[1:]+poly[:1]):
        x,y = q.x-p.x, q.y-p.y
        if x:
            directions.add((Q(1), y/x))
        else:
            directions.add((Q(0), Q(1)))
    for axis in sorted(directions):
        for sign in (1,-1):
            u,v = (sign*c for c in axis)
            projected = _project(poly,(u,v))
            edges = list(zip(projected, projected[1:]+projected[:1]))
            xs = sorted({x for (x,_),(xx,_) in edges if x == xx})
            ys = sorted({y for (_,y),(_,yy) in edges if y == yy})
            for left,right in combinations(xs,2):
                L = right-left
                if (L*(1-2*e))**2 <= 2*(u*u+v*v)*cutoff**2:
                    continue
                for bottom in ys:
                    frame = Frame((u,v),left,bottom,L)
                    hit = _frame_test(projected,frame,e,pad)
                    if hit is not None:
                        return hit
    return None
