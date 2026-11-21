"""Partial pairing on a fixed side wall and exact repair of four near contacts."""
from fractions import Fraction as F
import json

from polygon_orbit_census_v223 import (
    _evaluate, _on_seg, _solve, add, dot, mul, pt, sub, validate_jordan_polygon,
)
from half_turn_degree_v313 import ear_census, profile

SIDE = F(2, 5)
START = F(9, 10)


def norm_inf(v):
    return max(abs(v.x), abs(v.y))


def polygon(right, upper_x, lower_x, upper_scale, lower_scale):
    interior = profile()[1:-1]
    upper = [add(pt(upper_x, 1), mul(p, upper_scale)) for p in reversed(interior)]
    lower = [sub(pt(lower_x, -1), mul(p, lower_scale)) for p in reversed(interior)]
    return validate_jordan_polygon(
        tuple([pt(right, 1)] + upper + [pt(-4, 1), pt(-4, -1)] + lower + [pt(right, -1)])
    )


def bracket(critical_squared, depth):
    if depth < 0 or not START**2 < critical_squared < 1:
        raise ValueError("invalid bisection range")
    lo, hi = START, F(1)
    for _ in range(depth):
        mid = (lo + hi) / 2
        if mid * mid == critical_squared:
            raise ValueError("expected a strict irrational crossing")
        if mid * mid < critical_squared:
            lo = mid
        else:
            hi = mid
    return lo, hi


def displacement(z, p, q, eps):
    weight = max(F(0), min(F(1), 2 - norm_inf(sub(z, p)) / (2 * eps)))
    return mul(sub(q, p), weight)


def repair(curve, contacts, corners, eps):
    """Apply four disjoint PL translations, with the old edge parametrization."""
    if eps <= 0 or len(contacts) != 4 or len(corners) != 4:
        raise ValueError("need four marked points and positive tolerance")
    boundary = list(zip(curve, curve[1:] + curve[:1]))
    side = sub(corners[1], corners[0])
    turn = pt(-side.y, side.x)
    if (norm_inf(side) == 0 or corners[2] != add(corners[1], turn)
            or corners[3] != add(corners[0], turn)):
        raise ValueError("expected a counterclockwise square")
    for i, (p, q) in enumerate(zip(contacts, corners)):
        if norm_inf(sub(p, q)) > eps or not any(_on_seg(p, a, b) for a, b in boundary):
            raise ValueError("invalid contact or displacement")
        if any(norm_inf(sub(p, prev)) <= 8 * eps for prev in contacts[:i]):
            raise ValueError("overlapping contact supports")

    result = []
    axes = (pt(1, 0), pt(0, 1))
    diagonal_normals = (pt(1, 1), pt(1, -1))
    for a, b in boundary:
        direction = sub(b, a)
        ticks = {F(0)}
        for p in contacts:
            if (max(a.x, b.x) < p.x - 4 * eps or min(a.x, b.x) > p.x + 4 * eps
                    or max(a.y, b.y) < p.y - 4 * eps or min(a.y, b.y) > p.y + 4 * eps):
                continue
            # The max-norm plateau, support walls and changes of active coordinate.
            lines = [(axis, k * eps) for axis in axes for k in (-4, -2, 2, 4)]
            lines.extend((n, F(0)) for n in diagonal_normals)
            for n, offset in lines:
                slope = dot(n, direction)
                if slope:
                    t = (offset - dot(n, sub(a, p))) / slope
                    if 0 < t < 1:
                        ticks.add(t)
        for t in sorted(ticks):
            z = add(a, mul(direction, t))
            shifted = z
            for p, q in zip(contacts, corners):
                shifted = add(shifted, displacement(z, p, q, eps))
            result.append(shifted)
    repaired = validate_jordan_polygon(tuple(result))
    segments = list(zip(repaired, repaired[1:] + repaired[:1]))
    if not all(any(_on_seg(q, a, b) for a, b in segments) for q in corners):
        raise ArithmeticError("repaired polygon missed a square vertex")
    return repaired


def event_data():
    roots = ear_census()
    if len(roots) != 2:
        raise ArithmeticError("unexpected normalized census")
    smaller, larger = roots
    side_squared, ids, tags = larger
    critical = SIDE**2 / side_squared
    if not smaller[0] < SIDE**2 or not START**2 < critical < 1:
        raise ArithmeticError("cutoff isolation failed")
    edges = list(zip(profile(), profile()[1:]))
    support = [edges[i] for i in ids]
    root = _evaluate(support, tags)
    normals = []
    for j, (a, b) in enumerate(support):
        d = sub(b, a)
        normal = pt(-d.y, d.x)
        pulled = normal
        for _ in range(j):
            pulled = pt(pulled.y, -pulled.x)
        normals.append([normal.x, normal.y, pulled.x, pulled.y])
    if _solve(normals, [0]*4) is None or not all(0 < t < 1 for t in tags):
        raise ArithmeticError("center contact system is singular or at an edge end")
    return critical, root, side_squared


def replay(depth=8):
    critical, root, size = event_data()
    lo, hi = bracket(critical, depth)
    eps = hi - lo
    b = -1 + eps
    curve = polygon(6, 1, b, 1, lo)
    contacts = tuple(sub(pt(b, -1), mul(p, lo)) for p in root)
    corners = tuple(sub(pt(b, -1), mul(p, hi)) for p in root)
    repaired = repair(curve, contacts, corners, eps)
    direction = sub(corners[1], corners[0])
    if dot(direction, direction) != hi**2 * size or not lo**2 * size < SIDE**2 < hi**2 * size:
        raise ArithmeticError("failed to straddle fixed cutoff")
    # Both paired events occur above/below the flat parts of the limiting boundary.
    # Aperture is measured from the old ear center, including baseline contacts.
    aperture = F(7, 16) + F(17, 32)
    horizontal_gap = 2 - aperture - F(7, 16)
    if horizontal_gap <= SIDE:
        raise ArithmeticError("distant crossing separation failed")
    return {
        "critical_scale_squared": str(critical),
        "normalized_square_side_squared": str(size),
        "cutoff_side": str(SIDE),
        "depth": depth,
        "lower_scale": str(lo), "upper_scale": str(hi),
        "repair_budget_sup": str(eps),
        "lower_ear_center_x": str(b),
        "repaired_square_side_squared": str(dot(direction, direction)),
        "repaired_polygon_vertices": len(repaired),
        "horizontal_clearance_from_opposite_ear": str(horizontal_gap),
        "remote_crossing_distance_lower_bound": str(min(SIDE, horizontal_gap)),
    }


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2))
