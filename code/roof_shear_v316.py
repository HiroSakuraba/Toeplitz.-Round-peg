"""Mesh-independent square certificates from an upper arc and a base arc."""
from fractions import Fraction as Q
from itertools import product
import json

from polygon_orbit_census_v223 import (
    _affine_system, _box_vertices, _evaluate, _side_sq, _solve, pt,
)
from boundary_or_winding_v313 import lerp


def terrain_certificate(heights):
    """Unit grid columns, with a common base. Return an exact largest base square."""
    if not heights or any(type(h) is not int or h <= 0 for h in heights):
        raise ValueError("expected a nonempty list of positive integer column heights")
    width = len(heights)
    roof = (0, *heights, 0)

    def contains(side):
        run = 0
        for h in heights:
            run = run + 1 if h >= side else 0
            if run >= side:
                return True
        return False

    low, high = 0, min(width, max(heights)) + 1
    while high - low > 1:
        mid = (low + high) // 2
        if contains(mid):
            low = mid
        else:
            high = mid

    high_limit = min(width, max(heights))
    for side in range(high_limit, low - 1, -1):
        hits = {x for x in range(width + 1)
                if min(roof[x:x+2]) <= side <= max(roof[x:x+2])}
        start = next((x for x in sorted(hits) if x + side in hits), None)
        if start is not None:
            return {"contained_side": low, "boundary_side": side, "base_left": start}
    raise ArithmeticError(f"no square between {low} and {high_limit}")


def band_intervals(arc, height):
    """Closed horizontal projections of portions below or on a horizontal line."""
    intervals = []
    for p, q in zip(arc, arc[1:]):
        if min(p.y, q.y) > height:
            continue
        if max(p.y, q.y) > height:
            cut = lerp(p, q, (height - p.y) / (q.y - p.y))
            if p.y > height:
                p = cut
            else:
                q = cut
        intervals.append(tuple(sorted((p.x, q.x))))
    return intervals


def gate_certificate(upper, base, delta):
    if len(upper) < 3 or len(base) < 2 or delta < 0:
        raise ValueError("invalid arcs or strip width")
    h = upper[0].y
    if (h <= delta or upper[-1].y != h or upper[0].x >= upper[-1].x
            or any(p.y <= h for p in upper[1:-1])):
        raise ValueError("expected a left-to-right excursion above positive height")
    if any(abs(p.y) > delta for p in base):
        raise ValueError("base exits its strip")
    span = min(p.x for p in upper), max(p.x for p in upper)
    if not (base[0].x < span[0] - 2*delta and base[-1].x > span[1] + 2*delta):
        raise ValueError("insufficient clearance from base endpoints")

    intervals = sorted(band_intervals(upper, h + 2*delta))
    # Any valid split separates the complete low-band projection by a gap > h+delta.
    left_min = intervals[0][0]
    right_max = max(b for _, b in intervals)
    left_max = left_min
    for (_, end), (right_min, _) in zip(intervals, intervals[1:]):
        left_max = max(left_max, end)
        if (left_max - left_min < h - delta and right_max - right_min < h - delta
                and right_min - left_max > h + delta
                and upper[0].x <= left_max and upper[-1].x >= right_min):
            return {"height": str(h), "strip": str(delta), "side_floor": str(h-delta),
                    "gates": [[str(left_min), str(left_max)], [str(right_min), str(right_max)]]}
    raise ValueError("low-band projections fail the two-gate inequalities")


def digital_trapezoid(n):
    if type(n) is not int or n < 4:
        raise ValueError("need integer refinement at least four")
    x = y = 0
    path = [pt(0, 0)]
    for count, dx, dy in ((12*n, 1, 1), (2*n, -4, 1), (6*n, -1, -1), (2*n, 1, -4)):
        for _ in range(count):
            x += dx
            path.append(pt(Q(x+y, 2*n), Q(y-x, 2*n)))
            y += dy
            path.append(pt(Q(x+y, 2*n), Q(y-x, 2*n)))
    base = path[:24*n+1]
    roof = list(reversed(path[24*n:]))
    h = 5 - Q(2, n)
    upper = []
    for p, q in zip(roof, roof[1:]):
        if (p.y < h <= q.y) or (p.y > h >= q.y):
            upper.append(lerp(p, q, (h-p.y)/(q.y-p.y)))
        if q.y > h:
            upper.append(q)
    return upper, base



def folded_example():
    upper = [pt(1, 2), pt(2, 5), pt(7, Q(11, 2)), pt(9, 2)]
    base = [pt(-2, 0), pt(3, 0), pt(4, Q(1, 8)), pt(Q(15, 4), Q(1, 5)),
            pt(8, Q(1, 12)), pt(11, -Q(1, 8)), pt(13, 0)]
    cert = gate_certificate(upper, base, Q(1, 4))
    top_edges = list(zip(upper, upper[1:]))
    base_edges = list(zip(base, base[1:]))
    sides = []
    for chart in product(top_edges, top_edges, base_edges, base_edges):
        matrix, rhs = _affine_system(chart, -1)
        tags = _solve(matrix, rhs)
        if tags is None:
            if next(_box_vertices(matrix, rhs), None) is not None:
                raise ArithmeticError("feasible singular folded-example chart")
        elif all(0 <= t <= 1 for t in tags):
            if min(tags) == 0 or max(tags) == 1:
                raise ArithmeticError("folded-example contact on a vertex")
            sides.append(_side_sq(_evaluate(chart, tags)))
    if len(sides) % 2 != 1 or min(sides) < Q(cert["side_floor"])**2:
        raise ArithmeticError("contact census contradicts the gate certificate")
    return dict(cert, regular_marked_squares=len(sides),
                squared_sides=list(map(str, sorted(sides))))


def report():
    cases = []
    strict = 0
    for heights in product(range(1, 5), repeat=6):
        cert = terrain_certificate(heights)
        if cert["boundary_side"] > cert["contained_side"]:
            strict += 1
            if not cases:
                cases.append(dict(heights=list(heights), **cert))
    refinements = []
    for n in (4, 16, 64):
        upper, base = digital_trapezoid(n)
        cert = gate_certificate(upper, base, Q(1, 2*n))
        refinements.append(dict(mesh=f"1/{n}", physical_floor_squared=str(2*Q(cert["side_floor"])**2),
                                longest_grid_edge=str(Q(4, n)), **cert))
    return {"six_column_profiles": 4**6, "strict_improvements": strict, "first_strict": cases,
            "folded_base": folded_example(), "rotated_trapezoid": refinements}


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
