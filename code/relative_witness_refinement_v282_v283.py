"""Exact controls for refinement of a relative rectangle witness."""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from adversarial_rectangles_v278 import moving_notch
from liouville_primitive_experiments_v276 import block_vertices, closed_truncation
from relative_rectangle_v269_v273 import notch_square_certificate


ARTIFACT = Path(__file__).resolve().parents[1] / "results/v282-v283/refinement-witnesses.json"


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def difference(a, b):
    return tuple(x-y for x, y in zip(a, b))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def rotate(a):
    return (-a[1], a[0])


def on_segment(point, start, end):
    offset = difference(point, start)
    edge = difference(end, start)
    return dot(offset, rotate(edge)) == 0 and 0 <= dot(offset, edge) <= dot(edge, edge)


def diagonal_dot(corners):
    a, b, c, d = corners
    return dot(difference(c, a), difference(d, b))


def notch_rectangle(n, low, high):
    """One fixed labeling of the local notch triangle, including its flat edge."""
    certificate = notch_square_certificate(n)
    eps = Q(1, n)
    if not 0 <= low <= high <= eps:
        raise ValueError("heights must satisfy 0 <= low <= high <= 1/n")
    left, base = certificate["polygon"][6]
    right = left+eps
    return ((left, base+low), (right, base+low),
            (right, base+high), (left, base+high))


def notch_audit(n):
    eps = Q(1, n)
    square = notch_rectangle(n, 0, eps)
    flat = notch_rectangle(n, eps/2, eps/2)
    control = moving_notch(n)
    require(diagonal_dot(square) == 0
            and set(square) == set(notch_square_certificate(n)["selected_square"]),
            "triangle apex is not the selected notch square")
    require(diagonal_dot(flat) == -eps*eps,
            "the notch flat frame does not agree with the fixed labeling")
    return {
        "n": n,
        "flat_dot": str(diagonal_dot(flat)),
        "local_relative_index": 0,
        "upper_squared_radius_on_local_triangle": str(eps*eps/2),
        "uniform_coordinate_error_upper": control["uniform_coordinate_error_upper"],
        "other_square_side_squared": control["other_square_side_squared"],
    }


def support():
    arm = tuple(block_vertices(0, False))
    reflected = tuple((-a, -b) for a, b in arm)
    bottom = (Q(0), Q(-1))
    fixed = ((reflected[0], reflected[1]), (reflected[0], bottom), (bottom, arm[0]))
    outer_walk_vertices = arm + tuple(reversed(reflected[1:]))
    return fixed, outer_walk_vertices


def diamond_vertices():
    h = Q(1, 2)
    return ((h, 0), (0, h), (-h, 0), (0, -h))


def verify_spiral_domain():
    fixed, outer = support()
    candidates = outer + diamond_vertices()
    for a, b in candidates:
        offset = 7+7*a+8*b
        require(abs(a)+abs(b) <= 1 and b >= -Q(7, 8) and 0 <= offset <= 14
                and a+offset >= 0 and 1+a-b > 0 and 1+b-a >= 0,
                f"spiral walk leaves the proved domain at {(a, b)}")
        require((offset == 0) == ((a, b) == outer[-1]), "unexpected zero offset")
        require((1+b-a == 0) == ((a, b) == outer[0]), "unexpected zero side factor")
    require(all(p[0]-p[1] <= Q(7, 8) for p in fixed[0])
            and all(p[0]-p[1] == 1 for p in fixed[2]), "fixed support distance changed")
    return fixed, outer


def contact_corners(q, right):
    """Compute the corners after supplying the positive root of the quadratic."""
    a, b = q
    left = 15*right-(7+7*a+8*b)
    return ((a-left-right, b+left-right), (-left, left-1), (right, right-1), q)


def check_contacts(fixed, corners):
    require(all(on_segment(p, *edge) for p, edge in zip(corners, fixed)),
            "the first three corners left their fixed support")
    require(difference(corners[2], corners[1]) == difference(corners[3], corners[0]),
            "mismatched parallel sides")


def spiral_audit():
    fixed, outer = verify_spiral_domain()
    start = contact_corners(outer[0], Q(1))
    terminal = contact_corners(outer[-1], Q(0))
    check_contacts(fixed, start)
    check_contacts(fixed, terminal)
    require(start[0] == start[1] and start[2] == start[3]
            and terminal[0] == terminal[3] and terminal[1] == terminal[2],
            "flat types at the ends changed")
    signs = (diagonal_dot(start), diagonal_dot(terminal))
    require(signs[0] > 0 > signs[1], "the lift lacks opposite flat signs")

    u = Q(240, 241)
    a, b = outer[1]
    q = (1-u+u*a, u*b)
    corners = contact_corners(q, Q(225, 241))
    check_contacts(fixed, corners)
    require(on_segment(q, *outer[:2]), "square is outside the fixed first segment")
    side = difference(corners[2], corners[1])
    require(rotate(side) == difference(corners[3], corners[2]), "square residual changed")
    squared = dot(side, side)
    require(squared > 1, "fixed square lacks the strict unit-side lower bound")
    return {
        "fixed_support_segments": [[list(map(str, p)) for p in edge] for edge in fixed],
        "fourth_corner_outer_vertices": [list(map(str, p)) for p in outer],
        "inner_replacement_l1_radius_upper": "1/2",
        "squared_radius_lower_on_entire_witness": str(Q(1, 512)),
        "flat_dot_products": list(map(str, signs)),
        "relative_index_on_witness": 1,
        "one_shared_square_corners_ccw": [list(map(str, p)) for p in corners],
        "one_shared_square_squared_side": str(squared),
        "one_shared_square_fourth_corner_edge_parameter": str(u),
    }


def truncation_audit(blocks, counterwind):
    points = closed_truncation(blocks, counterwind)
    fixed, outer = support()
    require(tuple(points[:5]) == outer[:5], "outer arm changed")
    begin, end = points.index((Q(1, 2), 0)), points.index((-Q(1, 2), 0))
    require(begin < end and all(abs(a)+abs(b) <= Q(1, 2) for a, b in points[begin:end+1]),
            "inner arm replacement leaves its diamond")
    require(tuple(points[end:end+4]) == outer[5:], "opposite outer arm changed")
    require(tuple(points[-4:]) == (fixed[0][1], fixed[0][0], fixed[1][1], fixed[2][1]),
            "fixed supporting segments changed")
    return {"blocks": blocks, "counterwinding": counterwind, "changed_arc_vertices": end-begin+1}


def report():
    return {
        "notch": [notch_audit(n) for n in (4, 16, 128)],
        "spiral": spiral_audit(),
        "truncations": [truncation_audit(n, counterwind)
                        for counterwind in (False, True) for n in (1, 2, 3, 4)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = json.dumps(report(), indent=2, sort_keys=True) + "\n"
    if args.write:
        ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
        ARTIFACT.write_bytes(data.encode("utf-8"))
    else:
        require(ARTIFACT.read_text(encoding="utf-8") == data, "refinement artifact differs from replay")
    print("Verified the notch triangle and the common spiral witness.")


if __name__ == "__main__":
    main()
