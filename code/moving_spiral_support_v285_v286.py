"""Exact moving spiral supports and a transversality check at both flats."""

import argparse
from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path

from relative_witness_refinement_v282_v283 import diamond_vertices, dot, require


ARTIFACT = Path(__file__).resolve().parents[1] / "results/v284-v287/moving-spiral-supports.json"
CAP_DEPTHS = (Q(1), Q(101, 100))
ARM_HEIGHTS = (Q(7, 8), Q(177, 200))
RADIUS_FLOOR = Q(1, 25)


def outer_arm(k):
    return ((Q(1), Q(0)), (Q(0), k), (-Q(3, 4), Q(0)),
            (Q(0), -Q(5, 8)), (Q(1, 2), Q(0)))


def walking_vertices(k):
    arm = outer_arm(k)
    return arm+tuple((-a, -b) for a, b in reversed(arm[1:]))


def coefficients(h, k, q):
    a, b = q
    B = b+h
    L = k*(1+a)+b
    return (L*(h*B-a), 2*h*a-2*k*h*B+(1-h*h)*L, -2*h*(1-k*h))


def polynomial_value(poly, r):
    c0, c1, c2 = poly
    return c0+c1*r+c2*r*r


def comparisons(h, k, q):
    """Arguments and factorizations for the five contact and orientation bounds."""
    a, b = q
    B = h+b
    L = k*(1+a)+b
    D, M = h-k, h+k
    v = b+k
    return (
        (L/M, D*L*((1-k*h)*(h*a-b)+h**3-k)/M**2),
        ((L+D)/M, D*(B-h*a)*(-(1-h*k)*(k*a+b)-h*k*k-2*h-k)/M**2),
        ((h*a+b+k)/(2*h), D*(h*(a-v)**2+v*(2*h*D-(h-1)**2*a))/(2*h)),
        ((h*a+b+h)/(2*h), D*(B-h*a)*(h*b-a-1)/(2*h)),
        (B/h, -(B-h*a)*(2*h-k+b-k*a)/h),
    )


def check_moving_domain(h, k):
    vertices = walking_vertices(k)
    start, terminal = vertices[0], vertices[-1]
    for a, b in vertices+diamond_vertices():
        q = (a, b)
        B = h+b
        L = k*(1+a)+b
        require(abs(a)+abs(b) <= 1 and b >= -k and 0 <= L <= 2*k
                and h*a+b+k >= 0 and B-h*a >= 0
                and -1 < h*a-b and h*b-a < 1,
                f"contact inequality failed at {(h, k, a, b)}")
        require((L == 0) == (q == terminal) and (B-h*a == 0) == (q == start),
                "zero boundary changed")
        G = coefficients(h, k, q)
        for i, (r, value) in enumerate(comparisons(h, k, q)):
            require(r >= 0 and value == polynomial_value(G, r), "comparison identity failed")
            require(value >= 0 if i in (0, 2) else value <= 0, "comparison sign failed")

    start_coeffs = coefficients(h, k, start)
    terminal_coeffs = coefficients(h, k, terminal)
    require(polynomial_value(start_coeffs, Q(1)) == 0 and start_coeffs[1]+2*start_coeffs[2] < 0
            and terminal_coeffs[0] == 0 and terminal_coeffs[1] < 0,
            "flat endpoints lost their selected root")


def support_audit():
    hmin, hmax = CAP_DEPTHS
    kmin, kmax = ARM_HEIGHTS
    require(Q(1, 2) < Q(5, 8) < Q(3, 4) < kmin <= kmax < 1 <= hmin <= hmax,
            "outer support ranges changed")
    require(hmax*kmax < 1 and 2*hmin*(hmin-kmax) > (hmax-1)**2
            and 2*hmin-kmax > 1, "scalar comparison bounds failed")
    squared_floor = (hmin-kmax)**2/(4*(hmax*hmax+1))
    require(squared_floor > RADIUS_FLOOR**2, "uniform radius bound failed")
    for h, k in product(CAP_DEPTHS, ARM_HEIGHTS):
        check_moving_domain(h, k)
    return {
        "cap_depth_interval": list(map(str, CAP_DEPTHS)),
        "first_arm_height_interval": list(map(str, ARM_HEIGHTS)),
        "inner_replacement_l1_radius_upper": "1/2",
        "uniform_circumradius_strict_lower": str(RADIUS_FLOOR),
        "squared_circumradius_lower": str(squared_floor),
        "flat_diagonal_dot_signs": [1, -1],
        "relative_index_on_whole_labeled_path": 1,
    }


def endpoint_differential(h, k, q, q_direction, r):
    """Differentiate in the parameter that points from the flat into the path."""
    _, c1, c2 = coefficients(h, k, q)
    a, b = q
    da, db = q_direction
    L, dL = k*(1+a)+b, k*da+db
    c0_derivative = dL*(h*(b+h)-a)+L*(h*db-da)
    c1_derivative = 2*h*da-2*k*h*db+(1-h*h)*dL
    root_derivative = c1+2*c2*r
    require(root_derivative < 0, "flat lacks a simple larger root")
    dr = -(c0_derivative+c1_derivative*r)/root_derivative
    dl = ((h+k)*dr-dL)/(h-k)
    corners = ((da-dl-dr, db+h*(dl-dr)), (-dl, h*dl), (dr, h*dr), (da, db))
    return root_derivative, corners


def rounding_audit(h=CAP_DEPTHS[0], k=ARM_HEIGHTS[0]):
    ends = (
        ("horizontal_flat", (Q(1), Q(0)), (-Q(1), k), Q(1),
         ((0, 1), (3, 2)), ((1, 0), (-1, 0))),
        ("vertical_flat", (Q(0), -k), (Q(3, 4), k), Q(0),
         ((0, 3), (1, 2)), ((0, 1), (0, 1))),
    )
    expected = {
        "horizontal_flat": (-2*(h-k), [(h+k)/(h-k), 1], (0, -(h+k))),
        "vertical_flat": (-2*k*h*(h-k), [k, 7*h/8], (Q(7, 4), 0)),
    }
    results = {}
    for name, q, direction, root, pairs, inward in ends:
        root_derivative, p = endpoint_differential(h, k, q, direction, root)
        velocities = []
        short_sides = []
        for (i, j), axis in zip(pairs, inward):
            velocities.append((dot(p[i], axis)+dot(p[j], axis))/2)
            short = tuple(b-a for a, b in zip(p[i], p[j]))
            require(dot(short, axis) == 0 and dot(short, short) > 0,
                    "missing transverse short-side derivative")
            short_sides.append(short)
        require(short_sides[0] == short_sides[1] and min(velocities) > 0,
                "flat tracing lacks inward longitudinal monotonicity")
        require((root_derivative, velocities, short_sides[0]) == expected[name],
                "flat derivative formula changed")
        results[name] = {"root_derivative": str(root_derivative),
                         "inward_midpoint_derivatives": list(map(str, velocities)),
                         "short_side_derivative": list(map(str, short_sides[0]))}
    return results


def report():
    return {
        "moving_supports": support_audit(),
        "rounded_collar_input_at_reference": rounding_audit(),
        "rounded_collar_inputs_at_parameter_box_corners": [
            {"h": str(h), "k": str(k), "flat_inputs": rounding_audit(h, k)}
            for h, k in product(CAP_DEPTHS, ARM_HEIGHTS)
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = json.dumps(report(), indent=2, sort_keys=True) + "\n"
    if args.write:
        ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
        ARTIFACT.write_bytes(data.encode("utf-8"))
    else:
        require(ARTIFACT.read_text(encoding="utf-8") == data, "moving-support artifact differs")
    print("Verified the moving-support inequalities and both flat-collar derivatives.")


if __name__ == "__main__":
    main()
