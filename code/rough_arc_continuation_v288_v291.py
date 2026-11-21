"""Exact certificates for continuing the concave rectangle path through rough arcs."""

import argparse
from fractions import Fraction as Q
from itertools import product
import json
from pathlib import Path

from big_peg_multichart_v249 import certify_sign
from moving_polygon_events_v224 import (
    derivative, determinant_poly, evaluate, iszero, minus, plus, polynomial, times,
)
from protected_rectangle_band_v284 import affine_matrix, concave_edges, lerp
from relative_witness_refinement_v282_v283 import diagonal_dot, difference, dot, require


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results/v288-v291/rough-arc-continuation.json"
DEN = polynomial((1, 3, 1))
NUMS = tuple(map(polynomial, ((1, -Q(3, 4), -Q(1, 4)), (0, 4, 1),
                              (1, 0, -1), (0, Q(10, 3), Q(5, 3)))))
BOTTOM_WINDOW = (Q(1, 2), Q(3, 4))
OPPOSITE_WINDOW = (Q(9, 100), Q(11, 100))
ANGLE_WINDOW = (Q(2, 5), Q(3, 5))
HALFWIDTH = Q(1, 100)
EPSILON = Q(1, 1000)
INVERSE_BOUND = Q(2)
COIL_BOX = ((Q(2), Q(3)), (-Q(1, 8), Q(1, 8)))
OPPOSITE_BOX = ((Q(22, 25), Q(23, 25)), (Q(117, 100), Q(123, 100)))


def baseline_parameters(x):
    d = evaluate(DEN, x)
    return tuple(evaluate(n, x)/d for n in NUMS)


def baseline_corners(x):
    return tuple(lerp(e, u) for e, u in zip(concave_edges(), baseline_parameters(x)))


def linear_polynomial_matrix():
    edges = concave_edges()
    a0, c0 = affine_matrix(edges, Q(0))
    a1, c1 = affine_matrix(edges, Q(1))
    matrix = tuple(tuple(polynomial((p, q-p)) for p, q in zip(r, s)) for r, s in zip(a0, a1))
    constant = tuple(polynomial((p, q-p)) for p, q in zip(c0, c1))
    return matrix, constant


def reference_audit():
    A, constant = linear_polynomial_matrix()
    require(determinant_poly(A) == times((Q(12),), DEN), "reference determinant changed")
    for row, c in zip(A, constant):
        total = times(c, DEN)
        for coefficient, n in zip(row, NUMS):
            total = plus(total, times(coefficient, n))
        require(iszero(total), "baseline solution does not satisfy the rectangle equations")
    signs = (-1, 1, -1, 1)
    for n, sign in zip(NUMS, signs):
        numerator = minus(times(derivative(n), DEN), times(n, derivative(DEN)))
        require(certify_sign(times((Q(sign),), numerator), Q(0), Q(1), strict=True),
                "contact projection ceased to be strictly monotone")
    require(baseline_parameters(0) == (1, 0, 1, 0)
            and baseline_parameters(1) == (0, 1, 0, 1), "flat endpoint tags changed")
    square = baseline_corners(Q(1, 2))
    side = difference(square[1], square[0])
    require(dot(side, side) == Q(145, 121) and diagonal_dot(square) == 0,
            "the interior square changed")
    return {
        "source_coordinate_monotonicity": list(signs),
        "flat_diagonal_dot_products": [str(diagonal_dot(baseline_corners(x))) for x in (0, 1)],
        "reference_square_squared_side": str(dot(side, side)),
        "unaltered_circumradius_lower": "1/2",
    }


def contact_polynomial(edge, q):
    """Scalar perpendicularity polynomial after eliminating the two midpoint equations."""
    a, b = q
    if edge == 0:
        c, d = 4-a-3*b, 5+a
        return (c*d, -(5*c+d+9*a), Q(5))
    if edge == 2:
        u = (4-a)/3
        return (u*(b-u), (2*u-b)/3-a, -Q(1, 9))
    raise ValueError("the supported anchor edges are 0 and 2")


def lifted_corners(edge, q, root):
    a, b = q
    if edge == 0:
        z = root
        u = (4-a-z)/3
        return (q, (4-3*u, u), (z, 3-2*z), (0, b+(5+a-5*z)/3))
    if edge == 2:
        p = root
        u = (4-a-p)/3
        return ((p, 0), (4-3*u, u), q, (0, b-u))
    raise ValueError("the supported anchor edges are 0 and 2")


def bottom_audit():
    (alo, ahi), (blo, bhi) = COIL_BOX
    x0 = baseline_corners(Q(1, 4))[0][0]
    require(x0 < alo and (alo, ahi) == tuple(4*u for u in BOTTOM_WINDOW),
            "bottom joins or phase bound changed")
    cmin = 4-ahi-3*bhi
    dmin = 5+alo
    derivative_upper = 10-(5*cmin+dmin+9*alo)
    # F(1/2) is largest when b is lowest. The resulting polynomial decreases in a.
    coefficient = Q(7, 2)-3*blo
    max_at_half = polynomial((coefficient*Q(5, 2), coefficient-7, -1))
    require(all(evaluate(derivative(max_at_half), a) < 0 for a in (alo, ahi)),
            "upper comparison at z=1/2 must decrease")
    upper = evaluate(max_at_half, alo)
    require(cmin*dmin > 0 and upper < 0 and derivative_upper < 0,
            "bottom anchor lacks a unique interior root")
    edge1_min = (4-ahi-Q(1, 2))/3
    require(bhi < edge1_min and (4-alo)/3 < 1
            and blo+(dmin-Q(5, 2))/3 > 0 and bhi+(5+ahi)/3 < 3,
            "bottom lift left its support segments")
    return {"anchor_rectangle": [[str(lo), str(hi)] for lo, hi in COIL_BOX],
            "scalar_root_interval": ["0", "1/2"],
            "minimum_at_lower_wall": str(cmin*dmin),
            "maximum_at_upper_wall": str(upper),
            "strict_derivative_upper": str(derivative_upper),
            "other_contact_edge2_parameter_strict_lower": "1/2",
            "circumradius_lower_on_inserted_arc": str(alo/2)}


def opposite_audit():
    (alo, ahi), (blo, bhi) = OPPOSITE_BOX
    lo, hi = Q(3, 20), Q(1, 2)
    u_min_lo = (4-ahi-lo)/3
    u_max = (4-alo-lo)/3
    u_min = (4-ahi-hi)/3
    lower = u_min_lo*(blo-u_max)-ahi*lo
    upper = bhi*bhi/4-alo*hi
    derivative_upper = -alo+Q(2, 3)
    require(lower > 0 and upper < 0 and derivative_upper < 0,
            "opposite anchor lacks a unique root in its guarded box")
    require(0 < u_min < u_max < 1 and 0 < blo-u_max and bhi-u_min < 3,
            "opposite lift left its support segments")
    first, last = (lerp(concave_edges()[2], s) for s in OPPOSITE_WINDOW)
    require(all(lo_coord < v < hi_coord for q in (first, last)
                for v, (lo_coord, hi_coord) in zip(q, OPPOSITE_BOX)),
            "opposite joins left the enclosure")
    require(OPPOSITE_WINDOW[1] < baseline_parameters(ANGLE_WINDOW[1])[2]
            and hi/4 < baseline_parameters(Q(1, 2))[0]-HALFWIDTH
            and hi/4 < BOTTOM_WINDOW[0],
            "opposite arc lost phase isolation")
    return {"anchor_rectangle": [[str(lo_coord), str(hi_coord)] for lo_coord, hi_coord in OPPOSITE_BOX],
            "source_interval": list(map(str, OPPOSITE_WINDOW)),
            "scalar_root_interval": [str(lo), str(hi)],
            "minimum_at_lower_wall": str(lower), "maximum_at_upper_wall": str(upper),
            "strict_derivative_upper": str(derivative_upper),
            "other_contact_edge0_parameter_strict_upper": str(hi/4),
            "circumradius_lower_on_inserted_arc": str(blo/2)}


def inverse_audit(matrix, low, high):
    """Bound the row sums of the adjugate by checking every signed sum."""
    bound = times((12*INVERSE_BOUND,), DEN)
    for i in range(4):
        cofactors = []
        for j in range(4):
            minor = tuple(tuple(v for c, v in enumerate(row) if c != i)
                          for r, row in enumerate(matrix) if r != j)
            cofactors.append(times((Q((-1)**(i+j)),), determinant_poly(minor)))
        for signs in product((-1, 1), repeat=4):
            remaining = bound
            for p, s in zip(cofactors, signs):
                remaining = minus(remaining, times((Q(s),), p))
            require(certify_sign(remaining, low, high, strict=True), "matrix inverse bound failed")


def overlap_audit():
    low, high = ANGLE_WINDOW
    matrix, _ = linear_polynomial_matrix()
    inverse_audit(matrix, low, high)
    ends = (baseline_parameters(low), baseline_parameters(high))
    middle = baseline_parameters(Q(1, 2))
    require(all(HALFWIDTH < u < 1-HALFWIDTH for end in ends for u in end),
            "angle-band tube hits a polygon vertex")
    gap = min(abs(end[i]-middle[i]) for end in ends for i in (0, 2))
    require(gap > 2*HALFWIDTH and 4*INVERSE_BOUND*EPSILON < HALFWIDTH,
            "modified support can reach an end plate or lateral boundary")
    require(max(end[0] for end in ends)+HALFWIDTH < BOTTOM_WINDOW[0] and
            min(end[2] for end in ends)-HALFWIDTH > OPPOSITE_WINDOW[1] and
            middle[2]+HALFWIDTH < Q(1, 2),
            "wide-phase replacements can intersect the simultaneous-contact tube")
    return {"angle_interval": list(map(str, ANGLE_WINDOW)),
            "patched_edges": [0, 2],
            "patches_centered_at_square_with_halfwidth": str(HALFWIDTH),
            "moving_contact_box_halfwidth": str(HALFWIDTH),
            "minimum_end_parameter_gap_from_patch_center": str(gap),
            "arc_coordinate_error_upper": str(EPSILON),
            "matrix_inverse_infinity_upper": str(INVERSE_BOUND),
            "lateral_margin_after_preconditioning": str(HALFWIDTH-4*INVERSE_BOUND*EPSILON),
            "opposite_edge_coordinate_gap_lower": "1",
            "whole_witness_circumradius_lower": str(Q(1, 2)-EPSILON),
            "tagged_aspect_ratios_from_flats": "all positive ratios"}


def simultaneous_coil_audit():
    """Two late-coil copies meet the unique square fiber of the coupled tube."""
    first_index = 4
    d = Q(1, 2**first_index)
    deviation_over_displacement = 2*d*d/(1-d*d)
    edges = concave_edges()
    maximum_edge_length = 4
    for i in (0, 2):
        direction = difference(*edges[i])
        require(dot(direction, direction) <= maximum_edge_length**2, "coil scaling bound failed")
    K = maximum_edge_length*deviation_over_displacement
    require(K*HALFWIDTH < EPSILON and 4*INVERSE_BOUND*K < 1,
            "anchored coils cannot be placed within the overlap certificate")
    center = baseline_parameters(Q(1, 2))
    for i in (0, 2):
        a, b = (lerp(edges[i], center[i]+sign*HALFWIDTH) for sign in (-1, 1))
        require(EPSILON < min(a[0], b[0]) <= max(a[0], b[0]) < 1-EPSILON,
                "coil enclosure may hit an adjacent edge")
        if i == 2:
            require(min(a[1], b[1])-EPSILON > 1, "opposite coil enclosure is too low")
    return {"first_retained_coil_index": first_index,
            "copy_source_length": str(HALFWIDTH),
            "physical_pointwise_error_over_source_distance_upper": str(K),
            "preconditioned_square_fiber_error_factor_upper": str(4*INVERSE_BOUND*K),
            "unique_square_in_tube_is_reference_square": True,
            "same_square_meets_both_coil_accumulation_points": True}


def tent_edges():
    """Each patched edge, with a midpoint displaced off the supporting line."""
    e = concave_edges()
    u = baseline_parameters(Q(1, 2))
    result = {}
    for i, offset in ((0, (0, -EPSILON/2)), (2, (EPSILON/2, 0))):
        center = lerp(e[i], u[i])
        kink = tuple(v+d for v, d in zip(center, offset))
        result[i] = (lerp(e[i], u[i]-HALFWIDTH), kink, lerp(e[i], u[i]+HALFWIDTH))
    return result


def adversarial_audit():
    tents = tent_edges()
    q = tents[0][1]
    z0 = baseline_corners(Q(1, 2))[2][0]
    F = contact_polynomial(0, q)
    lower_wall, upper_wall = z0-HALFWIDTH, z0+HALFWIDTH
    fx = tuple(evaluate(F, t) for t in (lower_wall, upper_wall))
    require(fx[0] > 0 > fx[1] and evaluate(derivative(F), 1) < 0,
            "one-anchor trap lacks a unique straight-contact root in the removed arc")
    require(4*INVERSE_BOUND*(EPSILON/2)/HALFWIDTH < 1,
            "tent preconditioned equation lost its contraction bound")
    fold = json.loads((ROOT / "results/v280/moving-concave-atlas.json").read_text(encoding="utf-8"))
    counts = tuple(len(values)-1 for values in fold["framed_relative_intervals_negative_zero_positive"])
    require(counts == (0, 1, 2) and fold["moving_components"][1]["regular_time_total_parity"] == 0,
            "fold control changed")
    return {"sequential_straight_root_bracket": list(map(str, (lower_wall, upper_wall))),
            "sequential_straight_root_wall_signs": [1, -1],
            "straight_root_contact_deleted_on_edge": 2,
            "square_in_coupled_tube_meets_both_tents": True,
            "preconditioned_tent_contraction_upper": str(4*INVERSE_BOUND*(EPSILON/2)/HALFWIDTH),
            "fold_control_additional_interval_counts": list(counts),
            "fold_control_regular_time_parity": 0}


def flat_collar_audit():
    """Transversality for the v286 rounding at the concave path's four vertices."""
    reference = concave_edges()
    data = {}
    for x, orientation, pairs, axes in (
        (0, 1, ((0, 1), (3, 2)), ((-4, 3), (4, -3))),
        (1, -1, ((0, 3), (1, 2)), ((1, 1), (1, 1))),
    ):
        denominator = evaluate(DEN, x)
        derivatives = []
        for n, e in zip(NUMS, reference):
            numerator = minus(times(derivative(n), DEN), times(n, derivative(DEN)))
            velocity = orientation*evaluate(numerator, x)/denominator**2
            derivatives.append(tuple(velocity*v for v in difference(e[1], e[0])))
        shorts = [difference(derivatives[j], derivatives[i]) for i, j in pairs]
        longitudinal = [(dot(derivatives[i], axis)+dot(derivatives[j], axis))/2
                        for (i, j), axis in zip(pairs, axes)]
        require(min(longitudinal) > 0 and shorts[0] == shorts[1] and dot(*shorts) > 0
                and all(dot(shorts[0], axis) == 0 for axis in axes),
                "concave flat collar failed transversality")
        data[str(x)] = {"short_side_derivative": list(map(str, shorts[0])),
                        "inward_longitudinal_projections": list(map(str, longitudinal))}
    return data


def report():
    return {"reference": reference_audit(), "flat_collar_inputs": flat_collar_audit(),
            "v288_bottom_lift": bottom_audit(),
            "v289_opposite_lift": opposite_audit(), "v290_overlapping_arc_continuation": overlap_audit(),
            "v290_two_coil_square": simultaneous_coil_audit(),
            "v291_sequential_obstruction_and_fold": adversarial_audit()}


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
        require(ARTIFACT.read_text(encoding="utf-8") == data, "rough-arc artifact differs from replay")
    print("Verified the two lifts, overlapping continuation tube, and countercontrols.")


if __name__ == "__main__":
    main()
