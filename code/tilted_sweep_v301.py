"""Exact fixed-shear square census, mass obstruction and repair for v301."""

import argparse
from fractions import Fraction as Q
from itertools import combinations_with_replacement, permutations
import json
from math import comb
from pathlib import Path

from balanced_sweep_v300 import LinearMap, cumulative
from polygon_orbit_census_v223 import _affine_system, _evaluate, _solve, _square_residual, pt


ARTIFACT = Path(__file__).resolve().parents[1] / "results/v301/tilted-sweep.json"
WINNER = (1, 2, 3, 3)
CORNER_RANGES = ((0, Q(1, 256)), (Q(9, 17), Q(8, 15)),
                 (Q(15, 17), Q(8, 9)), (Q(7, 20), Q(71, 200)))
SAMPLE_TIMES = tuple(Q(i, 4) for i in range(5))
BERNSTEIN = [[comb(4, k) * t**k * (1-t)**(4-k) for k in range(5)] for t in SAMPLE_TIMES]
PERMUTATIONS = [(p, (-1)**sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4)))
                for p in permutations(range(4))]


def seed(drift):
    return tuple(pt(x, y + drift*x) for x, y in (
        (0, 0), (Q(1, 2), Q(-3, 8)), (1, 0), (Q(1, 2), Q(3, 4))))


def determinant(matrix):
    return sum(sign * matrix[0][p[0]] * matrix[1][p[1]] *
               matrix[2][p[2]] * matrix[3][p[3]] for p, sign in PERMUTATIONS)


def cramer(drift, chart):
    vertices = seed(drift)
    edges = [(vertices[i], vertices[(i+1) % 4]) for i in chart]
    matrix, rhs = _affine_system(edges, 1)
    nums = []
    for i in range(4):
        changed = [list(row) for row in matrix]
        for row, value in zip(changed, rhs):
            row[i] = value
        nums.append(determinant(changed))
    return (determinant(matrix), *nums)


def strict_sign(values):
    """Prove a sign on (0,1] for a polynomial of degree at most four."""
    c = _solve(BERNSTEIN, values)
    if min(c) >= 0 and c[-1] > 0:
        return 1
    if max(c) <= 0 and c[-1] < 0:
        return -1
    return 0


def audit_charts():
    records = []
    for chart in combinations_with_replacement(range(4), 4):
        if len(set(chart)) < 3:
            continue
        rows = [cramer(t/100, chart) for t in SAMPLE_TIMES]
        det = [r[0] for r in rows]
        orientation = strict_sign(det)
        if not orientation:
            raise ArithmeticError("unresolved determinant", chart)
        faces = []
        for i in range(4):
            numerator = [r[i+1] for r in rows]
            faces.append(strict_sign(numerator) * orientation)
            faces.append(strict_sign([d-n for d, n in zip(det, numerator)]) * orientation)
        if chart == WINNER:
            ordered = [r[4]-r[3] for r in rows]
            if min(faces) != 1 or strict_sign(ordered) != orientation:
                raise ArithmeticError("winner touches a contact wall")
            reason = "strict interior with ordered repeated-edge contacts"
        else:
            if -1 not in faces:
                raise ArithmeticError("chart lacks a uniform exclusion", chart)
            wall = faces.index(-1)
            reason = f"tag {wall//2} " + ("above one" if wall % 2 else "below zero")
        records.append({"edges": list(chart), "det_sign": orientation, "certificate": reason})
    return records


def winning_tags(d):
    D = 4*d*d + 12*d + 17
    return tuple(z/D for z in (4*d*d+12*d+1, 4, 5-4*d*d, 17-4*d*d))


def corner_x(d):
    a, b, c, e = winning_tags(d)
    return ((1-e)/2, (1+a)/2, 1-b/2, (1-c)/2)


def end_packet_cdf():
    return LinearMap([0, Q(1,16), Q(15,16), 1], [0, Q(9,32), Q(23,32), 1])


def red_fraction(H, xs):
    l, m, n, p = xs
    return (H(l) + H(m) + H(n) - H(p)) / 2


def mass_bounds(H):
    lo, hi = zip(*CORNER_RANGES)
    return (red_fraction(H, (*lo[:3], hi[3])), red_fraction(H, (*hi[:3], lo[3])))


def audit_ranges():
    # The denominator is positive and the cleared numerators are quadratic.
    for i, (lo, hi) in enumerate(CORNER_RANGES):
        lower, upper = [], []
        for t in SAMPLE_TIMES:
            d = t/100
            den = 4*d*d+12*d+17
            x = corner_x(d)[i]
            lower.append(den*(x-lo))
            upper.append(den*(hi-x))
        if strict_sign(lower) != 1 or strict_sign(upper) != 1:
            raise ArithmeticError("coordinate enclosure failed", i)


def audit_tube():
    # h is a source-edge parameter halfwidth, each edge has length < 1.
    # The worst drift is 1/100. The normalized estimates below are per unit drift.
    h = Q(1,8)
    inv = Q(400)
    ratio_change = Q(1, 20000)
    positional_change = Q(1, 100000)
    # The unperturbed first side has each component < 1 throughout the box.
    # At shape ratio < 2, each positional-error row costs at most 6 eps.
    used = inv * (ratio_change + 6*positional_change)
    if used >= h or Q(2,3) <= h:
        raise ArithmeticError("insufficient degree margin")
    max_tag_motion = h / 100
    mass_change = Q(9,2) * max_tag_motion  # each |dx| <= h*d/2
    e_hi = mass_bounds(end_packet_cdf())[1] + mass_change
    t_lo = mass_bounds(cumulative(1))[0] - mass_change
    if not e_hi < Q(1,2) < t_lo:
        raise ArithmeticError("the balancing mixture can exit [0,1]")
    return {"box_halfwidth_per_drift": str(h), "inverse_norm_upper": str(inv),
            "shape_band_halfwidth_per_drift": str(ratio_change),
            "independent_position_budget_per_drift": str(positional_change),
            "preconditioned_error_per_drift": str(used)}


def audit_annulus():
    eps = Q(1,24) + Q(1,100) + Q(1,1000)
    inner, outer = Q(29,100), Q(69,100)
    # Lower and upper supporting-line distances squared are 4/25 and 25/208.
    if not (inner+eps)**2 < min(Q(4,25), Q(25,208)) or Q(5,8)+eps >= outer:
        raise ArithmeticError("annulus clearance failed")
    if (outer / inner - 1)**2 >= 2:
        raise ArithmeticError("annulus too thick")
    side_lower = Q(41,100)
    if side_lower**2 >= 2*inner**2 or 2*eps >= inner:
        raise ArithmeticError("large square orientation clearance failed")
    if (side_lower - 2*eps)**2 <= Q(13,144) or side_lower <= 2*eps:
        raise ArithmeticError("source gap too short for the packet")
    return {"total_euclidean_displacement": str(eps), "inner_radius": str(inner),
            "outer_radius": str(outer), "square_side_squared_lower": str(2*inner**2),
            "source_gap_lower": "1/12", "sliding_density_cap": "13/2"}


def report():
    audit_ranges()
    # Clearing D leaves cubic contact residuals; five exact zeros suffice.
    for t in SAMPLE_TIMES:
        d = t / 100
        vertices = seed(d)
        edges = [(vertices[i], vertices[(i+1) % 4]) for i in WINNER]
        points = _evaluate(edges, winning_tags(d))
        if any(_square_residual(points, 1)):
            raise ArithmeticError("winning tag formula failed")

    bounds = {k: list(map(str, mass_bounds(H))) for k, H in
              (("old_wide", cumulative(0)), ("tall", cumulative(1)), ("end_packets", end_packet_cdf()))}
    xs = corner_x(Q(1,100))
    q_end = red_fraction(end_packet_cdf(), xs)
    q_tall = red_fraction(cumulative(1), xs)
    sample = {"drift": "1/100", "corner_x": list(map(str, xs)),
              "tall_weight_in_repaired_mixture": str((Q(1,2)-q_end)/(q_tall-q_end))}
    return {"positive_drift_charts": audit_charts(), "red_mass_bounds": bounds, "sample": sample,
            "tilted_tube": audit_tube(), "rough_family_annulus": audit_annulus()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = json.dumps(report(), indent=2, sort_keys=True) + "\n"
    if args.write:
        ARTIFACT.write_text(data, encoding="utf-8")
    elif ARTIFACT.read_text(encoding="utf-8") != data:
        raise SystemExit("v301 certificate differs from replay")
    print("Verified the tilted square census, weight bounds, degree margin and annulus.")
