"""Exact adversaries for the rectangle atlas, relative index and noescape gates."""

import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

from accumulating_notches_v220 import collision_chart_data, normalized_residual_at_collision
from big_peg_multichart_v249 import certify_sign
from big_peg_pair_death_v254 import QPOLY, ROOT_LEFT, ROOT_RIGHT
from moving_polygon_events_v224 import evaluate, iszero, pairdot, pairminus, plus, polynomial, times
from pentagon_complete_ledger_v263 import audit as pentagon_census
from polygon_orbit_census_v223 import (
    PositiveDimensionalSquares, census_polygon, pt, validate_jordan_polygon,
)
from rectangle_monodromy_v268 import (
    Rat, U, certify_circle_rectangle_loop, rectangle, rotate, unit_quarter,
)
from relative_rectangle_v269_v273 import (
    ellipse_relative_certificate, notch_square_certificate, relative_interval_parity,
)


MANIFEST = Path(__file__).resolve().parents[1] / "results" / "v278" / "adversarial-ledger.json"
END = Q(1, 100)
ZERO = polynomial((0,))
CANDIDATES = {
    "circle-continuum": ("G2", "Every positive-size rectangle component has finitely many square zeros."),
    "ellipse-interval": ("G3", "Every Jordan curve has an odd-monodromy geometric rectangle loop."),
    "false-closure": ("G1", "Equal center, radius and perpendicularity value certify rectangle-loop closure."),
    "bump-cell": ("G1", "Complete analysis in one four-edge cell exhausts the whole polygon."),
    "concave-vertex-pair": ("G3", "Every finite positive geometric-square orbit census is odd."),
    "convex-pair-death": ("G2", "A positive-size polygonal square pair can die only at matrix rank loss or a size cutoff."),
    "moving-notch": ("G5", "Uniform convergence of Jordan embeddings keeps every selected square above a size floor."),
    "fixed-collision-link": ("G2", "The bare normalized parameter-gap direction determines the normalized residual at collision."),
}


def _require(condition, message):
    if not condition:
        raise AssertionError(message)


def _points(points):
    return [[str(x), str(y)] for x, y in points]


def _dot(a, b):
    return a[0]*b[0] + a[1]*b[1]


def _sub(a, b):
    return (a[0]-b[0], a[1]-b[1])


def _square_side_squared(points):
    p, q, r, s = points
    u, v = _sub(q, p), _sub(s, p)
    _require(_sub(r, q) == v and v == (-u[1], u[0]) and _dot(u, u) > 0,
             "expected four geometric CCW corners of a genuine square")
    return _dot(u, u)


def _on_boundary(points, polygon):
    segments = tuple(zip(polygon, polygon[1:]+polygon[:1]))
    for p in points:
        _require(any(
            (p[0]-a[0])*(b[1]-a[1]) == (p[1]-a[1])*(b[0]-a[0])
            and _dot(_sub(p, a), _sub(p, b)) <= 0
            for a, b in segments
        ), f"point {p} has no segment contact")


def circle_continuum():
    """Verify an exact continuum and refuse its conversion into a finite count."""
    c, s = unit_quarter(U)
    corners = rotate(rectangle(1, 1), c, s)
    for i, j in zip(range(4), (1, 2, 3, 0)):
        side = tuple(corners[j][k]-corners[i][k] for k in range(2))
        _require(sum(q*q for q in side) == Rat.of(4), "circle side identity failed")
    for i in range(4):
        x, y = corners[i]
        _require(x*x+y*y == Rat.of(2), "circle membership identity failed")
    try:
        certify_circle_rectangle_loop((corners,), radius_squared=Q(2), strict_side_floor=Q(3))
    except ValueError as exc:
        _require("positive-dimensional interval" in str(exc), "unexpected finite loop rejection")
    else:
        raise AssertionError("circle square interval received a finite count")
    return {"square_zero_count": None, "square_side_squared": "4", "radius_squared": "2",
            "evidence": "rational rotation over a quarter turn, with every parameter a square"}


def ellipse_interval():
    cert = ellipse_relative_certificate(3, 2)
    _require(relative_interval_parity(*cert["flat_end_dot_signs"]) == 1,
             "ellipse lost relative endpoint parity")
    # v268 proves that its full geometric rectangle locus is an open interval.
    return {"semiaxes": ["3", "2"], "labelled_flat_signs": [-1, 1],
            "relative_interval_parity": 1, "whole_locus": "open interval",
            "proof_source": "research/v268-rectangle-monodromy-and-protected-square-criterion.md"}


def false_closure():
    """One arc has matching scalar endpoint data but different corner sets."""
    c, s = unit_quarter(U/2)
    corners = rotate(rectangle(3, 4), c, s)
    endpoints = tuple(tuple(tuple(x.at(t) for x in p) for p in corners) for t in (0, 1))

    def scalars(points):
        center = tuple((points[0][j]+points[2][j])/2 for j in range(2))
        d0, d1 = _sub(points[2], points[0]), _sub(points[3], points[1])
        return center, _dot(d0, d0)/4, _dot(d0, d1)

    _require(scalars(endpoints[0]) == scalars(endpoints[1]), "scalar endpoint data changed")
    _require(set(endpoints[0]) != set(endpoints[1]), "unexpected geometric closure")
    try:
        certify_circle_rectangle_loop((corners,), radius_squared=Q(25), strict_side_floor=Q(30))
    except ValueError as exc:
        _require("not same geometric rectangle" in str(exc), "unexpected loop rejection")
    else:
        raise AssertionError("different point sets were treated as a closed loop")
    center, radius, perpendicularity = scalars(endpoints[0])
    return {"both_centers": list(map(str, center)), "both_squared_radii": str(radius),
            "both_diagonal_dots": str(perpendicularity),
            "start": _points(endpoints[0]), "end": _points(endpoints[1])}


def _bump_polygon(h, eps):
    _require(0 < h < h+eps < 1, "the bump must lie strictly inside the top edge")
    return ((Q(0), Q(0)), (Q(1), Q(0)), (Q(1), Q(1)), (h+eps, Q(1)),
            (h+eps, 1+eps), (h, 1+eps), (h, Q(1)), (Q(0), Q(1)))


def bump_cell():
    """A square continuum in the selected chart, and a square outside it."""
    polygon = _bump_polygon(Q(3, 4), Q(1, 10))
    rational_polygon = validate_jordan_polygon(tuple(pt(*p) for p in polygon))
    t = polynomial((0, 1))
    one = polynomial((1,))
    rev = polynomial((1, -1))
    cont = ((one, t), (rev, one), (ZERO, rev), (t, ZERO))
    u, v = pairminus(cont[1], cont[0]), pairminus(cont[3], cont[0])
    _require(pairminus(cont[2], cont[1]) == v and (polynomial(-c for c in u[1]), u[0]) == v,
             "square-interval polynomial identity failed")
    exterior = tuple(polygon[i] for i in (0, 1, 2, 7))
    _on_boundary(exterior, polygon)
    side = _square_side_squared(exterior)
    try:
        census_polygon(rational_polygon)
    except PositiveDimensionalSquares:
        pass
    else:
        raise AssertionError("whole polygon square continuum received a finite census")
    return {"polygon": _points(polygon), "selected_cell_parameter_interval": ["2/5", "3/5"],
            "square_branch": "u=1-t", "square_side_squared_polynomial": list(map(str, pairdot(u, u))),
            "entire_square_orbit_count": None, "square_outside_selected_cell": _points(exterior),
            "outside_square_side_squared": str(side), "whole_rectangle_atlas_verified_here": False}


def concave_polygon(t):
    """Independent v274 fixture. The negatively turning tip is vertex 2."""
    t = Q(t)
    if abs(t) > END:
        raise ValueError("outside the certified perturbation interval")
    return ((Q(0), Q(0)), (Q(4), t), (Q(1), Q(1)), (Q(0), Q(3)))


def _poly_scale(p, c):
    return polynomial(c*x for x in p)


def _poly_sub(a, b):
    return plus(a, _poly_scale(b, -1))


def _strict_on_window(poly):
    # Existing Bernstein verifier expects a subset of [0,1].
    parameter = polynomial((-END, 2*END))
    reparameterized = ZERO
    for coefficient in reversed(poly):
        reparameterized = plus(times(reparameterized, parameter), polynomial((coefficient,)))
    return certify_sign(reparameterized, Q(0), Q(1), strict=True)


def _contact_checks(edge, point, denominator):
    """Check exact supporting-line membership and return endpoint inequalities."""
    x, y = point
    if edge == 0:
        line = _poly_sub(_poly_scale(y, 4), times(polynomial((0, 1)), x))
        box = (x, _poly_sub(_poly_scale(denominator, 4), x))
    elif edge == 2:
        line = _poly_sub(plus(y, _poly_scale(x, 2)), _poly_scale(denominator, 3))
        box = (x, _poly_sub(denominator, x))
    elif edge == 3:
        line = x
        box = (y, _poly_sub(_poly_scale(denominator, 3), y))
    else:
        raise ValueError("unsupported selected edge")
    _require(iszero(line), "incorrect polynomial contact identity")
    return box


def _concave_birth_charts():
    """Derive the two positive-radius germs by all-time polynomial signs."""
    def p(*coeffs):
        return polynomial(coeffs)

    return (
        (p(12, 1), ((ZERO, p(0, 3)), (p(12), p(0, 3)),
                    (p(12), p(12, 3)), (ZERO, p(12, 3))),
         (3, 0, 2, 3), ((0, 0), (2, 1))),
        (p(48, 4, 1), ((p(0, 12), p(0, 0, 3)), (p(48, 12), p(0, 12, 3)),
                       (p(48), p(48, 12, 3)), (ZERO, p(48, 0, 3))),
         (0, 0, 2, 3), ((0, 0), (2, 1))),
    )


def concave_vertex_pair():
    """Two branches coalesce at a two-polygon-vertex square on P(0)."""
    floors = Q(9, 10)
    # Edges 0 and 2 are separated in y by at least 1-END. Edges 1
    # and 3 are separated in x by at least 1. Check all successive turns.
    _require(0 < END < 1, "insufficient disjoint-edge gaps for P(t)")
    for turn_with_target_sign in ((12,), (4, -1), (5, 1), (3,)):
        _strict_on_window(polynomial(turn_with_target_sign))
    event_squares, positive_squares, branches = [], [], []
    for denominator, points, edges, critical in _concave_birth_charts():
        _strict_on_window(denominator)
        for i, (edge, point) in enumerate(zip(edges, points)):
            for end, inequality in enumerate(_contact_checks(edge, point, denominator)):
                if (i, end) in critical:
                    _require(inequality[0] == 0, "critical contact must vanish at the event")
                    _strict_on_window(polynomial(inequality[1:]))
                else:
                    _strict_on_window(inequality)
        u, v = pairminus(points[1], points[0]), pairminus(points[3], points[0])
        _require(pairminus(points[2], points[1]) == v
                 and v == (_poly_scale(u[1], -1), u[0]), "selected branch is not a square")
        _strict_on_window(_poly_sub(pairdot(u, u),
                                    _poly_scale(times(denominator, denominator), floors)))

        def get_points(t):
            d = evaluate(denominator, t)
            return tuple(tuple(evaluate(q, t)/d for q in point) for point in points)

        event_squares.append(get_points(Q(0)))
        positive_squares.append(get_points(END))
        branches.append({"denominator_in_ascending_order": list(map(str, denominator)),
                         "corner_numerators_in_ascending_order": [
                             [list(map(str, x)), list(map(str, y))] for x, y in points],
                         "edge_choices": edges, "critical_contact_indices": critical})
    _require(event_squares[0] == event_squares[1], "selected germs do not glue at common square")
    _require(set(positive_squares[0]) != set(positive_squares[1]), "born squares coincide")

    samples = []
    for t, count in ((-END, 1), (Q(0), 2), (END, 3)):
        polygon = concave_polygon(t)
        squares = census_polygon(validate_jordan_polygon(tuple(pt(*p) for p in polygon))).squares
        _require(len(squares) == count, "independent positive geometric orbit census disagrees")
        orbits = []
        for sq in squares:
            vertices = tuple((p.x, p.y) for p in sq.vertices)
            _on_boundary(vertices, polygon)
            _require(_square_side_squared(vertices) == sq.side_squared, "census side disagrees")
            orbits.append({"vertices_ccw": _points(vertices), "diagonal_pairing": [[0, 2], [1, 3]],
                           "size_squared": str(sq.side_squared),
                           "source_positions": list(map(str, sq.curve_positions)),
                           "source_order_sector": list(sq.cyclic_order_sector),
                           "support_labels": sq.labelings_checked,
                           "vertex_contact": sq.has_polygon_vertex,
                           "can_use_repeated_edge": sq.uses_repeated_edge})
        expected = () if t < 0 else (event_squares[0],) if t == 0 else positive_squares
        keys = {frozenset((p.x, p.y) for p in s.vertices) for s in squares}
        _require(all(frozenset(square) in keys for square in expected),
                 "a certified branch is missing from the independent census")
        samples.append({"time": str(t), "orbit_count": count, "geometric_orbits": orbits})
    return {"certified_branch_time_window": [str(-END), str(END)],
            "Jordan_nonadjacent_edge_coordinate_gaps": [str(1-END), "1"],
            "selected_branch_counts_negative_zero_positive": [0, 1, 2],
            "selected_side_squared_strict_lower": str(floors),
            "geometric_double_vertex_event_square": _points(event_squares[0]),
            "at_event_raw_orbit_count": samples[1]["orbit_count"],
            "at_event_total_support_labels": sum(o["support_labels"]
                                                 for o in samples[1]["geometric_orbits"]),
            "exact_selected_branches": branches, "independent_whole_polygon_slices": samples}


def convex_pair_death():
    cert = pentagon_census()
    _require(cert["all_real_time_exhaustiveness_proved"], "uncertified all-time census")
    counts = cert["all_real_time_positive_geometric_square_count"]
    order = ("before_algebraic_collision", "at_algebraic_collision", "after_algebraic_collision")
    _require(tuple(counts[stage] for stage in order) == (3, 2, 1), "pair death profile changed")
    _require(cert["all_real_time_every_positive_square_side_squared_exceeds_one"],
             "pair death was confused with size loss")
    return {"all_time_positive_orbits": counts, "event_polynomial": list(map(str, QPOLY)),
            "event_isolating_interval": [str(ROOT_LEFT), str(ROOT_RIGHT)],
            "strict_lower_for_all_side_squared": "1",
            "proof_source": "code/pentagon_complete_ledger_v263.py"}


def moving_notch(n=16):
    cert = notch_square_certificate(n)
    polygon = cert["polygon"]
    validate_jordan_polygon(tuple(pt(*p) for p in polygon))
    eps = Q(1, n)
    top = polygon[2:]
    length = Q(0)
    deviations = []
    nodes = []
    for i, (x, y) in enumerate(top):
        if i:
            a, b = top[i-1]
            length += abs(x-a)+abs(y-b)
        t = length/(1+2*eps)
        nodes.append([str(t), str(x), str(y)])
        deviations.append(max(abs(x-(1-t)), abs(y-1)))
    _require(length == 1+2*eps and max(deviations) <= 2*eps,
             "piecewise-linear uniform embedding bound failed")
    big_square = tuple(polygon[i] for i in (0, 1, 2, 7))
    return {"n": n, "top_edge_parameter_nodes": nodes,
            "uniform_coordinate_error_upper": str(2*eps),
            "selected_square_side_squared": str(_square_side_squared(cert["selected_square"])),
            "other_square_side_squared": str(_square_side_squared(big_square)),
            "other_square": _points(big_square), "common_limit": "unit-square Jordan embedding"}


def fixed_collision_link():
    checked_levels = (1, 4, 8)
    for n in checked_levels:
        square = collision_chart_data(n, "square")
        flat = collision_chart_data(n, "flat")
        _require(square[1:] == flat[1:], "gap directions or scales differ")
        square_residual = normalized_residual_at_collision(n, "square")
        flat_residual = normalized_residual_at_collision(n, "flat")
        _require(square_residual == (0,)*4 and flat_residual == (Q(37, 18), -Q(37, 18),
                                                                 Q(37, 6), -Q(37, 18)),
                 "normalized collision residual no longer matches proof")
    return {"gap_direction": ["1/3"]*3, "square_limit": list(map(str, square_residual)),
            "flat_control_limit": list(map(str, flat_residual)), "checked_levels": checked_levels,
            "curve_length": "37/6",
            "proof_source": "research/v220-fixed-rectifiable-collision-link-obstruction.md"}


def report():
    controls = {
        "circle-continuum": circle_continuum, "ellipse-interval": ellipse_interval,
        "false-closure": false_closure, "bump-cell": bump_cell,
        "concave-vertex-pair": concave_vertex_pair, "convex-pair-death": convex_pair_death,
        "moving-notch": moving_notch, "fixed-collision-link": fixed_collision_link,
    }
    _require(controls.keys() == CANDIDATES.keys(), "missing or unmapped adversarial control")
    return {"status": "PROVED_RESTRICTED", "specimens": [
        {"id": key, "gate": gate, "candidate": claim, "status": "DISPROVED_CANDIDATE",
         "evidence": controls[key]()}
        for key, (gate, claim) in CANDIDATES.items()
    ], "whole_curve_limit_claim": "OPEN", "unrestricted_toeplitz": "OPEN"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="recompute the checked-in ledger")
    mode.add_argument("--write", action="store_true", help="write the exact ledger")
    args = parser.parse_args()
    data = json.dumps(report(), indent=2, sort_keys=True) + "\n"
    if args.write:
        MANIFEST.parent.mkdir(parents=True, exist_ok=True)
        MANIFEST.write_bytes(data.encode("utf-8"))
    else:
        # Ignore newline conversion on Windows checkout.
        _require(MANIFEST.read_text(encoding="utf-8") == data, "adversarial ledger differs from replay")
    print(f"Verified {len(CANDIDATES)} exact adversarial controls.")


if __name__ == "__main__":
    main()
