"""Rational degree bounds for a connected band of perturbed rectangles."""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from polygon_orbit_census_v223 import _affine_system, _solve, pt
from relative_witness_refinement_v282_v283 import difference, dot, require, rotate, spiral_audit, support


ARTIFACT = Path(__file__).resolve().parents[1] / "results/v284-v287/protected-angle-bands.json"
HALF = Q(1, 2)


def norm_inf(values):
    return max(map(abs, values))


def lerp(edge, t):
    a, b = edge
    return tuple(x+t*(y-x) for x, y in zip(a, b))


def rectangle_residual(corners, x):
    a, b, c, d = corners
    midpoint = tuple(u+v-w-z for u, v, w, z in zip(a, c, b, d))
    side, turned = difference(d, a), rotate(difference(b, a))
    return midpoint+tuple(x*u-(1-x)*v for u, v in zip(side, turned))


def affine_matrix(edges, x):
    starts = tuple(a for a, b in edges)
    constant = rectangle_residual(starts, x)
    columns = []
    for i, (_, end) in enumerate(edges):
        moved = starts[:i]+(end,)+starts[i+1:]
        columns.append(difference(rectangle_residual(moved, x), constant))
    return tuple(zip(*columns)), constant


def certify_band(edges, halfwidth, angle_halfwidth, tolerance, radius_floor):
    """Certify fixed parameter boxes and a uniform preconditioned degree margin."""
    if len(edges) != 4 or not (0 < halfwidth < HALF and 0 < angle_halfwidth < HALF
                               and tolerance > 0 and radius_floor > 0):
        raise ValueError("four edges and positive nontrivial certificate bounds are required")
    A, constant = affine_matrix(edges, HALF)
    center = _solve(A, tuple(-v for v in constant))
    require(center is not None, "the central rectangle chart is singular")
    require(all(halfwidth < t < 1-halfwidth for t in center),
            "the contact box hits an edge endpoint")
    square = tuple(lerp(edge, t) for edge, t in zip(edges, center))
    require(rectangle_residual(square, HALF) == (0,)*4, "central residual disagrees")
    square_edges = tuple(tuple(pt(*p) for p in edge) for edge in edges)
    square_system, square_rhs = _affine_system(square_edges, 1)
    require(_solve(square_system, square_rhs) == center, "independent square solve disagrees")

    inverse_columns = tuple(_solve(A, tuple(int(i == j) for i in range(4))) for j in range(4))
    inverse_bound = max(sum(map(abs, row)) for row in zip(*inverse_columns))
    edge_norms = tuple(norm_inf(difference(b, a)) for a, b in edges)
    a, b, c, d = square
    transverse = tuple(u+v for u, v in zip(difference(d, a), rotate(difference(b, a))))
    angle_bound = norm_inf(transverse)+halfwidth*(2*edge_norms[0]+edge_norms[1]+edge_norms[3])
    # Evaluate the uniform perturbation estimate in the fixed preconditioned box.
    perturbation = inverse_bound*(angle_halfwidth*angle_bound+4*tolerance)
    require(perturbation < halfwidth, "homotopy may hit a side of the contact box")
    diagonal_lower = norm_inf(difference(c, a))-halfwidth*(edge_norms[0]+edge_norms[2])-2*tolerance
    require(diagonal_lower > 2*radius_floor, "insufficient opposite-corner separation")
    x0, x1 = HALF-angle_halfwidth, HALF+angle_halfwidth
    return {
        "reference_edges": [[list(map(str, p)) for p in edge] for edge in edges],
        "square_edge_parameters": list(map(str, center)),
        "square_corners_ccw": [list(map(str, p)) for p in square],
        "square_side_squared": str(dot(difference(b, a), difference(b, a))),
        "contact_box_halfwidth": str(halfwidth),
        "angle_coordinate_interval": [str(x0), str(x1)],
        "maximum_long_to_short_ratio": str(x1/x0),
        "coordinate_perturbation_infinity_upper": str(tolerance),
        "circumradius_strict_lower": str(radius_floor),
        "inverse_infinity_norm": str(inverse_bound),
        "angle_derivative_bound": str(angle_bound),
        "preconditioned_side_margin": str(halfwidth-perturbation),
        "opposite_corner_distance_lower": str(diagonal_lower),
        "end_diagonal_dot_signs": [1, -1],
        "preconditioned_slice_degree": 1,
    }


def diamond_edges():
    fixed, outer = support()
    return fixed+(outer[:2],)


def concave_edges():
    vertices = ((Q(0), Q(0)), (Q(4), Q(0)), (Q(1), Q(1)), (Q(0), Q(3)))
    return tuple(zip(vertices, vertices[1:]+vertices[:1]))


def coil_enclosure_audit(concave):
    """Guard the imported finite-coil arc and the four unchanged support windows."""
    # The nth disk sits at (2+2**(-n), 0) with radius 2**(-3*n), n >= 1.
    # It stays to the right of x=2 and to the left of x=2+1/2+1/8.
    x_right = Q(2)+Q(1, 2)+Q(1, 8)
    y_upper = Q(1, 8)
    edge1_min_height = (4-x_right)/3
    require(edge1_min_height > y_upper, "a coil disk can meet the second polygon edge")
    center = Q(concave["square_edge_parameters"][0])
    window_right = 4*(center+Q(concave["contact_box_halfwidth"]))
    require(window_right < 2, "a coil disk can enter the first contact window")
    return {
        "replacement_segment": [["2", "0"], ["3", "0"]],
        "coil_box_x": ["2", str(x_right)],
        "coil_box_y": [str(-y_upper), str(y_upper)],
        "second_edge_height_gap": str(edge1_min_height-y_upper),
        "gap_to_first_tagged_window": str(2-window_right),
        "parent_band": "concave_quadrilateral",
        "analytic_arc_source": "Asano and Ike, Liouville regular Jordan curves, Proposition 6.6",
    }


def report():
    diamond = certify_band(diamond_edges(), Q(1, 1000), Q(1, 20000), Q(1, 30000), Q(9, 10))
    require(diamond["square_corners_ccw"] == spiral_audit()["one_shared_square_corners_ccw"],
            "the band is centered at a different square from v283")
    concave = certify_band(concave_edges(), Q(1, 100), Q(1, 1000), Q(1, 2000), Q(7, 10))
    return {"diamond_spiral": diamond, "concave_quadrilateral": concave,
            "finite_coil_in_concave_polygon": coil_enclosure_audit(concave)}


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
        require(ARTIFACT.read_text(encoding="utf-8") == data, "angle-band certificate differs from replay")
    print("Verified two rational continuation boxes for continuous arc perturbations.")


if __name__ == "__main__":
    main()
