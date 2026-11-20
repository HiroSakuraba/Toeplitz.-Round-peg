"""Exact endpoint algebra for the rounded concave rectangle witness."""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from moving_polygon_events_v224 import determinant_poly
from relative_witness_refinement_v282_v283 import difference, dot, require
from rough_arc_continuation_v288_v291 import baseline_corners, flat_collar_audit


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results/v292-v295/moving-flat-endpoints.json"


def chord_hessian(length, f_second, g_second):
    """Hessian of half the squared chord length at two horizontal normals."""
    return ((1-length*f_second, -1), (-1, 1+length*g_second))


def binormal_jacobian(length, f_second, g_second):
    """Center-coordinate derivative of the divided-difference collar equations."""
    return ((f_second, -g_second), (length*f_second-1, 1))


def normalized_contact_jacobian(s, k, s_over_k, a, b, c, d, twist, length):
    """Rows for short-vector matching and long/short perpendicular similarity.

    The variables are the means and half differences of the two pairs of
    source-parameter variations. The two indicated rows have been divided by
    k = the derivative of the old collar parameter with respect to the new one.
    """
    return (
        (s*twist, a, -s*twist, -c),
        (1, k*b, -1, -k*d),
        (-s*k*a-length, -s*s*k*twist-length*k*b, s*k*c, s*s*k*twist),
        (-s*b+length*s*twist, -s_over_k+length*a, s*d, s_over_k),
    )


def mean_and_slope(arc, center, r):
    value, derivative = arc
    if r == 0:
        return value(center), derivative(center)
    low, high = value(center-r), value(center+r)
    return (high+low)/2, (high-low)/(2*r)


def fixed_aspect_collar_residual(left, right, u, v, r, ratio):
    """Append the short-to-long ratio in the fixed transverse coordinates."""
    mf, df = mean_and_slope(left, u, r)
    mg, dg = mean_and_slope(right, v, r)
    longitudinal = mg-mf
    return (df-dg, longitudinal*df+v-u, 2*r-ratio*longitudinal)


def rational_det(rows):
    return determinant_poly(tuple(tuple((Q(v),) for v in row) for row in rows))[0]


def reference_endpoint_audit():
    jets = flat_collar_audit()
    result = {}
    for x, other, b_axis_sign in ((0, 3, -1), (1, 1, 1)):
        points = baseline_corners(x)
        chord = difference(points[other], points[0])
        length_sq = dot(chord, chord)
        data = jets[str(x)]
        short = tuple(map(Q, data["short_side_derivative"]))
        require(dot(short, short) == length_sq and dot(short, chord) == 0,
                "wrong flat direction or short-length scale")
        a, c = (2*Q(p)/length_sq for p in data["inward_longitudinal_projections"])
        c *= b_axis_sign
        require(a != 0 and c != 0, "rounding transition loses longitude rank")
        linear_over_length = 2*(c-a)
        quadratic = -4*length_sq*a*c
        if x == 0:
            # L=5. The Hessian determinant is positive for every kappa >= 1.
            require(length_sq == 25 and quadratic > 0 and 5*linear_over_length+quadratic > 0,
                    "first binormal Hessian loses nondegeneracy")
            sign = 1
        else:
            require(linear_over_length < 0 and quadratic < 0,
                    "second binormal Hessian loses nondegeneracy")
            sign = -1
        result[str(x)] = {
            "squared_flat_chord_length": str(length_sq),
            "longitudinal_midpoint_derivatives_in_short_half_length": [str(a), str(c)],
            "hessian_linear_coefficient_divided_by_flat_chord_length": str(linear_over_length),
            "hessian_quadratic_coefficient": str(quadratic),
            "hessian_determinant_sign_for_curvature_parameter_at_least_one": sign,
            "normalized_transition_limit_determinant": str(length_sq*a*c),
        }
    return result


def parallel_tilt_control():
    left = (lambda t: Q(0), lambda t: Q(0))
    displacement = Q(1, 10000)
    right = (lambda t: 5+displacement*t, lambda t: displacement)
    for r in (Q(0), Q(1, 100), Q(1, 10)):
        require(fixed_aspect_collar_residual(left, right, Q(0), Q(0), r, Q(0))[0] == -displacement,
                "tilted parallel collar residual changed")
    return {"illustrative_tilt_slope": str(displacement),
            "first_collar_residual_for_any_centers_or_projected_half_length": str(-displacement),
            "binormals_between_the_tilted_straight_arcs": 0,
            "root_degree_after_tilt": 0}


def fold_control():
    atlas = json.loads((ROOT / "results/v280/moving-concave-atlas.json").read_text(encoding="utf-8"))
    model = atlas["moving_components"][1]
    counts = [len(values)-1 for values in atlas["framed_relative_intervals_negative_zero_positive"]]
    require(model["time_map"] == "t = abs(s)" and counts == [0, 1, 2]
            and model["regular_time_total_parity"] == 0, "v281 fold control differs")
    return {"geometric_interval_counts": counts, "transverse_fold_degree": 0,
            "positive_time_local_degrees": [-1, 1]}


def report():
    return {
        "v292_exact_endpoint_jets": reference_endpoint_audit(),
        "v293_local_equations": ["D_r f(u) - D_r g(v)",
                                 "(M_r g(v) - M_r f(u)) D_r f(u) + v - u",
                                 "2r - lambda (M_r g(v) - M_r f(u))"],
        "v294_rounded_concave_witness": {
            "endpoint_arc_topology_for_connected_witness": "C1",
            "endpoint_arc_topology_for_unique_local_branches": "C2",
            "interior_topology_with_compatible_source_parameters": "C0",
            "uniform_circumradius_lower": "1/3",
            "perturbation_tolerance": "positive, set by compactness",
        },
        "v295_parallel_tilt": parallel_tilt_control(),
        "v295_positive_radius_fold": fold_control(),
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
        require(ARTIFACT.read_text(encoding="utf-8") == data, "moving flat artifact differs from replay")
    print("Verified both flat jets and the tilt and fold controls.")


if __name__ == "__main__":
    main()
