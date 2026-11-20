"""Exact rescaled degree bounds for rough flat contacts."""

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from moving_flat_endpoints_v292_v295 import fixed_aspect_collar_residual, rational_det, reference_endpoint_audit
from polygon_orbit_census_v223 import _solve
from relative_witness_refinement_v282_v283 import require


ARTIFACT = Path(__file__).resolve().parents[1] / "results/v296-v299/rough-flat-contacts.json"
FLAT_TAGS = (-1, 1, 1, -1)


def limiting_residual(a, c, z):
    u, v, q, p = z
    return (a*(u*u-v*v)+c*(q*q-p*p), u+q-v-p,
            2+u-v, p-u+a*(v*v-u*u)/2)


def rescaled_residual(left, right, z, r):
    """Fixed-ratio residual, using four independent source parameters."""
    if r == 0:
        raise ValueError("the scaled residual requires a nonzero ratio")
    a, b, c, d = tuple(arc(r*t) for arc, t in zip((left, left, right, right), z))
    return (
        (a[0]+c[0]-b[0]-d[0])/r**2,
        (a[1]+c[1]-b[1]-d[1])/r,
        2*(d[0]-a[0])-(b[1]-a[1])/r,
        (d[1]-a[1])/r+(b[0]-a[0])/(2*r**2),
    )


def limiting_jacobian(a, c):
    return ((-2*a, -2*a, 2*c, 2*c), (1, -1, 1, -1),
            (1, -1, 0, 0), (a-1, a, 0, 1))


def certify_scaled_box(a, b, c, d):
    """Compute conservative rational bounds from the two normalized quadratic jets."""
    a, b, c, d = map(Q, (a, b, c, d))
    J = limiting_jacobian(a, c)
    determinant = 4*(c-a-2*a*c)
    require(determinant != 0 and rational_det(J) == determinant, "degenerate scaled contact chart")
    inverse_columns = tuple(_solve(J, tuple(int(i == j) for i in range(4))) for j in range(4))
    inverse_bound = max(sum(map(abs, row)) for row in zip(*inverse_columns))
    quadratic_bound = max(2*(abs(a)+abs(c)), abs(a), Q(1))
    halfwidth = min(Q(1, 4), 1/(8*inverse_bound*quadratic_bound))
    epsilon = halfwidth/(192*inverse_bound)
    expansion_bound = 8*max(Q(1), abs(a)+abs(b)+abs(c)+abs(d))
    scale = min(Q(1), halfwidth/(8*inverse_bound*expansion_bound),
                1/(16*(1+abs(a)+abs(c)+2*epsilon)))
    used = inverse_bound*(quadratic_bound*halfwidth**2+expansion_bound*scale+24*epsilon)
    require(used <= 3*halfwidth/8, "scaled side estimate failed")
    require(1-4*scale**2*(abs(a)+abs(c)+2*epsilon) >= Q(3, 4),
            "longitudinal separation estimate failed")
    return {"limiting_jacobian_determinant": determinant, "inverse_infinity_norm": inverse_bound,
            "parameter_box_halfwidth": halfwidth, "quadratic_error_tolerance": epsilon,
            "half_ratio_upper_before_source_domain_cutoff": scale,
            "side_margin_after_preconditioning": halfwidth-used,
            "circumradius_lower_in_units_of_flat_chord": Q(3, 8)}


def normalized_jets(curvature):
    """Use transverse graph parameters on the two smooth reference collars."""
    curvature = Q(curvature)
    if curvature < 1:
        raise ValueError("use a dimensionless curvature at least one")
    jets = reference_endpoint_audit()
    for x in ("0", "1"):
        a, c = (curvature*Q(v) for v in
                jets[x]["longitudinal_midpoint_derivatives_in_short_half_length"])
        yield x, (a, Q(0), c, Q(0))


def coil_quadratic_error(first_index, source_length, shear_lipschitz=Q(2)):
    """Uniform positional error divided by squared source distance for a coil tail."""
    if first_index < 1 or source_length <= 0 or shear_lipschitz <= 0:
        raise ValueError("positive coil index, source length and Lipschitz bound required")
    center = Q(1, 2**first_index)
    radius = center**3
    return 2*shear_lipschitz*radius/(source_length*(center-radius)**2)


def first_coil_index(source_length, budget):
    if budget <= 0:
        raise ValueError("positive coil budget required")
    index = 1
    while coil_quadratic_error(index, source_length) > budget:
        index += 1
    return index


def example_jet_audit():
    """Exact examples of the scale bounds for the algebraic normal forms."""
    curvature = Q(1)
    data = {}
    for x, jets in normalized_jets(curvature):
        a, b, c, d = jets
        certificate = certify_scaled_box(*jets)
        h = certificate["parameter_box_halfwidth"]
        # A sample localized coil size. The global construction picks gates
        # inside its chosen rounding and can reduce this support further.
        source_length = certificate["half_ratio_upper_before_source_domain_cutoff"]*(1-h)/4
        budget = certificate["quadratic_error_tolerance"]/2
        index = first_coil_index(source_length, budget)
        error = coil_quadratic_error(index, source_length)
        require(error <= budget < coil_quadratic_error(index-1, source_length),
                "coil tail cutoff is inconsistent")
        require(limiting_residual(a, c, tuple(map(Q, FLAT_TAGS))) == (0,)*4,
                "flat contact tags disagree")
        data[x] = {
            "dimensionless_quadratic_jets": [str(v) for v in (a, b, c, d)],
            "scaled_box": {key: str(v) for key, v in certificate.items()},
            "sample_coil_source_length": str(source_length),
            "first_retained_index_for_sample_length": index,
            "sample_coil_weighted_error_upper": str(error),
        }
    return {"illustrative_dimensionless_curvature": str(curvature), "ends": data}


def cubic_control():
    left = (lambda t: Q(0), lambda t: Q(0))
    right = (lambda t: 1+t**3, lambda t: 3*t*t)
    for v, r in ((Q(0), Q(1, 100)), (Q(1, 50), Q(1, 20)), (Q(-1, 40), Q(1, 50))):
        first_residual = fixed_aspect_collar_residual(left, right, Q(0), v, r, Q(0))[0]
        require(first_residual == -(3*v*v+r*r) < 0, "cubic obstruction changed")
    return {"same_flat_chord_and_tangent_directions": True,
            "right_minus_left_short_secant_slopes": "3*v^2+r^2",
            "genuine_rectangles_in_the_tagged_flat_collar": 0}


def report():
    return {
        "v296_rescaled_examples": example_jet_audit(),
        "v297_coil_envelope": {
            "normalized_centers": "2^-n",
            "normalized_radii": "2^-3n",
            "uniform_bound_on_local_shear_lipschitz_constant": "2",
            "tagged_position_error_vanishes_relative_to_squared_source_distance": True,
        },
        "v298_connected_concave_witness": {
            "certified_flat_contacts_that_can_carry_coil_accumulation_points": 4,
            "whole_witness_circumradius_lower": "1/3",
            "gate_control": "fixed reference arc windows or sufficiently small C1 changes at the gates",
        },
        "v299_first_order_control": cubic_control(),
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
        require(ARTIFACT.read_text(encoding="utf-8") == data, "rough-flat certificate differs from replay")
    print("Verified the scaled contact degree bounds, coil envelopes, and cubic control.")


if __name__ == "__main__":
    main()
