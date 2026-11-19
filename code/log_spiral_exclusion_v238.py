"""v238: exact rational audit of the a >= 1 open logarithmic-spiral theorem.

For S_a = {exp((-a+i)t):t>=0} union {0}, a >= 1 implies no
nondegenerate Euclidean square is inscribed in S_a.

This is a code audit of exact numerical inequalities, NOT an automated
verification of the geometric lemmas; see research/v238-log-spiral-threshold-one.md.
S_a is an arc, not a closed Jordan curve.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import cos, exp, sin

INNER_RADIUS = F(1, 4)
PROVED_THRESHOLD = F(1)


def positive_exponential_partial_sum(x, last_degree):
    """Exact strict lower bound on exp(x) when x > 0."""
    x = F(x)
    if x <= 0 or last_degree < 0:
        raise ValueError("x must be positive and degree nonnegative")
    term = F(1)
    out = term
    for k in range(1, last_degree + 1):
        term *= x / k
        out += term
    return out


@dataclass(frozen=True)
class Exclusion238:
    openness: F
    inner_cutoff: F
    triangle_width_floor_over_R: F
    strip_width_strict_upper_over_R: F
    tangent_half_span_strict_upper: F
    status: str

    @property
    def excluded(self):
        return self.status == "square-free-open-spiral"


def exclusion_certificate(a):
    """Audit exact sufficient inequalities for all rational a >= 1.

    A float is rejected: construct a Fraction or decimal string to preserve
    an intended exact boundary (e.g. '0.999999999999999999999').
    For 0<a<1, return NOT ESTABLISHED, never square-bearing.
    """
    if isinstance(a, (float, bool)):
        raise TypeError("use an integer, exact string, or Fraction, not float/bool")
    if not isinstance(a, (int, str, F)):
        raise TypeError("a must be an integer, exact string, or Fraction")
    a = F(a)
    if a <= 0:
        raise ValueError("a must be positive")

    eps = INNER_RADIUS
    triangle_floor = (1 - eps) / 2
    # For each outer point: y/R = exp(-a*theta)sin(theta)
    # <= theta exp(-a*theta) <= 1/(a*e) < 3/(8*a).
    strip_upper = F(3, 8) / a
    # Tangent angles lie in an interval of length delta=log(4)/a.
    # Rotate the graph axis to its midpoint: half-span
    # delta/2 < 3/(4*a), which is < pi/4 when a>=1.
    half_span_upper = F(3, 4) / a

    if a < PROVED_THRESHOLD:
        return Exclusion238(a, eps, triangle_floor, strip_upper,
                            half_span_upper, "not-established")

    # e > sum_{k=0}^3 1/k! = 8/3, so 1/e < 3/8.
    if positive_exponential_partial_sum(F(1), 3) != F(8, 3):
        raise ArithmeticError("failed exact exponential lower bound")
    # exp(3/2) > 67/16 > 4, so log(4) < 3/2.
    if positive_exponential_partial_sum(F(3, 2), 3) != F(67, 16):
        raise ArithmeticError("failed exact logarithm upper bound")
    if not (F(67, 16) > 4):
        raise ArithmeticError("log(4) bound failed")

    # Any two corners in the inner disk have distance <2eps*R.
    # A farthest and inner vertex imply s >= (1-eps)R/sqrt(2).
    if not (2 * (2 * eps) ** 2 < (1 - eps) ** 2):
        raise ArithmeticError("inner-pair separation bound failed")

    # Strict real strip width < strip_upper*R <= width_floor*R.
    if not (strip_upper <= triangle_floor):
        raise ArithmeticError("outer triangle strip bound failed")

    # Strict real tangent half-span < half_span_upper <= 3/4
    # and pi/4 > 3/4 by the elementary geometric fact pi>3.
    if not (half_span_upper <= F(3, 4)):
        raise ArithmeticError("rotated Lipschitz graph bound failed")

    return Exclusion238(a, eps, triangle_floor, strip_upper,
                        half_span_upper, "square-free-open-spiral")


def transverse_coordinate(a, theta):
    """Floating-point visual sanity check only, never used for certification."""
    return exp(-a * theta) * sin(theta)


def report():
    choices = [F(1), F(1001, 1000), F(4), F(100)]
    certs = [exclusion_certificate(a) for a in choices]
    return {
        "theorem": "For all real a>=1, S_a contains no nondegenerate square",
        "curve": "S_a={exp((-a+i)t):t>=0} union {0}",
        "scope": "compactified open arc only; not a closed Jordan curve",
        "proof_inputs": {
            "exp_1_lower_partial_sum": str(positive_exponential_partial_sum(1, 3)),
            "exp_3_over_2_lower_partial_sum":
                str(positive_exponential_partial_sum(F(3, 2), 3)),
            "pi_lower_from_inscribed_hexagon": "pi>3",
            "two_inner_squared_comparison": "1/2 < 9/16",
            "outer_triangle_width_floor_over_R": "3/8",
            "outer_strip_strict_upper_over_R_at_a1": "3/8",
            "rotated_tangent_half_span_strict_upper_at_a1": "3/4",
        },
        "certified_parameters": [str(c.openness) for c in certs if c.excluded],
        "below_one": "not established; v99 proves one special pitch has a square",
        "limits": "Neither arbitrary a in (0,1) nor Jordan closures are decided.",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(report(), indent=2, sort_keys=True))
