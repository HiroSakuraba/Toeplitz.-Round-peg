"""v239: exact-arithmetic inequalities for square-free open log spirals, a >= 9/10.

The geometric proof is in research/v239-log-spiral-nine-tenths.md.
This module certifies the rational strict margins, not the geometric lemmas.
The open arc is NOT a closed Jordan curve.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import exp, sin

THRESHOLD = F(9, 10)
INNER_CUTOFF = F(13, 50)
THREE_VERTEX_WIDTH = (1 - INNER_CUTOFF) / 2


def exp_lower(x, degree):
    """Strict lower bound on e**x when x>0, using positive Taylor terms."""
    x = F(x)
    if x <= 0 or degree < 0:
        raise ValueError("x must be positive and degree nonnegative")
    total = term = F(1)
    for n in range(1, degree + 1):
        term *= x / n
        total += term
    return total


def exact_margins():
    """Check each theorem-critical constant without floating point."""
    q = INNER_CUTOFF
    if not (2 * (2*q)**2 < (1-q)**2):
        raise ArithmeticError("two inner vertices are not excluded")

    # log(50/13)<27/20, since e**(27/20)>50/13.
    log_lower_exponential = exp_lower(F(27, 20), 6)
    if not log_lower_exponential > F(50, 13):
        raise ArithmeticError("angular-span bound not proven")

    # To lower bound arctan(10/9) by 4/5:
    # sin(4/5) < x-x^3/6+x^5/120; cos(4/5) > 1-x^2/2.
    # Both denominators are positive, and tan(4/5)<10/9.
    x = F(4, 5)
    sin_upper = x-x**3/6+x**5/120
    cos_lower = 1-x**2/2
    if not (cos_lower > 0 and sin_upper < F(10, 9)*cos_lower):
        raise ArithmeticError("arctan bound not proven")

    # sqrt(1+(9/10)^2) = sqrt(181)/10 > 4/3.
    if not (F(181, 100) > F(16, 9)):
        raise ArithmeticError("square-root bound not proven")

    # e**(18/25)>75/37, so (3/4)*e**(-18/25)<37/100.
    transverse_exponential = exp_lower(F(18, 25), 3)
    if not transverse_exponential > F(75, 37):
        raise ArithmeticError("transverse strip bound not proven")

    # delta <= log(50/13)/(9/10) < (27/20)/(9/10) = 3/2 < pi/2.
    # The final strict step is the elementary pi>3 inscribed-hexagon bound.
    if not F(27, 20)/THRESHOLD == F(3, 2):
        raise ArithmeticError("tangent cone bound not proven")

    return {
        "cutoff": q,
        "minimum_triangle_width_over_R": THREE_VERTEX_WIDTH,
        "two_inner_left_squared": 2*(2*q)**2,
        "two_inner_right_squared": (1-q)**2,
        "exp_27_20_lower": log_lower_exponential,
        "exp_27_20_target": F(50, 13),
        "tan_4_5_sine_upper": sin_upper,
        "tan_4_5_cosine_lower": cos_lower,
        "exp_18_25_lower": transverse_exponential,
        "exp_18_25_target": F(75, 37),
        "tangent_span_upper_at_threshold": F(3, 2),
        "strict_transverse_upper_over_R": THREE_VERTEX_WIDTH,
    }


@dataclass(frozen=True)
class ExclusionCertificate:
    a: F
    inner_cutoff: F
    minimum_triangle_width_over_R: F
    status: str

    @property
    def excluded(self):
        return self.status == "square-free-open-spiral"


def exclusion_certificate(a):
    """Return true only in the proven regime; no negative inference below it."""
    if isinstance(a, (float, bool)):
        raise TypeError("provide integer, Fraction or decimal/rational string, not float")
    if not isinstance(a, (int, str, F)):
        raise TypeError("a must be an integer, exact string or Fraction")
    a = F(a)
    if a <= 0:
        raise ValueError("a must be positive")
    if a < THRESHOLD:
        return ExclusionCertificate(a, INNER_CUTOFF, THREE_VERTEX_WIDTH,
                                    "not-established")
    exact_margins()
    return ExclusionCertificate(a, INNER_CUTOFF, THREE_VERTEX_WIDTH,
                                "square-free-open-spiral")


def transverse_coordinate(a, theta):
    """Floating-point sanity check only; not part of certificate."""
    return exp(-float(a)*theta)*sin(theta)


def report():
    m = exact_margins()
    samples = [F(9,10), F(901,1000), F(1), F(4)]
    return {
        "result": "For every real a >= 9/10, S_a has no nondegenerate inscribed square",
        "S_a": "{exp((-a+i)t): t>=0} union {0}",
        "status": "open-arc theorem, not a Jordan-curve result",
        "exact_margins": {k: str(v) for k,v in m.items()},
        "certified_samples": [str(a) for a in samples if exclusion_certificate(a).excluded],
        "below_threshold": "not established; do not infer existence or exclusion",
        "v99_control": "a special pitch near 0.1002747 has a square",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(report(), indent=2, sort_keys=True))
