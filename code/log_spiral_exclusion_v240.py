"""v240 exact rational audit: open logarithmic spirals are square-free for a>=17/20.

The complete geometric proof is research/v240-log-spiral-seventeen-twentieths.md.
These rational checks audit theorem constants; code execution is NOT a proof
of the minimum-width and tangent-graph lemmas. An open arc is not a Jordan curve.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import exp, sin

THRESHOLD = F(17, 20)
CUTOFF = F(4, 15)
WIDTH_FLOOR = (1-CUTOFF)/2


def exp_lower(x, degree):
    """Positive Taylor partial sum, strictly less than exp(x) for x>0."""
    x = F(x)
    if x <= 0 or degree < 0:
        raise ValueError("positive argument and nonnegative degree required")
    term = total = F(1)
    for k in range(1, degree+1):
        term *= x/k
        total += term
    return total


def exact_margins():
    """All numerical inequalities for the universal real-a theorem at a=17/20."""
    a0, q = THRESHOLD, CUTOFF
    if not (3*q < 1):
        raise ArithmeticError("adjacent-inner separation not established")

    # The half-perimeter of a unit-radius inscribed regular 24-gon is
    # 24 sin(pi/24) = 12 sqrt(2-sqrt(2+sqrt(3))) < pi.
    # Exact bounds on the nested square roots show pi > 25/8.
    sqrt3_upper = F(1733, 1000)
    nested_root_upper = F(17807, 9216)
    if not (sqrt3_upper**2 > 3 and
            nested_root_upper**2 > 2+sqrt3_upper and
            nested_root_upper < 2 and
            144*(2-nested_root_upper) == F(25, 8)**2):
        raise ArithmeticError("24-gon lower bound pi>25/8 not established")

    # exp(53/40)>15/4 => log(15/4)<53/40.
    exp_log_bound = exp_lower(F(53, 40), 5)
    if not exp_log_bound > F(15, 4):
        raise ArithmeticError("logarithmic angular bound not established")

    # delta=log(15/4)/a <= log(15/4)/a0 < 53/34
    # < 25/16 < pi/2.
    delta_upper = F(53, 40)/a0
    if not (delta_upper == F(53, 34) and
            delta_upper < F(25, 16)):
        raise ArithmeticError("rotated graph angle bound not established")

    # tan(108/125) < 20/17 implies atan(20/17) > 108/125.
    # Alternating Taylor bounds at x=108/125 in (0,1) are strict.
    x = F(108, 125)
    sin_upper = x-x**3/6+x**5/120
    cos_lower = 1-x**2/2+x**4/24-x**6/720
    if not (0 < x < 1 and cos_lower > 0 and
            17*sin_upper < 20*cos_lower):
        raise ArithmeticError("lower bound on arctangent not established")

    # sqrt(1+(17/20)^2)>131/100 by exact squared comparison.
    sqrt_lower = F(131, 100)
    if not sqrt_lower**2 < 1+a0*a0:
        raise ArithmeticError("positive square root bound not established")

    # Exact transverse envelope at a0:
    # exp(-a0*atan(20/17))/sqrt(1+a0^2)
    # < exp(-459/625)/(131/100) < (1-q)/2=11/30.
    exponent_lower = a0*x
    exp_transverse_lower = exp_lower(exponent_lower, 4)
    if not (exponent_lower == F(459, 625) and
            WIDTH_FLOOR == F(11, 30) and
            exp_transverse_lower > F(3000, 1441) and
            F(3000,1441) == 1/(sqrt_lower*WIDTH_FLOOR)):
        raise ArithmeticError("strict transverse envelope margin not established")

    return {
        "adjacent_inner_cutoff": q,
        "three_outer_triangle_width_floor_over_R": WIDTH_FLOOR,
        "sqrt3_strict_upper": sqrt3_upper,
        "nested_root_strict_upper": nested_root_upper,
        "pi_strict_lower": F(25, 8),
        "exp_log_bound_partial_sum": exp_log_bound,
        "log_ratio_strict_upper": F(53, 40),
        "tangent_span_strict_upper": delta_upper,
        "atan_strict_lower": x,
        "sin_taylor_strict_upper": sin_upper,
        "cos_taylor_strict_lower": cos_lower,
        "sqrt_denominator_strict_lower": sqrt_lower,
        "exp_transverse_partial_sum": exp_transverse_lower,
        "exp_transverse_required_lower": F(3000, 1441),
    }


@dataclass(frozen=True)
class Certificate:
    a: F
    cutoff: F
    triangle_width_floor_over_R: F
    status: str

    @property
    def excluded(self):
        return self.status == "square-free-open-spiral"


def exclusion_certificate(a):
    """Fail closed for 0<a<17/20; never accept rounded float input."""
    if isinstance(a, (bool, float)) or not isinstance(a, (int, str, F)):
        raise TypeError("supply an integer, exact string or Fraction")
    a = F(a)
    if a <= 0:
        raise ValueError("a must be positive")
    if a < THRESHOLD:
        return Certificate(a, CUTOFF, WIDTH_FLOOR, "not-established")
    exact_margins()
    return Certificate(a, CUTOFF, WIDTH_FLOOR, "square-free-open-spiral")


def y_over_R(a, theta):
    """Floating-point illustration only, never used in the certificate."""
    return exp(-float(a)*theta)*sin(theta)


def report():
    certs = [exclusion_certificate(a) for a in
             (F(17, 20), F(851, 1000), F(9,10), F(1), F(4))]
    return {
        "theorem": "No nondegenerate inscribed square in S_a for all real a>=17/20",
        "curve": "S_a={exp((-a+i)t): t>=0} union {0}",
        "proof_checks": {k: str(v) for k,v in exact_margins().items()},
        "excluded_samples": [str(c.a) for c in certs if c.excluded],
        "smaller_a": "not established; v99 certifies a square-bearing pitch near 0.1002747",
        "scope": "compactified open arc, NOT a closed Jordan curve or Toeplitz resolution",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(report(), indent=2, sort_keys=True))
