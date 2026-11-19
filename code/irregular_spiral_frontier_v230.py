"""v230 exact-rational *conditional* irregular power-spiral benchmark.

This module proves elementary inequalities and regularity thresholds; it does
NOT certify the transcendental Poincare-Miranda face signs needed for square
existence. Those belong to the restored v99/v101 verifiers. In particular,
claims about the square are explicitly CONDITIONAL on the v100/v101 theorem.

Gamma_p(t)=t exp(i/t**p) (0<t<=1). D_p is the union of Gamma_p, -Gamma_p,
and the outer unit semicircle joining their endpoints. For p>=2, the v100
Box A is c in [1,19/10], d in [9/20,18/25], q=1/(p+1).
"""
from dataclasses import dataclass
from fractions import Fraction as F


@dataclass(frozen=True)
class RegularityProfile:
    exponent: F
    rectifiable: bool
    finite_variation_below_two: bool
    liouville_regular: bool
    box_a_applies: bool


def profile(p):
    p = F(p)
    if p <= 0:
        raise ValueError('power exponent must be positive')
    return RegularityProfile(p, p < 1, p < 2, p < 2, p >= 2)


def finite_q_variation(p, q):
    """Exact sharp threshold for the compactified spiral in q>=1.

    For p<1, finite total length implies every q>=1 works. For p>=1,
    finite q-variation iff q>p. The same applies to the double spiral,
    since finitely many concatenated arms and an outer arc preserve finiteness.
    """
    p, q = F(p), F(q)
    if p <= 0 or q < 1:
        raise ValueError('p>0 and variation exponent q>=1 required')
    return p < 1 or q > p


@dataclass(frozen=True)
class UniformSquareBenchmark:
    """Rational bounds, CONDITIONALLY transferred from v100 Box A.

    Mathematically: every p>=2 Box-A Miranda square has minimum vertex
    radius STRICTLY > 1/5, maximum vertex radius STRICTLY < 1 (v100 section 9),
    and side length STRICTLY > 1/4. (The original v230 statement was > 1/8;
    the audit strengthened it.) No floating-point approximation or unverified
    root coordinate is used in these inequalities.
    """
    min_vertex_radius_strictly_greater_than: F
    side_length_strictly_greater_than: F
    cutoff_radius_that_preserves_all_vertices: F
    c_min: F
    c_max: F
    d_min: F
    d_max: F
    needs_external_face_certificate: bool
    max_vertex_radius_strictly_less_than: F = F(1)


def conditional_box_a_bounds(p):
    if F(p) < 2:
        raise ValueError('v230 uniform bounds only cover p>=2 (v100 Box A)')
    return UniformSquareBenchmark(
        F(1, 5), F(1, 4), F(1, 5),
        F(1), F(19, 10), F(9, 20), F(18, 25), True, F(1))


def protected_by_box_a(p, cutoff_radius):
    """Sufficient certificate for keeping the *known* v100 square in truncations.

    Returns False when outside the certified region; False is NOT evidence
    that the spiral fails to contain an inscribed square after truncation.
    """
    return F(p) >= 2 and F(0) < F(cutoff_radius) <= F(1, 5)


def proof_inequality_checks():
    """Validate all purely rational comparisons used in the theorem.

    Set q=1/(p+1), p>=2: 1-q>=2/3, c>=1, Delta<3*pi<10.
    exp(-x)<=1/(1+x), so 1-exp(-c(1-q))>=2/5;
    t0**p=(1-exp(-c(1-q)))/Delta>1/25.
    Thus t0>25**(-1/p)>=1/5. delta=d*q<=6/25<pi/2.

    Side (audit, replaces the original >1/8 step): with rho>=1 and
    cos(delta)>=1-delta^2/2>0,
        |1+rho e^(i delta)|^2 = 1+rho^2+2 rho cos(delta) >= 2+2cos(delta)
                               >= 4-delta^2 >= 4-(6/25)^2 = 2464/625 > 25/8,
    so side^2 = t0^2 |1+rho e^(i delta)|^2 / 2 > (1/25)(25/8)/2 = 1/16.

    Outer vertex (v100 section 9): t3^p = (exp(c(1-q))-1)/Delta with
    exp(c(1-q)) <= exp(19/10) < e^2 < (272/100)^2 < 37/5 and
    Delta = 3*pi - d*q > 9 - 6/25 = 219/25, so t3^p < (32/5)/(219/25) < 1.
    """
    assert F(18, 25) / 3 == F(6, 25)
    assert F(2, 3) / (1 + F(2, 3)) == F(2, 5)
    assert F(2, 5) / F(10) == F(1, 25)
    # pi>3 and 3*pi<10 are standard Archimedean estimates.
    assert F(6, 25) < F(3, 2)
    assert 3 * F(22, 7) < 10
    # side > 1/4
    assert 4 - F(6, 25) ** 2 == F(2464, 625)
    assert F(2464, 625) > F(25, 8)
    assert F(1, 25) * F(25, 8) / 2 == F(1, 16)
    # outer vertex radius < 1  (e < 2.72 is a standard estimate)
    assert F(272, 100) ** 2 < F(37, 5)
    assert 3 * 3 - F(6, 25) == F(219, 25)
    assert (F(37, 5) - 1) / F(219, 25) < 1
    return True


def report():
    proof_inequality_checks()
    return {
        'profiles': {str(p): {'rectifiable': profile(p).rectifiable,
                              'finite_q_variation_for_some_q_lt_2': profile(p).finite_variation_below_two,
                              'liouville_regular': profile(p).liouville_regular}
                     for p in (F(1, 2), F(1), F(3, 2), F(2), F(5, 2), F(3))},
        'non_liouville_family': 'p>=2',
        'conditional_uniform_radius_floor': '>1/5',
        'conditional_uniform_square_side_floor': '>1/4',
        'conditional_outer_vertex_radius_ceiling': '<1',
        'source_of_square_existence': 'v100 theorem, corrected/reproduced by v101 Arb (not computed by this script)',
        'scope': 'double spiral D_p and powers; NOT arbitrary Jordan curves',
    }


if __name__ == '__main__':
    import json
    print(json.dumps(report(), indent=2, sort_keys=True))
