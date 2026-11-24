"""v375: the comb law for two spiral arms, as a proved bound, checked on exact drawings.

Arms: Theta -> Theta^(-a) e^{i Theta} and Theta -> Theta^(-a) e^{i(Theta + pi + phi(r))}, Theta >= theta0.
Hypotheses: |phi| <= phibar < pi, |r phi'(r)| <= L, g = pi - phibar, R = theta0^(-a).
"""
from fractions import Fraction as Q
from math import atan, pi, sqrt


def spacing_bound(r, a=0.5, g=pi):
    """d_g(r): two points of the arc at radius >= r that are not neighbours along one arm are this far apart."""
    K = a * r ** (1 / a)
    return g * K * r / (1 + pi * K / 2)


def confinement_radius(t, a=0.5, g=pi, hi=10.0):
    """rho(t) = d_g^{-1}(t): every small square of side t on the arc has a corner within this radius."""
    lo = 0.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if spacing_bound(mid, a, g) < t:
            lo = mid
        else:
            hi = mid
    return hi


def hypotheses(t, rho, a=0.5, theta0=3.0, L=0.0):
    """The smallness conditions of Theorem 375.1 for a square of side t whose corners lie at radius >= rho."""
    R = theta0 ** (-a)
    KR = a * R ** (1 / a)
    eps1 = KR * L
    tau = atan(KR / (1 - eps1)) if eps1 < 1 else pi / 2
    eta = (pi / sqrt(2) + sqrt(2) * L) * t / rho
    ok = KR <= 0.25 and eps1 < 1 and eta <= pi / 4 and eta * (1 + eps1) + tau < pi / 2
    return dict(eta=eta, eps1=eps1, tau=tau, KR=KR, ok=ok)


def corner_radii(entry):
    return [sqrt(float(p[0]) ** 2 + float(p[1]) ** 2) for p in entry[1]]


def check_drawing(V, s, a=0.5, g=pi):
    """For every square of side below s: (side, min corner radius, rho(side))."""
    from wild_point_parity_v370 import sector_counts
    _, degenerate, squares = sector_counts(V, Q(s))
    if degenerate:
        raise ArithmeticError("degenerate small-square candidate")
    rows = []
    for e in squares.values():
        t = sqrt(float(e[0]))
        rows.append((t, min(corner_radii(e)), confinement_radius(t, a, g)))
    return rows
