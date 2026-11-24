"""v374: radially monotone wild points without symmetry.

Arm two is arm one turned, vertex by vertex, by pi + phi(r), where phi depends on the radius. The gap
path scales phi by lambda in [0, 1]: lambda = 0 is the half-turn-symmetric v370 drawing. Turned vertices keep
their radii and a positive angular gap, but the chords between them can cross once the turn spacing is below
the chord sagitta (the default drawing is not simple at depth 75), so parity counts check simplicity first.
"""
from fractions import Fraction as Q
from math import hypot, log, pi, sin, sqrt, tan

import numpy as np

from polygon_orbit_census_v223 import _intersect, pt
from wild_point_parity_v370 import DEN, _unit, arm, sector_counts


def assert_jordan(V):
    """Exact simplicity check: float bounding boxes prune pairs, `_intersect` decides the rest exactly.

    With a fixed number of vertices per turn, chords of neighbouring turns eventually cross once the turn
    spacing drops below the chord sagitta (for the defaults this happens before depth 75), so every
    parity computation checks simplicity first.
    """
    P = [pt(*p) for p in V]
    n = len(P)
    if len(set(P)) != n:
        raise ValueError('repeated vertex')
    A = np.array([[float(x), float(y)] for x, y in V])
    B = np.roll(A, -1, axis=0)
    lo, hi = np.minimum(A, B) - 1e-9, np.maximum(A, B) + 1e-9
    for i in range(n):
        m = np.all(lo[i + 2:] <= hi[i], axis=1) & np.all(hi[i + 2:] >= lo[i], axis=1)
        for j in (np.nonzero(m)[0] + i + 2):
            if i == 0 and j == n - 1:
                continue
            if _intersect(P[i], P[(i + 1) % n], P[j], P[(j + 1) % n]):
                raise ValueError(f'nonadjacent polygon edges {i} and {j} intersect')
    return V


def default_phi(r):
    """An oscillating, radius-dependent twist: no symmetry at any scale."""
    return 0.4 * sin(5 * log(1 / r))


def asymmetric_drawing(depth, lam=Q(1), phi=default_phi, a=0.5, theta0=3.0, per_turn=24, outer=Q(2)):
    lam = Q(lam)
    A = arm(a, theta0, theta0 + depth * pi, per_turn)
    B = []
    for x, y in reversed(A):
        r = hypot(float(x), float(y))
        u = Q(round(tan(float(lam) * phi(r) / 2) * DEN), DEN)
        c, sn = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
        B.append((-(c * x - sn * y), -(sn * x + c * y)))
    nrm = Q(round(sqrt(float(A[0][0] ** 2 + A[0][1] ** 2)) * DEN), DEN)
    far = lambda p: (p[0] * outer / nrm, p[1] * outer / nrm)
    ux, uy = _unit(theta0 + pi / 2)
    top = (ux * outer * 2, uy * outer * 2)
    return A[:-1] + [A[-1]] + B[:-1] + [B[-1], far(B[-1]), top, far(A[0])]


def large_graceful_parity(V, s):
    assert_jordan(V)
    counts, degenerate, squares = sector_counts(V, s)
    if degenerate:
        raise ArithmeticError("degenerate small-square candidate")
    return (1 + counts['aligned'] + counts['reversed']) % 2, counts, squares


def radius_over_side(entry):
    side_sq, pts = entry[0], entry[1]
    cx = sum(float(p[0]) for p in pts) / 4
    cy = sum(float(p[1]) for p in pts) / 4
    return hypot(cx, cy) / sqrt(float(side_sq))


def comb_radius(s, a=0.5, phi=None, lam=1.0):
    """Radius where adjacent arm passes of r = theta^(-a) are s apart.

    One turn of an arm moves the radius by about 2 pi a r^(1 + 1/a); the two arms split a turn in the
    angular gaps pi + lam phi(r) and pi - lam phi(r), so the closer pair of passes is a g r^(1 + 1/a) apart,
    g = pi - |lam phi(r)|. Small squares of side below s need passes at most about s apart.
    """
    r = (s / (a * pi)) ** (a / (1 + a))
    if phi is not None:
        for _ in range(60):
            r = (s / (a * (pi - abs(lam * phi(r))))) ** (a / (1 + a))
    return r


def locate_event(depth, s, lo, hi, steps=12, **kw):
    p_lo, _, sq_lo = large_graceful_parity(asymmetric_drawing(depth, lo, **kw), s)
    p_hi, _, sq_hi = large_graceful_parity(asymmetric_drawing(depth, hi, **kw), s)
    if p_lo == p_hi:
        raise ValueError("no parity change on this bracket")
    for _ in range(steps):
        mid = (lo + hi) / 2
        p_mid, _, sq_mid = large_graceful_parity(asymmetric_drawing(depth, mid, **kw), s)
        if p_mid == p_lo:
            lo, sq_lo = mid, sq_mid
        else:
            hi, sq_hi = mid, sq_mid
    return (lo, hi), sq_lo, sq_hi
