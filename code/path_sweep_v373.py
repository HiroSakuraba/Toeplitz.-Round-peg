"""v373: the path sweep, checked on an exact twist path.

The twist path turns arm two of the v370 double spiral by pi + phi, where u = tan(phi / 2) is a
rational parameter, so every drawing on the path is exact. On a Jordan polygon the large graceful
parity at cutoff s is 1 + (small aligned) + (small reversed) (Lemma 370.1).
"""
from fractions import Fraction as Q
from math import pi, sqrt

from wild_point_parity_v370 import DEN, _unit, arm, sector_counts


def twisted_drawing(depth, u, a=0.5, theta0=3.0, per_turn=24, outer=Q(2)):
    """v370 drawing whose arm two is arm one turned by pi + phi, with tan(phi/2) = u (rational)."""
    u = Q(u)
    c, sn = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
    A = arm(a, theta0, theta0 + depth * pi, per_turn)
    B = [(-(c * x - sn * y), -(sn * x + c * y)) for x, y in reversed(A)]
    nrm = Q(round(sqrt(float(A[0][0] ** 2 + A[0][1] ** 2)) * DEN), DEN)
    far = lambda p: (p[0] * outer / nrm, p[1] * outer / nrm)
    ux, uy = _unit(theta0 + pi / 2)
    top = (ux * outer * 2, uy * outer * 2)
    return A[:-1] + [A[-1]] + B[:-1] + [B[-1], far(B[-1]), top, far(A[0])]


def large_graceful_parity(V, s):
    counts, degenerate, squares = sector_counts(V, s)
    if degenerate:
        raise ArithmeticError("degenerate small-square candidate")
    return (1 + counts['aligned'] + counts['reversed']) % 2, squares


def locate_event(depth, s, u_lo, u_hi, steps=12, **kw):
    """Bisect the twist parameter between two drawings of different large graceful parity.
    Returns the final bracket and the graceful squares of both endpoint drawings."""
    p_lo, sq_lo = large_graceful_parity(twisted_drawing(depth, u_lo, **kw), s)
    p_hi, sq_hi = large_graceful_parity(twisted_drawing(depth, u_hi, **kw), s)
    if p_lo == p_hi:
        raise ValueError("no parity change on this bracket")
    for _ in range(steps):
        mid = (u_lo + u_hi) / 2
        p_mid, sq_mid = large_graceful_parity(twisted_drawing(depth, mid, **kw), s)
        if p_mid == p_lo:
            u_lo, sq_lo = mid, sq_mid
        else:
            u_hi, sq_hi = mid, sq_mid
    return (u_lo, u_hi), sq_lo, sq_hi


def nearest_to_cutoff(squares, s):
    """The small square whose side is closest to s (side^2 < s^2 by construction)."""
    return max(squares.values(), key=lambda e: e[0]) if squares else None
