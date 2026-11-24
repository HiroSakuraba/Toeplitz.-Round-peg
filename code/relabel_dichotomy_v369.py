"""v369: exact checks for the relabeling dichotomy.

Rows are lifted times t0 <= t1 <= t2 <= t3 <= t0 + 1 modulo a common integer shift. The cyclic
relabeling is T(t0, t1, t2, t3) = (t1, t2, t3, t0 + 1). The square residual of four points is
(r1, r2) = (p2 - p1 - J v, p3 - p0 - J v) with v = p1 - p0 and J the counterclockwise quarter turn.
"""
from fractions import Fraction as Q
from itertools import combinations

from sweep_search_v366 import _seg_dist_sq


def J(v, sigma=1):
    return (-sigma * v[1], sigma * v[0])


def residual(points, sigma=1):
    p0, p1, p2, p3 = points
    v = (p1[0] - p0[0], p1[1] - p0[1])
    Jv = J(v, sigma)
    return ((p2[0] - p1[0] - Jv[0], p2[1] - p1[1] - Jv[1]),
            (p3[0] - p0[0] - Jv[0], p3[1] - p0[1] - Jv[1]))


def L(r, sigma=1):
    """The linear map with residual(shifted points) = L(residual(points)) (Lemma 369.3)."""
    r1, r2 = r
    Jr1 = J(r1, sigma)
    return ((r2[0] - r1[0] - Jr1[0], r2[1] - r1[1] - Jr1[1]), (-Jr1[0], -Jr1[1]))


def shift(points):
    return tuple(points[1:]) + (points[0],)


def T(row):
    t0, t1, t2, t3 = row
    return (t1, t2, t3, t0 + 1)


def gaps(row):
    t0, t1, t2, t3 = row
    return (t1 - t0, t2 - t1, t3 - t2, t0 + 1 - t3)


def is_row(row):
    return all(g >= 0 for g in gaps(row))


def same_row(x, y):
    """Equality modulo a common integer shift."""
    d = y[0] - x[0]
    return d == int(d) and all(b - a == d for a, b in zip(x, y))


def big_gap_label(row):
    """Index of the unique gap exceeding 1/2, or None."""
    big = [i for i, g in enumerate(gaps(row)) if g > Q(1, 2)]
    return big[0] if big else None


def circular(s, t):
    d = (t - s) % 1
    return min(d, 1 - d)


def separated_pair(row):
    """For a row with every gap at most 1/2, a pair of times at circular distance >= 1/3."""
    if any(g > Q(1, 2) for g in gaps(row)):
        raise ValueError("a gap exceeds one half")
    best = max(combinations(range(4), 2), key=lambda ij: circular(row[ij[0]], row[ij[1]]))
    return best, circular(row[best[0]], row[best[1]])


def small_square_scale_sq(poly):
    """tau^2, the smallest squared distance between non-adjacent edges (polygons with n >= 4).
    A polygon has no inscribed square of side s with 2 s^2 < tau^2 (Lemma 369.6)."""
    n = len(poly.edges)
    if n < 4:
        raise ValueError("need at least four edges")
    best = None
    for i, j in combinations(range(n), 2):
        if (j - i) % n in (1, n - 1):
            continue
        (a, b), (c, d) = poly.edges[i], poly.edges[j]
        dd = _seg_dist_sq(a, b, c, d)
        best = dd if best is None else min(best, dd)
    return best
