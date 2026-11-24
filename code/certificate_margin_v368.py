"""v368: exact lower and upper bounds for the redraw margin of a v366 certificate.

For the box B = A x B x E of a v366 certificate put

    mu = min over the boundary of B of  max( |F1| / 4 , |rho(D)| / (2 + sqrt 2) ),

where F1 = gamma(b) + s J(gamma(b) - gamma(a)) - gamma(e) and D = gamma(a) + s J(gamma(b) - gamma(a)).
Every registered Jordan redraw within uniform distance delta < mu keeps the degree (Lemma 368.1).

The lower bound is rigorous: each face of B is cut into cells on which every corner is affine,
and on each cell |F1|_2 is bounded below through its value at the centre (Euclidean or max-norm) minus
an l1 bound on its variation, dist(D, curve) >= dist(D(center), curve) - |D - D(center)|_1, with
square roots bounded from below by integer square roots. Cells are refined best-first.
"""
import heapq
from fractions import Fraction as Q
from itertools import count
from math import isqrt

from sweep_search_v366 import _J, _add, _mul, _sub, locus

SQRT2_LO = Q(141421, 100000)
SQRT2_HI = Q(141422, 100000)
TWO_PLUS_SQRT2_HI = 2 + SQRT2_HI     # dividing by this bounds |rho|/(2+sqrt2) from below
TWO_PLUS_SQRT2_LO = 2 + SQRT2_LO


def sqrt_lo(q, scale=10**9):
    q = Q(q)
    if q <= 0:
        return Q(0)
    return Q(isqrt(q.numerator * scale * scale // q.denominator), scale)


def sqrt_hi(q, scale=10**9):
    q = Q(q)
    if q <= 0:
        return Q(0)
    r = sqrt_lo(q, scale)
    return r if r * r == q else r + Q(1, scale)


def dist_sq(poly, p):
    return poly.clearance_sq(p)


def _affine_point(poly, role_pieces, which):
    """Point of one role as (const, coef_u, coef_w), where `which` is ('fixed', point) or ('u'|'w', edge)."""
    kind, data = which
    zero = (Q(0), Q(0))
    if kind == 'fixed':
        return (data, zero, zero)
    a, b = poly.edges[data]
    d = _sub(b, a)
    return (a, d, zero) if kind == 'u' else (a, zero, d)


def _comb(terms):
    """Linear combination of affine points: terms = [(scalar_or_matrix_fn, affine)]."""
    out = [(Q(0), Q(0))] * 3
    for f, aff in terms:
        for k in range(3):
            out[k] = _add(out[k], f(aff[k]))
    return tuple(out)


def face_cells(poly, arcs, orient):
    """Cells covering the boundary of the box. Each cell: (F1 affine, D affine, (u0,u1), (w0,w1))."""
    pieces = [poly.pieces(*arc) for arc in arcs]
    ident = lambda v: v
    neg = lambda v: _mul(v, -1)
    J = lambda v: _J(v, orient)
    negJ = lambda v: _mul(_J(v, orient), -1)
    cells = []
    for fixed in range(3):
        for end in (0, 1):
            i, u0, u1 = pieces[fixed][0] if end == 0 else pieces[fixed][-1]
            a, b = poly.edges[i]
            pt = _add(a, _mul(_sub(b, a), u0 if end == 0 else u1))
            free = [r for r in range(3) if r != fixed]
            for (eu, su0, su1) in pieces[free[0]]:
                for (ew, sw0, sw1) in pieces[free[1]]:
                    which = {fixed: ('fixed', pt), free[0]: ('u', eu), free[1]: ('w', ew)}
                    Pa, Pb, Pe = (_affine_point(poly, None, which[r]) for r in range(3))
                    F1 = _comb([(ident, Pb), (J, Pb), (negJ, Pa), (neg, Pe)])
                    D = _comb([(ident, Pa), (J, Pb), (negJ, Pa)])
                    cells.append((F1, D, (su0, su1), (sw0, sw1)))
    return cells


def _eval(aff, u, w):
    return _add(_add(aff[0], _mul(aff[1], u)), _mul(aff[2], w))


def _cell_lower(poly, F1, D, ur, wr):
    um, wm = (ur[0] + ur[1]) / 2, (wr[0] + wr[1]) / 2
    ru, rw = (ur[1] - ur[0]) / 2, (wr[1] - wr[0]) / 2
    f = _eval(F1, um, wm)
    spreadF = (abs(F1[1][0]) + abs(F1[1][1])) * ru + (abs(F1[2][0]) + abs(F1[2][1])) * rw
    lbF = max(abs(f[0]) - abs(F1[1][0]) * ru - abs(F1[2][0]) * rw,
              abs(f[1]) - abs(F1[1][1]) * ru - abs(F1[2][1]) * rw,
              sqrt_lo(f[0] ** 2 + f[1] ** 2) - spreadF, Q(0))
    d = _eval(D, um, wm)
    spread = (abs(D[1][0]) + abs(D[1][1])) * ru + (abs(D[2][0]) + abs(D[2][1])) * rw
    lbR = max(sqrt_lo(dist_sq(poly, d)) - spread, Q(0))
    return max(lbF / 4, lbR / TWO_PLUS_SQRT2_HI)


def margin_lower(poly, arcs, orient, iterations=4000):
    """Rigorous lower bound for mu (0 means nothing could be certified)."""
    heap = []
    tick = count()
    for F1, D, ur, wr in face_cells(poly, arcs, orient):
        heapq.heappush(heap, (_cell_lower(poly, F1, D, ur, wr), next(tick), F1, D, ur, wr))
    for _ in range(iterations):
        lb, _, F1, D, ur, wr = heapq.heappop(heap)
        su = (abs(F1[1][0]) + abs(F1[1][1]) + abs(D[1][0]) + abs(D[1][1])) * (ur[1] - ur[0])
        sw = (abs(F1[2][0]) + abs(F1[2][1]) + abs(D[2][0]) + abs(D[2][1])) * (wr[1] - wr[0])
        if su >= sw:
            m = (ur[0] + ur[1]) / 2
            parts = [((ur[0], m), wr), ((m, ur[1]), wr)]
        else:
            m = (wr[0] + wr[1]) / 2
            parts = [(ur, (wr[0], m)), (ur, (m, wr[1]))]
        for pu, pw in parts:
            heapq.heappush(heap, (_cell_lower(poly, F1, D, pu, pw), next(tick), F1, D, pu, pw))
    return heap[0][0]


def margin_upper(poly, arcs, orient):
    """Upper bound for mu from the outer locus endpoints, where F1 = 0 (None if there are none)."""
    best = None
    for e0, e1 in locus(poly, *arcs, orient):
        for P, Qp, R, D, outer in (e0, e1):
            if outer:
                v = sqrt_hi(dist_sq(poly, D)) / TWO_PLUS_SQRT2_LO
                best = v if best is None else min(best, v)
    return best


def registered_distance_sq(p1, p2):
    """Exact squared uniform distance between two registered polygons at equal times."""
    times = sorted(set(p1.times) | set(p2.times))
    best = Q(0)
    for t in times[:-1]:
        x, y = p1.point(t), p2.point(t)
        best = max(best, (x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2)
    return best


def _value_hi(poly, F1, D, u, w):
    """Upper bound for max(|F1|_2/4, dist(D)/(2+sqrt2)) at one parameter point."""
    f = _eval(F1, u, w)
    d = _eval(D, u, w)
    return max(sqrt_hi(f[0] ** 2 + f[1] ** 2) / 4, sqrt_hi(dist_sq(poly, d)) / TWO_PLUS_SQRT2_LO)


def margin_bracket(poly, arcs, orient, iterations=4000):
    """(lower, upper, where): rigorous lower bound, and an upper bound from the value at the centre
    of the cell holding the lower bound and at the outer locus endpoints."""
    heap = []
    tick = count()
    for F1, D, ur, wr in face_cells(poly, arcs, orient):
        heapq.heappush(heap, (_cell_lower(poly, F1, D, ur, wr), next(tick), F1, D, ur, wr))
    upper = margin_upper(poly, arcs, orient)
    where = None
    for _ in range(iterations):
        lb, _, F1, D, ur, wr = heapq.heappop(heap)
        um, wm = (ur[0] + ur[1]) / 2, (wr[0] + wr[1]) / 2
        v = _value_hi(poly, F1, D, um, wm)
        if upper is None or v < upper:
            upper, where = v, (_eval(D, um, wm), _eval(F1, um, wm))
        su = (abs(F1[1][0]) + abs(F1[1][1]) + abs(D[1][0]) + abs(D[1][1])) * (ur[1] - ur[0])
        sw = (abs(F1[2][0]) + abs(F1[2][1]) + abs(D[2][0]) + abs(D[2][1])) * (wr[1] - wr[0])
        if su >= sw:
            m = um
            parts = [((ur[0], m), wr), ((m, ur[1]), wr)]
        else:
            m = wm
            parts = [(ur, (wr[0], m)), (ur, (m, wr[1]))]
        for pu, pw in parts:
            heapq.heappush(heap, (_cell_lower(poly, F1, D, pu, pw), next(tick), F1, D, pu, pw))
    return heap[0][0], upper, where
