"""v366: candidate sweeps read off an unmarked polygon.

Slide P along one stretch and Q along another, complete PQ to a square P, Q, R, D, and keep the
configurations with R on a third stretch of the boundary. These configurations form a piecewise
linear locus L in the box B of the three traversal parameters. The fourth corner D is affine
along each piece of L.

Certificate. Let F = (R - Gamma(tau), rho(D)) on B, with rho the signed distance to the whole
boundary. If no endpoint of L on the outer faces of B has D on the boundary, then F has no zero on
the boundary of B, and its mod-two degree is the parity of the outer endpoints of L with D strictly
inside. An odd degree gives a square with P, Q, R in the three stretches on every sufficiently close
registered redraw. All tests below are exact rational arithmetic.
"""
from fractions import Fraction as Q
from itertools import product

from polygon_orbit_census_v223 import pt, validate_jordan_polygon


def _J(v, orient):
    return (-orient * v[1], orient * v[0])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _mul(a, k):
    return (a[0] * k, a[1] * k)


def _det(a, b):
    return a[0] * b[1] - a[1] * b[0]


class RegisteredPolygon:
    """A rational polygon with registration times t_0 = 0 < ... < t_n = 1 at its vertices."""

    def __init__(self, vertices, times=None):
        poly = validate_jordan_polygon([pt(*map(Q, v)) for v in vertices])
        self.V = [(p.x, p.y) for p in poly]
        n = len(self.V)
        self.times = [Q(i, n) for i in range(n + 1)] if times is None else [Q(t) for t in times]
        if len(self.times) != n + 1 or self.times[0] != 0 or self.times[-1] != 1 or any(
                a >= b for a, b in zip(self.times, self.times[1:])):
            raise ValueError("expected increasing registration times from 0 to 1")
        self.edges = [(self.V[i], self.V[(i + 1) % n]) for i in range(n)]
        self._pt_edges = list(zip(poly, poly[1:] + poly[:1]))

    def point(self, t):
        t = Q(t) % 1
        for i in range(len(self.edges)):
            a, b = self.times[i], self.times[i + 1]
            if a <= t <= b:
                u = (t - a) / (b - a)
                p, q = self.edges[i]
                return _add(p, _mul(_sub(q, p), u))
        raise AssertionError

    def pieces(self, start, length):
        """Edge pieces (edge, u0, u1) covering the closed arc [start, start+length], in order."""
        key = (Q(start) % 1, Q(length))
        cache = self.__dict__.setdefault('_pieces', {})
        if key not in cache:
            cache[key] = self._pieces_uncached(*key)
        return cache[key]

    def _pieces_uncached(self, start, length):
        start, length = Q(start) % 1, Q(length)
        if not 0 < length < 1:
            raise ValueError("arc length must be in (0,1)")
        out = []
        end = start + length
        for shift in (0, 1):
            for i in range(len(self.edges)):
                a, b = self.times[i] + shift, self.times[i + 1] + shift
                lo, hi = max(a, start), min(b, end)
                if lo < hi:
                    out.append((i, (lo - a) / (b - a), (hi - a) / (b - a)))
        def key(piece):
            i, u0, _ = piece
            t = self.times[i] + u0 * (self.times[i + 1] - self.times[i])
            return t if t >= start else t + 1
        out.sort(key=key)
        return out

    def side_of(self, point):
        """True inside, False outside, None on the boundary (exact ray parity)."""
        x, y = point
        inside = False
        for (ax, ay), (bx, by) in self.edges:
            dx, dy = bx - ax, by - ay
            cross = (x - ax) * dy - (y - ay) * dx
            if cross == 0 and min(ax, bx) <= x <= max(ax, bx) and min(ay, by) <= y <= max(ay, by):
                return None
            if (ay > y) != (by > y) and x < ax + (y - ay) * dx / dy:
                inside = not inside
        return inside

    def clearance_sq(self, point):
        x, y = point
        best = None
        for (ax, ay), (bx, by) in self.edges:
            dx, dy = bx - ax, by - ay
            t = max(Q(0), min(Q(1), ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy)))
            d = (x - ax - t * dx) ** 2 + (y - ay - t * dy) ** 2
            best = d if best is None else min(best, d)
        return best


def _segment_hits(poly, p0, p1):
    """Parameters lam in [0,1] where p0 + lam (p1 - p0) meets the boundary; None if collinear overlap."""
    d = _sub(p1, p0)
    out = set()
    for a, b in poly.edges:
        e = _sub(b, a)
        den = _det(d, e)
        w = _sub(a, p0)
        if den == 0:
            if d != (0, 0) and _det(w, d) == 0:
                # collinear: overlap only matters if the two segments share a point
                dd = d[0] * d[0] + d[1] * d[1]
                u0 = (w[0] * d[0] + w[1] * d[1]) / dd
                u1 = ((b[0] - p0[0]) * d[0] + (b[1] - p0[1]) * d[1]) / dd
                if max(min(u0, u1), 0) <= min(max(u0, u1), 1):
                    return None  # overlap along a line: not a transverse crossing
            continue
        lam = _det(w, e) / den
        mu = _det(w, d) / den
        if 0 <= lam <= 1 and 0 <= mu <= 1:
            out.add(lam)
    return sorted(out)


def _chart_line(poly, ia, ib, ie, orient):
    """Solution line x + lam n of the 2x3 system in local edge parameters, cached per polygon."""
    cache = poly.__dict__.setdefault('_lines', {})
    key = (ia, ib, ie, orient)
    if key in cache:
        return cache[key]
    A0, A1 = poly.edges[ia]
    B0, B1 = poly.edges[ib]
    E0, E1 = poly.edges[ie]
    dA, dB, dE = _sub(A1, A0), _sub(B1, B0), _sub(E1, E0)
    m1 = _mul(_J(dA, orient), -1)
    m2 = _add(dB, _J(dB, orient))
    m3 = _mul(dE, -1)
    c = _sub(_add(B0, _J(_sub(B0, A0), orient)), E0)
    n = (_det(m2, m3), _det(m3, m1), _det(m1, m2))
    if n == (0, 0, 0):
        cache[key] = None
        return None
    k = max(range(3), key=lambda i: abs(n[i]))
    cols = [m1, m2, m3]
    i, j = [t for t in range(3) if t != k]
    det = _det(cols[i], cols[j])
    rhs = _mul(c, -1)
    x = [Q(0)] * 3
    x[i] = _det(rhs, cols[j]) / det
    x[j] = _det(cols[i], rhs) / det
    cache[key] = (x, n, (A0, dA, B0, dB, E0, dE))
    return cache[key]


def locus(poly, arcA, arcB, arcE, orient):
    """Pieces of L for stretches given as (start, length). Returns a list of segments
    (end0, end1) where each end is (P, Q, R, D, outer) and `outer` says whether the end lies on an
    outer face of B. Raises ArithmeticError on a degenerate (rank-deficient) chart that meets B."""
    pa, pb, pe = (poly.pieces(*arc) for arc in (arcA, arcB, arcE))
    segments = []
    for (ia, a0, a1), (ib, b0, b1), (ie, e0, e1) in product(pa, pb, pe):
        line = _chart_line(poly, ia, ib, ie, orient)
        if line is None:
            raise ArithmeticError("rank-deficient chart")
        x, n, (A0, dA, B0, dB, E0, dE) = line
        lo = [a0, b0, e0]
        hi = [a1, b1, e1]
        lam_lo, lam_hi = None, None
        ok = True
        for t in range(3):
            if n[t] == 0:
                if not lo[t] <= x[t] <= hi[t]:
                    ok = False
                    break
                continue
            u, v = (lo[t] - x[t]) / n[t], (hi[t] - x[t]) / n[t]
            if u > v:
                u, v = v, u
            lam_lo = u if lam_lo is None else max(lam_lo, u)
            lam_hi = v if lam_hi is None else min(lam_hi, v)
        if not ok or lam_lo is None or lam_lo > lam_hi:
            continue
        ends = []
        for lam in (lam_lo, lam_hi):
            s, r, tau = (x[t] + lam * n[t] for t in range(3))
            P = _add(A0, _mul(dA, s))
            Qp = _add(B0, _mul(dB, r))
            R = _add(E0, _mul(dE, tau))
            D = _add(P, _J(_sub(Qp, P), orient))
            faces = set()
            for name, val, idx, piece_list in (('A', s, ia, pa), ('B', r, ib, pb), ('E', tau, ie, pe)):
                if (idx, val) == (piece_list[0][0], piece_list[0][1]):
                    faces.add(name + '0')
                if (idx, val) == (piece_list[-1][0], piece_list[-1][2]):
                    faces.add(name + '1')
            outer = frozenset(faces)
            ends.append((P, Qp, R, D, outer))
        segments.append(tuple(ends))
    return segments


def certificate(poly, arcA, arcB, arcE, orient):
    """Evaluate one stretch triple. Returns a dict, or None if the boundary of B carries a zero or a
    chart is degenerate. 'degree' is the parity of inside outer endpoints; 'squares' lists exact
    squares found as boundary crossings of D along L (an independent count)."""
    try:
        segs = locus(poly, arcA, arcB, arcE, orient)
    except ArithmeticError:
        return None
    inside_ends = 0
    clearance = None
    outer_count = 0
    squares = []
    for e0, e1 in segs:
        if e0[:3] != e1[:3] and e0[4] & e1[4]:
            return None       # a piece of L runs inside an outer face: not transverse to the boundary
        for P, Qp, R, D, outer in (e0, e1):
            if not outer:
                continue
            side = poly.side_of(D)
            if side is None:
                return None   # a zero on the boundary of B
            c = poly.clearance_sq(D)
            clearance = c if clearance is None else min(clearance, c)
            if e0[:3] == e1[:3]:
                continue      # a chart met at a single point contributes twice: nothing mod two
            # every piece end on an outer face counts once; where L touches an outer face and
            # returns, the two piece ends cancel
            inside_ends += side
            outer_count += 1
        if e0[:3] == e1[:3]:
            continue
        hits = _segment_hits(poly, e0[3], e1[3])
        if hits is None:
            return None
        for lam in hits:
            if 0 < lam < 1 or (lam == 0 and not e0[4]) or (lam == 1 and not e1[4]):
                P = _add(e0[0], _mul(_sub(e1[0], e0[0]), lam))
                Qp = _add(e0[1], _mul(_sub(e1[1], e0[1]), lam))
                R = _add(e0[2], _mul(_sub(e1[2], e0[2]), lam))
                D = _add(e0[3], _mul(_sub(e1[3], e0[3]), lam))
                squares.append((P, Qp, R, D))
    uniq = []
    for sq in squares:
        if sq not in uniq:
            uniq.append(sq)
    if len(uniq) % 2 != inside_ends % 2:
        return None           # the independent crossing count disagrees: treat as inconclusive
    return {'degree': inside_ends % 2, 'outer_ends': outer_count, 'clearance_sq': clearance,
            'squares': uniq}


def _seg_dist_sq(p, q, a, b):
    def pd(x, a, b):
        d = _sub(b, a)
        t = max(Q(0), min(Q(1), ((x[0] - a[0]) * d[0] + (x[1] - a[1]) * d[1]) / (d[0] ** 2 + d[1] ** 2)))
        return (x[0] - a[0] - t * d[0]) ** 2 + (x[1] - a[1] - t * d[1]) ** 2
    d1, d2 = _sub(q, p), _sub(b, a)
    if _det(d1, d2) != 0:
        den = _det(d1, d2)
        w = _sub(a, p)
        u, v = _det(w, d2) / den, _det(w, d1) / den
        if 0 <= u <= 1 and 0 <= v <= 1:
            return Q(0)
    return min(pd(p, a, b), pd(q, a, b), pd(a, p, q), pd(b, p, q))


def arc_distance_sq(poly, arc1, arc2):
    best = None
    for i, u0, u1 in poly.pieces(*arc1):
        p0, p1 = poly.edges[i]
        a = _add(p0, _mul(_sub(p1, p0), u0)); b = _add(p0, _mul(_sub(p1, p0), u1))
        for j, w0, w1 in poly.pieces(*arc2):
            q0, q1 = poly.edges[j]
            c = _add(q0, _mul(_sub(q1, q0), w0)); d = _add(q0, _mul(_sub(q1, q0), w1))
            dd = _seg_dist_sq(a, b, c, d)
            best = dd if best is None else min(best, dd)
    return best


def search(poly, bins=16, offset=Q(0), spans=(1,), gap=1, orients=(1, -1)):
    """Every certificate on stretch triples made of whole bins, with at least `gap` empty bins
    between consecutive stretches in boundary order. Returns a list sorted by the squared side floor
    max(d(A,B), d(B,E), d(A,E)/2), where d is the squared stretch distance."""
    out = []
    for span in spans:
        L = Q(span, bins)
        for a, b, e in product(range(bins), repeat=3):
            arcs = [(offset + Q(x, bins), L) for x in (a, b, e)]
            starts = sorted((a, 'A'), (b, 'B'), (e, 'E')) if False else sorted([a, b, e])
            # pairwise cyclic gaps of at least `gap` bins
            gaps_ok = True
            order = sorted([a, b, e])
            for x, y in zip(order, order[1:] + [order[0] + bins]):
                if y - x < span + gap:
                    gaps_ok = False
            if not gaps_ok:
                continue
            for orient in orients:
                cert = certificate(poly, *arcs, orient)
                if cert and cert['degree'] == 1:
                    dAB = arc_distance_sq(poly, arcs[0], arcs[1])
                    dBE = arc_distance_sq(poly, arcs[1], arcs[2])
                    dAE = arc_distance_sq(poly, arcs[0], arcs[2])
                    floor = max(dAB, dBE, dAE / 2)
                    out.append({'arcs': (a, b, e), 'span': span, 'orient': orient,
                                'side_sq_floor': floor, **cert})
    out.sort(key=lambda c: -c['side_sq_floor'])
    return out


def rational_circle(n, denominator=4096):
    """n points on the unit circle with rational coordinates (half-angle tangent rounded)."""
    from math import pi, tan
    pts = []
    for k in range(n):
        theta = 2 * pi * k / n
        if abs(theta - pi) < 1e-12:
            pts.append((Q(-1), Q(0)))
            continue
        u = Q(round(tan(theta / 2) * denominator), denominator) if abs(theta - pi) > 1e-9 else None
        pts.append(((1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)))
    return pts


def symmetric_circle(m, denominator=4096, bump=Q(0)):
    """4m rational points on the unit circle, exactly invariant under the quarter turn.
    `bump` pushes the first vertex radially outward by the given factor, breaking the symmetry."""
    from math import pi, tan
    quad = []
    for k in range(m):
        u = Q(round(tan(pi / 4 * k / m) * denominator), denominator)
        quad.append(((1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)))
    pts = []
    for r in range(4):
        for x, y in quad:
            for _ in range(r):
                x, y = -y, x
            pts.append((x, y))
    if bump:
        x, y = pts[0]
        pts[0] = (x * (1 + bump), y * (1 + bump))
    return pts
