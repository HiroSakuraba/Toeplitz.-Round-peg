"""v358: rectangles with any finite number of rectangular roof notches.

Exact rational tools:
  * an exact decision of sqrt(2) r <= side, where r is the inradius (ties included);
  * closed-form robust square components (bottom-supported and wall-anchored squares, and
    full-height or floor-level slides with an interleaving parity rule);
  * a three-level certifier: general theorems, closed-form components, and an exhaustive
    transverse-square search followed by v341 component tracing.

Standard library only (an optional numpy search accelerator is used when available).
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
from math import lcm

from tall_strip_v355 import roof_disk_fits


# ---------------------------------------------------------------------------
# Shapes

@dataclass(frozen=True)
class NotchedRectangle:
    W: Q
    H: Q
    cuts: tuple  # sorted tuples (a, b, k): bite (a, b) x (k, H]

    @staticmethod
    def make(W, H, cuts):
        W, H = Q(W), Q(H)
        cs = tuple((Q(a), Q(b), Q(k)) for a, b, k in cuts)
        previous = Q(0)
        if not cs:
            raise ValueError("expected at least one roof cut")
        for a, b, k in cs:
            if not (previous < a < b < W and 0 < k < H):
                raise ValueError("cuts must be sorted, disjoint, proper and leave positive roof")
            previous = b
        return NotchedRectangle(W, H, cs)

    def polygon(self):
        P = [(Q(0), Q(0)), (self.W, Q(0)), (self.W, self.H)]
        for a, b, k in reversed(self.cuts):
            P += [(b, self.H), (b, k), (a, k), (a, self.H)]
        P.append((Q(0), self.H))
        return P

    def mirror(self):
        return NotchedRectangle(self.W, self.H, tuple((self.W - b, self.W - a, k) for a, b, k in reversed(self.cuts)))

    def towers(self):
        xs = [Q(0)] + [c for a, b, k in self.cuts for c in (a, b)] + [self.W]
        return [(xs[2 * t], xs[2 * t + 1]) for t in range(len(self.cuts) + 1)]

    def cut_walls(self):
        """(abscissa, floor height) of every notch wall; the wall spans (floor, H)."""
        return [(c, k) for a, b, k in self.cuts for c in (a, b)]

    def in_open_tower(self, x):
        return any(s < x < e for s, e in self.towers())

    def as_args(self):
        return self.W, self.H, [tuple(c) for c in self.cuts]


# ---------------------------------------------------------------------------
# Exact inradius comparison

def _sign_lin(a, b, q):
    """Sign of a + b sqrt(q) for rationals a, b and q >= 0."""
    if b == 0 or q == 0:
        return (a > 0) - (a < 0)
    sb = 1 if b > 0 else -1
    if a == 0:
        return sb
    sa = 1 if a > 0 else -1
    if sa == sb:
        return sa
    d = a * a - b * b * q
    return sa if d > 0 else (-sa if d < 0 else 0)


class QuadraticNumber:
    """(p + s sqrt(q)) / m with rational p, s, q >= 0 and m > 0."""

    def __init__(self, p, s, q, m):
        p, s, q, m = Q(p), Q(s), Q(q), Q(m)
        if m < 0:
            p, s, m = -p, -s, -m
        self.p, self.s, self.q, self.m = p, s, q, m

    def compare(self, x):
        return _sign_lin(self.p - self.m * Q(x), self.s, self.q)

    def compare_square(self, R):
        return _sign_lin(self.p ** 2 + self.s ** 2 * self.q - Q(R) * self.m ** 2, 2 * self.p * self.s, self.q)


def _roots(A, B, C):
    if A == 0:
        return [QuadraticNumber(-C, 0, 0, B)] if B else []
    disc = B * B - 4 * A * C
    return [] if disc < 0 else [QuadraticNumber(-B, sg, disc, 2 * A) for sg in (1, -1)]


def inradius_candidates(shape):
    """A finite list containing the inradius.

    The grounded disk centre must lie at least ell(s) from each active obstacle, with ell = s at a wall
    or a low notch and ell = sqrt(2ks - k^2) at a notch corner. The inradius is H/2, half a floor
    height (a notch becoming active over the centre), or a tight pair ell_1 + ell_2 = D. Each pair gives
    a linear or quadratic equation.
    """
    W, H = shape.W, shape.H
    out = [QuadraticNumber(H, 0, 0, 2), QuadraticNumber(W, 0, 0, 2)]
    out += [QuadraticNumber(k, 0, 0, 2) for a, b, k in shape.cuts]
    lowers = [(Q(0), None)] + [(b, k) for a, b, k in shape.cuts]
    uppers = [(W, None)] + [(a, k) for a, b, k in shape.cuts]
    for lc, lk in lowers:
        for uc, uk in uppers:
            D = uc - lc
            for t1 in ['s'] + (['r'] if lk is not None else []):
                for t2 in ['s'] + (['r'] if uk is not None else []):
                    if t1 == t2 == 's':
                        out.append(QuadraticNumber(D, 0, 0, 2))
                    elif t1 == t2 == 'r':
                        a1, b1, a2, b2 = 2 * lk, lk * lk, 2 * uk, uk * uk
                        g, d = a1 + a2, b1 + b2
                        out += _roots(g * g - 4 * a1 * a2, -2 * g * (D * D + d) + 4 * (a1 * b2 + a2 * b1),
                                      (D * D + d) ** 2 - 4 * b1 * b2)
                    else:
                        k = lk if t1 == 'r' else uk
                        out += _roots(Q(1), -2 * (D + k), D * D + k * k)
    return [c for c in out if c.compare(0) > 0]


def radius_bracket(shape, iterations=60):
    lo, hi = Q(0), min(shape.W, shape.H) / 2
    if roof_disk_fits(*shape.as_args(), hi):
        return hi, hi
    for _ in range(iterations):
        mid = (lo + hi) / 2
        if roof_disk_fits(*shape.as_args(), mid):
            lo = mid
        else:
            hi = mid
    return lo, hi


def reaches_disk_target(shape, side_sq, max_iter=400):
    """Exactly decide side_sq >= 2 r^2. Returns True or False (None only if undecided)."""
    F = Q(side_sq)
    args = shape.as_args()
    lo, hi = Q(0), min(shape.W, shape.H) / 2
    if roof_disk_fits(*args, hi):
        return 2 * hi * hi <= F
    cands = inradius_candidates(shape)
    for _ in range(max_iter):
        if 2 * hi * hi <= F:
            return True
        if 2 * lo * lo > F:
            return False
        inside = [c for c in cands if c.compare(lo) >= 0 and c.compare(hi) <= 0]
        if inside:
            le = [c.compare_square(F / 2) <= 0 for c in inside]
            if all(le):
                return True
            if not any(le):
                return False
        mid = (lo + hi) / 2
        if roof_disk_fits(*args, mid):
            lo = mid
        else:
            hi = mid
    return None


# ---------------------------------------------------------------------------
# Level 0: general theorems

def largest_axis_square(shape):
    xs = sorted({Q(0), shape.W} | {c for a, b, k in shape.cuts for c in (a, b)})
    best = Q(0)
    for i in range(len(xs)):
        for j in range(i + 1, len(xs)):
            x0, x1 = xs[i], xs[j]
            m = min([shape.H] + [k for a, b, k in shape.cuts if a < x1 and b > x0])
            best = max(best, min(x1 - x0, m))
    return best


def covering_theorem(shape):
    """Name of a previously proved theorem that gives the robust disk bound, or None."""
    W, H = shape.W, shape.H
    if H > W:
        return 'tall-walls-v340'
    if H == W:
        return 'square-box-v353'
    T = shape.towers()
    if T[0][1] - T[0][0] > H or T[-1][1] - T[-1][0] > H:
        return 'end-tower-v351'
    if T[0][1] - T[0][0] == H or T[-1][1] - T[-1][0] == H:
        return 'end-plateau-v351'
    i = largest_axis_square(shape)
    low = [(a, b) for a, b, k in shape.cuts if k < i]
    for (l, r) in zip([Q(0)] + [b for a, b in low], [a for a, b in low] + [W]):
        if r - l >= i and l > 0 and r < W:
            return 'interior-slot-v338'
    return None


# ---------------------------------------------------------------------------
# Level 1: closed-form components

def square_from(p0, e):
    (x, y), (a, b) = p0, e
    return ((x, y), (x + a, y + b), (x + a - b, y + b + a), (x - b, y + a))


def _interleave(points, A, B):
    labels = ''.join(lab for _, lab in sorted([(p, 'A') for p in A] + [(p, 'B') for p in B]))
    return labels in ('ABAB', 'BABA')


def slide_parity(point, turn, Al, Ar, Bl, Br):
    """Mod-2 degree of a maximal slide of axis squares on one horizontal level.

    A is the level piece under the top-left vertex, B the sheared piece under the top-right vertex.
    Al, Ar, Bl, Br say whether A or B leaves the level at the left or right end of the overlap.
    The boundary leaves a roof piece downward and a floor piece upward; A leaves vertically and B
    along the image of a vertical wall under (x, y) -> (x - y, y).
    """
    if turn == 'down':
        pos = ({'w': 0, 'Lv': 1, 'Rv': 1, 'Ld': 2, 'Rd': 2, 'e': 3} if point else
               {'w': 0, 'Lv': 1, 'Ld': 2, 'Rv': 3, 'Rd': 4, 'e': 5})
    else:
        pos = ({'e': 0, 'Rv': 1, 'Lv': 1, 'Rd': 2, 'Ld': 2, 'w': 3} if point else
               {'e': 0, 'Rv': 1, 'Rd': 2, 'Lv': 3, 'Ld': 4, 'w': 5})
    A = [pos['Lv' if Al else 'w'], pos['Rv' if Ar else 'e']]
    B = [pos['Ld' if Bl else 'w'], pos['Rd' if Br else 'e']]
    return _interleave(pos, A, B)


@dataclass(frozen=True)
class Component:
    kind: str
    side_sq: Q
    seed: tuple          # one square of the component
    odd: bool            # local mod-2 degree
    detail: tuple = ()


def level_slides(shape):
    W, H = shape.W, shape.H
    wallx = {c for c, k in shape.cut_walls()}
    levels = {H: ('down', shape.towers())}
    for a, b, k in shape.cuts:
        levels.setdefault(k, ('up', []))[1].append((a, b))
    for c, (turn, pieces) in levels.items():
        for (s1, e1), (s2, e2) in product(pieces, repeat=2):
            x0, x1 = max(s1, s2 - c), min(e1, e2 - c)
            if x0 > x1:
                continue
            if (x0 == 0 and c in wallx) or (x1 + c == W and W - c in wallx):
                continue  # a box-corner end that can turn along a wall: outside this rule
            odd = slide_parity(x0 == x1, turn, x0 == s1, x1 == e1, x0 + c == s2, x1 + c == e2)
            name = 'roof-slide' if turn == 'down' else 'floor-slide'
            yield Component(name, c * c, square_from((x0, Q(0)), (c, Q(0))), odd, (c, x0, x1))


def bottom_supported(shape):
    """Squares with two vertices in the open bottom edge: wall pairs, wall-floor squares and slides."""
    H = shape.H
    walls = shape.cut_walls()
    for (w1, k1), (w2, k2) in product(walls, repeat=2):
        s = w2 - w1
        if s > 0 and max(k1, k2) < s < H:
            yield Component('wall-pair', s * s, square_from((w1, Q(0)), (s, Q(0))), True, (w1, w2))
    for (w, kw), (a, b, k) in product(walls, shape.cuts):
        if kw < k and a < w + k < b:
            yield Component('wall-floor', k * k, square_from((w, Q(0)), (k, Q(0))), True, (w, k))
        if kw < k and a < w - k < b:
            yield Component('floor-wall', k * k, square_from((w - k, Q(0)), (k, Q(0))), True, (w, k))
    yield from level_slides(shape)


def _right_facing_anchored(shape):
    """Squares (c+u,0), (c+t,u), (c+v,t), (c,v), t = u+v, anchored on the bottom and on a wall x = c
    whose region lies to its right (the outer left wall or the right wall of a notch)."""
    H, W = shape.H, shape.W
    anchors = [(Q(0), Q(0))] + [(b, k) for a, b, k in shape.cuts]
    verticals = [(c, k) for c, k in shape.cut_walls()]
    p1_verticals = verticals + [(W, Q(0))]  # p1 may reach the outer right wall; p2 never can (v < t)
    for c, kc in anchors:
        def ok(u, v, t):
            return u > 0 and kc < v < H and 0 < c + u < W
        for (a, b, k) in shape.cuts:            # p1 on a floor: u = k
            u = k
            for w, kw in verticals:              # p2 on a wall: c + v = w
                v = w - c
                t = u + v
                if ok(u, v, t) and a < c + t < b and kw < t < H:
                    yield u, v, c, 'floor-wall'
            for (a2, b2, k2) in shape.cuts:      # p2 on a floor: t = k2
                t = k2
                v = t - u
                if ok(u, v, t) and a < c + t < b and a2 < c + v < b2:
                    yield u, v, c, 'floor-floor'
            t = H                                # p2 on the roof
            v = t - u
            if ok(u, v, t) and a < c + t < b and shape.in_open_tower(c + v):
                yield u, v, c, 'floor-roof'
        for (w, kw), (w2, kw2) in product(p1_verticals, verticals):   # p1, p2 on walls
            t, v = w - c, w2 - c
            u = t - v
            if ok(u, v, t) and kw < u < H and kw2 < t < H:
                yield u, v, c, 'wall-wall'


def wall_anchored(shape):
    for u, v, c, name in _right_facing_anchored(shape):
        sq = ((c + u, Q(0)), (c + u + v, u), (c + v, u + v), (c, v))
        yield Component('anchored-' + name, u * u + v * v, sq, True, (c, u, v))
    W = shape.W
    for u, v, c, name in _right_facing_anchored(shape.mirror()):
        sq = ((W - c - u, Q(0)), (W - c, v), (W - c - v, u + v), (W - c - u - v, u))
        yield Component('anchored-' + name + '-mirror', u * u + v * v, sq, True, (W - c, u, v))


def closed_form_components(shape):
    yield from bottom_supported(shape)
    yield from wall_anchored(shape)


# ---------------------------------------------------------------------------
# Level 2: exhaustive transverse search, then traced components

_ROWS = {
    (0, 'x'): (1, 0, 0, 0), (0, 'y'): (0, 1, 0, 0), (1, 'x'): (1, 0, 1, 0), (1, 'y'): (0, 1, 0, 1),
    (2, 'x'): (1, 0, 1, -1), (2, 'y'): (0, 1, 1, 1), (3, 'x'): (1, 0, 0, -1), (3, 'y'): (0, 1, 1, 0),
}


def _det(m):
    if len(m) == 1:
        return m[0][0]
    return sum((-1) ** j * m[0][j] * _det([r[:j] + r[j + 1:] for r in m[1:]]) for j in range(len(m)) if m[0][j])


def _patterns():
    out = []
    for pat in product('xy', repeat=4):
        M = [list(_ROWS[(i, pat[i])]) for i in range(4)]
        d = _det(M)
        if not d:
            continue
        A = [[(-1) ** (i + j) * _det([r[:i] + r[i + 1:] for kk, r in enumerate(M) if kk != j])
              for j in range(4)] for i in range(4)]
        free = [_ROWS[(i, 'y' if pat[i] == 'x' else 'x')] for i in range(4)]
        F = [[sum(free[i][r] * A[r][j] for r in range(4)) for j in range(4)] for i in range(4)]
        out.append((pat, d, F, A))
    return out


_PATTERNS = _patterns()


def _scaled_edges(polygon):
    den = 1
    for x, y in polygon:
        den = lcm(den, x.denominator, y.denominator)
    pts = [(int(x * den), int(y * den)) for x, y in polygon]
    E = {'x': [], 'y': []}
    for i, ((x0, y0), (x1, y1)) in enumerate(zip(pts, pts[1:] + pts[:1])):
        if y0 == y1:
            E['y'].append((y0, min(x0, x1), max(x0, x1)))
        else:
            E['x'].append((x0, min(y0, y1), max(y0, y1)))
    return den, E


def _candidate_tuples(E, pat, d, F, A, den, min_side_sq, strict):
    """Edge tuples that may carry a square of this contact pattern (a superset; rechecked exactly)."""
    lists = [E[c] for c in pat]
    try:
        import numpy as np
    except ImportError:
        yield from product(*lists)
        return
    arrs = [np.array(lst, dtype=np.int64) for lst in lists]
    shape = [len(a) for a in arrs]
    C = [arrs[i][:, 0].reshape([-1 if j == i else 1 for j in range(4)]) for i in range(4)]
    a = sum(A[2][j] * C[j] for j in range(4))
    b = sum(A[3][j] * C[j] for j in range(4))
    cut = float(min_side_sq) * d * d * den * den
    ok = np.broadcast_to(((a * a + b * b).astype(float) >= cut * (1 - 1e-9) - 1) & ((a != 0) | (b != 0)), shape).copy()
    sg = 1 if d > 0 else -1
    for i in range(4):
        f = sum(F[i][j] * C[j] for j in range(4)) * sg
        lo = arrs[i][:, 1].reshape([-1 if j == i else 1 for j in range(4)]) * d * sg
        hi = arrs[i][:, 2].reshape([-1 if j == i else 1 for j in range(4)]) * d * sg
        ok &= (np.minimum(lo, hi) <= f) & (f <= np.maximum(lo, hi))
    for ids in np.argwhere(ok):
        yield tuple(lists[i][int(ids[i])] for i in range(4))


def contact_squares(polygon, min_side_sq=Q(0), strict=True):
    """Every square whose vertices lie on edges of an axis-parallel polygon with a nonsingular contact
    system (each vertex fixes one coordinate). strict=True requires open edges: such a square is a
    regular zero of the four-point residual, hence has local degree +-1. strict=False allows corners and
    is used only to seed component traces. numpy, when present, only prunes the enumeration; every
    reported square is recomputed and checked in exact arithmetic."""
    den, E = _scaled_edges(polygon)
    out = {}
    for pat, d, F, A in _PATTERNS:
        sg = 1 if d > 0 else -1
        for es in _candidate_tuples(E, pat, d, F, A, den, min_side_sq, strict):
            c = (es[0][0], es[1][0], es[2][0], es[3][0])
            a = sum(A[2][j] * c[j] for j in range(4))
            b = sum(A[3][j] * c[j] for j in range(4))
            if a == 0 and b == 0:
                continue
            if Q(a * a + b * b, d * d * den * den) < min_side_sq:
                continue
            good = True
            for i in range(4):
                f = sum(F[i][j] * c[j] for j in range(4)) * sg
                lo, hi = es[i][1] * d * sg, es[i][2] * d * sg
                lo, hi = min(lo, hi), max(lo, hi)
                if (strict and not lo < f < hi) or (not strict and not lo <= f <= hi):
                    good = False
                    break
            if good:
                x = Q(sum(A[0][j] * c[j] for j in range(4)), d * den)
                y = Q(sum(A[1][j] * c[j] for j in range(4)), d * den)
                sq = square_from((x, y), (Q(a, d * den), Q(b, d * den)))
                key = min(sq[j:] + sq[:j] for j in range(4))
                out[key] = Q(a * a + b * b, d * d * den * den)
    return sorted(((s, k) for k, s in out.items()), reverse=True)


def traced_component(polygon, seed):
    """(minimum squared side, absolute degree, number of charts) of the v341 component through seed."""
    from polygon_component_v341 import PolygonComponents, degree, segment_floor
    from polygon_orbit_census_v223 import pt
    engine = PolygonComponents([tuple(p) for p in polygon])
    cells = engine.trace(tuple(pt(*p) for p in seed))
    floor = min(segment_floor(cell.squares) for cell in cells.values())
    return floor, abs(sum(s for _, s in degree(cells)[1])), len(cells)


@dataclass(frozen=True)
class Certificate:
    level: str
    kind: str
    side_sq: Q
    seed: tuple = None


def certify(shape, use_traces=True):
    """Return a Certificate for the robust disk bound, or None when nothing was found."""
    theorem = covering_theorem(shape)
    if theorem:
        return Certificate('theorem', theorem, None)
    best = None
    for comp in closed_form_components(shape):
        if comp.odd and (best is None or comp.side_sq > best.side_sq):
            best = comp
    if best is not None and reaches_disk_target(shape, best.side_sq):
        return Certificate('closed-form', best.kind, best.side_sq, best.seed)
    polygon = shape.polygon()
    lo, hi = radius_bracket(shape)
    found = contact_squares(polygon, 2 * lo * lo, strict=True)
    if found and reaches_disk_target(shape, found[0][0]):
        return Certificate('transverse', 'regular-square', found[0][0], found[0][1])
    if not use_traces:
        return None
    from polygon_component_v341 import UnsupportedComponent
    seen = set()
    for side_sq, seed in contact_squares(polygon, 2 * lo * lo, strict=False):
        if seed in seen:
            continue
        try:
            floor, deg, _ = traced_component(polygon, seed)
        except (UnsupportedComponent, ArithmeticError, ValueError):
            continue
        seen.add(seed)
        if deg % 2 == 1 and reaches_disk_target(shape, floor):
            return Certificate('traced', 'odd-component', floor, seed)
    return None


# ---------------------------------------------------------------------------
# Finite redraw tolerance for a transverse square

def _edge_of(polygon, p):
    for (a, b) in zip(polygon, polygon[1:] + polygon[:1]):
        if a[1] == b[1] == p[1] and min(a[0], b[0]) < p[0] < max(a[0], b[0]):
            return 'h', min(a[0], b[0]), max(a[0], b[0])
        if a[0] == b[0] == p[0] and min(a[1], b[1]) < p[1] < max(a[1], b[1]):
            return 'v', min(a[1], b[1]), max(a[1], b[1])
    raise ValueError("vertex is not interior to an edge")


def _inverse(M):
    n = len(M)
    m = [list(row) + [Q(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for i in range(n):
        piv = next(j for j in range(i, n) if m[j][i])
        m[i], m[piv] = m[piv], m[i]
        v = m[i][i]
        m[i] = [x / v for x in m[i]]
        for j in range(n):
            if j != i and m[j][i]:
                f = m[j][i]
                m[j] = [x - f * y for x, y in zip(m[j], m[i])]
    return [row[n:] for row in m]


def redraw_tolerance(polygon, square, e):
    """For a transverse square with every vertex interior to an edge, return a rational lower bound for
    the side of a square on every registered Jordan redraw of Euclidean error at most e, or None if the
    sufficient test fails. Each residual coordinate changes by at most 4e; the reference residual is
    M (z - z*), so zeros stay in the box |z - z*| <= 4e ||M^-1|| and the box of radius rho (the
    smallest edge clearance) has a zero-free boundary when 4e ||M^-1|| < rho."""
    edges = [_edge_of(polygon, p) for p in square]
    dirs = [(Q(1), Q(0)) if kind == 'h' else (Q(0), Q(1)) for kind, lo, hi in edges]
    rho = min(min((p[0] if kind == 'h' else p[1]) - lo, hi - (p[0] if kind == 'h' else p[1]))
              for p, (kind, lo, hi) in zip(square, edges))

    def residual(P):
        p0, p1, p2, p3 = P
        ex, ey = p1[0] - p0[0], p1[1] - p0[1]
        return (p2[0] - p1[0] + ey, p2[1] - p1[1] - ex, p3[0] - p0[0] + ey, p3[1] - p0[1] - ex)

    base = residual(square)
    cols = []
    for i in range(4):
        moved = list(square)
        moved[i] = (square[i][0] + dirs[i][0], square[i][1] + dirs[i][1])
        cols.append([x - y for x, y in zip(residual(moved), base)])
    M = [[cols[j][i] for j in range(4)] for i in range(4)]
    inv_norm = max(sum(abs(x) for x in row) for row in _inverse(M))
    shift = 4 * Q(e) * inv_norm
    if not shift < rho:
        return None
    side_sq = (square[1][0] - square[0][0]) ** 2 + (square[1][1] - square[0][1]) ** 2
    # each of the first two vertices moves at most `shift` along its edge, then at most e off it
    loss = 2 * shift + 2 * Q(e)
    return side_sq, loss, rho, inv_norm
