"""v370: exact small-square parity near one spiral wild point.

A double spiral: arm one is r(theta) e^{i theta}, arm two is its half-turn image, both for theta
from theta0 to theta_n. At depth n the inner ends are joined by a polygonal connector inside the disk
of radius r(theta_n), and the outer ends are closed by an arc outside radius r(theta0). Every vertex
is rational (radii rounded, directions from rational points of the unit circle).

The counter finds every inscribed square of side below s whose corners are transverse (each corner in
the interior of an edge, nonsingular 4x4 contact system), optionally only those touching a given set
of edges. numpy prunes; every reported square is recomputed and checked in exact arithmetic.
"""
from fractions import Fraction as Q
from math import cos, pi, sin, tan

import numpy as np

DEN = 2 ** 24


def _unit(theta):
    """A rational point on the unit circle close to angle theta."""
    t = theta / 2
    k = round(t / pi)
    t -= k * pi
    if abs(abs(t) - pi / 2) < 1e-9:
        return (Q(-1), Q(0))
    u = Q(round(tan(t) * DEN), DEN)
    return ((1 - u * u) / (1 + u * u), 2 * u / (1 + u * u))


def radius(theta, a, scale=1.0):
    return Q(round(scale * theta ** (-a) * DEN), DEN)


def arm(a, theta0, theta1, per_turn, scale=1.0):
    """Arm-one vertices from theta0 inward to theta1 (theta increasing)."""
    steps = max(1, round((theta1 - theta0) * per_turn / (2 * pi)))
    pts = []
    for k in range(steps + 1):
        th = theta0 + (theta1 - theta0) * k / steps
        r = radius(th, a, scale)
        ux, uy = _unit(th)
        pts.append((r * ux, r * uy))
    return pts


def connector(P, kind):
    """Polygonal arc from P to -P strictly inside the disk |z| < |P| (except its ends)."""
    x, y = P
    if kind == 'straight':
        return [P, (-x, -y)]
    if kind == 'bent':          # through a point off the diameter
        return [P, (-y / 3, x / 3), (-x, -y)]
    if kind == 'zigzag':        # three kinks
        return [P, (x / 2 - y / 4, y / 2 + x / 4), (y / 4, -x / 4),
                (-x / 2 - y / 4, -y / 2 + x / 4), (-x, -y)]
    if kind == 'hook':          # loops around one side of the centre
        return [P, (x / 2 + y / 2, y / 2 - x / 2), (-x / 2 + y / 2, -y / 2 - x / 2),
                (-x * Q(3, 4) - y / 4, -y * Q(3, 4) + x / 4), (-x, -y)]
    raise ValueError(kind)


def drawing(a, theta0, half_turns, per_turn, kind='straight', scale=1.0, outer=Q(2), twist=0.0):
    """Vertex list of the depth-`half_turns` drawing, counterclockwise or clockwise as it falls.
    Returns (vertices, groups) where groups names index ranges: 'arm1', 'conn', 'arm2', 'outer'."""
    theta1 = theta0 + half_turns * pi
    A = arm(a, theta0, theta1, per_turn, scale)
    inner = A[-1]
    if twist:
        # arm two is arm one turned by pi + twist (exact rational rotation); the connector is a chord
        c, sn = _unit(pi + twist)
        B = [(c * x - sn * y, sn * x + c * y) for x, y in reversed(A)]
        conn = [inner, B[0]]
    else:
        conn = connector(inner, kind)
        B = [(-x, -y) for x, y in reversed(A)]              # arm two, from the centre outward
    z0 = A[0]
    r0 = (z0[0] ** 2 + z0[1] ** 2)
    R = outer
    # radial-ish exits and a polygonal outer arc on the side of angle theta0 + pi/2
    import math
    nrm = Q(round(math.sqrt(float(r0)) * DEN), DEN)
    far = lambda p: (p[0] * R / nrm, p[1] * R / nrm)
    ux, uy = _unit(theta0 + pi / 2)
    top = (-uy * 0 + ux * R * 2, uy * R * 2)
    V = A[:-1] + conn[:-1] + B[:-1] + [B[-1], far(B[-1]), top, far(A[0])]
    n1 = len(A) - 1
    n2 = n1 + len(conn) - 1
    n3 = n2 + len(B)
    groups = {'arm1': range(0, n1), 'conn': range(n1, n2), 'arm2': range(n2, n3 - 1),
              'outer': range(n3 - 1, len(V))}
    return V, groups


# ---------------------------------------------------------------------------
# Exact transverse square counting

def _slot_matrices(d, sigma):
    """Residual columns contributed by an edge direction d in each of the four slots."""
    Jd = (-sigma * d[1], sigma * d[0])
    return [((Jd[0], Jd[1]), (-d[0] + Jd[0], -d[1] + Jd[1])),     # slot 0: p0
            ((-d[0] - Jd[0], -d[1] - Jd[1]), (-Jd[0], -Jd[1])),     # slot 1: p1
            ((d[0], d[1]), (0, 0)),                                  # slot 2: p2
            ((0, 0), (d[0], d[1]))]                                  # slot 3: p3


def _residual(P, sigma):
    p0, p1, p2, p3 = P
    v = (p1[0] - p0[0], p1[1] - p0[1])
    Jv = (-sigma * v[1], sigma * v[0])
    return (p2[0] - p1[0] - Jv[0], p2[1] - p1[1] - Jv[1], p3[0] - p0[0] - Jv[0], p3[1] - p0[1] - Jv[1])


def _exact_square(E, ids, sigma):
    A = [E[i][0] for i in ids]
    D = [(E[i][1][0] - E[i][0][0], E[i][1][1] - E[i][0][1]) for i in ids]
    c = _residual(A, sigma)
    cols = []
    for j in range(4):
        moved = list(A)
        moved[j] = E[ids[j]][1]
        r = _residual(moved, sigma)
        cols.append([x - y for x, y in zip(r, c)])
    M = [[cols[j][i] for j in range(4)] for i in range(4)]
    rhs = [-x for x in c]
    # exact Gaussian elimination
    m = [row[:] + [b] for row, b in zip(M, rhs)]
    for i in range(4):
        piv = next((r for r in range(i, 4) if m[r][i] != 0), None)
        if piv is None:
            return None
        m[i], m[piv] = m[piv], m[i]
        p = m[i][i]
        m[i] = [x / p for x in m[i]]
        for r in range(4):
            if r != i and m[r][i] != 0:
                f = m[r][i]
                m[r] = [x - f * y for x, y in zip(m[r], m[i])]
    u = [m[i][4] for i in range(4)]
    pts = tuple((a[0] + t * d[0], a[1] + t * d[1]) for a, d, t in zip(A, D, u))
    return u, pts


def small_squares(V, s, touching=None, tol=1e-9):
    """All transverse inscribed squares of side < s (exact). If `touching` is a set of edge indices,
    only squares with a corner on one of those edges. Returns a dict keyed by the corner set with
    (side_sq, corners in square order, curve positions) and a list of degenerate candidates."""
    n = len(V)
    E = [(V[i], V[(i + 1) % n]) for i in range(n)]
    Ef = np.array([[float(a[0]), float(a[1]), float(b[0]), float(b[1])] for a, b in E])
    s_f = float(s)
    lo = np.minimum(Ef[:, :2], Ef[:, 2:])
    hi = np.maximum(Ef[:, :2], Ef[:, 2:])
    reach = s_f * 1.4143 + 1e-12
    seeds = range(n) if touching is None else sorted(touching)
    out, degenerate = {}, []
    s2 = Q(s) ** 2
    for e in seeds:
        gap = np.maximum(0, np.maximum(lo - hi[e], lo[e] - hi))
        near = np.nonzero((gap[:, 0] <= reach) & (gap[:, 1] <= reach))[0]
        k = len(near)
        if k == 0:
            continue
        A = Ef[near, :2]
        D = Ef[near, 2:] - Ef[near, :2]
        Ae, De = Ef[e, :2], Ef[e, 2:] - Ef[e, :2]
        for sigma in (1, -1):
            for slot in range(4):
                others = [j for j in range(4) if j != slot]
                # grids over the three other slots
                idx = np.stack(np.meshgrid(np.arange(k), np.arange(k), np.arange(k), indexing='ij'), -1).reshape(-1, 3)
                T = idx.shape[0]
                Apts = np.zeros((T, 4, 2))
                Dirs = np.zeros((T, 4, 2))
                Apts[:, slot] = Ae
                Dirs[:, slot] = De
                for col, j in enumerate(others):
                    Apts[:, j] = A[idx[:, col]]
                    Dirs[:, j] = D[idx[:, col]]
                p0, p1, p2, p3 = (Apts[:, j] for j in range(4))
                v = p1 - p0
                Jv = np.stack([-sigma * v[:, 1], sigma * v[:, 0]], -1)
                c = np.concatenate([p2 - p1 - Jv, p3 - p0 - Jv], -1)
                M = np.zeros((T, 4, 4))
                for j in range(4):
                    d = Dirs[:, j]
                    Jd = np.stack([-sigma * d[:, 1], sigma * d[:, 0]], -1)
                    if j == 0:
                        M[:, 0:2, 0] = Jd
                        M[:, 2:4, 0] = -d + Jd
                    elif j == 1:
                        M[:, 0:2, 1] = -d - Jd
                        M[:, 2:4, 1] = -Jd
                    elif j == 2:
                        M[:, 0:2, 2] = d
                    else:
                        M[:, 2:4, 3] = d
                det = np.linalg.det(M)
                # scale-relative test: det is a product of four edge-vector columns, so an absolute cutoff
                # would discard every candidate on short edges (edges of length 1e-5 give det ~ 1e-20)
                scale = np.prod(np.linalg.norm(Dirs, axis=2), axis=1)
                good = np.abs(det) > 1e-12 * scale
                u = np.full((T, 4), np.nan)
                if good.any():
                    u[good] = np.linalg.solve(M[good], -c[good][..., None])[..., 0]
                inside = good & np.all((u > -tol) & (u < 1 + tol), axis=1)
                if not inside.any():
                    continue
                q0 = p0 + u[:, 0:1] * Dirs[:, 0]
                q1 = p1 + u[:, 1:2] * Dirs[:, 1]
                side2 = np.sum((q1 - q0) ** 2, axis=1)
                cand = np.nonzero(inside & (side2 < s_f * s_f * (1 + 1e-9)) & (side2 > 0))[0]
                for t in cand:
                    ids = [0] * 4
                    ids[slot] = e
                    for col, j in enumerate(others):
                        ids[j] = int(near[idx[t, col]])
                    res = _exact_square(E, ids, sigma)
                    if res is None:
                        degenerate.append(('singular', tuple(ids), sigma))
                        continue
                    uu, pts = res
                    side_sq = (pts[1][0] - pts[0][0]) ** 2 + (pts[1][1] - pts[0][1]) ** 2
                    if not (0 < side_sq < s2) or not all(0 <= x <= 1 for x in uu):
                        continue
                    if any(x == 0 or x == 1 for x in uu):
                        degenerate.append(('vertex', tuple(ids), sigma))
                        continue
                    key = frozenset(pts)
                    if len(key) < 4:
                        continue
                    pos = tuple(Q(i) + x for i, x in zip(ids, uu))
                    out.setdefault(key, (side_sq, pts, pos, sigma))
    return out, degenerate


def sector(entry, n, orientation):
    """'aligned', 'reversed' or 'crossed': whether the corners occur along the curve in the square's
    cyclic order, with the curve's orientation (aligned) or against it (reversed)."""
    side_sq, pts, pos, sigma = entry
    order = sorted(range(4), key=lambda i: pos[i])
    k = order.index(0)
    cyc = tuple(order[k:] + order[:k])
    if cyc == (0, 1, 2, 3):
        geometric = sigma          # labels follow the square counterclockwise when sigma = 1
    elif cyc == (0, 3, 2, 1):
        geometric = -sigma
    else:
        return 'crossed'
    return 'aligned' if geometric == orientation else 'reversed'


def orientation(V):
    area = sum(V[i][0] * V[(i + 1) % len(V)][1] - V[(i + 1) % len(V)][0] * V[i][1] for i in range(len(V)))
    return 1 if area > 0 else -1


def sector_counts(V, s, touching=None):
    sq, degenerate = small_squares(V, s, touching)
    o = orientation(V)
    counts = {'aligned': 0, 'reversed': 0, 'crossed': 0}
    for entry in sq.values():
        counts[sector(entry, len(V), o)] += 1
    return counts, degenerate, sq
