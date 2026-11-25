"""v378: the local surgery rule for the large graceful parity, and near-square lattices at spiral points.

Surgery rule (Theorem 378.1): if two Jordan curves agree outside a ball of radius below s/2, their large graceful
parities at cutoff s differ by the number of graceful near-squares (squares of side s with three corners on the common
part and the fourth corner in the ball) whose free corner lies inside exactly one of the two curves.

This module finds near-squares on polygons in closed form and tests the rule with small bumps on exact drawings.
"""
from fractions import Fraction as Q
from math import hypot

import numpy as np

from asymmetric_spiral_v374 import large_graceful_parity
from wild_point_parity_v370 import orientation, sector


def _edges(V):
    n = len(V)
    P = np.array([[float(x), float(y)] for x, y in V])
    return P, np.roll(P, -1, axis=0)


def near_squares(V, s, center, radius, skip=(), tol=1e-12):
    """Squares A, B, C, D (counterclockwise, side s) with B, C, D on edges of V and A in the closed ball.

    Returns dicts with corners, edge indices and parameters. Edges in `skip` are not used for B, C, D."""
    P0, P1 = _edges(V)
    n = len(V)
    W = P1 - P0
    cx, cy = center
    # edges that can carry B or D (within s + radius of the center) or C (within s sqrt2 + radius)
    def dist_to_edges(R):
        d = P0 - np.array([cx, cy])
        t = np.clip(-(d * W).sum(1) / (W * W).sum(1), 0, 1)
        closest = P0 + t[:, None] * W - np.array([cx, cy])
        return np.hypot(closest[:, 0], closest[:, 1])
    dist = dist_to_edges(None)
    near_bd = [k for k in range(n) if dist[k] <= s + radius + 1e-9 and k not in skip]
    near_c = [k for k in range(n) if dist[k] <= s * 2 ** 0.5 + radius + 1e-9 and k not in skip]
    J = np.array([[0.0, -1.0], [1.0, 0.0]])
    M = {'B': np.eye(2), 'D': J, 'C': np.eye(2) + J}
    out = []
    for eb in near_bd:
        for ed in near_bd:
            for ec in near_c:
                rows, rhs = [], []
                for X, e in (('B', eb), ('C', ec), ('D', ed)):
                    w = W[e]
                    nrm = np.array([-w[1], w[0]])
                    rows.append(np.r_[nrm, s * nrm @ M[X]])
                    rhs.append(nrm @ P0[e])
                R, b = np.array(rows), np.array(rhs)
                u, sv, vt = np.linalg.svd(R)
                if sv[-1] < 1e-14:
                    continue
                kern = vt[-1]
                part = np.linalg.lstsq(R, b, rcond=None)[0]
                cp, ck = part[2:], kern[2:]
                qa, qb, qc = ck @ ck, 2 * cp @ ck, cp @ cp - 1
                disc = qb * qb - 4 * qa * qc
                if qa < 1e-30 or disc < 0:
                    continue
                for tau in ((-qb + disc ** 0.5) / (2 * qa), (-qb - disc ** 0.5) / (2 * qa)):
                    x = part + tau * kern
                    A, c = x[:2], x[2:]
                    if hypot(A[0] - cx, A[1] - cy) > radius:
                        continue
                    pts, params, ok = {}, {}, True
                    for X, e in (('B', eb), ('C', ec), ('D', ed)):
                        Z = A + s * (M[X] @ c)
                        tpar = ((Z - P0[e]) @ W[e]) / (W[e] @ W[e])
                        if not (tol < tpar < 1 - tol):
                            ok = False
                            break
                        pts[X], params[X] = Z, (e, tpar)
                    if ok:
                        out.append(dict(A=A, B=pts['B'], C=pts['C'], D=pts['D'], edges=params))
    # deduplicate (a root can repeat through different but touching edges only at vertices, excluded above)
    uniq = {}
    for q in out:
        key = tuple(np.round(np.r_[q['A'], q['B']], 10))
        uniq.setdefault(key, q)
    return list(uniq.values())


def graceful(q, a_pos, V):
    """Sector of the completed square when the free corner sits at curve position a_pos."""
    pos = (a_pos,) + tuple(q['edges'][X][0] + q['edges'][X][1] for X in 'BCD')
    entry = (None, None, pos, 1)
    return sector(entry, len(V), orientation(V))


def inside(V, z):
    """Even-odd point-in-polygon test (float)."""
    x, y = z
    P, Pn = _edges(V)
    c = False
    for (x0, y0), (x1, y1) in zip(P, Pn):
        if (y0 > y) != (y1 > y):
            xi = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if xi > x:
                c = not c
    return c


def bump(V, k, t1, t2, apex):
    """Replace the piece of edge k between parameters t1 < t2 by two segments through `apex` (exact)."""
    a, b = V[k], V[(k + 1) % len(V)]
    P1 = (a[0] + t1 * (b[0] - a[0]), a[1] + t1 * (b[1] - a[1]))
    P2 = (a[0] + t2 * (b[0] - a[0]), a[1] + t2 * (b[1] - a[1]))
    return V[:k + 1] + [P1, apex, P2] + V[k + 1:]


def surgery_trial(V, s, k, t1, t2, apex):
    """Exact parity change vs the rule's prediction for one bump on edge k."""
    t1, t2 = Q(t1), Q(t2)
    W = bump(V, k, t1, t2, apex)
    a, b = V[k], V[(k + 1) % len(V)]
    pts = [(float(a[0] + t * (b[0] - a[0])), float(a[1] + t * (b[1] - a[1]))) for t in (t1, t2)]
    tri = [pts[0], (float(apex[0]), float(apex[1])), pts[1]]
    cen = tuple(sum(p[i] for p in tri) / 3 for i in range(2))
    rad = max(hypot(p[0] - cen[0], p[1] - cen[1]) for p in tri) + 1e-12
    if rad >= float(s) / 2:
        raise ValueError('bump does not fit in a ball of radius s/2')
    qs = near_squares(V, float(s), cen, rad)
    a_pos = k + float((t1 + t2) / 2)
    swept = [q for q in qs if inside(V, q['A']) != inside(W, q['A'])]
    pred = sum(1 for q in swept if graceful(q, a_pos, V) != 'crossed') % 2
    p0 = large_graceful_parity(V, Q(s))[0]
    p1 = large_graceful_parity(W, Q(s))[0]
    return dict(edge=k, swept=len(swept), swept_graceful=sum(1 for q in swept if graceful(q, a_pos, V) != 'crossed'),
                predicted_flip=pred, actual_flip=(p0 + p1) % 2, near_in_ball=len(qs))


# ---------------------------------------------------------------- near-square lattice at a spiral point
def _h(x, y, s, which, arm, lam, a=0.5):
    """Comb phase of corner `which` (frame coordinates, theta removed), shifted so the arm-`arm` condition is 2 pi Z."""
    X, Y = {'B': (s + x, y), 'D': (x, s + y), 'C': (s + x, s + y)}[which]
    r = np.hypot(X, Y)
    val = r ** (-1 / a) - np.arctan2(Y, X)
    if arm == 2:
        val = val + np.pi + lam * 0.4 * np.sin(5 * np.log(1 / r))
    return val


def lattice_map(X, Y, s, arms, lam, a=0.5):
    """P(X, Y) = (h_B - h_D, h_B - h_C) with A = Ks (X, Y) in the rotating frame; near-squares are P^{-1}(2 pi Z^2)."""
    Ks = a * s ** (1 + 1 / a)
    x, y = Ks * X, Ks * Y
    hB = _h(x, y, s, 'B', arms[0], lam, a)
    hD = _h(x, y, s, 'D', arms[1], lam, a)
    hC = _h(x, y, s, 'C', arms[2], lam, a)
    return np.array([hB - hD, hB - hC]), hB


def lattice_count(s, M, lam, a=0.5):
    """Near-squares of side s with |A| <= M K(s) s, enumerated through the lattice map, for all 8 arm triples."""
    c = 2 ** (-1 / (2 * a) - 1)
    Lp = np.array([[-1.0, 1.0], [-(1 - c), c]])
    Li = np.linalg.inv(Lp)
    found = []
    for arms in [(i, j, k) for i in (1, 2) for j in (1, 2) for k in (1, 2)]:
        P0, _ = lattice_map(0.0, 0.0, s, arms, lam, a)
        # lattice points whose affine preimage lies within M + 1
        span = np.abs(Lp).sum(1) * (M + 1)
        lo = np.floor((P0 - span) / (2 * np.pi)).astype(int)
        hi = np.ceil((P0 + span) / (2 * np.pi)).astype(int)
        for m1 in range(lo[0], hi[0] + 1):
            for m2 in range(lo[1], hi[1] + 1):
                tgt = 2 * np.pi * np.array([m1, m2])
                z = Li @ (tgt - P0)
                if z @ z > (M + 1.5) ** 2:
                    continue
                for _ in range(30):
                    Pz, _ = lattice_map(z[0], z[1], s, arms, lam, a)
                    F = Pz - tgt
                    if np.max(np.abs(F)) < 1e-11:
                        break
                    h = 1e-6
                    Jm = np.column_stack([(lattice_map(z[0] + h, z[1], s, arms, lam, a)[0] - Pz) / h,
                                          (lattice_map(z[0], z[1] + h, s, arms, lam, a)[0] - Pz) / h])
                    z = z - np.linalg.solve(Jm, F)
                if np.max(np.abs(lattice_map(z[0], z[1], s, arms, lam, a)[0] - tgt)) < 1e-9 and z @ z <= M * M:
                    found.append((arms, z))
    return found


def predicted_count(M, a=0.5):
    c = 2 ** (-1 / (2 * a) - 1)
    return 2 * abs(1 - 2 * c) * M * M / np.pi


def lattice_sector(arms, z, s, lam, a=0.5):
    """Sector of a lattice near-square at the wild point: along the curve, arm one runs inward (Theta increasing),
    the wild point comes next, then arm two runs outward. The free corner A sits at the wild point."""
    Ks = a * s ** (1 + 1 / a)
    x, y = Ks * z[0], Ks * z[1]
    # (block, position): arm one with Theta increasing, then the wild point, then arm two with Theta decreasing.
    # (Float keys big - Theta round to big and tie whenever two corners are on arm two.)
    keys = [(1, 0.0)]
    for which, arm in zip('BCD', (arms[0], arms[2], arms[1])):
        X, Y = {'B': (s + x, y), 'D': (x, s + y), 'C': (s + x, s + y)}[which]
        Th = np.hypot(X, Y) ** (-1 / a)
        keys.append((0, Th) if arm == 1 else (2, -Th))
    order = sorted(range(4), key=lambda i: keys[i])
    k = order.index(0)
    cyc = tuple(order[k:] + order[:k])
    return 'graceful' if cyc in ((0, 1, 2, 3), (0, 3, 2, 1)) else 'crossed'
