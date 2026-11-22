"""v331: complete small-square census on spiral polygons, filtered by cap action residue.

At the square angle the upper spectral value is the action of some inscribed square
(spectrality). A small square whose action residue lies closer to A Z than the spectral
gap G = A - pi R^2/2 cannot be that carrier. This module lists every inscribed square of
half-diagonal < 1 whose residue survives that test.

Method:
1. Vectorized floating search: for each ordered edge pair (two adjacent corners, both
   orientations), solve for the third corner on every edge, then intersect the fourth
   corner's segment with every edge. Tolerance-widened, so it over-reports.
2. Each candidate's size and residue are computed in floating point. Only squares with
   q^2 < 1 and residue distance >= G - margin are kept.
3. Kept candidates are certified exactly: the two line conditions are solved in rationals,
   all four contacts are checked on their edges, and the residue is recomputed exactly.

Validated against the exhaustive v223 census (tests). No Floer homology is computed.
"""
from fractions import Fraction as Q
from multiprocessing import Pool
import json
import sys

import numpy as np

from angular_windows_v328_v330 import area, pi_enclosure, primitive, radius_squared, square_with_gadgets


def midside_fixture(N, eps=Q(1, 2), height=Q(2)):
    """[0,4]^2 with one v324 gadget of N turns at scale eps attached at mid-left side.

    Attached there, the gadget leaves the diameter at 4 sqrt 2 (R^2 = 8) for eps <= 1/2,
    so the spectral gap is positive already at this coarse scale."""
    return square_with_gadgets([(N, Q(eps), Q(height))])


# ---------------------------------------------------------------- floating search

def _search_rows(P, rows, qmax2, keep_dist, area_f, count_only=False, collect=False, tol=1e-9):
    n = len(P)
    A = np.array([[float(x), float(y)] for x, y in P])
    D = np.roll(A, -1, axis=0) - A
    dd = np.sum(D*D, axis=1)
    crossAD = A[:, 0]*D[:, 1]-A[:, 1]*D[:, 0]
    Fcum = np.concatenate([[0.0], np.cumsum(0.5*crossAD)])[:n]
    smax2 = 2*qmax2
    kept, small, best, degenerate = set(), 0, 0.0, 0
    best_count = [0]
    flagged = []

    def clip(c0, c1, lo, hi):          # c0 + c1*lam >= 0
        with np.errstate(divide='ignore', invalid='ignore'):
            r = -c0/c1
        pos, neg = c1 > 1e-15, c1 < -1e-15
        lo = np.where(pos, np.maximum(lo, r), lo)
        hi = np.where(neg, np.minimum(hi, r), hi)
        return lo, np.where(~(pos | neg) & (c0 < 0), -np.inf, hi)

    for sigma in (1, -1):
        J = np.array([[0, -1], [1, 0]])*sigma
        for i in rows:
            ai, di = A[i], D[i]
            for j in range(n):
                aj, dj = A[j], D[j]
                gap = [max(0.0, max(min(ai[c], ai[c]+di[c]), min(aj[c], aj[c]+dj[c]))
                           - min(max(ai[c], ai[c]+di[c]), max(aj[c], aj[c]+dj[c]))) for c in (0, 1)]
                if gap[0]**2+gap[1]**2 >= smax2:
                    continue
                w0 = aj-ai
                P3, A3, B3 = ai+J@w0, di-J@di, J@dj            # fourth corner v3 = P3 + s A3 + t B3
                P2, A2, B2 = aj+J@w0, -(J@di), dj+J@dj         # third corner  v2 = P2 + s A2 + t B2
                alpha = D[:, 0]*A3[1]-D[:, 1]*A3[0]
                beta = D[:, 0]*B3[1]-D[:, 1]*B3[0]
                gamma = -(D[:, 0]*(P3[1]-A[:, 1])-D[:, 1]*(P3[0]-A[:, 0]))
                nrm = np.abs(alpha)+np.abs(beta)
                deg = nrm < 1e-14
                degenerate += int(np.sum(deg & (np.abs(gamma) < 1e-12)))
                if collect:
                    flagged.extend((i, j, int(k), None, sigma) for k in np.nonzero(deg & (np.abs(gamma) < 1e-12))[0])
                ok = ~deg
                a2b2 = np.where(ok, alpha**2+beta**2, 1.0)
                p0s, p0t = gamma*alpha/a2b2, gamma*beta/a2b2
                ds, dt = beta/np.where(ok, nrm, 1.0), -alpha/np.where(ok, nrm, 1.0)
                lo, hi = np.full(n, -np.inf), np.full(n, np.inf)
                for c0, c1 in ((p0s, ds), (1-p0s, -ds), (p0t, dt), (1-p0t, -dt)):
                    lo, hi = clip(c0+tol, c1, lo, hi)
                u0 = ((P3[0]+p0s*A3[0]+p0t*B3[0]-A[:, 0])*D[:, 0]+(P3[1]+p0s*A3[1]+p0t*B3[1]-A[:, 1])*D[:, 1])/dd
                u1 = ((ds*A3[0]+dt*B3[0])*D[:, 0]+(ds*A3[1]+dt*B3[1])*D[:, 1])/dd
                lo, hi = clip(u0+tol, u1, lo, hi)
                lo, hi = clip(1-u0+tol, -u1, lo, hi)
                ks = np.nonzero(ok & (lo <= hi))[0]
                if len(ks) == 0:
                    continue
                s0, t0 = p0s[ks]+lo[ks]*ds[ks], p0t[ks]+lo[ks]*dt[ks]
                s1, t1 = p0s[ks]+hi[ks]*ds[ks], p0t[ks]+hi[ks]*dt[ks]
                Q0 = P2+np.outer(s0, A2)+np.outer(t0, B2)
                E = P2+np.outer(s1, A2)+np.outer(t1, B2)-Q0
                den = E[:, None, 0]*D[None, :, 1]-E[:, None, 1]*D[None, :, 0]
                rx, ry = A[None, :, 0]-Q0[:, None, 0], A[None, :, 1]-Q0[:, None, 1]
                with np.errstate(divide='ignore', invalid='ignore'):
                    lam = (rx*D[None, :, 1]-ry*D[None, :, 0])/den
                    mu = (rx*E[:, None, 1]-ry*E[:, None, 0])/den
                good = (np.abs(den) > 1e-15)&(lam >= -tol)&(lam <= 1+tol)&(mu >= -tol)&(mu <= 1+tol)
                point = np.sum(np.abs(E), axis=1) < 1e-14
                cr = D[None, :, 0]*(Q0[:, None, 1]-A[None, :, 1])-D[None, :, 1]*(Q0[:, None, 0]-A[None, :, 0])
                if point.any():
                    mu2 = ((Q0[:, None, 0]-A[None, :, 0])*D[None, :, 0]+(Q0[:, None, 1]-A[None, :, 1])*D[None, :, 1])/dd[None, :]
                    good |= point[:, None]&(np.abs(cr) < 1e-12)&(mu2 >= -tol)&(mu2 <= 1+tol)
                    lam = np.where(point[:, None], 0.0, lam)
                    mu = np.where(point[:, None], mu2, mu)
                collinear = (np.abs(den) <= 1e-15)&~point[:, None]&(np.abs(cr) < 1e-12)
                degenerate += int(np.sum(collinear))
                if collect:
                    flagged.extend((i, j, int(ks[a]), int(l), sigma) for a, l in zip(*np.nonzero(collinear)))
                aa, ll = np.nonzero(good)
                if len(aa) == 0:
                    continue
                la, kk = np.clip(lam[aa, ll], 0, 1), ks[aa]
                s, t = s0[aa]+la*(s1[aa]-s0[aa]), t0[aa]+la*(t1[aa]-t0[aa])
                w = aj+np.outer(t, dj)-ai-np.outer(s, di)
                side2 = np.sum(w*w, axis=1)
                size_ok = (side2/2 < qmax2)&(side2 > 1e-18)
                small += int(np.sum(size_ok))
                v3 = P3+np.outer(s, A3)+np.outer(t, B3)
                u3 = np.sum((v3-A[kk])*D[kk], axis=1)/dd[kk]
                F = lambda m, u: Fcum[m]+0.5*u*crossAD[m]
                F0, F1, F3, F2 = F(i, s), F(j, t), F(kk, u3), F(ll, mu[aa, ll])
                Jsum = (F0-F1)+(F2-F3) if sigma == 1 else (F0-F3)+(F2-F1)
                r = np.mod(-Jsum, area_f)
                dist = np.minimum(r, area_f-r)
                if size_ok.any():
                    best = max(best, float(np.max(np.where(size_ok, dist, 0.0))))
                passing = size_ok & (dist >= keep_dist)
                if count_only:
                    best_count[0] += int(np.sum(passing))
                    continue
                for a, l in zip(aa[passing], ll[passing]):
                    kept.add((i, j, int(ks[a]), int(l), sigma))
    if collect:
        return flagged
    if count_only:
        return best_count[0], small, best, degenerate
    return kept, small, best, degenerate


def _job(args):
    return _search_rows(*args)


def _collect_job(args):
    P, rows, qmax2 = args
    return _search_rows(P, rows, qmax2, -1.0, 1.0, collect=True)


def float_search(P, qmax2=1.0, keep_dist=-1.0, area_f=None, workers=1, chunks=16, count_only=False):
    """Floating candidate search. With count_only, return the number of raw edge-quadruple hits
    passing the filters instead of the set (each square is hit several times)."""
    n = len(P)
    if area_f is None:
        area_f = float(area(P))
    jobs = [(P, list(range(r, n, chunks)), qmax2, keep_dist, area_f, count_only) for r in range(chunks)]
    if workers > 1:
        with Pool(workers) as pool:
            parts = pool.map(_job, jobs)
    else:
        parts = [_job(job) for job in jobs]
    first = sum(p[0] for p in parts) if count_only else set().union(*[p[0] for p in parts])
    return (first, sum(p[1] for p in parts),
            max(p[2] for p in parts), sum(p[3] for p in parts))


# ---------------------------------------------------------------- exact certification

def certify(P, quadruples, qmax2=Q(1), F=None, A=None):
    """Exact squares behind candidate edge quadruples.

    Returns ({frozenset(vertices): (ccw vertex tuple, residue distance)}, singular count).
    The residue uses (323.1) with the primitive evaluated on the edges the contacts lie on;
    at a shared vertex both edges give the same value."""
    n = len(P)
    edges = [(P[m], P[(m+1) % n]) for m in range(n)]
    F = F or primitive(P+[P[0]])
    A = A or area(P)
    Fon = lambda m, p: F[m]+(edges[m][0][0]*p[1]-edges[m][0][1]*p[0])/2
    found, singular = {}, 0
    for i, j, k, l, sigma in quadruples:
        (ai, bi), (aj, bj), (ak, bk), (al, bl) = edges[i], edges[j], edges[k], edges[l]
        vec = lambda a, b: (b[0]-a[0], b[1]-a[1])
        rot = lambda v: (-sigma*v[1], sigma*v[0])
        di, dj, dk, dl = vec(ai, bi), vec(aj, bj), vec(ak, bk), vec(al, bl)
        w0, Jdi, Jdj = rot(vec(ai, aj)), rot(di), rot(dj)
        P3, A3, B3 = (ai[0]+w0[0], ai[1]+w0[1]), (di[0]-Jdi[0], di[1]-Jdi[1]), Jdj
        P2, A2, B2 = (aj[0]+w0[0], aj[1]+w0[1]), (-Jdi[0], -Jdi[1]), (dj[0]+Jdj[0], dj[1]+Jdj[1])
        row = lambda d, a, P0, As, Bs: (d[0]*As[1]-d[1]*As[0], d[0]*Bs[1]-d[1]*Bs[0],
                                        -(d[0]*(P0[1]-a[1])-d[1]*(P0[0]-a[0])))
        r1, r2 = row(dk, ak, P3, A3, B3), row(dl, al, P2, A2, B2)
        det = r1[0]*r2[1]-r1[1]*r2[0]
        if det == 0:
            singular += 1
            continue
        s = (r1[2]*r2[1]-r1[1]*r2[2])/det
        t = (r1[0]*r2[2]-r1[2]*r2[0])/det
        if not (0 <= s <= 1 and 0 <= t <= 1):
            continue
        v0, v1 = (ai[0]+s*di[0], ai[1]+s*di[1]), (aj[0]+t*dj[0], aj[1]+t*dj[1])
        w = rot(vec(v0, v1))
        v3, v2 = (v0[0]+w[0], v0[1]+w[1]), (v1[0]+w[0], v1[1]+w[1])
        on = lambda p, a, d: (d[0]*(p[1]-a[1])-d[1]*(p[0]-a[0])) == 0 and \
            0 <= (p[0]-a[0])*d[0]+(p[1]-a[1])*d[1] <= d[0]**2+d[1]**2
        side2 = (v1[0]-v0[0])**2+(v1[1]-v0[1])**2
        if side2 == 0 or side2/2 >= qmax2 or not (on(v3, ak, dk) and on(v2, al, dl)):
            continue
        key = frozenset((v0, v1, v2, v3))
        if key in found:
            continue
        if sigma == 1:
            verts, F_ccw = (v0, v1, v2, v3), (Fon(i, v0), Fon(j, v1), Fon(l, v2), Fon(k, v3))
        else:
            verts, F_ccw = (v0, v3, v2, v1), (Fon(i, v0), Fon(k, v3), Fon(l, v2), Fon(j, v1))
        r = (-((F_ccw[0]-F_ccw[1])+(F_ccw[2]-F_ccw[3]))) % A
        found[key] = (verts, min(r, A-r))
    return found, singular


def _certify_job(args):
    return certify(*args)


def certify_parallel(P, quadruples, workers=1, chunks=32):
    quads = sorted(quadruples)
    F, A = primitive(P+[P[0]]), area(P)
    jobs = [(P, quads[r::chunks], Q(1), F, A) for r in range(chunks)]
    if workers > 1:
        with Pool(workers) as pool:
            parts = pool.map(_certify_job, jobs)
    else:
        parts = [_certify_job(job) for job in jobs]
    found = {}
    for part, _ in parts:
        found.update(part)
    return found, sum(sing for _, sing in parts)


def residue_distance(P, square, F=None, A=None):
    """Exact distance from A Z of the cap residue of a counterclockwise square (c0,c1,c2,c3).

    Rotation by +pi/2 about the centre sends c0 -> c1 and c2 -> c3; the return arcs run
    c1 -> c0 and c3 -> c2 along the curve, so by (323.1) the residue is
    -[(F(c0)-F(c1)) + (F(c2)-F(c3))] modulo A."""
    n = len(P)
    F = F or primitive(P+[P[0]])
    A = A or area(P)

    def at(p):
        for m in range(n):
            a, b = P[m], P[(m+1) % n]
            d = (b[0]-a[0], b[1]-a[1])
            if d[0]*(p[1]-a[1])-d[1]*(p[0]-a[0]) == 0 and 0 <= (p[0]-a[0])*d[0]+(p[1]-a[1])*d[1] <= d[0]**2+d[1]**2:
                return F[m]+(a[0]*p[1]-a[1]*p[0])/2
        raise ValueError("point not on polygon")
    c0, c1, c2, c3 = square
    r = (-((at(c0)-at(c1))+(at(c2)-at(c3)))) % A
    return min(r, A-r)


def _feasible_segment(P, case):
    """Exact (s,t)-set where all four contacts of a flagged edge tuple hold, as a list of extreme
    points (empty, one point, or two endpoints of a segment)."""
    i, j, k, l, sigma = case
    n = len(P)
    E = lambda m: (P[m], P[(m+1) % n])
    vec = lambda a, b: (b[0]-a[0], b[1]-a[1])
    rot = lambda v: (-sigma*v[1], sigma*v[0])
    (ai, bi), (aj, bj) = E(i), E(j)
    di, dj = vec(ai, bi), vec(aj, bj)
    w0, Jdi, Jdj = rot(vec(ai, aj)), rot(di), rot(dj)
    v3 = ((ai[0]+w0[0], ai[1]+w0[1]), (di[0]-Jdi[0], di[1]-Jdi[1]), Jdj)
    v2 = ((aj[0]+w0[0], aj[1]+w0[1]), (-Jdi[0], -Jdi[1]), (dj[0]+Jdj[0], dj[1]+Jdj[1]))
    eqs, ineqs = [], [(1, 0, 0), (-1, 0, 1), (0, 1, 0), (0, -1, 1)]     # a s + b t + c >= 0
    for edge, (P0, As, Bs) in (([k], v3), ([l] if l is not None else [], v2)):
        for m in edge:
            a, b = E(m)
            d = vec(a, b)
            dd = d[0]**2+d[1]**2
            eqs.append((d[0]*As[1]-d[1]*As[0], d[0]*Bs[1]-d[1]*Bs[0], d[0]*(P0[1]-a[1])-d[1]*(P0[0]-a[0])))
            u = ((As[0]*d[0]+As[1]*d[1]), (Bs[0]*d[0]+Bs[1]*d[1]), (P0[0]-a[0])*d[0]+(P0[1]-a[1])*d[1])
            ineqs += [u, (-u[0], -u[1], dd-u[2])]
    if l is None:   # type A: the third corner still has to lie on some edge; handled by the caller
        return None
    eqs = [e for e in eqs if (e[0], e[1]) != (0, 0) or e[2] != 0]
    if any((e[0], e[1]) == (0, 0) for e in eqs):
        return []
    lines = eqs[:1]
    if len(eqs) == 2 and eqs[0][0]*eqs[1][1]-eqs[0][1]*eqs[1][0] != 0:
        return "isolated"
    if not lines:
        return "two-dimensional"
    a, b, c = lines[0]
    # parametrize a s + b t + c = 0 as (s,t) = p + lam*(b,-a)
    p0 = (-c*a/(a*a+b*b), -c*b/(a*a+b*b))
    lo, hi = None, None
    for x, y, z in ineqs:
        c0, c1 = x*p0[0]+y*p0[1]+z, x*b-y*a
        if c1 == 0:
            if c0 < 0:
                return []
            continue
        r = -c0/c1
        if c1 > 0:
            lo = r if lo is None else max(lo, r)
        else:
            hi = r if hi is None else min(hi, r)
    if lo is None or hi is None or lo > hi:
        return []
    return [(p0[0]+lam*b, p0[1]-lam*a) for lam in ((lo,) if lo == hi else (lo, hi))]


def degenerate_audit(P, gap, qmax2=Q(1), workers=1, chunks=16):
    """Classify every collinear configuration the floating search skipped, exactly.

    Returns (counts, families, point_survivors): families are positive-size continua with a
    square of half-diagonal below the cutoff, each with its residue distances at the two ends and
    the midpoint; point_survivors are isolated squares at range boundaries whose residue reaches
    the gap (these must already be among the certified survivors)."""
    n = len(P)
    jobs = [(P, list(range(r, n, chunks)), float(qmax2)) for r in range(chunks)]
    if workers > 1:
        with Pool(workers) as pool:
            parts = pool.map(_collect_job, jobs)
    else:
        parts = [_collect_job(job) for job in jobs]
    flagged = [c for part in parts for c in part]
    counts = {"flagged": len(flagged), "empty": 0, "single_point": 0, "segment": 0, "other": 0}
    families, point_survivors, seen = [], [], set()
    for case in flagged:
        res = _feasible_segment(P, case)
        if res is None or isinstance(res, str):
            counts["other"] += 1
            continue
        if not res:
            counts["empty"] += 1
            continue
        counts["single_point" if len(res) == 1 else "segment"] += 1
        i, j, k, l, sigma = case
        a, b, c, d = P[i], P[(i+1) % n], P[j], P[(j+1) % n]
        rot = lambda v: (-sigma*v[1], sigma*v[0])

        def square(st):
            s_, t_ = st
            v0 = (a[0]+s_*(b[0]-a[0]), a[1]+s_*(b[1]-a[1]))
            v1 = (c[0]+t_*(d[0]-c[0]), c[1]+t_*(d[1]-c[1]))
            w = rot((v1[0]-v0[0], v1[1]-v0[1]))
            v3, v2 = (v0[0]+w[0], v0[1]+w[1]), (v1[0]+w[0], v1[1]+w[1])
            return (v0, v1, v2, v3) if sigma == 1 else (v0, v3, v2, v1)
        points = res if len(res) == 1 else [res[0], ((res[0][0]+res[1][0])/2, (res[0][1]+res[1][1])/2), res[1]]
        squares = [square(st) for st in points]
        side2 = [(q[1][0]-q[0][0])**2+(q[1][1]-q[0][1])**2 for q in squares]
        smallest = min(side2)
        if len(res) == 2:
            # The side vector is affine along the segment, so side^2 is a convex quadratic in the
            # segment parameter mu in [0,1]; minimize it exactly (the minimum may be interior).
            w0 = (squares[0][1][0]-squares[0][0][0], squares[0][1][1]-squares[0][0][1])
            w1 = (squares[-1][1][0]-squares[-1][0][0], squares[-1][1][1]-squares[-1][0][1])
            dw = (w1[0]-w0[0], w1[1]-w0[1])
            dd = dw[0]**2+dw[1]**2
            if dd:
                mu = min(Q(1), max(Q(0), -(w0[0]*dw[0]+w0[1]*dw[1])/dd))
                smallest = (w0[0]+mu*dw[0])**2+(w0[1]+mu*dw[1])**2
        if max(side2) == 0 or smallest/2 >= qmax2:
            continue
        key = frozenset(frozenset(q) for q in squares)
        if key in seen:
            continue
        seen.add(key)
        dists = [residue_distance(P, q) if s2 > 0 else None for q, s2 in zip(squares, side2)]
        if len(res) == 1:
            if dists[0] is not None and dists[0] >= gap:
                point_survivors.append((dists[0], squares[0]))
        else:
            families.append({"edges": case, "sides": [float(x)**0.5 for x in side2],
                             "residue_distances": dists})
    return counts, families, point_survivors


def census(N, workers=4, margin=0.05):
    P = midside_fixture(N)
    A, R2 = area(P), radius_squared(P)
    p = pi_enclosure(30)
    gap = min(A/2, A-p[1]*R2/2)
    quads, small, best, degenerate = float_search(P, keep_dist=float(gap)-margin, area_f=float(A), workers=workers)
    found, singular = certify_parallel(P, quads, workers=workers)
    scored = sorted(((d, s) for s, d in found.values()), key=lambda x: (-x[0], [tuple(v) for v in x[1]]))
    survivors = [(d, s) for d, s in scored if d >= gap]
    tagged = (Q(19*N, 24)-Q(1, 384*N)-Q(1, 2))/4      # v324 marked square, eps = 1/2
    return {"N": N, "vertices": len(P), "area": A, "radius_sq": R2, "spectral_gap_lower": gap,
            "marked_square_residue": min(tagged % A, A-tagged % A), "small_square_hits": small,
            "max_float_residue_distance": best, "candidates_certified": len(found),
            "survivors": survivors, "degenerate_flags": degenerate, "singular_systems": singular}


def summary(row, top=5):
    out = {k: (float(v) if isinstance(v, Q) else v) for k, v in row.items() if k != "survivors"}
    out["max_float_residue_distance"] = round(row["max_float_residue_distance"], 6)   # floating search; rounded for replay
    out["survivor_count"] = len(row["survivors"])
    out["marked_square_survives"] = bool(row["marked_square_residue"] >= row["spectral_gap_lower"])
    out["top_survivors"] = [{"residue_distance": float(d),
                             "side": float(((s[1][0]-s[0][0])**2+(s[1][1]-s[0][1])**2))**0.5,
                             "vertices": [[str(x), str(y)] for x, y in s]} for d, s in row["survivors"][:top]]
    return out


def large_runs(workers=4):
    """Exact censuses at N = 24 and 30 and a floating count at the retuned N = 40.

    At N = 40 exact certification of every survivor needs more memory than a 16 GB machine
    has, so only raw floating hits are counted; a square with contacts in edge interiors is
    hit 8 times (4 starting corners x 2 orientations), the ratio measured at N = 30."""
    rows = [summary(census(N, workers=workers)) for N in (24, 30)]
    P = midside_fixture(40)
    A, R2 = area(P), radius_squared(P)
    p = pi_enclosure(30)
    gap = min(A/2, A-p[1]*R2/2)
    hits, small, best, _ = float_search(P, keep_dist=float(gap), area_f=float(A), workers=workers, count_only=True)
    tagged = (Q(19*40, 24)-Q(1, 384*40)-Q(1, 2))/4
    rows.append({"N": 40, "vertices": len(P), "spectral_gap_lower": float(gap),
                 "marked_square_residue": float(min(tagged % A, A-tagged % A)),
                 "small_square_hits": small, "raw_hits_at_or_above_gap": hits,
                 "estimated_survivors": hits//8, "max_float_residue_distance": round(best, 6),
                 "certified": False})
    return rows


def _survivors_of(P, workers=4):
    A, R2 = area(P), radius_squared(P)
    p = pi_enclosure(30)
    gap = min(A/2, A-p[1]*R2/2)
    quads, *_ = float_search(P, keep_dist=float(gap)-0.05, area_f=float(A), workers=workers)
    found, singular = certify_parallel(P, quads, workers=workers)
    return gap, [(d, s) for s, d in found.values() if d >= gap], singular


def perturbed(P, delta, seed):
    """Move every vertex by a seeded rational offset in [-delta, delta]^2."""
    import random
    rng = random.Random(seed)
    return [(x+Q(rng.randint(-1000, 1000), 1000)*delta, y+Q(rng.randint(-1000, 1000), 1000)*delta) for x, y in P]


def chamfered(P, eta):
    """Cut every corner: replace vertex b by the points at max-norm distance eta along both edges."""
    out, n = [], len(P)
    for i in range(n):
        b = P[i]
        for u in (P[i-1], P[(i+1) % n]):
            d = (u[0]-b[0], u[1]-b[1])
            f = eta/max(abs(d[0]), abs(d[1]))
            out.append((b[0]+f*d[0], b[1]+f*d[1]))
    return out


def genericity(N=24, size=Q(1, 10000), radius=2e-3, workers=4):
    """Rerun the census on perturbed and corner-cut versions of the N-turn fixture.

    Survivors are matched by centre within `radius`; since spiral arms are about 1.3e-3 apart,
    a match may be a neighbouring square, so only the survivor counts are matching-free."""
    from polygon_orbit_census_v223 import validate_jordan_polygon, pt
    base = midside_fixture(N)
    centre = lambda sq: (float(sum(v[0] for v in sq)/4), float(sum(v[1] for v in sq)/4))
    _, base_surv, _ = _survivors_of(base, workers)
    base_c = [centre(sq) for _, sq in base_surv]
    near = lambda c, pool: any(abs(c[0]-x) <= radius and abs(c[1]-y) <= radius for x, y in pool)
    rows = [{"polygon": "base", "vertices": len(base), "survivors": len(base_surv)}]
    for label, P in (("perturbed, seed 1", perturbed(base, size, 1)), ("perturbed, seed 2", perturbed(base, size, 2)),
                     ("corners cut", chamfered(base, size))):
        validate_jordan_polygon([pt(*v) for v in P])
        gap, surv, singular = _survivors_of(P, workers)
        c = [centre(sq) for _, sq in surv]
        rows.append({"polygon": label, "vertices": len(P), "spectral_gap_lower": float(gap), "survivors": len(surv),
                     "singular_systems": singular, "base_survivors_matched": sum(near(x, c) for x in base_c),
                     "new_survivors": sum(not near(x, base_c) for x in c)})
    return {"N": N, "displacement": str(size), "match_radius": radius, "rows": rows}


if __name__ == "__main__":
    if sys.argv[1:] == ["--genericity"]:
        print(json.dumps(genericity(), indent=2))
    elif sys.argv[1:] == ["--large"]:
        print(json.dumps(large_runs(), indent=2))
    else:
        turns = [int(x) for x in sys.argv[1:]] or [4, 8, 16, 20]
        print(json.dumps([summary(census(N)) for N in turns], indent=2))
