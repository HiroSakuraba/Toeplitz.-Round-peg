"""v255: the sharp sqrt(2) Big Peg bound holds for convex curves.

Theorem (proved in research/v255-convex-sharp-big-peg.md). If a compact
convex set K contains a closed disk of radius r, then EVERY square with all
four vertices on the boundary of K has side >= sqrt(2) r, with equality only
for a square inscribed in the disk whose four vertices are contact points.

It follows from a curve-free lemma about one square and one disk:

    Hidden-corner lemma. If no vertex of a square lies in the open disk D of
    radius r and the side is < sqrt(2) r, then some vertex lies in the
    interior of conv(D U {other three vertices}).

This module supplies computational support, not the proof.

Exact part (standard library, fractions.Fraction):
  * an exact decision procedure for "corner hidden" (closed form, below);
  * an exact grid census of the lemma on rational squares;
  * exact complete square censuses (v223) of rational convex polygons
    containing the unit disk: every square has side^2 >= 2;
  * Pochuev's ten-gon (v244 data): its unique small square has both
    diagonal lines missing the disk and exactly one hidden corner, which
    sits on the edge leaving the reflex vertex v2; the polygon's kernel
    misses the unit disk (exact rational kernel).

Numeric part (optional: numpy, scipy, cma; run with --numeric):
  random convex polygons, an optimizer trying to break the theorem,
  twin-spike regions, and searches over polygons star-shaped about the
  disk (kernel contains D) or about its centre only.

Closed form used for the exact hidden test. For a corner V of a square with
centre p, put w = V - p (|w| = h). The supporting directions of the square
at V form the closed cone of unit n within 45 degrees of w. V is NOT hidden
iff some n in that cone has n.V >= r. The largest n.V over the cone is |V|
when V is within 45 degrees of w, and otherwise (w.V + |Jw.V|)/(sqrt2 h).
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import argparse
import json
import math
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from polygon_orbit_census_v223 import (  # noqa: E402
    pt, sub, add, mul, dot, cross, validate_jordan_polygon, census_polygon,
)
from big_peg_exact_v243 import (  # noqa: E402
    read_polygon, maximum_square, contained_disk, _distance_sq_to_segment,
)

POCHUEV_JSON = ROOT / "research" / "v244-pochuev-original-ten-gon.json"
CHECKS_JSON = ROOT / "research" / "v255-convex-big-peg-checks.json"
FIGURE = ROOT / "research" / "figures" / "v255-convex-big-peg.png"


# ---------------------------------------------------------------------------
# Exact hidden-corner decision
# ---------------------------------------------------------------------------

def J(v):
    return pt(-v.y, v.x)


def square_from(p, u):
    """Square p+u, p+Ju, p-u, p-Ju (cyclic order). p, u exact Points."""
    return (add(p, u), add(p, J(u)), sub(p, u), sub(p, J(u)))


def corner_hidden(V, p, r2):
    """Exact: is corner V of the square with centre p strictly inside
    conv(D U other corners), where D is the closed disk |x|^2 <= r2 at 0?"""
    w = sub(V, p)
    h2, wv, vv = dot(w, w), dot(w, V), dot(V, V)
    if wv >= 0 and 2 * wv * wv >= vv * h2:      # V within 45 degrees of w: max n.V = |V|
        return vv < r2
    s = wv + abs(dot(J(w), V))                  # = sqrt(2) * h * (max n.V over the cone)
    return not (s >= 0 and s * s >= 2 * h2 * r2)


def hidden_corners(square, r2):
    p = mul(add(square[0], square[2]), Q(1, 2))
    return [k for k, V in enumerate(square) if corner_hidden(V, p, r2)]


def triangle_certificate(square, k, r2=Q(1), grid=10):
    """Independent exact witness that corner k is hidden: a rational z with
    |z|^2 < r2 and two other corners X, Y such that corner k lies strictly
    inside triangle (z, X, Y). Returns None if the grid search finds none."""
    H = square[k]
    others = [square[j] for j in range(4) if j != k]
    for X, Y in ((others[0], others[1]), (others[0], others[2]), (others[1], others[2])):
        for i, j in product(range(-grid, grid + 1), repeat=2):
            z = pt(Q(i, grid), Q(j, grid))
            if dot(z, z) >= r2:
                continue
            d = cross(sub(X, z), sub(Y, z))
            if d == 0:
                continue
            l1 = cross(sub(X, H), sub(Y, H)) / d
            l2 = cross(sub(Y, H), sub(z, H)) / d
            l3 = 1 - l1 - l2
            if l1 > 0 and l2 > 0 and l3 > 0:
                return {"z": [str(z.x), str(z.y)], "barycentric": [str(l1), str(l2), str(l3)]}
    return None


def lemma_grid_census(r2=Q(1)):
    """Exact lemma check on a rational grid of squares.

    Every square with no corner in the open disk and side^2 < 2 r2 must have a
    hidden corner. Squares with side^2 == 2 r2 and no hidden corner must be
    centred at the origin (equality case)."""
    small = small_hidden = equal = equal_free = equal_free_off_centre = 0
    us = [pt(Q(a, 8), Q(b, 8)) for a in range(1, 9) for b in range(0, 9)]
    unit_us = [pt(1, 0), pt(Q(3, 5), Q(4, 5)), pt(Q(4, 5), Q(3, 5)),
               pt(Q(5, 13), Q(12, 13)), pt(Q(12, 13), Q(5, 13))]
    centres = [pt(Q(i, 3), Q(j, 3)) for i in range(-9, 10) for j in range(-9, 10)]
    for p in centres:
        for u in us + unit_us:
            sq = square_from(p, u)
            if any(dot(V, V) < r2 for V in sq):
                continue
            h2 = dot(u, u)
            if h2 < r2:
                small += 1
                small_hidden += bool(hidden_corners(sq, r2))
            elif h2 == r2:
                equal += 1
                if not hidden_corners(sq, r2):
                    equal_free += 1
                    equal_free_off_centre += (p.x != 0 or p.y != 0)
    return {"small_squares": small, "small_with_hidden_corner": small_hidden,
            "equality_squares": equal, "equality_without_hidden_corner": equal_free,
            "equality_free_off_centre": equal_free_off_centre}


def lemma_float_sweep(n=100000, seed=255):
    """Seeded float Monte Carlo of the lemma with the same closed form."""
    rng = random.Random(seed)
    tested = failures = 0
    while tested < n:
        px, py = (rng.uniform(-4, 4) * rng.choice((0.05, 0.3, 1.0)) for _ in range(2))
        h = rng.uniform(0.01, 0.999999) if rng.random() < 0.7 else rng.uniform(0.95, 0.999999)
        a = rng.uniform(0, math.pi / 2)
        ux, uy = h * math.cos(a), h * math.sin(a)
        V = [(px + ux, py + uy), (px - uy, py + ux), (px - ux, py - uy), (px + uy, py - ux)]
        if any(x * x + y * y < 1 for x, y in V):
            continue
        tested += 1
        hidden = False
        for (vx, vy) in V:
            wx, wy = vx - px, vy - py
            vv = math.hypot(vx, vy)
            ang = math.atan2(wx * vy - wy * vx, wx * vx + wy * vy)   # angle from w to V
            best = vv if abs(ang) <= math.pi / 4 else vv * math.cos(abs(ang) - math.pi / 4)
            if best < 1 - 1e-12:
                hidden = True
                break
        failures += not hidden
    return {"tested": tested, "without_hidden_corner": failures}


# ---------------------------------------------------------------------------
# Exact convex censuses
# ---------------------------------------------------------------------------

def _turns(P):
    n = len(P)
    return [cross(sub(P[i], P[i - 1]), sub(P[(i + 1) % n], P[i])) for i in range(n)]


def is_strictly_convex_ccw(P):
    return all(t > 0 for t in _turns(P))


def convex_test_polygons():
    """Fixed rational convex polygons containing the closed unit disk."""
    polys = {
        "triangle": [(-3, -1), (3, -1), (0, 4)],
        "asymmetric_quad": [(0, -3), (Q(5, 2), 0), (Q(1, 2), 2), (-2, Q(1, 4))],
        "lopsided_pentagon": [(-1, -2), (4, -1), (3, 2), (-1, Q(5, 4)), (-2, 0)],
        "near_square_hexagon": [(Q(-11, 10), -1), (1, Q(-11, 10)), (Q(11, 10), 1),
                                (Q(-1, 5), Q(13, 10)), (Q(-6, 5), Q(6, 5)), (Q(-13, 10), 0)],
        "long_heptagon": [(-6, -1), (-1, Q(-13, 10)), (5, -1), (6, 0), (5, Q(11, 10)),
                          (-1, Q(13, 10)), (-6, 1)],
        # tangent-ish edges at 0, 90.3, 180.2, 270.1 degrees (support 1.0002):
        # nearly the equality configuration, smallest square side^2 ~ 2.001
        "octagon_near_equality": [(Q(5001, 5000), Q(1141, 2500)), (Q(227, 500), Q(5013, 5000)),
                                  (Q(-4873, 10000), Q(9977, 10000)), (Q(-10019, 10000), Q(483, 1000)),
                                  (Q(-4993, 5000), Q(-2361, 5000)), (Q(-2349, 5000), Q(-1001, 1000)),
                                  (Q(4997, 10000), Q(-9993, 10000)), (Q(5001, 5000), Q(-4989, 10000))],
    }
    rng = random.Random(2550)
    k = 0
    while k < 4:
        n = rng.randint(4, 6)
        ang = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
        gaps = [(ang[(i + 1) % n] - ang[i]) % (2 * math.pi) for i in range(n)]
        if max(gaps) >= 0.95 * math.pi:
            continue
        hs = [1 + rng.choice((0.05, 0.4, 2.0)) * rng.random() for _ in range(n)]
        V = []
        for i in range(n):
            a1, a2 = ang[i], ang[(i + 1) % n]
            det = math.sin(a2 - a1)
            x = (hs[i] * math.sin(a2) - hs[(i + 1) % n] * math.sin(a1)) / det
            y = (hs[(i + 1) % n] * math.cos(a1) - hs[i] * math.cos(a2)) / det
            V.append((Q(round(x * 200), 200), Q(round(y * 200), 200)))
        P = tuple(pt(*v) for v in V)
        try:
            validate_jordan_polygon(P)
        except Exception:
            continue
        if is_strictly_convex_ccw(P) and contained_disk(P, pt(0, 0), Q(1)):
            polys[f"random_convex_{k}"] = V
            k += 1
    return {name: tuple(pt(*v) for v in V) for name, V in polys.items()}


def exact_convex_census():
    out = {}
    for name, P in convex_test_polygons().items():
        assert is_strictly_convex_ccw(P), name
        assert contained_disk(P, pt(0, 0), Q(1)), name
        c = census_polygon(P)
        sides = [s.side_squared for s in c.squares]
        out[name] = {"vertices": len(P), "squares": len(sides),
                     "min_side_squared": str(min(sides)), "min_side_squared_float": float(min(sides)),
                     "max_side_squared_float": float(max(sides)),
                     "all_side_squared_at_least_2": all(s >= 2 for s in sides)}
    return out


# ---------------------------------------------------------------------------
# Pochuev's ten-gon (non-convex counterexample to sharp sqrt 2)
# ---------------------------------------------------------------------------

def _clip(poly, a, b):
    """Keep the part of a convex polygon left of (or on) the directed line a->b."""
    out = []
    n = len(poly)
    side = lambda x: cross(sub(b, a), sub(x, a))
    for i in range(n):
        cur, nxt = poly[i], poly[(i + 1) % n]
        sc, sn = side(cur), side(nxt)
        if sc >= 0:
            out.append(cur)
        if (sc > 0 and sn < 0) or (sc < 0 and sn > 0):
            t = sc / (sc - sn)
            out.append(add(cur, mul(sub(nxt, cur), t)))
    return out


def exact_kernel(P):
    """Exact kernel (set of points seeing the whole polygon) of a CCW polygon."""
    K = [pt(-1000, -1000), pt(1000, -1000), pt(1000, 1000), pt(-1000, 1000)]
    n = len(P)
    for i in range(n):
        K = _clip(K, P[i], P[(i + 1) % n])
        if not K:
            return []
    return K


def pochuev_certificate():
    P, _, _ = read_polygon(POCHUEV_JSON)
    if sum(cross(P[i], P[(i + 1) % len(P)]) for i in range(len(P))) < 0:
        raise AssertionError("expected counter-clockwise data")
    n = len(P)
    turns = _turns(P)
    reflex = [i for i, t in enumerate(turns) if t < 0]
    star_about_centre = all(cross(P[i], P[(i + 1) % n]) > 0 for i in range(n))
    assert contained_disk(P, pt(0, 0), Q(1))
    best = maximum_square(P)
    sq = best.witness                                   # p, q, r, s in cyclic order
    side2 = best.side_squared
    p = mul(add(sq[0], sq[2]), Q(1, 2))
    # squared distances from the disk centre to the two diagonal lines
    def line_d2(a, b):
        e = sub(b, a)
        return cross(a, e) ** 2 / dot(e, e)
    diag_d2 = [line_d2(sq[0], sq[2]), line_d2(sq[1], sq[3])]
    hid = hidden_corners(sq, Q(1))
    assert len(hid) == 1
    Hc = sq[hid[0]]
    # which edge carries the hidden corner, and how far from vertex v2
    on = [i for i in range(n)
          if cross(sub(P[(i + 1) % n], P[i]), sub(Hc, P[i])) == 0
          and 0 <= dot(sub(Hc, P[i]), sub(P[(i + 1) % n], P[i])) <= dot(sub(P[(i + 1) % n], P[i]), sub(P[(i + 1) % n], P[i]))]
    d2_v2 = dot(sub(Hc, P[2]), sub(Hc, P[2]))
    # independent barycentric certificate: hidden corner strictly inside a triangle
    # (z, X, Y) with z rational in the OPEN unit disk and X, Y other corners
    bary = triangle_certificate(sq, hid[0])
    assert bary is not None
    K = exact_kernel(P)
    kd2 = min(_distance_sq_to_segment(pt(0, 0), K[i], K[(i + 1) % len(K)]) for i in range(len(K)))
    return {
        "max_side_squared": str(side2), "max_side_squared_float": float(side2),
        "side_squared_below_2": side2 < 2,
        "half_diagonal_squared_float": float(side2 / 2),
        "centre": [float(p.x), float(p.y)], "centre_distance_float": math.hypot(float(p.x), float(p.y)),
        "diagonal_line_distance_squared_float": [float(d) for d in diag_d2],
        "both_diagonal_lines_miss_unit_disk": all(d > 1 for d in diag_d2),
        "hidden_corner_indices": hid,
        "hidden_corner": [str(Hc.x), str(Hc.y)], "hidden_corner_float": [float(Hc.x), float(Hc.y)],
        "hidden_corner_on_edges": on,
        "hidden_corner_distance_to_v2_float": math.sqrt(float(d2_v2)),
        "barycentric_certificate": bary,
        "reflex_vertices": reflex,
        "star_shaped_about_disk_centre": star_about_centre,
        "kernel_vertices_float": [[float(v.x), float(v.y)] for v in K],
        "kernel_distance_squared": str(kd2), "kernel_distance_float": math.sqrt(float(kd2)),
        "kernel_misses_unit_disk": kd2 > 1,
    }


def exact_report():
    return {"lemma_grid": lemma_grid_census(),
            "lemma_float_sweep": lemma_float_sweep(),
            "convex_census": exact_convex_census(),
            "pochuev": pochuev_certificate()}


# ---------------------------------------------------------------------------
# Numeric part (numpy / scipy / cma)
# ---------------------------------------------------------------------------

def _np():
    import numpy as np
    return np


def all_squares_float(P, tol=1e-9, return_u=False):
    """All inscribed squares of a float polygon via the 4x4 corner/edge-line
    systems (generic assignments; singular systems skipped).
    Returns (centres, half_diagonals)."""
    np = _np()
    P = np.asarray(P, float)
    n = len(P)
    E = np.roll(P, -1, 0) - P
    N = np.stack([E[:, 1], -E[:, 0]], 1)
    b = np.einsum('ij,ij->i', N, P)
    Jm = np.array([[0.0, -1.0], [1.0, 0.0]])
    corner = [(1, np.eye(2)), (1, Jm), (-1, np.eye(2)), (-1, Jm)]
    idx = np.array(list(product(range(n), repeat=4)))
    M = np.zeros((len(idx), 4, 4))
    rhs = np.zeros((len(idx), 4))
    for k, (s, Mk) in enumerate(corner):
        nk = N[idx[:, k]]
        M[:, k, 0:2] = nk
        M[:, k, 2:4] = s * nk @ Mk
        rhs[:, k] = b[idx[:, k]]
    scale = np.abs(M).max(axis=(1, 2)) ** 4
    ok = np.abs(np.linalg.det(M)) > 1e-12 * scale
    sol = np.full((len(idx), 4), np.nan)
    sol[ok] = np.linalg.solve(M[ok], rhs[ok][..., None])[..., 0]
    c, u = sol[:, :2], sol[:, 2:]
    good = ok & (np.linalg.norm(u, axis=1) > 1e-7)
    for k, (s, Mk) in enumerate(corner):
        q = c + s * u @ Mk.T
        e = idx[:, k]
        t = np.einsum('ij,ij->i', q - P[e], E[e]) / np.einsum('ij,ij->i', E[e], E[e])
        good &= (t > -tol) & (t < 1 + tol)
    if return_u:
        return c[good], u[good]
    return c[good], np.linalg.norm(u[good], axis=1)


def inradius_float(P):
    np = _np()
    from scipy.optimize import linprog
    P = np.asarray(P, float)
    E = np.roll(P, -1, 0) - P
    N = np.stack([E[:, 1], -E[:, 0]], 1)
    N = N / np.linalg.norm(N, axis=1)[:, None]
    res = linprog([0, 0, -1], A_ub=np.c_[N, np.ones(len(P))], b_ub=np.einsum('ij,ij->i', N, P),
                  bounds=[(None, None)] * 3)
    return -res.fun


def _support_polygon(np, ang, hs):
    n = len(ang)
    V = []
    for i in range(n):
        a1, a2 = ang[i], ang[(i + 1) % n]
        V.append(np.linalg.solve([[np.cos(a1), np.sin(a1)], [np.cos(a2), np.sin(a2)]], [hs[i], hs[(i + 1) % n]]))
    V = np.array(V)
    a = np.roll(V, -1, 0) - V
    bb = np.roll(a, -1, 0)
    if np.any(a[:, 0] * bb[:, 1] - a[:, 1] * bb[:, 0] <= 1e-9):
        return None
    return V


def random_convex_check(count=3000, seed=0):
    """Every inscribed square of random convex polygons with all support
    values >= 1: report min half-diagonal / inradius (theorem: >= 1)."""
    np = _np()
    rng = np.random.default_rng(seed)
    done = nsq = 0
    worst = math.inf
    while done < count:
        n = int(rng.integers(3, 13))
        ang = np.sort(rng.uniform(0, 2 * np.pi, n))
        if np.max(np.diff(np.r_[ang, ang[0] + 2 * np.pi])) >= np.pi:
            continue
        hs = 1 + rng.exponential(rng.choice([0.05, 0.5, 3.0]), n)
        P = _support_polygon(np, ang, hs)
        if P is None:
            continue
        c, h = all_squares_float(P)
        if len(h) == 0:
            continue
        done += 1
        nsq += len(h)
        worst = min(worst, h.min() / inradius_float(P))
    return {"polygons": done, "squares": int(nsq), "min_half_diagonal_over_inradius": float(worst)}


def convex_adversary(n, seed, restarts=4, maxfevals=2500):
    """CMA-ES minimising the SMALLEST inscribed half-diagonal over convex
    n-gons containing the unit disk (support values >= 1)."""
    np = _np()
    import cma
    rng = np.random.default_rng(seed)

    def f(x):
        g = np.exp(x[:n])
        g = 2 * np.pi * g / g.sum()
        if g.max() >= 0.98 * np.pi:
            return 5.0
        P = _support_polygon(np, np.cumsum(g), 1 + np.exp(x[n:]))
        if P is None:
            return 5.0
        c, h = all_squares_float(P)
        return float(h.min()) if len(h) else 5.0
    best = math.inf
    for rep in range(restarts):
        es = cma.CMAEvolutionStrategy(np.r_[rng.normal(0, .3, n), rng.normal(-1, 1, n)], 0.6,
                                      {'seed': 100 * n + rep + 1, 'maxfevals': maxfevals, 'verbose': -9,
                                       'tolflatfitness': 50})
        es.optimize(f)
        best = min(best, es.result.fbest)
    return float(best)


def twin_spike(gap, L=6.0, m=24):
    """Union of conv(D U T1) and conv(D U T2), T = (L, +-gap/2): star-shaped
    with respect to every point of the unit disk D."""
    np = _np()
    T1, T2 = np.array([L, gap / 2]), np.array([L, -gap / 2])

    def tang(T):
        R = np.linalg.norm(T)
        b = np.arctan2(T[1], T[0])
        d = np.arccos(1 / R)
        return b + d, b - d
    a1u, a1l = tang(T1)
    a2u, a2l = tang(T2)
    p1 = np.array([np.cos(a1l), np.sin(a1l)])
    p2 = np.array([np.cos(a2u), np.sin(a2u)])
    s, _ = np.linalg.solve(np.c_[T1 - p1, -(T2 - p2)], p2 - p1)
    notch = p1 + s * (T1 - p1)
    arc = np.linspace(a1u, a2l + 2 * np.pi, m)
    return np.array([T1] + [np.array([np.cos(a), np.sin(a)]) for a in arc] + [T2, notch])


def twin_spike_check():
    np = _np()
    out = {}
    for gap in (0.3, 0.6, 1.0, 1.6):
        P = twin_spike(gap)
        c, h = all_squares_float(P)
        out[str(gap)] = {"squares": int(len(h)), "min_half_diagonal": float(h.min()),
                         "max_half_diagonal": float(h.max()),
                         "small_squares_far_from_disk": int(np.sum((h < 1) & (np.linalg.norm(c, axis=1) > 2)))}
    return out


def _polar_polygon(np, x, n):
    g = np.exp(x[:n])
    g = 2 * np.pi * g / g.sum()
    phi = np.cumsum(g)
    R = 1.0 + np.exp(x[n:2 * n])
    return np.c_[R * np.cos(phi), R * np.sin(phi)], g


def _edge_line_dist(np, P):
    Qp = np.roll(P, -1, 0)
    return (P[:, 0] * Qp[:, 1] - P[:, 1] * Qp[:, 0]) / np.linalg.norm(Qp - P, axis=1)


def _seg_dist(np, P):
    Qp = np.roll(P, -1, 0)
    d = Qp - P
    t = np.clip(-(P * d).sum(1) / (d * d).sum(1), 0, 1)
    return np.linalg.norm(P + t[:, None] * d, axis=1)


def star_search(n, seed, mode, restarts=3, maxfevals=3000):
    """CMA-ES minimising the LARGEST inscribed half-diagonal over polygons
    star-shaped about the origin (vertices in increasing polar angle) with
      mode='kernel': every edge LINE at distance >= 1 (unit disk in kernel);
      mode='centre': every edge SEGMENT at distance >= 1 (unit disk inside)."""
    np = _np()
    import cma
    rng = np.random.default_rng(seed)
    dist = (lambda P: _edge_line_dist(np, P)) if mode == 'kernel' else (lambda P: _seg_dist(np, P))

    def f(x):
        P, g = _polar_polygon(np, x, n)
        pen = 50 * np.sum(np.maximum(0, 1 - dist(P))) + 50 * np.sum(np.maximum(0, g - 0.95 * np.pi))
        if pen > 0:
            return 2 + float(pen)
        c, h = all_squares_float(P)
        return float(h.max()) if len(h) else 5.0
    best = (math.inf, None)
    for rep in range(restarts):
        x0 = np.r_[rng.normal(0, 0.3, n), np.log(rng.uniform(0.2, 2.0, n))]
        es = cma.CMAEvolutionStrategy(x0, 0.5, {'seed': seed * 10 + rep + 1, 'maxfevals': maxfevals,
                                               'verbose': -9, 'tolflatfitness': 50})
        es.optimize(f)
        if es.result.fbest < best[0]:
            best = (float(es.result.fbest), es.result.xbest)
    P, _ = _polar_polygon(np, best[1], n)
    return {"n": n, "seed": seed, "mode": mode, "best_max_half_diagonal": best[0],
            "min_edge_line_distance": float(_edge_line_dist(np, P).min()),
            "min_segment_distance": float(_seg_dist(np, P).min()), "polygon": P.tolist()}


def _simple_float(P):
    n = len(P)

    def o(p, q, r):
        v = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
        return int(v > 0) - int(v < 0)
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            a, b, c, d = P[i], P[(i + 1) % n], P[j], P[(j + 1) % n]
            if o(a, b, c) != o(a, b, d) and o(c, d, a) != o(c, d, b):
                return False
    return True


def pochuev_push(seed, lams=(0.3, 1, 3, 10, 30, 100), maxfevals=2500):
    """Start at Pochuev's ten-gon and penalise failure to be star-shaped about
    the disk centre while keeping the unit disk inside; minimise the largest
    inscribed half-diagonal."""
    np = _np()
    import cma
    P0, _, _ = read_polygon(POCHUEV_JSON)
    x = np.array([[float(v.x), float(v.y)] for v in P0]).ravel()

    def viol(P):
        Qp = np.roll(P, -1, 0)
        return float(np.sum(np.maximum(0, -(P[:, 0] * Qp[:, 1] - P[:, 1] * Qp[:, 0]))))

    def f(y, lam):
        P = y.reshape(-1, 2)
        if not _simple_float(P):
            return 10.0
        Qp = np.roll(P, -1, 0)
        star = np.sum(np.maximum(0, 1e-3 - (P[:, 0] * Qp[:, 1] - P[:, 1] * Qp[:, 0])))
        pen = 50 * np.sum(np.maximum(0, 1.0 - _seg_dist(np, P))) + lam * star
        c, h = all_squares_float(P)
        return (float(h.max()) if len(h) else 5.0) + float(pen)
    trail = [{"lambda": 0, "max_half_diagonal": float(all_squares_float(x.reshape(-1, 2))[1].max()),
              "star_violation": viol(x.reshape(-1, 2))}]
    for lam in lams:
        es = cma.CMAEvolutionStrategy(x, 0.03, {'seed': seed + int(lam * 10) + 1, 'maxfevals': maxfevals,
                                               'verbose': -9, 'tolflatfitness': 50})
        es.optimize(lambda y: f(y, lam))
        x = es.result.xbest
        P = x.reshape(-1, 2)
        trail.append({"lambda": lam, "max_half_diagonal": float(all_squares_float(P)[1].max()),
                      "star_violation": viol(P), "min_segment_distance": float(_seg_dist(np, P).min())})
    return {"seed": seed, "trail": trail}


def _job(args):
    kind, a, b = args
    if kind == 'adversary':
        return ('adversary', a, convex_adversary(a, a))
    if kind in ('kernel', 'centre'):
        return (kind, (a, b), star_search(a, b, kind, maxfevals=3000 if kind == 'kernel' else 4000))
    if kind == 'push':
        return ('push', a, pochuev_push(a))
    raise ValueError(kind)


def numeric_report(processes=None):
    from multiprocessing import Pool
    jobs = ([('adversary', n, None) for n in (4, 5, 6, 8)]
            + [('kernel', n, s) for n in (5, 6, 7, 8, 10) for s in (1, 2)]
            + [('centre', n, s) for n in (6, 8, 10, 12) for s in (1, 2, 3)]
            + [('push', s, None) for s in (0, 7)])
    with Pool(processes) as pool:
        results = pool.map(_job, jobs)
    adv = {str(a): v for k, a, v in results if k == 'adversary'}
    ker = [v for k, a, v in results if k == 'kernel']
    cen = [v for k, a, v in results if k == 'centre']
    push = [v for k, a, v in results if k == 'push']
    strip = lambda rs: [{kk: vv for kk, vv in r.items() if kk != 'polygon'} for r in rs]
    return {
        "random_convex": random_convex_check(),
        "convex_adversary_min_half_diagonal": adv,
        "twin_spike": twin_spike_check(),
        "kernel_contains_disk_search": {"runs": strip(ker),
                                        "overall_min_max_half_diagonal": min(r["best_max_half_diagonal"] for r in ker)},
        "star_about_centre_search": {"runs": strip(cen),
                                     "overall_min_max_half_diagonal": min(r["best_max_half_diagonal"] for r in cen)},
        "pochuev_push_to_centre_star": push,
    }


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(path=FIGURE):
    np = _np()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MPoly, Circle

    fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
    t = np.linspace(0, 2 * np.pi, 400)

    # (a) convex polygon: every inscribed square is big
    P = np.array([[float(v.x), float(v.y)] for v in convex_test_polygons()["octagon_near_equality"]])
    ax = axs[0]
    ax.add_patch(MPoly(P, closed=True, fill=False, lw=2, color="k"))
    ax.plot(np.cos(t), np.sin(t), color="tab:blue", lw=1.5)
    polys = _square_polys(np, P)
    for k, (V, hh) in enumerate(polys):
        ax.add_patch(MPoly(V, closed=True, fill=False, lw=1.2, color=plt.cm.viridis(k / max(1, len(polys)))))
    ax.set_title(f"(a) convex octagon near the equality case: its {len(polys)}\n"
                 f"inscribed squares have half-diagonals {', '.join(f'{h:.4f}' for _, h in sorted(polys, key=lambda t: t[1]))}",
                 fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)

    # (b) Pochuev: the small square's hidden corner
    rep = pochuev_certificate()
    P0, _, _ = read_polygon(POCHUEV_JSON)
    P = np.array([[float(v.x), float(v.y)] for v in P0])
    best = maximum_square(P0)
    S = np.array([[float(v.x), float(v.y)] for v in best.witness])
    k = rep["hidden_corner_indices"][0]
    others = np.delete(S, k, 0)
    pts = np.vstack([np.c_[np.cos(t), np.sin(t)], others])
    from scipy.spatial import ConvexHull
    hull = pts[ConvexHull(pts).vertices]
    ax = axs[1]
    ax.add_patch(MPoly(hull, closed=True, color="tab:orange", alpha=0.18, lw=0))
    ax.add_patch(MPoly(P, closed=True, fill=False, lw=2, color="k"))
    ax.plot(np.cos(t), np.sin(t), color="tab:blue", lw=1.5)
    ax.add_patch(MPoly(S, closed=True, fill=False, lw=1.8, color="tab:red"))
    ax.plot(*S[k], "o", color="tab:red", ms=7)
    ax.plot(*P[2], "kx", ms=8)
    ax.annotate("hidden corner\n(next to reflex v₂)", S[k], xytext=(0.2, 1.9), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_title("(b) Pochuev's ten-gon: its only square (side 1.4104)\n"
                 "has a corner inside hull(disk + other corners)", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-3.4, 3.8)
    ax.set_ylim(-1.7, 2.6)

    # (c) twin spike: star-shaped about the whole disk, tiny squares exist
    P = twin_spike(0.6)
    ax = axs[2]
    ax.add_patch(MPoly(P, closed=True, fill=False, lw=2, color="k"))
    ax.plot(np.cos(t), np.sin(t), color="tab:blue", lw=1.5)
    polys = _square_polys(np, P)
    hs = [hh for _, hh in polys]
    for V, hh in (polys[int(np.argmin(hs))], polys[int(np.argmax(hs))]):
        ax.add_patch(MPoly(V, closed=True, fill=False, lw=1.6, color="tab:red" if hh < 1 else "tab:green"))
    ax.set_title(f"(c) twin spike (every disk point sees the boundary):\n"
                 f"smallest half-diagonal {min(hs):.2f}, largest {max(hs):.2f}", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-1.4, 6.4)
    ax.set_ylim(-2.2, 2.2)
    for ax in axs:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    return path


def _square_polys(np, P):
    """Geometric squares (deduplicated by vertex set) as (vertex array, half-diagonal)."""
    c, u = all_squares_float(P, return_u=True)
    out = []
    for ci, ui in zip(c, u):
        Ju = np.array([-ui[1], ui[0]])
        V = np.array([ci + ui, ci + Ju, ci - ui, ci - Ju])
        key = np.sort(np.round(V, 7).view([('x', float), ('y', float)]).ravel())
        if not any(np.array_equal(key, k) for k, _, _ in out):
            out.append((key, V, float(np.linalg.norm(ui))))
    return [(V, h) for _, V, h in out]


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--numeric", action="store_true", help="also run the numpy/scipy/cma campaign")
    ap.add_argument("--write", action="store_true", help="write research/v255-convex-big-peg-checks.json")
    ap.add_argument("--figure", action="store_true", help="write the v255 figure")
    a = ap.parse_args(argv)
    report = {"exact": exact_report()}
    ex = report["exact"]
    assert ex["lemma_grid"]["small_with_hidden_corner"] == ex["lemma_grid"]["small_squares"]
    assert ex["lemma_grid"]["equality_free_off_centre"] == 0
    assert ex["lemma_float_sweep"]["without_hidden_corner"] == 0
    assert all(v["all_side_squared_at_least_2"] for v in ex["convex_census"].values())
    po = ex["pochuev"]
    assert po["side_squared_below_2"] and po["both_diagonal_lines_miss_unit_disk"]
    assert po["hidden_corner_on_edges"] == [2] and po["reflex_vertices"] == [2, 9]
    assert po["kernel_misses_unit_disk"] and not po["star_shaped_about_disk_centre"]
    if a.numeric:
        report["numeric"] = numeric_report()
    if a.write:
        CHECKS_JSON.write_text(json.dumps(report, indent=1) + "\n")
    if a.figure:
        make_figure()
    print(json.dumps(report, indent=1)[:6000])
    return report


if __name__ == "__main__":
    main()
