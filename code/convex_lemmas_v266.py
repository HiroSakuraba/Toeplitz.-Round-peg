"""v266: a two-corner proof of the hidden-corner lemma, a stability lemma,
an automatic noncollapse cutoff for moving convex polygons, and a
graceful-or-off-disk lemma.  Proofs are in
research/v266-two-corner-stability-and-convex-cutoff.md; this module gives
exact (and optional numeric) support.

Lemma 0 (two near corners). Square: centre p, half-diagonal h, unit
diagonal directions u, v. Its NEAR corners are A = p - h sgn(p.u) u and
B = p - h sgn(p.v) v (the corner on the disk-centre side of each
diagonal). If neither near corner is hidden behind the unit disk D (i.e. in
the interior of conv(D U its two neighbours)), then h >= 1, with equality
only for p = 0.  Proof: with X = |p.u|, Y = |p.v| and supporting directions
inside the 45-degree cones,
   (h-X) cos f + Y sin f >= 1,   (h-Y) cos g + X sin g >= 1,
and adding gives 2 <= h(cos f + cos g) - X(cos f - sin g) - Y(cos g - sin f)
<= 2h because |f|, |g| <= pi/4.

Lemma 2 (stability). Same hypotheses with h <= 1 + eps, eps < 3 - 2 sqrt 2:
   X + Y <= T(eps) = 2 eps (1 + eps) / (1 - eps - 2 sqrt(eps)).
For convex curves (no corner hidden) and eps <= 1/100 this sharpens to
|p| <= sqrt2 eps (1 + 7 eps); |p| = sqrt2 eps is attained.

Lemma 1 (automatic cutoff). If a moving convex polygon contains a disk of
radius r0 at every time, every inscribed square at every time has
side^2 >= 2 r0^2 (v255). `certify_moving_disk` proves the containment for
v264-style linear motions P(t) = (1-t)P0 + t P1, exactly, with an affine
moving centre and Bernstein sign certificates on dyadic time slabs.

Lemma 3 (graceful or off the disk) is checked numerically on twin spikes.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb, sqrt
from pathlib import Path
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from polygon_orbit_census_v223 import pt, add, sub, mul, dot, cross  # noqa: E402
from big_peg_exact_v243 import read_polygon, maximum_square  # noqa: E402
import convex_big_peg_v255 as v255  # noqa: E402

CHECKS_JSON = ROOT / "research" / "v266-convex-lemmas-checks.json"


# ---------------------------------------------------------------------------
# Lemma 0: two near corners
# ---------------------------------------------------------------------------

def near_corner_indices(square):
    """Indices of the two near corners of a square given in cyclic order
    (p+u, p+Ju, p-u, p-Ju as in v255.square_from)."""
    p = mul(add(square[0], square[2]), Q(1, 2))
    u = sub(square[0], p)
    w = sub(square[1], p)                  # = Ju
    a = 2 if dot(p, u) >= 0 else 0         # p - u if p.u >= 0 else p + u
    b = 3 if dot(p, w) >= 0 else 1
    return a, b


def near_corners_hidden(square, r2=Q(1)):
    p = mul(add(square[0], square[2]), Q(1, 2))
    return [k for k in near_corner_indices(square) if v255.corner_hidden(square[k], p, r2)]


def _grid_squares(centres, us):
    for p in centres:
        for u in us:
            yield p, u, v255.square_from(p, u)


def two_corner_grid_census():
    """Every rational grid square with h < 1 has a hidden NEAR corner
    (no restriction on where the corners are)."""
    us = [pt(Q(a, 8), Q(b, 8)) for a in range(1, 9) for b in range(0, 9)]
    centres = [pt(Q(i, 3), Q(j, 3)) for i in range(-9, 10) for j in range(-9, 10)]
    small = ok = 0
    for p, u, sq in _grid_squares(centres, us):
        if dot(u, u) < 1:
            small += 1
            ok += bool(near_corners_hidden(sq))
    return {"small_squares": small, "with_hidden_near_corner": ok}


# ---------------------------------------------------------------------------
# Lemma 2: stability constants (exact when sqrt(eps) is rational)
# ---------------------------------------------------------------------------

def stability_constants(root_eps):
    """eps = root_eps^2 with root_eps rational. Returns exact c, s, m, T."""
    e = root_eps * root_eps
    c = (1 - e) / (1 + e)                  # lower bound for cos f, cos g
    s = 2 * root_eps / (1 + e)             # upper bound for |sin f|, |sin g|
    assert c * c + s * s == 1
    m = c - s                              # lower bound for cos f - sin g
    if m <= 0:
        raise ValueError("needs eps < 3 - 2 sqrt 2")
    T = 2 * e / m
    return {"eps": e, "cos_lower": c, "sin_upper": s, "margin": m, "T": T}


def _pythagorean_dirs():
    return [pt(1, 0), pt(Q(3, 5), Q(4, 5)), pt(Q(4, 5), Q(3, 5)), pt(Q(5, 13), Q(12, 13)),
            pt(Q(12, 13), Q(5, 13)), pt(Q(8, 17), Q(15, 17)), pt(Q(20, 29), Q(21, 29)),
            pt(Q(7, 25), Q(24, 25)), pt(Q(9, 41), Q(40, 41))]


def stability_grid_census(root_eps=Q(1, 10)):
    """Exact grid check of Lemma 2 (two-corner form) and of the convex form.

    Squares with rational centre p, rational unit diagonal direction and
    rational h in [1, 1 + eps]."""
    K = stability_constants(root_eps)
    e, T = K["eps"], K["T"]
    hs = [1 + e * Q(k, 8) for k in range(0, 9)]
    step = e / 4
    centres = [pt(i * step, j * step) for i in range(-12, 13) for j in range(-12, 13)]
    tested = two_ok = conv = conv_ok = 0
    worst_ratio = Q(0)
    for p in centres:
        for d in _pythagorean_dirs():
            for h in hs:
                u = mul(d, h)
                sq = v255.square_from(p, u)
                if near_corners_hidden(sq):
                    continue
                tested += 1
                X, Y = abs(dot(p, d)), abs(dot(p, pt(-d.y, d.x)))
                two_ok += (X + Y <= T)
                if not v255.hidden_corners(sq, Q(1)):
                    conv += 1
                    pp = dot(p, p)
                    conv_ok += (pp <= 2 * e * e * (1 + 7 * e) ** 2)
                    if e > 0:
                        worst_ratio = max(worst_ratio, pp / (2 * e * e))
    return {"eps": str(e), "T": str(T), "T_over_eps": float(T / e),
            "unhidden_near_corner_squares": tested, "satisfy_X_plus_Y_le_T": two_ok,
            "fully_unhidden_squares": conv, "satisfy_convex_centre_bound": conv_ok,
            "max_centre_sq_over_2eps2": float(worst_ratio)}


def sharpness_example(e=Q(1, 100)):
    """p = (eps, eps), u = (1+eps, 0): no hidden corner, |p| = sqrt2 eps."""
    p = pt(e, e)
    sq = v255.square_from(p, pt(1 + e, 0))
    return {"eps": str(e), "hidden_corners": v255.hidden_corners(sq, Q(1)),
            "corners_outside_open_disk": all(dot(V, V) >= 1 for V in sq),
            "centre_distance_squared_over_eps_squared": str(dot(p, p) / (e * e))}


def pochuev_near_corner():
    P, _, _ = read_polygon(v255.POCHUEV_JSON)
    sq = maximum_square(P).witness
    hid = v255.hidden_corners(sq, Q(1))
    return {"hidden": hid, "near": list(near_corner_indices(sq)),
            "hidden_corner_is_near": all(k in near_corner_indices(sq) for k in hid)}


# ---------------------------------------------------------------------------
# Lemma 1: exact moving-disk certificate for linear convex motions
# ---------------------------------------------------------------------------

def _padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def _pmul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def _pscale(a, c):
    return [c * x for x in a]


def _affine(p0, p1):
    """Coordinates of (1-t) p0 + t p1 as polynomials in t."""
    return ([Q(p0.x), Q(p1.x - p0.x)], [Q(p0.y), Q(p1.y - p0.y)])


def _bernstein(poly, lo, hi):
    """Bernstein coefficients of poly(t) on [lo, hi] (exact)."""
    n = len(poly) - 1
    # substitute t = lo + (hi - lo) s
    sub_ = [Q(0)]
    for k in reversed(range(n + 1)):
        sub_ = _padd(_pmul(sub_, [lo, hi - lo]), [poly[k]])
    sub_ += [Q(0)] * (n + 1 - len(sub_))
    return [sum(Q(comb(j, k), comb(n, k)) * sub_[k] for k in range(j + 1)) for j in range(n + 1)]


def _positive_on(poly, lo, hi, depth):
    """Exact certificate that poly > 0 on [lo, hi] by Bernstein subdivision."""
    b = _bernstein(poly, lo, hi)
    if all(x > 0 for x in b):
        return 1
    if depth == 0 or b[0] <= 0 or b[-1] <= 0:
        return 0
    mid = (lo + hi) / 2
    left = _positive_on(poly, lo, mid, depth - 1)
    if not left:
        return 0
    right = _positive_on(poly, mid, hi, depth - 1)
    return left + right if right else 0


def _chebyshev_centre_float(P):
    """Largest inscribed disk of a convex CCW polygon (float, small LP by
    bisection on the radius with half-plane intersection)."""
    n = len(P)
    lines = []
    for i in range(n):
        a, b = P[i], P[(i + 1) % n]
        ex, ey = b[0] - a[0], b[1] - a[1]
        L = sqrt(ex * ex + ey * ey)
        lines.append((-ey / L, ex / L, (-ey * a[0] + ex * a[1]) / L))   # inner normal . x >= c
    def feasible(r):
        # clip a big box by inner half-planes shifted by r; return a point or None
        poly = [(-1e6, -1e6), (1e6, -1e6), (1e6, 1e6), (-1e6, 1e6)]
        for nx, ny, c in lines:
            out = []
            m = len(poly)
            for k in range(m):
                p, q = poly[k], poly[(k + 1) % m]
                sp, sq = nx * p[0] + ny * p[1] - c - r, nx * q[0] + ny * q[1] - c - r
                if sp >= 0:
                    out.append(p)
                if sp * sq < 0:
                    t = sp / (sp - sq)
                    out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
            poly = out
            if not poly:
                return None
        return (sum(x for x, _ in poly) / len(poly), sum(y for _, y in poly) / len(poly))
    lo, hi = 0.0, 1e3
    best = feasible(0.0)
    for _ in range(80):
        mid = (lo + hi) / 2
        c = feasible(mid)
        if c is None:
            hi = mid
        else:
            lo, best = mid, c
    return best, lo


def certify_moving_disk(P0, P1, c0, c1, r2, lo=Q(0), hi=Q(1), depth=10):
    """Exact proof that the closed disk of squared radius r2 centred at
    (1-t)c0 + t c1 lies in the CCW polygon (1-t)P0 + t P1 for all t in
    [lo, hi]. For each edge line: L(t) = cross(E(t), c(t) - A(t)) > 0 and
    L(t)^2 - r2 |E(t)|^2 > 0. (Convexity for all t must be certified
    separately, e.g. by v264's certify_strict_convex_motion.)"""
    n = len(P0)
    cx, cy = _affine(c0, c1)
    cells = 0
    for i in range(n):
        ax, ay = _affine(P0[i], P1[i])
        bx, by = _affine(P0[(i + 1) % n], P1[(i + 1) % n])
        ex, ey = _padd(bx, _pscale(ax, -1)), _padd(by, _pscale(ay, -1))
        dx, dy = _padd(cx, _pscale(ax, -1)), _padd(cy, _pscale(ay, -1))
        L = _padd(_pmul(ex, dy), _pscale(_pmul(ey, dx), -1))
        E2 = _padd(_pmul(ex, ex), _pmul(ey, ey))
        g = _padd(_pmul(L, L), _pscale(E2, -r2))
        a = _positive_on(L, lo, hi, depth)
        b = _positive_on(g, lo, hi, depth)
        if not (a and b):
            return None
        cells += a + b
    return cells


def auto_disk_cutoff(P0, P1, slabs=1, shrink=Q(999, 1000), depth=10, max_slabs=64):
    """Find and exactly certify a moving disk on [0, 1] (dyadic slabs if
    needed). Returns the certified r0^2 lower bound (min over slabs), the
    cutoff 2 r0^2 for squared square sides, and per-slab data."""
    from convex_square_event_atlas_v264 import certify_strict_convex_motion
    certify_strict_convex_motion(P0, P1)
    while slabs <= max_slabs:
        out = []
        good = True
        for k in range(slabs):
            lo, hi = Q(k, slabs), Q(k + 1, slabs)
            at = lambda t: [(float((1 - t) * a.x + t * b.x), float((1 - t) * a.y + t * b.y)) for a, b in zip(P0, P1)]
            (c0f, r0f), (c1f, r1f) = _chebyshev_centre_float(at(lo)), _chebyshev_centre_float(at(hi))
            c0 = pt(Q(c0f[0]).limit_denominator(10 ** 6), Q(c0f[1]).limit_denominator(10 ** 6))
            c1 = pt(Q(c1f[0]).limit_denominator(10 ** 6), Q(c1f[1]).limit_denominator(10 ** 6))
            # the endpoint centres lie in the convex polygon at the endpoint times; express
            # the slab's centre path in global t
            def endpoint(c_lo, c_hi):
                # centre(t) = c_lo + (t - lo)/(hi - lo) (c_hi - c_lo) = (1-t) C0 + t C1
                d = sub(c_hi, c_lo)
                C1 = add(c_lo, mul(d, (1 - lo) / (hi - lo)))
                C0 = sub(c_lo, mul(d, lo / (hi - lo)))
                return C0, C1
            C0, C1 = endpoint(c0, c1)
            r = Q(min(r0f, r1f)).limit_denominator(10 ** 6) * shrink
            cert = None
            for _ in range(12):
                cert = certify_moving_disk(P0, P1, C0, C1, r * r, lo, hi, depth)
                if cert:
                    break
                r *= Q(9, 10)
            if not cert:
                good = False
                break
            out.append({"slab": [str(lo), str(hi)], "centre_at_0": [str(C0.x), str(C0.y)],
                        "centre_at_1": [str(C1.x), str(C1.y)], "r": str(r), "bernstein_cells": cert})
        if good:
            r2 = min(Q(s["r"]) ** 2 for s in out)
            return {"slabs": out, "r0_squared": r2, "side_squared_cutoff": 2 * r2,
                    "side_squared_cutoff_float": float(2 * r2)}
        slabs *= 2
    raise ValueError("no moving-disk certificate found")


def v264_families():
    """The moving convex families used by v264's tests."""
    from big_peg_pair_death_v254 import family_vertices, ROOT_LEFT, ROOT_RIGHT, START, END
    from big_peg_multichart_v249 import ppoint_at
    F = family_vertices()
    win = (tuple(ppoint_at(v, ROOT_LEFT) for v in F), tuple(ppoint_at(v, ROOT_RIGHT) for v in F))
    full = (tuple(ppoint_at(v, START) for v in F), tuple(ppoint_at(v, END) for v in F))
    quad = tuple(pt(x, y) for x, y in ((0, 0), (4, 0), (3, 2), (-1, 3)))
    return {
        "v254_collision_window": win,
        "v254_full_interval_3/5_to_13/20": full,
        "v264_translated_pentagon": (win[0], tuple(add(v, pt(Q(1, 7), -Q(2, 11))) for v in win[0])),
        "v264_scaled_pentagon": (win[1], tuple(mul(v, Q(3, 2)) for v in win[1])),
        "v264_translated_quadrilateral": (quad, tuple(add(v, pt(Q(1, 10), -Q(1, 7))) for v in quad)),
    }


def cutoff_comparison():
    from convex_square_event_atlas_v264 import certify_uniform_side_floor
    out = {}
    for name, (A, B) in v264_families().items():
        d = auto_disk_cutoff(A, B)
        proj = certify_uniform_side_floor(A, B)
        out[name] = {"disk_side_squared_floor": str(d["side_squared_cutoff"]),
                     "disk_side_squared_floor_float": d["side_squared_cutoff_float"],
                     "v264_projection_side_squared_floor_float": float(proj["every_positive_square_side_squared_lower_bound"]),
                     "improvement_factor": d["side_squared_cutoff_float"] / float(proj["every_positive_square_side_squared_lower_bound"]),
                     "slabs": len(d["slabs"])}
    return out


def exact_report():
    return {"two_corner_grid": two_corner_grid_census(),
            "stability_constants_eps_1_100": {k: str(v) for k, v in stability_constants(Q(1, 10)).items()},
            "stability_grid_eps_1_100": stability_grid_census(Q(1, 10)),
            "stability_grid_eps_1_400": stability_grid_census(Q(1, 20)),
            "sharpness_example": sharpness_example(),
            "pochuev_near_corner": pochuev_near_corner(),
            "cutoff_comparison": cutoff_comparison()}


# ---------------------------------------------------------------------------
# Lemma 3 (numeric): graceful or off the disk, on twin spikes
# ---------------------------------------------------------------------------

def _boundary_param(np, P, q):
    E = np.roll(P, -1, 0) - P
    t = np.clip(((q - P) * E).sum(1) / (E * E).sum(1), 0, 1)
    d = np.linalg.norm(P + t[:, None] * E - q, axis=1)
    i = int(np.argmin(d))
    return i + float(t[i])


def sector_check():
    import numpy as np
    out = {}
    for gap in (0.3, 0.6, 1.0, 1.6):
        P = v255.twin_spike(gap)
        rows = []
        for V, h in v255._square_polys(np, P):
            p = V.mean(axis=0)
            cr = lambda a, b: a[0] * b[1] - a[1] * b[0]
            d1 = abs(cr(V[2] - V[0], -V[0])) / np.linalg.norm(V[2] - V[0])
            d2 = abs(cr(V[3] - V[1], -V[1])) / np.linalg.norm(V[3] - V[1])
            s = [_boundary_param(np, P, q) for q in V]
            order = np.argsort(s)
            graceful = any(list(np.roll(order, k)) in ([0, 1, 2, 3], [3, 2, 1, 0]) for k in range(4))
            rows.append({"half_diagonal": float(h), "diag_line_distances": [float(d1), float(d2)],
                         # tolerance: some diagonals lie along a spike edge TANGENT to the circle
                         "a_diagonal_crosses_open_disk": bool(min(d1, d2) < 1 - 1e-9),
                         "a_diagonal_tangent_to_disk": bool(abs(min(d1, d2) - 1) <= 1e-9),
                         "graceful": bool(graceful)})
        out[str(gap)] = {
            "squares": rows,
            "small_squares_all_ungraceful": all(not r["graceful"] for r in rows if r["half_diagonal"] < 1),
            "lemma_holds": all((not r["a_diagonal_crosses_open_disk"]) or (r["graceful"] and r["half_diagonal"] >= 1)
                               for r in rows),
            "small_or_ungraceful_have_both_diagonals_off_disk": all(
                not r["a_diagonal_crosses_open_disk"] for r in rows if r["half_diagonal"] < 1 or not r["graceful"]),
        }
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--numeric", action="store_true")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args(argv)
    rep = {"exact": exact_report()}
    ex = rep["exact"]
    g = ex["two_corner_grid"]
    assert g["with_hidden_near_corner"] == g["small_squares"] > 10000
    for key in ("stability_grid_eps_1_100", "stability_grid_eps_1_400"):
        s = ex[key]
        assert s["satisfy_X_plus_Y_le_T"] == s["unhidden_near_corner_squares"] > 0
        assert s["satisfy_convex_centre_bound"] == s["fully_unhidden_squares"] > 0
    assert ex["sharpness_example"]["hidden_corners"] == []
    assert ex["pochuev_near_corner"]["hidden_corner_is_near"]
    if a.numeric:
        rep["numeric"] = {"twin_spike_sectors": sector_check()}
        assert all(v["lemma_holds"] for v in rep["numeric"]["twin_spike_sectors"].values())
    if a.write:
        CHECKS_JSON.write_text(json.dumps(rep, indent=1, default=str) + "\n")
    print(json.dumps(rep, indent=1, default=str)[:5000])
    return rep


if __name__ == "__main__":
    main()
