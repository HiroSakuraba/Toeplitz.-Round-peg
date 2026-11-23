"""v337: toward Conjecture 336.1 for single terrains.

Two tools.

1. Interior excursions (Theorem 337.1). For a grid terrain with interior side i, choose a largest
   interior square [a, a+i] x [0, i]. The v316 excursion at height h just below i leaves the roof at
   the wall x = l, where l - 1 is the last column left of the square with height < i, and returns at
   x = r, where r is the first column right of it with height < i. If no such column exists, l = 0
   or r = W and the excursion starts or ends on an end wall. When 0 < l and r < W, v316's bent-base
   theorem applies to every registered perturbation of size eps < delta, with L the perturbed base
   and h = i - 3 delta, and it gives a square of side at least i - 4 delta. `interior_case` finds such
   a square position; `bent_base_check` verifies the hypotheses exactly on a perturbed polygon.

2. Cutoff parity. For a generic polygon, N_sigma is the number of inscribed squares of side > sigma.
   If sigma is not the side of any square of the exact terrain, N_sigma mod 2 is the same for all
   generic polygons close enough to it, so N_sigma odd forces a robust square of side > sigma.
   `exact_side_ranges` and `side_gaps` find the admissible sigma; `certified_cutoff` reads the parity
   at each gap from a few perturbations; `parity_profile` records the sigma where N_sigma is odd.
"""
from fractions import Fraction as Q
from multiprocessing import Pool
import json
import random
import sys

from residue_census_v331 import _collect_job, _feasible_segment, certify, float_search
from terrain_gluing_v334 import inscribed_squares, terrain_i
from robust_squares_v336 import perturb, terrain_polygon


# ---------------------------------------------------------------- interior excursions

def excursion_walls(heights, a):
    """Walls (l, r) where the v316 excursion of the interior square starting at column a leaves and
    returns, for h just below i."""
    i = terrain_i(heights)
    low = [x for x in range(a) if heights[x] < i]
    high = [x for x in range(a+i, len(heights)) if heights[x] < i]
    return (low[-1]+1 if low else 0), (high[0] if high else len(heights))


def largest_square_positions(heights):
    i = terrain_i(heights)
    return [a for a in range(len(heights)-i+1) if min(heights[a:a+i]) >= i]


def interior_case(heights):
    """A position a of a largest interior square whose excursion avoids both end walls, or None."""
    W = len(heights)
    for a in largest_square_positions(heights):
        l, r = excursion_walls(heights, a)
        if 0 < l and r < W:
            return a
    return None


def _level_crossings(p, q, h):
    """Parameter u in [0, 1] where the segment p -> q meets height h (None if it does not, or lies on it)."""
    if p[1] == q[1]:
        return None
    u = (h-p[1])/(q[1]-p[1])
    return u if 0 <= u <= 1 else None


def bent_base_check(heights, Pp, delta):
    """Exact check of the v316 bent-base hypotheses on the perturbed polygon Pp (registered to
    terrain_polygon(heights)), with h = i - 3 delta and L the perturbed base edge.

    Returns (side_floor, A_x_range, low_band_intervals) or raises AssertionError."""
    P = terrain_polygon(heights)
    assert len(P) == len(Pp)
    i, W = terrain_i(heights), len(heights)
    a = interior_case(heights)
    assert a is not None, "excursion meets an end wall"
    l, r = excursion_walls(heights, a)
    h = i-3*delta
    n = len(P)
    # the base is the edge (0,0) -> (W,0); the roof runs counterclockwise from (W,0) back to (0,0)
    k0 = P.index((Q(0), Q(0)))
    assert P[(k0+1) % n] == (Q(W), Q(0))
    Lp = (Pp[k0], Pp[(k0+1) % n])
    assert all(abs(v[1]) <= delta for v in Lp)
    roof = [Pp[(k0+1+m) % n] for m in range(n)]           # (W,0) ... (0,0) on the perturbed curve
    exact = [P[(k0+1+m) % n] for m in range(n)]
    # registered walls x = r (entered first, going right to left) and x = l
    def wall_index(x):
        for m in range(len(exact)-1):
            p, q = exact[m], exact[m+1]
            if p[0] == q[0] == x and min(p[1], q[1]) < i <= max(p[1], q[1]):
                return m
        raise AssertionError("wall not found")
    mr, ml = wall_index(r), wall_index(l)
    ur, ul = _level_crossings(roof[mr], roof[mr+1], h), _level_crossings(roof[ml], roof[ml+1], h)
    assert ur is not None and ul is not None
    lerp = lambda p, q, u: (p[0]+u*(q[0]-p[0]), p[1]+u*(q[1]-p[1]))
    start, end = lerp(roof[mr], roof[mr+1], ur), lerp(roof[ml], roof[ml+1], ul)
    A = [start]+roof[mr+1:ml+1]+[end]                       # right to left; v316 orients either way
    assert all(v[1] > h for v in A[1:-1]), "excursion dips to height h"
    # low band h <= y <= h + 2 delta: only on the first and last segments (vertices are above it)
    assert all(v[1] > h+2*delta for v in A[1:-1])
    def band_x(p, q):
        us = [u for u in (_level_crossings(p, q, h), _level_crossings(p, q, h+2*delta)) if u is not None]
        us += [u for u, v in ((0, p), (1, q)) if h <= v[1] <= h+2*delta]
        xs = [lerp(p, q, u)[0] for u in us]
        return min(xs), max(xs)
    right_band, left_band = band_x(A[0], A[1]), band_x(A[-2], A[-1])
    assert right_band[1]-right_band[0] < h-delta and left_band[1]-left_band[0] < h-delta
    assert right_band[0]-left_band[1] > h+delta
    xs = [v[0] for v in A]
    Lx = sorted(v[0] for v in Lp)
    assert Lx[0] < min(xs)-2*delta and Lx[1] > max(xs)+2*delta
    return h-delta, (min(xs), max(xs)), (left_band, right_band)


# ---------------------------------------------------------------- cutoff parity

def exact_side_ranges(P):
    """Closed intervals [lo^2, hi^2] (exact squared sides) covering every inscribed square of P:
    a point for each isolated square, and for each positive-dimensional family the exact range of
    its squared side (the side vector is affine along the family, so the minimum may be interior)."""
    quads, *_ = float_search(P, qmax2=1e12)
    found, _ = certify(P, quads, qmax2=Q(10**12))
    out = []
    for sq, *_ in found.values():
        s2 = (sq[1][0]-sq[0][0])**2+(sq[1][1]-sq[0][1])**2
        out.append((s2, s2))
    n = len(P)
    for case in _collect_job((P, list(range(n)), 1e12)):
        res = _feasible_segment(P, case)
        if not res or isinstance(res, str):
            continue
        i, j = case[0], case[1]
        a, b, c, d = P[i], P[(i+1) % n], P[j], P[(j+1) % n]
        vecs = []
        for s_, t_ in res:
            v0 = (a[0]+s_*(b[0]-a[0]), a[1]+s_*(b[1]-a[1]))
            v1 = (c[0]+t_*(d[0]-c[0]), c[1]+t_*(d[1]-c[1]))
            vecs.append((v1[0]-v0[0], v1[1]-v0[1]))
        sq = [u[0]**2+u[1]**2 for u in vecs]
        lo, hi = min(sq), max(sq)
        if len(vecs) == 2:
            Da, Db = vecs
            E = (Db[0]-Da[0], Db[1]-Da[1])
            ee = E[0]**2+E[1]**2
            if ee:
                lam = -(Da[0]*E[0]+Da[1]*E[1])/ee
                if 0 < lam < 1:
                    lo = min(lo, (Da[0]+lam*E[0])**2+(Da[1]+lam*E[1])**2)
        if hi:
            out.append((lo, hi))
    return out


def side_gaps(P, min_width=0.02):
    """Open intervals of sides containing no exact square side of P (width >= min_width), as floats."""
    rngs = sorted((float(lo)**0.5, float(hi)**0.5) for lo, hi in exact_side_ranges(P))
    merged = []
    for lo, hi in rngs:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    ends = [0.0]+[x for iv in merged for x in iv]
    gaps = [(ends[k], ends[k+1]) for k in range(0, len(ends)-1, 2)]
    return [g for g in gaps if g[1]-g[0] >= min_width], merged


def certified_cutoff(P, lists):
    """Parity of N_sigma at the midpoint of every gap, for each perturbation's side list. Returns
    (rho, gaps): rho is the upper end of the highest gap whose midpoint has odd parity for every
    perturbation (a robust square has side at least rho, given the cutoff parity lemma), and gaps
    lists (gap, parities)."""
    gaps, _ = side_gaps(P)
    table = []
    rho = 0.0
    for lo, hi in gaps:
        mid = (lo+hi)/2
        par = [sum(x > mid for x in s) % 2 for s in lists]
        table.append(((lo, hi), par))
        if all(par):
            rho = max(rho, hi)
    return rho, table


def square_sides(P):
    """Sides of the inscribed squares of P, one entry per square (descending)."""
    seen = {}
    for s2, sq in inscribed_squares(P):
        seen[frozenset(sq)] = float(s2)**0.5
    return sorted(seen.values(), reverse=True)


def odd_set(sides):
    """The sigma > 0 with an odd number of sides above sigma: [s2, s1) u [s4, s3) u ..."""
    s = list(sides)+[0.0]
    return [(s[k+1], s[k]) for k in range(0, len(s)-1, 2)]


def intersect(sets):
    out = sets[0]
    for other in sets[1:]:
        out = [(max(a, c), min(b, d)) for a, b in out for c, d in other if max(a, c) < min(b, d)]
    return out


def parity_profile(heights, seeds, eps=Q(1, 1000)):
    P = terrain_polygon(heights)
    lists = [square_sides(perturb(P, s, eps)) for s in seeds]
    common = intersect([odd_set(s) for s in lists])
    return {"counts": [len(s) for s in lists], "robust": min(s[0] for s in lists), "odd_everywhere": common}


# ---------------------------------------------------------------- experiment

def _job(seed, trials=6):
    rng = random.Random(seed)
    heights = [rng.randint(1, 8) for _ in range(rng.randint(2, 8))]
    i = terrain_i(heights)
    if i < 2:
        return None
    P = terrain_polygon(heights)
    lists = [square_sides(perturb(P, 1000*seed+k, Q(1, 1000))) for k in range(trials)]
    rho, table = certified_cutoff(P, lists)
    return {"heights": heights, "i": i, "case": "interior" if interior_case(heights) is not None else "end wall",
            "robust_ratio": min(s[0] for s in lists)/i, "total_counts_odd": all(len(s) % 2 for s in lists),
            "certified_ratio": rho/i, "parity_consistent": all(len(set(par)) == 1 for _, par in table),
            "gaps": len(table)}


def experiment(n=300, workers=4):
    with Pool(workers) as pool:
        rows = [r for r in pool.map(_job, range(n)) if r]
    out = {"eps": "1/1000", "perturbations_per_shape": 6, "shapes": len(rows)}
    for case in ("interior", "end wall"):
        sub = [r for r in rows if r["case"] == case]
        out[case] = {"shapes": len(sub), "min_robust_ratio": min(r["robust_ratio"] for r in sub),
                     "all_total_counts_odd": all(r["total_counts_odd"] for r in sub),
                     "parity_consistent_across_perturbations": sum(r["parity_consistent"] for r in sub),
                     "min_certified_ratio": min(r["certified_ratio"] for r in sub),
                     "certified_below_one_over_sqrt2": [r for r in sub if r["certified_ratio"] < 2**-0.5-1e-9],
                     "worst_certified": sorted(sub, key=lambda r: r["certified_ratio"])[:5]}
    return out


if __name__ == "__main__":
    print(json.dumps(experiment(), indent=2))
