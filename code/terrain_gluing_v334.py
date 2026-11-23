"""v334: inscribed squares on two grid terrains glued along a common base segment.

An upper terrain (column heights above y = 0) and a lower terrain (column depths below) whose base
intervals overlap form a vertically convex polyomino P. The overlap, the seam, lies in the interior.

* i(P): largest side of an open axis-aligned square in the interior, exactly
  max over column intervals I of min(|I|, min_I f - max_I g); it may straddle the seam.
* M(P): largest side of a square with all four vertices on the boundary, computed exactly over
  every inscribed square: isolated squares from the v331 search with exact certification, plus the
  endpoints of every positive-dimensional family (side^2 is convex along a family).
* For each terrain alone, the v316 base-supported square and whether gluing pushes its base corners
  into the open seam.

The searches are seeded and reproducible but heuristic; a lower ratio may exist.
"""
from fractions import Fraction as Q
from multiprocessing import Pool
import json
import random
import sys

from residue_census_v331 import _collect_job, _feasible_segment, certify, float_search


# ---------------------------------------------------------------- geometry

def column_maps(top, bot, top_start, bot_start):
    """Upper heights f >= 0 and lower depths g <= 0 per unit column."""
    if not max(top_start, bot_start) < min(top_start+len(top), bot_start+len(bot)):
        raise ValueError("the two bases must overlap in at least one column")
    if min(list(top)+list(bot)) < 1:
        raise ValueError("heights and depths must be positive integers")
    L, R = min(top_start, bot_start), max(top_start+len(top), bot_start+len(bot))
    f = {x: (top[x-top_start] if top_start <= x < top_start+len(top) else 0) for x in range(L, R)}
    g = {x: (-bot[x-bot_start] if bot_start <= x < bot_start+len(bot) else 0) for x in range(L, R)}
    return L, R, f, g


def polygon(L, R, f, g):
    """Counterclockwise boundary with collinear vertices removed."""
    pts = []
    for x in range(L, R):
        pts += [(x, g[x]), (x+1, g[x])]
    for x in range(R-1, L-1, -1):
        pts += [(x+1, f[x]), (x, f[x])]
    out = []
    for p in pts:
        if not out or out[-1] != p:
            out.append(p)
    if out[0] == out[-1]:
        out.pop()
    changed = True
    while changed:
        changed = False
        for k in range(len(out)):
            a, b, c = out[k-1], out[k], out[(k+1) % len(out)]
            if (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]) == 0:
                out.pop(k)
                changed = True
                break
    return [(Q(x), Q(y)) for x, y in out]


def interior_square(L, R, f, g):
    """(side, (first column, bottom)) of a largest open axis-aligned interior square."""
    best, where = 0, None
    for a in range(L, R):
        lo, hi = -10**12, 10**12
        for b in range(a, R):
            lo, hi = max(lo, g[b]), min(hi, f[b])
            s = min(b-a+1, hi-lo)
            if s > best:
                best, where = s, (a, lo)
    return best, where


def inscribed_squares(P):
    """Exact (side^2, counterclockwise vertices) for every isolated square and both ends of every
    positive-dimensional family of squares with vertices on the polygon."""
    quads, *_ = float_search(P, qmax2=1e12)
    found, _ = certify(P, quads, qmax2=Q(10**12))
    out = [((s[1][0]-s[0][0])**2+(s[1][1]-s[0][1])**2, s) for s, *_ in found.values()]
    n = len(P)
    for case in _collect_job((P, list(range(n)), 1e12)):
        res = _feasible_segment(P, case)
        if not res or isinstance(res, str):
            continue
        i, j, _, _, sigma = case
        a, b, c, d = P[i], P[(i+1) % n], P[j], P[(j+1) % n]
        for s_, t_ in res:
            v0 = (a[0]+s_*(b[0]-a[0]), a[1]+s_*(b[1]-a[1]))
            v1 = (c[0]+t_*(d[0]-c[0]), c[1]+t_*(d[1]-c[1]))
            w = (-sigma*(v1[1]-v0[1]), sigma*(v1[0]-v0[0]))
            v3, v2 = (v0[0]+w[0], v0[1]+w[1]), (v1[0]+w[0], v1[1]+w[1])
            side2 = (v1[0]-v0[0])**2+(v1[1]-v0[1])**2
            if side2:
                out.append((side2, (v0, v1, v2, v3) if sigma == 1 else (v0, v3, v2, v1)))
    return out


def largest_square(P):
    squares = inscribed_squares(P)
    return max(squares, key=lambda t: t[0]) if squares else (Q(0), None)


def edges_through(P, p):
    n = len(P)
    out = []
    for m in range(n):
        a, b = P[m], P[(m+1) % n]
        d = (b[0]-a[0], b[1]-a[1])
        if d[0]*(p[1]-a[1])-d[1]*(p[0]-a[0]) == 0 and 0 <= (p[0]-a[0])*d[0]+(p[1]-a[1])*d[1] <= d[0]**2+d[1]**2:
            out.append(m)
    return out


def base_square(heights, start):
    """v316: (side, left abscissa) of a largest lattice square with lower side on the base."""
    W = len(heights)
    H = lambda x: heights[x-start] if start <= x < start+W else 0
    roof = lambda X: (min(H(X-1), H(X)), max(H(X-1), H(X)))
    best = None
    for x in range(start, start+W+1):
        for s in range(1, start+W-x+1):
            r1, r2 = roof(x), roof(x+s)
            if r1[0] <= s <= r1[1] and r2[0] <= s <= r2[1] and (best is None or s > best[0]):
                best = (s, x)
    return best


def terrain_i(heights):
    return interior_square(0, len(heights), {x: h for x, h in enumerate(heights)},
                           {x: 0 for x in range(len(heights))})[0]


def analyse(top, bot, bot_start, top_start=0):
    L, R, f, g = column_maps(top, bot, top_start, bot_start)
    P = polygon(L, R, f, g)
    i, where = interior_square(L, R, f, g)
    side2, sq = largest_square(P)
    report = {"top": list(top), "bot": list(bot), "top_start": top_start, "bot_start": bot_start,
              "polygon": [[int(x), int(y)] for x, y in P], "i": i, "interior_square_at": list(where),
              "i_upper": terrain_i(top), "i_lower": terrain_i(bot),
              "M_squared": str(side2), "M_over_i": float(side2)**0.5/i,
              "largest_square": [[str(x), str(y)] for x, y in sq] if sq else None,
              "supporting_edges": [[[int(c) for c in P[m]] + [int(c) for c in P[(m+1) % len(P)]]
                                    for m in edges_through(P, v)] for v in sq] if sq else None,
              "base_squares": {}}
    for name, heights, start, sign in (("upper", top, top_start, 1), ("lower", bot, bot_start, -1)):
        found = base_square(heights, start)
        if found:
            s, x = found
            corners = [(Q(x), Q(0)), (Q(x+s), Q(0)), (Q(x), Q(sign*s)), (Q(x+s), Q(sign*s))]
            report["base_squares"][name] = {"side": s, "left": x, "corners_lost_to_seam":
                                            [[int(c[0]), int(c[1])] for c in corners if not edges_through(P, c)]}
    return report


# ---------------------------------------------------------------- seeded searches

def _ratio(state):
    top, bot, bs = state
    try:
        L, R, f, g = column_maps(top, bot, 0, bs)
    except ValueError:
        return None
    i, _ = interior_square(L, R, f, g)
    if i < 2:
        return None
    side2, _ = largest_square(polygon(L, R, f, g))
    return float(side2)**0.5/i


def _random_job(seed, hmax=6, wmax=7):
    rng = random.Random(seed)
    top = tuple(rng.randint(1, hmax) for _ in range(rng.randint(1, wmax)))
    bot = tuple(rng.randint(1, hmax) for _ in range(rng.randint(1, wmax)))
    bs = rng.randint(-len(bot)+1, len(top)-1)
    return (_ratio((top, bot, bs)), (top, bot, bs))


def _climb_job(seed, steps=700, hmax=12, wmax=12):
    rng = random.Random(seed)
    if seed == 0:
        state = ((6, 6, 6, 6, 5, 5), (5, 1, 1), 3)
    else:
        top = tuple(rng.randint(1, hmax) for _ in range(rng.randint(3, 9)))
        bot = tuple(rng.randint(1, hmax) for _ in range(rng.randint(2, 7)))
        state = (top, bot, rng.randint(-len(bot)+1, len(top)-1))
    cur = _ratio(state) or 9.0
    best = (cur, state)
    for _ in range(steps):
        top, bot, bs = list(state[0]), list(state[1]), state[2]
        k = rng.randrange(7)
        if k == 0:
            j = rng.randrange(len(top))
            top[j] = max(1, min(hmax, top[j]+rng.choice((-2, -1, 1, 2))))
        elif k == 1:
            j = rng.randrange(len(bot))
            bot[j] = max(1, min(hmax, bot[j]+rng.choice((-2, -1, 1, 2))))
        elif k == 2 and len(top) < wmax:
            top.insert(rng.randrange(len(top)+1), rng.randint(1, hmax))
        elif k == 3 and len(bot) < wmax:
            bot.insert(rng.randrange(len(bot)+1), rng.randint(1, hmax))
        elif k == 4 and len(top) > 1:
            top.pop(rng.randrange(len(top)))
        elif k == 5 and len(bot) > 1:
            bot.pop(rng.randrange(len(bot)))
        else:
            bs += rng.choice((-1, 1))
        nb = (tuple(top), tuple(bot), bs)
        r = _ratio(nb) or 9.0
        if r <= cur or rng.random() < 0.02:
            state, cur = nb, r
            if r < best[0]:
                best = (r, nb)
    return best


def searches(workers=4):
    with Pool(workers) as pool:
        sampled = [r for r in pool.map(_random_job, range(4000)) if r[0] is not None]
        climbs = sorted(pool.map(_climb_job, range(24)))
    sampled.sort()
    return {"random_pairs": len(sampled), "random_min_ratio": sampled[0][0],
            "random_below_one": sum(r < 1 for r, _ in sampled),
            "random_worst": [list(map(list, s[:2]))+[s[2]] for _, s in sampled[:3]],
            "climb_best": [[r, list(map(list, s[:2]))+[s[2]]] for r, s in climbs[:6]]}


WITNESSES = {"ratio_2_over_sqrt5": ([6, 6, 6, 6, 5, 5], [5, 1, 1], 3),
             "ratio_5_over_7": ([8, 8, 8, 7, 7, 7, 7, 3], [1, 1], 6)}


if __name__ == "__main__":
    if sys.argv[1:] == ["--search"]:
        print(json.dumps(searches(), indent=2))
    else:
        print(json.dumps({name: analyse(*w) for name, w in WITNESSES.items()}, indent=2))
