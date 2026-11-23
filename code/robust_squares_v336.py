"""v336: how large are the robust inscribed squares of grid terrains?

v335 showed that the side-i square of a terrain can rest on corner contacts and vanish under a
perturbation of size 1e-3. Here the robust size is estimated: the largest inscribed square of the
exact polygon and of several random perturbations (every vertex moved by a seeded rational offset of
max-norm at most eps), each perturbed polygon checked to be simple, with squares computed exactly by
the v334 machinery. The minimum over perturbations is an upper estimate of the robust size; an
adversarial perturbation could do worse.
"""
from fractions import Fraction as Q
from multiprocessing import Pool
import json
import random
import sys

from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from terrain_gluing_v334 import column_maps, interior_square, largest_square, polygon, terrain_i


def terrain_polygon(heights):
    return polygon(0, len(heights), {x: h for x, h in enumerate(heights)}, {x: 0 for x in range(len(heights))})


def perturb(P, seed, eps):
    rng = random.Random(seed)
    for _ in range(20):
        Pp = [(x+Q(rng.randint(-1000, 1000), 1000)*eps, y+Q(rng.randint(-1000, 1000), 1000)*eps) for x, y in P]
        try:
            validate_jordan_polygon([pt(*p) for p in Pp])
            return Pp
        except ValueError:
            continue
    raise ArithmeticError("no simple perturbation found")


def robust_estimate(P, seeds, eps=Q(1, 1000)):
    exact = float(largest_square(P)[0])**0.5
    perturbed = [float(largest_square(perturb(P, s, eps))[0])**0.5 for s in seeds]
    return exact, min(perturbed), perturbed


def _terrain_job(seed, trials=12):
    rng = random.Random(seed)
    heights = [rng.randint(1, 8) for _ in range(rng.randint(2, 8))]
    i = terrain_i(heights)
    if i < 2:
        return None
    exact, robust, _ = robust_estimate(terrain_polygon(heights), [1000*seed+k for k in range(trials)])
    return {"heights": heights, "i": i, "M": exact, "robust": robust, "ratio": robust/i, "exact_ratio": exact/i}


def _union_job(seed, trials=12):
    rng = random.Random(seed)
    top = [rng.randint(1, 8) for _ in range(rng.randint(2, 7))]
    bot = [rng.randint(1, 8) for _ in range(rng.randint(1, 6))]
    bs = rng.randint(-len(bot)+1, len(top)-1)
    L, R, f, g = column_maps(top, bot, 0, bs)
    i, _ = interior_square(L, R, f, g)
    if i < 2:
        return None
    exact, robust, _ = robust_estimate(polygon(L, R, f, g), [1000*seed+k for k in range(trials)])
    return {"top": top, "bot": bot, "bot_start": bs, "i": i, "M": exact, "robust": robust, "ratio": robust/i, "exact_ratio": exact/i}


def summarize(rows):
    rows = sorted((r for r in rows if r), key=lambda r: r["ratio"])
    ratios = [r["ratio"] for r in rows]
    q = lambda p: ratios[min(len(ratios)-1, int(p*len(ratios)))]
    return {"shapes": len(rows), "min_robust_ratio": ratios[0], "quantiles_10_50_90": [q(0.1), q(0.5), q(0.9)],
            "fraction_robust_below_exact": sum(r["robust"] < r["M"]-1e-6 for r in rows)/len(rows),
            "min_exact_ratio": min(r["exact_ratio"] for r in rows), "worst": rows[:6]}


def experiment(n_terrains=300, n_unions=200, workers=4):
    with Pool(workers) as pool:
        terrains = pool.map(_terrain_job, range(n_terrains))
        unions = pool.map(_union_job, range(n_unions))
    return {"eps": "1/1000", "perturbations_per_shape": 12,
            "terrains": summarize(terrains), "two_terrain_unions": summarize(unions)}


if __name__ == "__main__":
    print(json.dumps(experiment(), indent=2))


# ---------------------------------------------------------------- the diagonal degree and its corner zeros

def corner_square_zero(heights):
    """The fragile corner square at the left base corner: (0,0), (c,0), (c,c), (0,c) with (c,c) on the
    roof. It is a zero of the base-to-roof diagonal problem on the boundary b = (0,0), which is why that
    degree is undefined on the naive parameter domain. Returns c (the roof crosses y = x at x = c)."""
    for x, h in enumerate(heights):
        lo, hi = Q(x), Q(x+1)
        if lo <= h <= hi:            # horizontal roof piece at height h meets y = x
            return Q(h)
        if h < lo:                   # the roof dropped below the diagonal at the wall x = lo
            return lo
        prev = heights[x+1] if x+1 < len(heights) else 0
        if h > hi and prev <= hi:    # the wall at x+1 descends through height x+1
            return hi
    return Q(len(heights))
