"""v341: the incircle reading of 1/sqrt2 and a test of the disk version.

Incircle reading. The largest interior axis-aligned square of side i has an incircle of radius i/2;
the square inscribed in that circle has side i/sqrt2, the constant of Conjecture 336.1.

Disk version (Conjecture 341.1). With r the inradius (largest inscribed disk), every Jordan curve
carries a robust inscribed square of side at least sqrt2 * r. Since a disk of radius r contains an
axis-aligned square of side sqrt2 * r, and every interior square of side i contains a disk of radius
i/2, we have i/sqrt2 <= sqrt2 * r <= i, so the disk version implies Conjecture 336.1 and is implied
on the v338 interior case by Theorem 338.1.

The test compares the v336 robust estimate (minimum largest square over seeded perturbations, an
upper estimate of the robust size) with sqrt2 * r, where r is computed numerically (grid search
refined by a pattern search; the distance is exact for each point).
"""
from fractions import Fraction as Q
from multiprocessing import Pool
import json
import math
import random
import sys

import numpy as np

from robust_squares_v336 import robust_estimate, terrain_polygon
from terrain_gluing_v334 import column_maps, interior_square, polygon, terrain_i


# ---------------------------------------------------------------- inradius

def _dist_inside(V, pts):
    """Distance to the boundary for points inside the polygon V (float array), -1 outside."""
    A, B = V, np.roll(V, -1, axis=0)
    d = B-A
    rel = pts[:, None, :]-A[None]
    t = np.clip((rel*d[None]).sum(-1)/(d*d).sum(-1)[None], 0, 1)
    dist = np.sqrt(((rel-t[..., None]*d[None])**2).sum(-1)).min(1)
    x, y = pts[:, :1], pts[:, 1:]
    cross = ((A[None, :, 1] > y) != (B[None, :, 1] > y))
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = A[None, :, 0]+(y-A[None, :, 1])*(B[None, :, 0]-A[None, :, 0])/(B[None, :, 1]-A[None, :, 1])
    inside = (cross & (x < xi)).sum(1) % 2 == 1
    return np.where(inside, dist, -1.0)


def inradius(P, grid=160):
    """(r, centre) of the largest disk inside the polygon P, to about 1e-9 relative accuracy."""
    V = np.array([[float(x), float(y)] for x, y in P])
    lo, hi = V.min(0), V.max(0)
    h = float((hi-lo).max())/grid
    xs = np.arange(lo[0]+h/2, hi[0], h)
    ys = np.arange(lo[1]+h/2, hi[1], h)
    G = np.array([(x, y) for x in xs for y in ys])
    D = _dist_inside(V, G)
    best = (-1.0, None)
    for k in np.argsort(-D)[:12]:
        c, r, step = G[k].copy(), D[k], h
        while step > 1e-10*h*grid:
            moves = c+step*np.array([[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [1, -1], [-1, 1], [-1, -1]])
            dm = _dist_inside(V, moves)
            j = int(np.argmax(dm))
            if dm[j] > r:
                c, r = moves[j], dm[j]
            else:
                step /= 2
        if r > best[0]:
            best = (float(r), (float(c[0]), float(c[1])))
    return best


# ---------------------------------------------------------------- shapes

def union_polygon(top, bot, bs):
    L, R, f, g = column_maps(top, bot, 0, bs)
    return polygon(L, R, f, g), interior_square(L, R, f, g)[0]


def rational_polygon(points, den=1000):
    return [(Q(round(x*den), den), Q(round(y*den), den)) for x, y in points]


def regular_polygon(n, R=4.0, phase=0.1):
    return rational_polygon([(R*math.cos(phase+2*math.pi*k/n), R*math.sin(phase+2*math.pi*k/n)) for k in range(n)])


def ellipse_polygon(a, b, n=24, phase=0.05):
    return rational_polygon([(a*math.cos(phase+2*math.pi*k/n), b*math.sin(phase+2*math.pi*k/n)) for k in range(n)])


def bulb_on_stick(n=16, R=2.0, stick=6.0, w=0.6):
    """A regular n-gon of circumradius R with a thin rectangular stick of width w attached on the right."""
    pts = []
    y0 = w/2
    a0 = math.asin(y0/R)
    for k in range(n+1):
        t = a0+(2*math.pi-2*a0)*k/n
        pts.append((R*math.cos(t), R*math.sin(t)))
    pts += [(R*math.cos(a0)+stick, -y0), (R*math.cos(a0)+stick, y0)]
    return rational_polygon(pts)


def notched_square(s=4, notch=1, depth=3):
    """An s x s square with a slot of width `notch` cut from the top, `depth` deep."""
    a = (s-notch)/2
    return rational_polygon([(0, 0), (s, 0), (s, s), (a+notch, s), (a+notch, s-depth), (a, s-depth), (a, s), (0, s)])


def random_star(seed, n=20):
    """A star-shaped polygon: n vertices at jittered angles with radii uniform in [2, 4]."""
    rng = random.Random(seed)
    angles = sorted(2*math.pi*(k+rng.uniform(-0.35, 0.35))/n for k in range(n))
    radii = [rng.uniform(2, 4) for _ in angles]
    return rational_polygon([(q*math.cos(t), q*math.sin(t)) for q, t in zip(radii, angles)])


def slotted_circle(seed, n=32, R=4.0):
    """A regular n-gon of circumradius R with 3 to 6 thin inward slots at random positions: each slot
    replaces one edge's midpoint by a notch of width about R/10 and random depth up to R/2."""
    rng = random.Random(seed)
    slots = set(rng.sample(range(n), rng.randint(3, 6)))
    pts = []
    for k in range(n):
        t0, t1 = 2*math.pi*k/n, 2*math.pi*(k+1)/n
        p0 = (R*math.cos(t0), R*math.sin(t0))
        pts.append(p0)
        if k in slots:
            depth = rng.uniform(0.2, 0.5)*R
            tm, half = (t0+t1)/2, (t1-t0)/4
            for t, rad in ((tm-half, R*math.cos(half)), (tm-half, R-depth), (tm+half, R-depth), (tm+half, R*math.cos(half))):
                pts.append((rad*math.cos(t), rad*math.sin(t)))
    return rational_polygon(pts)


def named_shapes():
    out = {"terrain (2,7)": terrain_polygon([2, 7]), "v334 glued witness A":
           [(Q(x), Q(y)) for x, y in ((0, 0), (6, 0), (6, -1), (8, -1), (8, 3), (7, 3), (7, 7), (3, 7), (3, 8), (0, 8))]}
    for n in (8, 12, 16, 20):
        out[f"regular {n}-gon"] = regular_polygon(n)
    for a, b in ((3, 2), (4, 1.5), (5, 1)):
        out[f"24-gon ellipse {a}x{b}"] = ellipse_polygon(a, b)
    out["16-gon bulb on a stick"] = bulb_on_stick()
    out["square with slot 1 x 3"] = notched_square()
    out["square with slot 2 x 3"] = notched_square(notch=2)
    return out


# ---------------------------------------------------------------- experiment

def row(P, i=None, seeds=range(12)):
    r, c = inradius(P)
    exact, robust, _ = robust_estimate(P, list(seeds))
    out = {"r": r, "centre": c, "M": exact, "robust": robust, "disk_ratio": robust/(math.sqrt(2)*r)}
    if i is not None:
        out.update({"i": i, "i_over_r": i/r})
    return out


def _terrain_job(seed):
    rng = random.Random(seed)
    heights = [rng.randint(1, 8) for _ in range(rng.randint(2, 8))]
    i = terrain_i(heights)
    if i < 2:
        return None
    return {"heights": heights, **row(terrain_polygon(heights), i, [1000*seed+k for k in range(12)])}


def _union_job(seed):
    rng = random.Random(seed)
    top = [rng.randint(1, 8) for _ in range(rng.randint(2, 7))]
    bot = [rng.randint(1, 8) for _ in range(rng.randint(1, 6))]
    bs = rng.randint(-len(bot)+1, len(top)-1)
    P, i = union_polygon(top, bot, bs)
    if i < 2:
        return None
    return {"top": top, "bot": bot, "bot_start": bs, **row(P, i, [1000*seed+k for k in range(12)])}


def _named_job(item):
    name, P = item
    return {"shape": name, **row(P)}


def _star_job(seed):
    try:
        return {"star_seed": seed, **row(random_star(seed), seeds=[1000*seed+k for k in range(12)])}
    except (ArithmeticError, ValueError):     # the sampled polygon is not simple
        return None


def _slot_job(seed):
    return {"slot_seed": seed, **row(slotted_circle(seed), seeds=[1000*seed+k for k in range(12)])}


def summarize(rows):
    rows = sorted((r for r in rows if r), key=lambda r: r["disk_ratio"])
    return {"shapes": len(rows), "min_disk_ratio": rows[0]["disk_ratio"],
            "below_one_minus_1e-2": sum(r["disk_ratio"] < 0.99 for r in rows),
            "i_over_r_range": [min(r["i_over_r"] for r in rows), max(r["i_over_r"] for r in rows)] if "i" in rows[0] else None,
            "worst": rows[:6]}


def experiment(n_terrains=300, n_unions=200, n_stars=60, n_slots=40, workers=4):
    with Pool(workers) as pool:
        named = pool.map(_named_job, list(named_shapes().items()))
        terrains = pool.map(_terrain_job, range(n_terrains))
        unions = pool.map(_union_job, range(n_unions))
        stars = pool.map(_star_job, range(n_stars))
        slots = pool.map(_slot_job, range(n_slots))
    return {"eps": "1/1000", "perturbations_per_shape": 12, "named": sorted(named, key=lambda r: r["disk_ratio"]),
            "terrains": summarize(terrains), "two_terrain_unions": summarize(unions),
            "random_stars": summarize(stars), "slotted_circles": summarize(slots)}


if __name__ == "__main__":
    print(json.dumps(experiment(), indent=2))
