"""v335: a counterexample to the localized half-square proposal, and the fragility of v316 squares.

Part A. For x, y on the curve put S1 = g(x) + J(g(y) - g(x)), S2 = g(y) + J(g(y) - g(x)) and
h(x, y) = (f(S1), f(S2)), with f the signed distance to the curve (negative inside) and J the
quarter turn. The proposal claimed that if h has no zeros on the band y - x in [a, 1 - a], every
core loop there has odd winding. On the unit circle with 1/4 < a < 1/2 the band is zero-free, yet
the antipodal loop has h constant = (sqrt5 - 1, sqrt5 - 1), so winding zero. The same holds on a
generic ellipse. Here h is computed with exact closed forms on the circle and with dense sampling of
the signed distance on the ellipse; windings are integers and are rounded only after checking.

Part B. Exact largest inscribed squares (v334) on deformations of the 5/7 witness: a bump of depth t
under the base, and random vertex perturbations of the terrain alone.
"""
from fractions import Fraction as Q
import math
import random

import numpy as np


# ---------------------------------------------------------------- part A

J = np.array([[0.0, -1.0], [1.0, 0.0]])


def ellipse(a, b, samples=20000):
    """Parametrization by turns and a sampled signed distance (negative inside)."""
    g = lambda t: np.stack([a*np.cos(2*np.pi*t), b*np.sin(2*np.pi*t)], -1)
    C = g(np.linspace(0, 1, samples, endpoint=False))

    def f(P):
        d = np.sqrt(((P[:, None, :]-C[None, :, :])**2).sum(-1)).min(1)
        return np.where((P[:, 0]/a)**2+(P[:, 1]/b)**2 < 1, -d, d)
    return g, f


def core_loop_winding(a, b, c, points=1500):
    """Winding of h along the core loop y = x + c and the minimum of |h|_inf along it."""
    g, f = ellipse(a, b)
    xs = np.linspace(0, 1, points, endpoint=False)
    X, Y = g(xs), g(xs+c)
    W = (Y-X) @ J.T
    vals = np.stack([f(X+W), f(Y+W)], -1)
    ang = np.angle(vals[:, 0]+1j*vals[:, 1])
    d = (np.diff(np.concatenate([ang, ang[:1]]))+math.pi) % (2*math.pi)-math.pi
    w = d.sum()/(2*math.pi)
    if abs(w-round(w)) > 1e-6:
        raise ArithmeticError("winding is not an integer; refine the sampling")
    return round(w), round(float(np.abs(vals).max(1).min()), 6)


def circle_antipodal_corners():
    """On the unit circle, y = x + 1/2 gives S1 = g(x)(1 - 2i), S2 = -g(x)(1 + 2i): both at radius sqrt 5."""
    return abs(complex(1, -2)), abs(complex(1, 2))


# ---------------------------------------------------------------- part B

WITNESS_TERRAIN = [(0, 0), (8, 0), (8, 3), (7, 3), (7, 7), (3, 7), (3, 8), (0, 8)]


def bumped(t, c=Q(6)):
    """The 5/7 witness with a bump of depth t under the base on [c, 8] (t = 0: the terrain alone)."""
    if t == 0:
        return [(Q(x), Q(y)) for x, y in WITNESS_TERRAIN]
    pts = [(0, 0), (c, 0), (c, -t), (8, -t), (8, 3), (7, 3), (7, 7), (3, 7), (3, 8), (0, 8)]
    return [(Q(x), Q(y)) for x, y in pts]


def perturbed_terrain(seed, size=Q(1, 1000)):
    rng = random.Random(seed)
    return [(Q(x)+rng.randint(-1000, 1000)*size/1000, Q(y)+rng.randint(-1000, 1000)*size/1000) for x, y in WITNESS_TERRAIN]


def report():
    from terrain_gluing_v334 import largest_square
    from polygon_orbit_census_v223 import pt, validate_jordan_polygon
    part_a = {f"{a}x{b}": [[c, *core_loop_winding(a, b, c)] for c in (0.3, 0.4, 0.5, 0.6, 0.7)]
              for a, b in ((1.0, 1.0), (1.05, 1.0))}
    depths = {str(t): str(largest_square(bumped(t))[0]) for t in (Q(0), Q(1, 1000), Q(1, 100), Q(1, 10), Q(1, 2), Q(1))}
    ends = {str(c): str(largest_square(bumped(Q(1), c))[0]) for c in (Q(7), Q(69, 10), Q(13, 2), Q(6), Q(4), Q(1))}
    perturbed = []
    for seed in range(12):
        P = perturbed_terrain(seed)
        validate_jordan_polygon([pt(*p) for p in P])
        perturbed.append(round(float(largest_square(P)[0])**0.5, 4))
    return {"core_loop_windings_and_margins": part_a,
            "antipodal_corner_radii_on_circle": [round(r, 12) for r in circle_antipodal_corners()],
            "claimed_side_at_a_3_10": round(2*math.sin(0.3*math.pi), 6), "unit_circle_square_side": round(math.sqrt(2), 6),
            "bump_depth_to_M_squared": depths, "bump_left_end_to_M_squared": ends,
            "terrain_perturbed_1e-3_M": perturbed}


if __name__ == "__main__":
    import json
    print(json.dumps(report(), indent=2))
