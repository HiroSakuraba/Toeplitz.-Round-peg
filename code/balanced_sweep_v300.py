"""Exact piecewise-linear checks for the v300 boundary-mass sweep."""

import argparse
from bisect import bisect_left, bisect_right
from fractions import Fraction as Q
import json
from pathlib import Path


ARTIFACT = Path(__file__).resolve().parents[1] / "results/v300/balanced-sweep.json"
WIDE_PACKETS = ((Q(0), Q(1, 16)), (Q(15, 32), Q(17, 32)), (Q(15, 16), Q(1)))
TALL_PACKET = (Q(7, 16), Q(9, 16))


class LinearMap:
    def __init__(self, knots, values):
        self.knots = tuple(map(Q, knots))
        self.values = tuple(map(Q, values))
        if len(self.knots) != len(self.values) or len(self.knots) < 2:
            raise ValueError("one value per knot, with at least two knots")
        if any(a >= b for a, b in zip(self.knots, self.knots[1:])):
            raise ValueError("knots must be strictly increasing")

    def __call__(self, x):
        x = Q(x)
        if not self.knots[0] <= x <= self.knots[-1]:
            raise ValueError("outside the interpolation domain")
        i = min(bisect_right(self.knots, x) - 1, len(self.knots) - 2)
        a, b = self.knots[i:i + 2]
        u, v = self.values[i:i + 2]
        return u + (v - u) * (x - a) / (b - a)

    def inverse(self):
        return LinearMap(self.values, self.knots)


def cumulative(weight):
    """Integral of the prescribed convex mixture of the wide and tall densities."""
    weight = Q(weight)
    if not 0 <= weight <= 1:
        raise ValueError("mixture outside [0,1]")
    knots = sorted({p for interval in (*WIDE_PACKETS, TALL_PACKET) for p in interval})
    values = [Q(0)]
    for a, b in zip(knots, knots[1:]):
        middle = (a + b) / 2
        wide = Q(1, 4) + 4 * any(l < middle < r for l, r in WIDE_PACKETS)
        tall = Q(1, 2) + 4 * (TALL_PACKET[0] < middle < TALL_PACKET[1])
        values.append(values[-1] + (b - a) * ((1 - weight) * wide + weight * tall))
    if values[-1] != 1:
        raise ArithmeticError("unnormalized density")
    return LinearMap(knots, values)


def profile_value(depth, x):
    if depth < 0 or not 0 <= x <= 1:
        raise ValueError("invalid refinement depth or source")
    x = Q(x)
    y = Q(9, 8) * min(x, 1 - x)
    multiplier = 1 if x < Q(1, 2) else 2
    for n in range(1, depth + 1):
        phase = (4**n * x) % 1
        y += multiplier * 2 * min(phase, 1 - phase) / 2**(n + 6)
    return y


def profile(depth):
    if depth < 0:
        raise ValueError("negative depth")
    divisions = 2 * 4**depth
    knots = [Q(i, divisions) for i in range(divisions + 1)]
    return LinearMap(knots, [profile_value(depth, x) for x in knots])


def zero_intervals(knots, residual):
    """Complete zero set of a map affine between the supplied sorted knots."""
    roots = []
    values = [residual(x) for x in knots]
    for a, b, u, v in zip(knots, knots[1:], values, values[1:]):
        if u == 0 and v == 0:
            lo, hi = a, b
        elif u == 0:
            lo = hi = a
        elif v == 0:
            lo = hi = b
        elif u * v < 0:
            lo = hi = a - u * (b - a) / (v - u)
        else:
            continue
        if roots and roots[-1][1] == lo:
            roots[-1] = (roots[-1][0], hi)
        else:
            roots.append((lo, hi))
    return roots


def balanced_slice(height, weight):
    """Isolate ALL equal-level half-mass intervals for one rational weight.

    The coordinate p is the mass preceding the left endpoint, between 0 and 1/2.
    The output gives maximal zero intervals and exact imbalance extrema over them.
    """
    mass = cumulative(weight)
    source = mass.inverse()
    half = Q(1, 2)
    breaks = {mass(x) for x in (*height.knots, *mass.knots)}
    cuts = sorted({p for b in breaks for p in (b, b - half) if 0 <= p <= half})

    def residual(p):
        return height(source(p)) - height(source(p + half))

    roots = zero_intervals(cuts, residual)
    imbalances = []
    for lo, hi in roots:
        for p in {lo, hi, *cuts[bisect_left(cuts, lo):bisect_right(cuts, hi)]}:
            s, t = source(p), source(p + half)
            imbalances.append(t - s - 2 * height(s))
    return roots, (min(imbalances), max(imbalances))


def physical_squares(height):
    """Enumerate the equal-level square locus before imposing a weight.

    Split t=s+2h(s) whenever it crosses a profile knot, also at t=1.
    Return positive squares as closed intervals of left coordinates.
    """
    roots = []
    for a, b in zip(height.knots, height.knots[1:]):
        ta, tb = a + 2 * height(a), b + 2 * height(b)
        cuts = {a, b}
        if ta != tb:
            lower, upper = sorted((ta, tb))
            for k in height.knots[bisect_right(height.knots, lower):bisect_left(height.knots, upper)]:
                cuts.add(a + (k - ta) * (b - a) / (tb - ta))
        ordered = sorted(cuts)
        for left, right in zip(ordered, ordered[1:]):
            middle = (left + right) / 2
            if middle + 2 * height(middle) > 1:
                continue

            def gap(s):
                return height(s + 2 * height(s)) - height(s)

            for lo, hi in zero_intervals([left, right], gap):
                if lo == hi and height(lo) == 0:
                    continue
                if height(lo) == 0 or height(hi) == 0:
                    raise ValueError("an interval of squares meets a degenerate endpoint")
                if roots and roots[-1][1] == lo:
                    roots[-1] = (roots[-1][0], hi)
                else:
                    roots.append((lo, hi))
    return roots


def balancing_weight(s, t):
    q0, q1 = (mass(t) - mass(s) for mass in (cumulative(0), cumulative(1)))
    if q0 == q1:
        return (Q(0), Q(1)) if q0 == Q(1, 2) else None
    weight = (Q(1, 2) - q0) / (q1 - q0)
    return (weight, weight) if 0 <= weight <= 1 else None


def atlas_cost(parameters):
    """Four normalized source points in increasing cyclic order."""
    a, b, c, d = map(Q, parameters)
    if not 0 <= a < b < c < d < 1:
        raise ValueError("expected four ordered source points on the unit circle")
    red = b - a + d - c
    return red, 1 / (2 * min(red, 1 - red))


def atlas_band_red(x):
    return x * (169 + 71 * x) / (48 * (1 + 3 * x + x*x))


def report():
    atlas = {}
    for name, tags in [
        ("unit", [Q(0), Q(1, 16), Q(1, 2), Q(11, 12)]),
        ("persistent", [Q(9, 176), Q(5, 11), Q(25, 44), Q(31, 33)]),
    ]:
        red, cap = atlas_cost(tags)
        atlas[name] = {"red_fraction": str(red), "minimum_density_cap": str(cap)}
    rows = []
    for depth in range(5):
        height = profile(depth)
        wide, wide_range = balanced_slice(height, 0)
        tall, tall_range = balanced_slice(height, 1)
        if wide_range[0] < Q(5, 32) or tall_range[1] > -Q(1, 4):
            raise ArithmeticError("endpoint aspect bound failed")
        squares = physical_squares(height)
        examples = []
        for lo, hi in squares:
            # These fixtures have isolated roots. The enumerator itself also handles intervals.
            if lo != hi:
                raise ArithmeticError("fixture has a continuum of level squares")
            t = lo + 2 * height(lo)
            weights = balancing_weight(lo, t)
            if weights is not None:
                examples.append({"s": str(lo), "t": str(t), "weight": str(weights[0])})
        if not examples:
            raise ArithmeticError("missing balanced square")
        rows.append({"depth": depth, "wide_zero_components": len(wide),
                     "wide_min_width_minus_height": str(wide_range[0]),
                     "tall_zero_components": len(tall),
                     "tall_max_width_minus_height": str(tall_range[1]),
                     "level_square_components": len(squares), "balanced_square_samples": examples})
    return {"density_cap": "9/2", "square_side_floor": "1/9",
            "atlas_controls": atlas, "atlas_band_cap": str(1 / (2 * (1 - atlas_band_red(Q(3, 5))))),
            "refinements": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    data = json.dumps(report(), indent=2, sort_keys=True) + "\n"
    if args.write:
        ARTIFACT.write_text(data, encoding="utf-8")
    elif ARTIFACT.read_text(encoding="utf-8") != data:
        raise SystemExit("v300 certificate differs from replay")
    print("Verified the boundary-mass sweep and the exact refinement records.")


if __name__ == "__main__":
    main()
