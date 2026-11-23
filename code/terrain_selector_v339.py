"""v339 strict contact selector and v340 component certificates for grid terrains."""
import csv
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
import random
import sys


@dataclass(frozen=True)
class Certificate:
    method: str
    side_squared: int | Fraction
    contacts: tuple  # (left, roof, right) for crossings, (plateau end,) for plateaus


def contact_crossings(heights, *, permit_one_corner=False):
    """Base and three vertical supports, optionally crossing a single wall corner."""
    walls = {}
    for x, (before, after) in enumerate(zip((0, *heights), (*heights, 0))):
        if before != after:
            walls[x] = before, after

    def interior(x, y):
        before, after = walls[x]
        return min(before, after) < y < max(before, after)

    def crosses_corner(x, y, slope):
        # Traverse the roof left to right. At a wall endpoint the two
        # incident rays are one horizontal and one vertical. Compare the
        # signs of X + slope*Y minus its value at the corner on these rays.
        before, after = walls[x]
        vertical = slope * (after-before)
        if y == before and x > 0:
            return vertical > 0  # incoming horizontal ray points left
        if y == after and x < len(heights):
            return vertical > 0  # outgoing horizontal ray points right
        return False

    for left, roof, right in combinations(walls, 3):
        u, v = roof-left, right-roof
        if not interior(roof, u+v):
            continue
        left_open, right_open = interior(left, u), interior(right, v)
        if left_open and right_open:
            method = "strict"
        elif permit_one_corner and (
            (right_open and crosses_corner(left, u, +1))
            or (left_open and crosses_corner(right, v, -1))
        ):
            method = "one corner"
        else:
            continue
        yield Certificate(method, u*u+v*v, (left, roof, right))


def end_component(heights):
    width = len(heights)
    if min(heights[0], heights[-1]) > width:
        return Certificate("tall walls", width**2, ())
    for orientation, h in (("left", heights), ("right", heights[::-1])):
        if h[0] != width or h[-1] <= width:
            continue
        a = next(x for x, height in enumerate(h) if height != width)
        if h[a] > width:
            u = min(a, Fraction(width, 2))
            return Certificate(f"{orientation} plateau", u*u+(width-u)**2,
                               (a if orientation == "left" else width-a,))
    return None


def sample():
    """Replay v338's seed law and original end-wall classifier, preserving repeats."""
    from robust_parity_v338 import interior_case, terrain_i

    for seed in range(300):
        rng = random.Random(seed)
        heights = tuple(rng.randint(1, 8) for _ in range(rng.randint(2, 8)))
        interior_side = terrain_i(heights)
        if interior_side < 2 or interior_case(heights) is not None:
            continue
        crossings = list(contact_crossings(heights, permit_one_corner=True))
        best = max((c for c in crossings if c.method == "strict"),
                   key=lambda c: c.side_squared, default=None)
        strict_ok = best is not None and 2*best.side_squared >= interior_side**2
        # Give the strict selector first credit so added coverage is visible.
        cert = best if strict_ok else end_component(heights)
        if cert is None:
            cert = max((c for c in crossings if c.method == "one corner"
                        and 2*c.side_squared >= interior_side**2),
                       key=lambda c: c.side_squared, default=None)
        yield (seed, heights, interior_side, best, cert)


def write_report(stream):
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("seed", "heights", "area", "interior_side", "best_strict_side_squared",
                     "method", "certified_side_squared", "support_abscissas"))
    for seed, heights, i, best, cert in sample():
        writer.writerow((seed, " ".join(map(str, heights)), sum(heights), i,
                         best.side_squared if best else "", cert.method if cert else "unresolved",
                         cert.side_squared if cert else "",
                         " ".join(map(str, cert.contacts)) if cert else ""))


if __name__ == "__main__":
    write_report(sys.stdout)
