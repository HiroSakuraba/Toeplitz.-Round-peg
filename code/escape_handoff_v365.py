"""v365: chaining comparisons through an escaping crossing, and changing the third stretch.

Exact data reuse the v363 triangle and the v364 sloping-roof fixture. Only the installation
schedule is new: the crossing now passes the seam between two comparison stages while its
fourth corner is still inside the target.
"""
from fractions import Fraction as Q
from itertools import product

from cutoff_interactions_v363 import box_size_intervals
from three_patch_sweep_v364 import completion_boxes

ROOT = Q(13, 14)          # full-mask installed amount at which the square has the cutoff side
KNOTS = (Q(0), Q(1, 3), Q(2, 3), Q(1))
SEAM = Q(2, 3)            # stage one is t in [0, 2/3], stage two is t in [2/3, 1]


def levels(lam):
    """Installed amounts at the four knots. Both interior levels rise through ROOT."""
    lam = Q(lam)
    if not 0 <= lam <= 1:
        raise ValueError("invalid sweep parameter")
    return (Q(0), ROOT - Q(3, 56) + lam / 8, ROOT - Q(1, 56) + lam / 14, Q(1))


def crossings(lam):
    """Edit times t at which the full-mask square has exactly the cutoff side."""
    z = levels(lam)
    times = set()
    for i in range(3):
        a, b = z[i], z[i + 1]
        if a != b and min(a, b) <= ROOT <= max(a, b):
            times.add(KNOTS[i] + (ROOT - a) / (b - a) * (KNOTS[i + 1] - KNOTS[i]))
        elif a == b == ROOT:
            raise ArithmeticError("a whole segment sits on the cutoff")
    return tuple(sorted(times))


def seam_events():
    """Sweep positions at which the crossing passes an interior knot."""
    z1 = lambda lam: levels(lam)[1]
    z2 = lambda lam: levels(lam)[2]
    # both interior levels are affine in lam
    out = {}
    for name, f in (('t=1/3', z1), ('t=2/3', z2)):
        a, b = f(0), f(1)
        out[name] = (ROOT - a) / (b - a)
    return out


def stage_counts(lam):
    """Number of marked crossings at a fixed sweep position in each stage (seam counted in both)."""
    ts = crossings(lam)
    return (sum(1 for t in ts if t <= SEAM), sum(1 for t in ts if t >= SEAM))


# ---------------------------------------------------------------------------
# Changing the third stretch: Boolean bookkeeping and the separation test

def cube_difference(p, bits, background=frozenset()):
    """XOR of p over all masks background | S with S a subset of bits."""
    total = 0
    bits = tuple(bits)
    for chosen in product((0, 1), repeat=len(bits)):
        total ^= p(frozenset(background) | {b for b, c in zip(bits, chosen) if c})
    return total


def union_bit_difference(p, first, second, third, extra):
    """Cubic difference when one bit switches the edits in `third` and `extra` together."""
    total = 0
    for chosen in product((0, 1), repeat=3):
        mask = {b for b, c in zip((first, second, 'joint'), chosen) if c}
        if 'joint' in mask:
            mask = (mask - {'joint'}) | {third, extra}
        total ^= p(frozenset(mask))
    return total


def _box_distance(a, b):
    dx = max(0, a[0] - b[1], b[0] - a[1])
    dy = max(0, a[2] - b[3], b[2] - a[3])
    return dx, dy


def union_separation(k1, k2, k3, k3p, side):
    """Squared lower bound for the distance from every possible fourth corner to all four patches,
    over both marked systems (k1,k2,k3) and (k1,k2,k3p), or None if it may vanish.

    A side-`side` square with corners in k1, k2 and the union k3 | k3p has its fourth corner in one
    of the completion boxes; the bound is the smallest box-to-patch distance."""
    patches = [tuple(map(Q, b)) for b in (k1, k2, k3, k3p)]
    best = None
    for third in (patches[2], patches[3]):
        for box in completion_boxes([patches[0], patches[1], third], side):
            for patch in patches:
                dx, dy = _box_distance(box, patch)
                d2 = dx * dx + dy * dy
                if d2 == 0:
                    return None
                best = d2 if best is None else min(best, d2)
    return best


def four_patch_hit(boxes, side):
    """True when some side-`side` square meets all four boxes at distinct corners."""
    s2 = Q(side) ** 2
    return any(lo <= s2 <= hi for lo, hi in box_size_intervals(boxes))
