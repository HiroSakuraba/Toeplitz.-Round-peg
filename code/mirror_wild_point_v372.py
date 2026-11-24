"""v372: mirror-symmetric wild points.

The mirror is the y-axis and the wild point is the origin. Arm one lies in x > 0, arm two is its mirror
image. A drawing at depth n keeps the outer part of the arm and replaces the inner part by the segment
from the origin to the arm's last exit point, so the local arc stays mirror-symmetric.

Squares fixed by the reflection have their mirror line as a midline: corners (+-u, y), (+-u, y + 2u)
with (u, y) and (u, y + 2u) on arm one. These are the chords of arm one joining a point P to the
sheared point S(P) = P + (0, 2 x_P); they are counted exactly below.
"""
from fractions import Fraction as Q
from math import pi

from wild_point_parity_v370 import _unit


def curl(center_x, size, turns, per_turn, inner=Q(3, 10)):
    """A finite polygonal double spiral of outer radius `size` about (center_x, 0), entering at angle
    pi and leaving at angle 0, with `turns` turns on each arm (all rational)."""
    steps = turns * per_turn
    armA = []
    for k in range(steps + 1):
        t = Q(k, steps)
        r = size * (1 - (1 - inner) * t)
        ux, uy = _unit(pi + 2 * pi * turns * float(t))
        armA.append((center_x + r * ux, r * uy))
    # arm B is the half-turn image of arm A about the centre, traversed outward
    armB = [(2 * center_x - x, -y) for x, y in reversed(armA)]
    return armA + armB[1:]


def mirror_arm(depth, D=Q(1), turns=2, per_turn=8, zig=None):
    """Arm one from the origin outward at a given depth: a segment from the origin to the entry of
    curl `depth`, then curls depth, ..., 1 at x = D 2^{-k} of size D 2^{-k}/4, then out to x = 2D.
    `zig`, if given, is a list of extra outer vertices inserted after the last curl (for controls)."""
    pts = [(Q(0), Q(0))]
    for k in range(depth, 0, -1):
        c = D / 2 ** k
        pts += curl(c, c / 4, turns, per_turn)
    if zig:
        pts += [tuple(map(Q, z)) for z in zig]
    pts.append((2 * D, Q(0)))
    return pts


GENERIC_TOP = [(Q(2), Q(7, 10)), (Q(1, 3), Q(13, 10)), (Q(-2), Q(9, 10))]


def mirror_drawing(depth, top=None, **kw):
    """Closed polygon: arm two (mirror of arm one) inward, arm one outward, then an outer closure.
    The default closure is generic (not mirror-symmetric), so only the local arc is symmetric."""
    A = mirror_arm(depth, **kw)
    left = [(-x, y) for x, y in reversed(A)]
    return left[:-1] + A + list(top or GENERIC_TOP)


def in_cone(points):
    """Every point other than the origin satisfies |y| < x (the open 45-degree cone)."""
    return all(x > abs(y) for x, y in points if (x, y) != (0, 0))


def _singular_meets_box(m11, m12, r1, m21, m22, r2):
    """For a singular 2x2 system in (t, w), whether it has a solution with 0 <= t, w <= 1."""
    rows = [(m11, m12, r1), (m21, m22, r2)]
    live = [r for r in rows if r[0] != 0 or r[1] != 0]
    if any(r[0] == 0 and r[1] == 0 and r[2] != 0 for r in rows):
        return False
    if not live:
        return True
    a, b, c = live[0]
    for a2, b2, c2 in live[1:]:
        # rows are proportional (det = 0); consistent only if right-hand sides match the ratio
        if (a2 * c - a * c2) != 0 or (b2 * c - b * c2) != 0:
            return False
    values = [a * t + b * w for t in (0, 1) for w in (0, 1)]
    return min(values) <= c <= max(values)


def symmetric_chords(arm, max_u=None):
    """Exact chords of the polyline `arm` (in x > 0) of the form P, P + (0, 2 x_P), with u = x_P > 0.
    Returns sorted list of (u, P). Raises ArithmeticError when an edge and a sheared edge overlap,
    where chords may form a continuum; absence is certified only when nothing is raised."""
    segs = list(zip(arm, arm[1:]))
    out = set()
    for a, b in segs:
        for c, d in segs:
            # P = a + t (b - a);  S(P) = (Px, Py + 2 Px) = c + w (d - c)
            dx, dy = b[0] - a[0], b[1] - a[1]
            ex, ey = d[0] - c[0], d[1] - c[1]
            # Px = c_x + w e_x ;  Py + 2 Px = c_y + w e_y
            # a_x + t dx = c_x + w ex ;  a_y + t dy + 2 (a_x + t dx) = c_y + w ey
            m11, m12, r1 = dx, -ex, c[0] - a[0]
            m21, m22, r2 = dy + 2 * dx, -ey, c[1] - a[1] - 2 * a[0]
            det = m11 * m22 - m12 * m21
            if det == 0:
                if _singular_meets_box(m11, m12, r1, m21, m22, r2):
                    raise ArithmeticError("degenerate chord family: an edge and a sheared edge overlap")
                continue
            t = (r1 * m22 - m12 * r2) / det
            w = (m11 * r2 - r1 * m21) / det
            if 0 <= t <= 1 and 0 <= w <= 1:
                P = (a[0] + t * dx, a[1] + t * dy)
                if P[0] > 0 and (max_u is None or P[0] < max_u):
                    out.add((P[0], P))
    return sorted(out)


def reflect(points):
    return frozenset((-x, y) for x, y in points)
