"""Select the anchored square family on an inward-tilted L, in rational arithmetic."""
from fractions import Fraction as Q


def section(frontier, t, axis):
    values = []
    for a, b in zip(frontier, frontier[1:]):
        x, y = a[axis], a[1-axis]
        dx, dy = b[axis]-x, b[1-axis]-y
        if min(x, x+dx) <= t <= max(x, x+dx):
            if dx:
                values.append(y+(t-x)*dy/dx)
            else:
                values.extend((y, y+dy))
    return min(values), max(values)


def anchor_interval(frontier):
    """Return (span, least u, greatest u) for a coordinate-monotone frontier."""
    frontier = [tuple(map(Q, p)) for p in frontier]
    m = min(frontier[0][0], frontier[-1][1])
    knots = sorted({x for p in frontier for x in p if 0 <= x <= m} | {m})
    previous = None
    for t in knots:
        a, A = section(frontier, t, 0)
        b, B = section(frontier, t, 1)
        if a+b <= t <= A+B:
            return t, max(a, t-B), min(A, t-b)
        # Inside each open knot interval, both sections are affine and single-valued.
        if previous is not None and A+B < t:
            x, gap = previous
            root = x+(t-x)*gap/(gap-(A+B-t))
            u, _ = section(frontier, root, 0)
            return root, u, u
        previous = t, a+b-t
    raise ValueError("no anchor crossing in the common span")


def witness(W, H, a, b, c):
    W, H, a, b, c = map(Q, (W, H, a, b, c))
    if not (0 < c < a < W and 0 < b < H):
        raise ValueError("require a nondegenerate inward tilt")
    if b >= W:
        d = min(W-a, W/2) if b == W else 0
        return ((0, 0), (W, 0), (W, W), (0, W)), d*d+(W-d)**2
    if c >= H:
        return ((0, 0), (H, 0), (H, H), (0, H)), H*H
    frontier = ((W, 0), (W, b), (a, b), (c, H), (0, H))
    t, lo, hi = anchor_interval(frontier)
    u = max(lo, min(t/2, hi))
    v = t-u
    return ((u, 0), (t, u), (v, t), (0, v)), u*u+v*v
