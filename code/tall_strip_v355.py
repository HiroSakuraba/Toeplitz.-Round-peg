"""Exact two-notch square witness and grounded-disk feasibility checks."""
from fractions import Fraction as Q


def tall_strip_witness(W, H, p, q, g, A, B, h):
    W, H, p, q, g, A, B, h = map(Q, (W, H, p, q, g, A, B, h))
    if not (0 < p < q < A < B < W and 0 < min(g, h) <= max(g, h) < H):
        raise ValueError("require two disjoint proper roof cuts")
    t = max(A, g, h+p)
    lo, hi = h, t-p
    if t > A:
        lo = hi = h
    if t > g:
        lo = hi = t-p
    if W-B > A or t >= min(H, B) or t-lo >= q:
        return None
    u = max(lo, min(t/2, hi))
    v = t-u
    return ((u, 0), (t, u), (v, t), (0, v)), u*u+v*v


def roof_disk_fits(W, H, cuts, radius):
    """Decide closed-disk feasibility on a rational rectangle with sorted roof cuts."""
    W, H, r = map(Q, (W, H, radius))
    if W <= 0 or H <= 0 or r < 0:
        raise ValueError("invalid box or radius")
    active = []
    previous = 0
    for cut in cuts:
        a, b, h = map(Q, cut)
        if not (previous < a < b < W and 0 < h < H):
            raise ValueError("expected sorted disjoint proper roof cuts")
        previous = b
        if h < 2*r:
            reach2 = r*r-max(0, h-r)**2
            active.append((a, b, reach2))
    if 2*r > min(W, H):
        return False

    def separated(lower, upper):
        x, s = lower
        y, t = upper
        gap = y-x
        excess = gap*gap-s-t
        return gap >= 0 and excess >= 0 and excess*excess >= 4*s*t

    for split in range(len(active)+1):
        lower = [(0, r*r)] + [(b, s) for _, b, s in active[:split]]
        upper = [(W, r*r)] + [(a, s) for a, _, s in active[split:]]
        if all(separated(left, right) for left in lower for right in upper):
            return True
    return False
