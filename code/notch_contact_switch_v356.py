"""Exact square coordinates for the roof switch and the far-wall replacement."""
from fractions import Fraction as Q


def switch_witnesses(W, H, p, q, g, A, B, h):
    W, H, p, q, g, A, B, h = map(Q, (W, H, p, q, g, A, B, h))
    if not (0 < p < q < A < B < W and 0 < min(g, h) <= max(g, h) < H):
        raise ValueError("require two disjoint proper roof cuts")
    result = []
    if W-B > A:
        return result
    if max(A, g) < h+p < B and H > A:
        t = min(H, h+p)
        v = t-h
        result.append(("roof switch", ((h, 0), (t, h), (v, t), (0, v)), h*h+v*v))
    side = max(g, B-q)
    x = B-side
    if side >= A and h < side < H and p < x:
        result.append(("floor square", ((x, 0), (B, 0), (B, side), (x, side)), side*side))
    return result
