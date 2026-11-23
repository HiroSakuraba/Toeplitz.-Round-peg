"""Replacement components at a simultaneous loss of two marked square pairs."""
from fractions import Fraction as Q


def double_loss_witness(W, H, p, q, g, A, B, h):
    W, H, p, q, g, A, B, h = map(Q, (W, H, p, q, g, A, B, h))
    if not (0 < p < q < A < B < W and 0 < min(g, h) <= max(g, h) < H):
        raise ValueError("require two disjoint proper roof cuts")
    if W-B > A or not (A < h+q < min(B, H) and g < B):
        return None
    t = min(H, h+A)
    if t < B:
        u = h
    elif H > B:
        t, u = B, B-q
    else:
        t = B
        u = max(h, min(B/2, B-q))
    v = t-u
    return ((u, 0), (t, u), (v, t), (0, v)), u*u+v*v
