"""Certified squares when a second roof notch covers the turning endpoint."""
from fractions import Fraction as Q


def endpoint_witness(A, B, d, h, p, q, g):
    """Return a marked square and the minimum squared side on its component."""
    A, B, d, h, p, q, g = map(Q, (A, B, d, h, p, q, g))
    if not (0 < d < p and 0 < h and p+h < A < B and
            2*h <= A and p < B-h <= q < A and 0 < g < B):
        raise ValueError("dimensions outside the endpoint theorem")
    t = max(A, g)
    if g < A:
        u = A-p
    elif g > A:
        u = h
    else:
        u = max(h, min(t/2, A-p))
    v = t-u
    return ((u, 0), (t, u), (v, t), (0, v)), u*u+v*v
