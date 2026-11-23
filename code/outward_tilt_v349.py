"""Rational witnesses for an outward-tilted L with a thin lower arm."""

from fractions import Fraction as Q


def witnesses(W, H, a, b, c):
    """Return (method, marked square, robust squared floor) for each witness."""
    W, H, a, b, c = map(Q, (W, H, a, b, c))
    if not (0 < a < c < W and 0 < b < min(H, W)):
        raise ValueError("require outward tilt and a thin lower arm")
    if H <= a:
        return [("strip", ((0, 0), (H, 0), (H, H), (0, H)), H**2)]

    k = (c-a)/(H-b)
    m = min(W, H)
    t = min(m, b+a/(1-k)) if k < 1 else m
    f = a+k*(t-b)
    if t < W:
        u = b
    elif t < H:
        u = t-f
    else:
        u = max(t-c, min(t/2, b))
    v = t-u
    row = ("anchors", ((u, 0), (t, u), (v, t), (0, v)), u*u+v*v)
    if t == m:
        return [row]

    d = c/(1+k+k*k)
    y = b+(d-a)/k
    upper = ((d, y), ((1+k)*d, y+d), (k*d, H), (0, y+k*d))
    return [row, ("upper", upper, d*d*(1+k*k))]
