"""Exact sufficient certificates for a rectangle with one roof notch."""
from fractions import Fraction as Q

from convex_l_disk_v344_v345 import OrthogonalL


def notch_certificate(W, H, left, right, h):
    """Return (wider comparison L, method, squared robust floor), or None floor."""
    W, H, left, right, h = map(Q, (W, H, left, right, h))
    n = W-left-right
    if min(n, left, right, h) <= 0 or h >= H:
        raise ValueError("require two positive towers and a proper roof notch")
    A = max(left, right)
    model = OrthogonalL(W, H, A, h)
    if H > W:
        return model, "full width", W*W
    if W > H and W-H < min(left, right):
        return model, "full height", H*H
    cert = model.certificate()
    if A >= H or (cert.span is not None and cert.span < A+n):
        return model, cert.method, cert.side_squared
    return model, "unresolved", None
