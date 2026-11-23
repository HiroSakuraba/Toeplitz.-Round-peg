"""Complete the v351 sufficient tests at both exceptional box dimensions."""
from fractions import Fraction as Q
from rectangular_notch_v351 import notch_certificate


def complete_notch_certificate(W, H, left, right, h):
    model, method, floor = notch_certificate(W, H, left, right, h)
    if floor is not None:
        return model, method, floor
    W, H, left, right = map(Q, (W, H, left, right))
    if W == H:
        return model, "rotating square in square box", H**2/2
    assert H == W-min(left, right)
    return model, "sliding and rotating at notch wall", H**2/2
