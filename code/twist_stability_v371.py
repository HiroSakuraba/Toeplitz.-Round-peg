"""v371: stability of the odd large graceful degree, and twisted double spirals.

By Lemma 370.1 a Jordan polygon has aligned degree 1 and reversed degree 0, so the parity of its
large graceful squares (side above s, aligned plus reversed) is 1 + small aligned + small reversed,
computable from the small squares alone.
"""
from wild_point_parity_v370 import drawing, sector_counts


def large_graceful_parity(V, s):
    """Mod-two degree of aligned plus reversed squares of side above s (no square of side exactly s
    and all small squares transverse are checked by the caller through the returned counts)."""
    counts, degenerate, squares = sector_counts(V, s)
    if degenerate:
        raise ArithmeticError("degenerate small square candidate")
    return (1 + counts['aligned'] + counts['reversed']) % 2, counts


def twisted(depth, twist, per_turn=24, a=0.5, theta0=3.0):
    return drawing(a, theta0, depth, per_turn, twist=twist)
