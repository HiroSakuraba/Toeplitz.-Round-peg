"""Exact normalized scalar negative controls; not a Floer computation."""
from fractions import Fraction as Q
from itertools import combinations


def plateau(b, x, mixing=0):
    b, x, mixing = map(Q, (b, x, mixing))
    if not 0 < b <= 1 or not 0 <= x <= 1 or not 0 <= mixing <= 1:
        raise ValueError("outside normalized domain")
    a = b / 2
    return (1-mixing)*(min(x, a)+max(0, x-1+a))+mixing*b*x


def subdivision_vertices(b):
    """All vertices of the closed triangle cut by x,y,x+y break lines."""
    b = Q(b)
    if not 0 < b <= 1:
        raise ValueError("b must be in (0,1]")
    a = b / 2
    lines = {(1,0,0), (0,1,0), (1,1,1)}
    for t in (a, 1-a):
        lines.update(((1,0,t), (0,1,t), (1,1,t)))
    vertices = set()
    for (p,q,c),(r,s,d) in combinations(lines,2):
        det = p*s-q*r
        if det:
            x, y = Q(c*s-q*d, det), Q(p*d-c*r, det)
            if x >= 0 and y >= 0 and x+y <= 1:
                vertices.add((x,y))
    return sorted(vertices)


def audit(b, mixing=0):
    """Exact checks on all affine cells, not sampling a continuous inequality."""
    b, mixing = Q(b), Q(mixing)
    a = b / 2
    knots = sorted({Q(0),a,1-a,Q(1)})
    f = lambda x: plateau(b,x,mixing)
    slopes = [(f(y)-f(x))/(y-x) for x,y in zip(knots,knots[1:])]
    # Reflection error is affine between these symmetric knots.
    assert f(0) == 0 and f(1) == b
    assert all(0 <= s <= 1 for s in slopes)
    assert all(f(x)+f(1-x) == b for x in knots)
    vertices = subdivision_vertices(b)
    # Every affine cell is a bounded convex polygon; check all its vertices.
    slack = min(f(x)+f(y)-f(x+y) for x,y in vertices)
    assert slack >= 0
    return {"b":b, "mixing":mixing, "vertices":len(vertices),
            "min_subadditivity_slack":slack, "slopes":slopes,
            "central_slope":mixing*b if b < 1 else Q(1),
            "plateau_width":1-b if mixing == 0 else Q(0)}


def angular_floor(b, h):
    """Normalized central mean slope bound, h = half-window in units of pi."""
    b,h = Q(b),Q(h)
    if not 0 < b <= 1 or not 0 < h <= Q(1,2):
        raise ValueError("outside normalized domain")
    return max(Q(0), 1-(1-b)/(2*h))


def central_mean(b,h,mixing=0):
    h = Q(h)
    if not 0 < h <= Q(1,2):
        raise ValueError("outside normalized domain")
    return (plateau(b,Q(1,2)+h,mixing)-
            plateau(b,Q(1,2)-h,mixing))/(2*h)


if __name__ == "__main__":
    for b in (Q(1,100),Q(1,2),Q(9,10),Q(999,1000),Q(1)):
        for mixing in (Q(0),Q(1,1000),Q(1)):
            print(audit(b,mixing))
