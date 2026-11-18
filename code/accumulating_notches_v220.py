"""v220: fixed rectifiable Jordan curve with arbitrarily small exact square pegs.

The limit curve is a 2x1 rectangle with countably many pairwise-disjoint
outward rectangular notches, of widths/depths s_n=4^(-n-1), accumulating
at the origin from the positive x-axis. Coordinates, lengths, areas,
parameter gaps and residuals are exact fractions.

This is a rigorous obstruction to a *zero-free normalized collision-link*
for arbitrary rectifiable Jordan curves, not a counterexample to Toeplitz.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations
from typing import Tuple


@dataclass(frozen=True)
class Point:
    x: Q
    y: Q


def point(x, y):
    return Point(Q(x), Q(y))


def notch(n: int) -> Tuple[Point, Point, Point, Point]:
    if type(n) is not int or n < 1:
        raise ValueError("notch index must be a positive integer")
    a, s = Q(1, 2**n), Q(1, 4**(n+1))
    return (point(a, 0), point(a, -s), point(a+s, -s), point(a+s, 0))


def flat_control(n: int) -> Tuple[Point, Point, Point, Point]:
    """Four equally spaced points on the flat bottom just before notch n.

    The interval is disjoint from all smaller notches by the inequality
    a_n - 3 s_n > a_(n+1) + s_(n+1), for every n >= 1.
    """
    p = notch(n)
    a = p[0].x
    s = p[3].x - a
    if not a - 3*s > notch(n+1)[3].x:
        raise AssertionError("flat reference interval meets a smaller notch")
    return tuple(point(a - (3-j)*s, 0) for j in range(4))


def polygon(N: int) -> Tuple[Point, ...]:
    """Cyclic simple polygon with notches N,...,1 in x-increasing order."""
    if type(N) is not int or N < 0:
        raise ValueError("N must be nonnegative integer")
    p = [point(0,0)]
    for n in range(N,0,-1):
        p.extend(notch(n))
    p.extend((point(2,0), point(2,1), point(0,1)))
    return tuple(p)


def square_residual(p: Tuple[Point, ...], orientation=1) -> Tuple[Q,Q,Q,Q]:
    """v218's four real equations for an oriented (possibly collapsed) square."""
    if len(p) != 4 or orientation not in (-1,1):
        raise ValueError("four points and orientation +/-1 required")
    a,b,c,d=p
    ex,ey=b.x-a.x,b.y-a.y
    return (c.x-b.x+orientation*ey, c.y-b.y-orientation*ex,
            d.x-a.x+orientation*ey, d.y-a.y-orientation*ex)


def side_squared(p: Tuple[Point, ...]) -> Q:
    return (p[1].x-p[0].x)**2 + (p[1].y-p[0].y)**2


def squared_diameter(p: Tuple[Point,...]) -> Q:
    return max((a.x-b.x)**2+(a.y-b.y)**2 for a,b in combinations(p,2))


def signed_area(poly: Tuple[Point, ...]) -> Q:
    return sum(a.x*b.y - a.y*b.x for a,b in zip(poly,poly[1:]+poly[:1]))/2


def rectilinear_perimeter(poly: Tuple[Point, ...]) -> Q:
    total=Q(0)
    for a,b in zip(poly,poly[1:]+poly[:1]):
        if a.x!=b.x and a.y!=b.y:
            raise ValueError("only axis-aligned polygon edges supported")
        total+=abs(a.x-b.x)+abs(a.y-b.y)
    return total


def predicted_finite_area(N: int) -> Q:
    return Q(2) + Q(1,240)*(1-Q(1,16**N))


def predicted_finite_length(N: int) -> Q:
    return Q(6) + Q(1,6)*(1-Q(1,4**N))


def limiting_length() -> Q:
    """Sum of all edge lengths of a genuine rectifiable Jordan curve."""
    return Q(37,6)


def limiting_area() -> Q:
    return Q(2)+Q(1,240)


def omitted_tail_vertical(N: int) -> Q:
    """Upper bound on Hausdorff distance to the infinite-notch curve."""
    if type(N) is not int or N<0:
        raise ValueError("nonnegative N required")
    return Q(1,4**(N+2))


def _on_segment(a,b,p):
    return ((b.x-a.x)*(p.y-a.y)-(b.y-a.y)*(p.x-a.x)==0
            and min(a.x,b.x)<=p.x<=max(a.x,b.x)
            and min(a.y,b.y)<=p.y<=max(a.y,b.y))


def _cross(a,b,c):
    return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x)


def _intersect(a,b,c,d):
    o1,o2,o3,o4=_cross(a,b,c),_cross(a,b,d),_cross(c,d,a),_cross(c,d,b)
    if o1*o2<0 and o3*o4<0:
        return True
    return any((z==0 and _on_segment(u,v,p)) for z,u,v,p in
               ((o1,a,b,c),(o2,a,b,d),(o3,c,d,a),(o4,c,d,b)))


def assert_simple_polygon(poly: Tuple[Point, ...]):
    """Independent exact nonadjacent-segment intersection test."""
    n=len(poly)
    if n<3 or len(set(poly))!=n:
        raise ValueError("repeated/insufficient polygon vertices")
    for i in range(n):
        if poly[i]==poly[(i+1)%n]:
            raise ValueError("zero edge")
    for i in range(n):
        for j in range(i+1,n):
            if j==i+1 or (i==0 and j==n-1):
                continue
            if _intersect(poly[i],poly[(i+1)%n],poly[j],poly[(j+1)%n]):
                raise ValueError("nonadjacent polygon edges intersect")
    return True


def parameter_span(n:int) -> Q:
    """Span t3-t0 under normalized arclength on the INFINITE curve."""
    return Q(3)*Q(1,4**(n+1))/limiting_length()


def collision_chart_data(n:int, kind="square") -> tuple:
    """(starting normalized arclength, (a,b,c), scale epsilon).

    Both square and flat configurations have gaps (1/3,1/3,1/3) and
    starting parameter tending to zero; they reach the same collision
    chart point, even though their scaled square residuals disagree.
    """
    a,s = Q(1,2**n),Q(1,4**(n+1))
    # Sum of notch widths over k>n equals s_n/3.
    earlier_extra = 2*s/3
    if kind=="square":
        t0=(a+earlier_extra)/limiting_length()
    elif kind=="flat":
        t0=(a-3*s+earlier_extra)/limiting_length()
    else:
        raise ValueError("kind must be square or flat")
    return t0,(Q(1,3),)*3,parameter_span(n)


def normalized_residual_at_collision(n:int,kind="square") -> tuple:
    pts=notch(n) if kind=="square" else flat_control(n)
    eps=parameter_span(n)
    return tuple(x/eps for x in square_residual(pts,1))


def report(N=8):
    poly=polygon(N)
    assert_simple_polygon(poly)
    return {
        "N":N,
        "polygon_vertices":len(poly),
        "exact_area":str(signed_area(poly)),
        "exact_perimeter":str(rectilinear_perimeter(poly)),
        "limiting_Jordan_curve_perimeter":str(limiting_length()),
        "limiting_area":str(limiting_area()),
        "smallest_explicit_square_side":str(Q(1,4**(N+1))),
        "smallest_square_residual":tuple(str(x) for x in square_residual(notch(N))),
        "normalized_residual_at_notch":tuple(str(x) for x in normalized_residual_at_collision(N,"square")),
        "normalized_residual_on_same_gap_flat_sequence":tuple(str(x) for x in normalized_residual_at_collision(N,"flat")),
        "normalized_gap_direction":tuple(str(x) for x in collision_chart_data(N)[1]),
        "upper_Hausdorff_limit_error":str(omitted_tail_vertical(N)),
        "scope":"a fixed rectifiable Jordan curve with infinitely many shrinking squares; not a Toeplitz counterexample",
    }


if __name__=="__main__":
    import json
    print(json.dumps(report(),sort_keys=True,indent=2))