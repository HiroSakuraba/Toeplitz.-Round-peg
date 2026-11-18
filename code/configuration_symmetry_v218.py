"""v218: exact cyclic-equivariance, collision-stratum and sector audit.

Four labeled points on a Jordan curve, with cyclic relabeling T.
No claim of a universal square, an all-Jordan global degree, or Floer theory.
Arithmetic for algebraic checks and finite polygon examples is exact Fraction.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import permutations, product
from typing import Tuple


@dataclass(frozen=True)
class Point:
    x: Q
    y: Q


def pt(x, y):
    return Point(Q(x), Q(y))


def add(a, b):
    return pt(a.x + b.x, a.y + b.y)


def sub(a, b):
    return pt(a.x - b.x, a.y - b.y)


def square_residual(p, orientation=1):
    """2 complex equations (4 real scalar conditions)."""
    if len(p) != 4 or orientation not in (-1, 1):
        raise ValueError("Four points and orientation +/-1 required")
    a, b, c, d = p
    ex, ey = b.x-a.x, b.y-a.y
    jx, jy = -orientation*ey, orientation*ex
    return (c.x-b.x-jx, c.y-b.y-jy,
            d.x-a.x-jx, d.y-a.y-jy)


def residual_cyclic_matrix(orientation=1):
    """F_sigma(T p)=B_sigma F_sigma(p), T=(1,2,3,0)."""
    if orientation == 1:
        return ((-1, 1, 1, 0),
                (-1,-1, 0, 1),
                ( 0, 1, 0, 0),
                (-1, 0, 0, 0))
    if orientation == -1:
        return ((-1,-1, 1, 0),
                ( 1,-1, 0, 1),
                ( 0,-1, 0, 0),
                ( 1, 0, 0, 0))
    raise ValueError("orientation +/-1 required")


def matmul(M, vec):
    return tuple(sum(Q(z)*Q(v) for z,v in zip(row,vec)) for row in M)


def determinant(M):
    """Exact Gaussian-elimination determinant for rational square matrix."""
    n=len(M)
    if n == 0 or any(len(row)!=n for row in M):
        raise ValueError("square matrix required")
    A=[list(map(Q,row)) for row in M]
    det=Q(1)
    for i in range(n):
        piv=next((j for j in range(i,n) if A[j][i]),None)
        if piv is None:return Q(0)
        if piv!=i:
            A[i],A[piv]=A[piv],A[i]
            det=-det
        z=A[i][i]
        det*=z
        for j in range(i+1,n):
            fac=A[j][i]/z
            for k in range(i+1,n):A[j][k]-=fac*A[i][k]
            A[j][i]=0
    return det


# Coordinate order: t0, g1, g2, g3, with g4=1-g1-g2-g3.
# T takes (t0,g1,g2,g3) -> (t0+g1,g2,g3,g4), locally mod 1 in t0.
CYCLIC_CHART_JACOBIAN = ((1, 1, 0, 0),
                         (0, 0, 1, 0),
                         (0, 0, 0, 1),
                         (0,-1,-1,-1))


def positive_gap_chart(t0, g1, g2, g3):
    g=(Q(g1),Q(g2),Q(g3))
    if any(x<=0 for x in g) or sum(g)>=1:
        raise ValueError("four strictly positive cyclic gaps required")
    return (Q(t0)%1,*g)


def cyclic_chart_shift(chart):
    t,g1,g2,g3=chart
    return ((t+g1)%1,g2,g3,1-g1-g2-g3)


def orbit_local_signs(sign=1):
    """For transverse zeros: sgn(DT)=-1 and det(B)=+1."""
    if sign not in (-1,1):
        raise ValueError("sign +/-1")
    domain=determinant(CYCLIC_CHART_JACOBIAN)
    target=determinant(residual_cyclic_matrix(1))
    assert domain==-1 and target==1
    return tuple(sign * (target*domain)**k for k in range(4))


def order_sectors():
    """Six connected cyclic order components of Conf_4(S^1) with labels.

    Each sector is recorded starting with label 0 around the parameter circle.
    'Graceful' means the square labels appear in their geometric cyclic
    order, either direction, on the source Jordan curve.
    """
    return tuple((0,)+p for p in permutations((1,2,3)))


def graceful_sector(s):
    return tuple(s) in ((0,1,2,3),(0,3,2,1))


def closed_gap_chart(g1,g2,g3,g4):
    gaps=tuple(map(Q,(g1,g2,g3,g4)))
    if any(x<0 for x in gaps) or sum(gaps)!=1:
        raise ValueError("closed four-gap simplex required")
    return gaps


def total_collision_vertex(gaps):
    """All parameters coincide modulo 1 iff gaps are a simplex vertex."""
    if len(gaps)!=4 or any(g<0 for g in gaps) or sum(gaps)!=1:
        raise ValueError("closed simplex gaps required")
    return sum(x==1 for x in gaps)==1


def polygon_point(vertices, t):
    """Cyclic PL curve with equal parameter time per polygon edge."""
    p=tuple(vertices)
    if len(p)<3:raise ValueError("3+ vertices required")
    t=Q(t)%1
    v=t*len(p);i=int(v);s=v-i
    a,b=p[i],p[(i+1)%len(p)]
    return pt(a.x+s*(b.x-a.x),a.y+s*(b.y-a.y))


def chart_points(vertices,t,gaps):
    if len(gaps)!=4:raise ValueError("four gaps required")
    closed_gap_chart(*gaps)
    a=Q(t)
    return tuple(polygon_point(vertices, a+sum(gaps[:i])) for i in range(4))


def is_nondegenerate_square(points, orientation):
    a,b=points[:2]
    return square_residual(points,orientation)==(0,0,0,0) and (
        (a.x-b.x)**2+(a.y-b.y)**2)>0


def count_orbit_signs(signs):
    return sum(signs)


def _root_on_segments(chosen, orientation):
    """Exact transversal square on four chosen directed segments, or None.

    This is a *diagnostic* for four distinct open edges, not an all-squares
    decision: squares on edge endpoints or occupying repeated edges are not
    included. Nontransverse segment systems also require another method.
    """
    starts=tuple(a for a,b in chosen)
    base=square_residual(starts,orientation)
    columns=[]
    for i in range(4):
        shifted=list(starts)
        shifted[i]=chosen[i][1]
        next_res=square_residual(shifted,orientation)
        columns.append(tuple(x-y for x,y in zip(next_res,base)))
    A=[[columns[j][i] for j in range(4)] for i in range(4)]
    b=[-v for v in base]
    if determinant(A)==0:return None
    # Exact Gauss-Jordan solve.
    M=[list(map(Q,row))+[Q(rhs)] for row,rhs in zip(A,b)]
    for i in range(4):
        p=next((j for j in range(i,4) if M[j][i]),None)
        assert p is not None
        M[i],M[p]=M[p],M[i]
        d=M[i][i];M[i]=[z/d for z in M[i]]
        for j in range(4):
            if i!=j:
                c=M[j][i]
                M[j]=[a-c*z for a,z in zip(M[j],M[i])]
    u=tuple(row[-1] for row in M)
    if not all(Q(0)<x<Q(1) for x in u):return None
    pts=tuple(pt(a.x+(b.x-a.x)*t, a.y+(b.y-a.y)*t)
              for (a,b),t in zip(chosen,u))
    if not is_nondegenerate_square(pts,orientation):return None
    return pts,determinant(A)


def isolated_distinct_edge_squares(vertices):
    """Exact finite audit on 4 DISTINCT OPEN polygon edges.

    Returns geometric square orbits with all eight orientation/rotation
    labelings represented when the curve's four square vertices lie in
    distinct open edges. Absence here is *inconclusive* for general polygons.
    """
    p=tuple(vertices);n=len(p)
    if n<4:raise ValueError("four polygon edges required")
    edges=tuple((p[i],p[(i+1)%n]) for i in range(n))
    found={}
    for indices in permutations(range(n),4):
        choice=tuple(edges[j] for j in indices)
        for sigma in (1,-1):
            root=_root_on_segments(choice,sigma)
            if root is None:continue
            pts,det=root
            key=tuple(sorted((q.x,q.y) for q in pts))
            found.setdefault(key,[]).append((indices,sigma,1 if det>0 else -1))
    return found


def run_audit():
    square=[pt(0,0),pt(2,0),pt(3,3),pt(0,2)]
    orbits=isolated_distinct_edge_squares(square)
    return {
        'domain_shift_det':str(determinant(CYCLIC_CHART_JACOBIAN)),
        'target_shift_det_plus':str(determinant(residual_cyclic_matrix(1))),
        'target_shift_det_minus':str(determinant(residual_cyclic_matrix(-1))),
        'orbit_local_signs':tuple(map(int,orbit_local_signs())),
        'orbit_sign_sum':int(count_orbit_signs(orbit_local_signs())),
        'cyclic_sectors':len(order_sectors()),
        'graceful_sectors':sum(map(graceful_sector,order_sectors())),
        'distinct_open_edge_square_orbits':len(orbits),
        'distinct_edge_labelings':tuple(sorted(len(records) for records in orbits.values())),
        'not_a_global_square_existence_test':True,
    }


if __name__=='__main__':
    import json
    print(json.dumps(run_audit(),indent=2,sort_keys=True))