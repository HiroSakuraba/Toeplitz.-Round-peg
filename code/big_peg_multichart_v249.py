"""v249: an exact finite-chart square continuation certificate across a vertex event.

A one-parameter family of genuinely changing convex Jordan polygons has an
inscribed nondegenerate square of side > 27/10 FOR EVERY lambda in [0,1].
At lambda=1/2 a supporting square corner meets a *moving subdivision
vertex*. The first feasible four-edge chart ceases to apply, but a second
chart takes over. Exact polynomial Bernstein bounds, including all
endpoints, prove full coverage, and exact Cramer roots agree at the seam.

This is a controlled PL family with a reparametrization/subdivision event;
it does NOT prove that arbitrary square birth/death or singular orbits
persist, or give a universal Big Peg bound for all Jordan curves.
Standard library only; no float proof decisions or interval libraries.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
from math import comb
import json

from polygon_orbit_census_v223 import pt, validate_jordan_polygon, sub, cross, dot
from big_peg_persistence_v248 import OCTAGON, CORNER_COEFFS

Z=(F(0),)
O=(F(1),)
HALF=F(1,2)
ALPHA=F(39906,73777)
S2_LIMIT=F(27,10)**2


def pol(*coefficients):
    p=tuple(F(x) for x in coefficients)
    while len(p)>1 and p[-1]==0:
        p=p[:-1]
    return p or Z


def padd(a,b):
    return pol(*( (a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                  for i in range(max(len(a),len(b))) ))


def pneg(a):
    return tuple(-x for x in a)


def psub(a,b):
    return padd(a,pneg(b))


def pmul(a,b):
    v=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            v[i+j]+=x*y
    return pol(*v)


def pscale(a,x):
    return pol(*(F(x)*v for v in a))


def peval(p,t):
    z=F(0)
    for x in reversed(p):
        z=z*t+x
    return z


def bcoeff(p,a,b):
    """Exact degree-n Bernstein coefficients of p(a+(b-a)t), 0<=t<=1."""
    n=len(p)-1
    power=[]
    for k in range(n+1):
        power.append(sum(p[j]*comb(j,k)*a**(j-k)*(b-a)**k
                         for j in range(k,n+1)))
    return tuple(sum(power[j]*F(comb(i,j),comb(n,j))
                     for j in range(i+1)) for i in range(n+1))


def certify_sign(p,a,b,strict=False,max_depth=10):
    """Prove p>=0 (or >0) throughout [a,b] by Bernstein subdivision.

    If a polynomial crosses the negative axis, or the chosen depth cannot
    settle a sign, this procedure FAILS CLOSED with ValueError.
    """
    if not (F(0)<=a<b<=F(1)):
        raise ValueError("nonempty rational time interval inside [0,1] needed")
    if type(max_depth) is not int or max_depth<0:
        raise ValueError("nonnegative integer max_depth required")
    stack=[(a,b,0)]
    checked=0
    while stack:
        low,high,depth=stack.pop()
        checked+=1
        coeff=bcoeff(p,low,high)
        if all(x>0 for x in coeff) if strict else all(x>=0 for x in coeff):
            continue
        if depth>=max_depth:
            raise ValueError("polynomial Bernstein sign not certified")
        mid=(low+high)/2
        stack.extend(((low,mid,depth+1),(mid,high,depth+1)))
    return checked


def ppoint(x,y):
    return (pol(F(x)),pol(F(y)))


def ppoint_at(P,t):
    return pt(peval(P[0],t),peval(P[1],t))


def pvadd(a,b):
    return (padd(a[0],b[0]),padd(a[1],b[1]))


def pvsub(a,b):
    return (psub(a[0],b[0]),psub(a[1],b[1]))


def pvdot(a,b):
    return padd(pmul(a[0],b[0]),pmul(a[1],b[1]))


def pvcross(a,b):
    return psub(pmul(a[0],b[1]),pmul(a[1],b[0]))


def family_vertices():
    """Nine vertices: octagon with a moving split point in old edge 0.

    Old vertex 3 additionally moves in direction (1/20,-1/30), so the
    underlying polygon changes geometry, not merely parametrization.
    """
    base=[ppoint(p.x,p.y) for p in OCTAGON]
    moved=(padd(base[3][0],pol(0,F(1,20))),
           psub(base[3][1],pol(0,F(1,30))))
    base[3]=moved
    first=base[0]
    second=base[1]
    q=pol(ALPHA-F(1,4),HALF)
    split=pvadd(first,(pmul(q,psub(second[0],first[0])),
                         pmul(q,psub(second[1],first[1]))))
    return tuple((base[0],split,*base[1:]))


def signed_det4(rows):
    """Exact determinant of 4x4 polynomial matrix, no symbolic CAS."""
    answer=Z
    for perm in permutations(range(4)):
        inversions=sum(perm[i]>perm[j] for i in range(4)
                       for j in range(i+1,4))
        product=O
        for i,j in enumerate(perm):
            product=pmul(product,rows[i][j])
        answer=padd(answer,product if inversions%2==0 else pneg(product))
    return answer


def chart_polynomials(P,assignment):
    """Square line-incidence M(t) z(t) = b(t); Cramer rational functions."""
    if len(assignment)!=4 or len(set(assignment))!=4 or any(
        type(i) is not int or i<0 or i>=len(P) for i in assignment
    ):
        raise ValueError("four distinct ordered edges required")
    M=[];rhs=[]
    for k,i in enumerate(assignment):
        a=P[i]
        e=pvsub(P[(i+1)%len(P)],a)
        normal=(pneg(e[1]),e[0])
        x,y=CORNER_COEFFS[k]
        M.append(tuple(padd(pscale(normal[0],x[j]),
                             pscale(normal[1],y[j])) for j in range(4)))
        rhs.append(pvdot(normal,a))
    D=signed_det4(M)
    N=[]
    for j in range(4):
        C=tuple(tuple(rhs[i] if k==j else M[i][k]
                       for k in range(4)) for i in range(4))
        N.append(signed_det4(C))
    return D,tuple(N)


def root_at(D,N,t):
    divisor=peval(D,t)
    if divisor==0:
        raise ValueError("singular chart at the seam")
    return tuple(peval(n,t)/divisor for n in N)


def corner_numerators(N,k):
    x,y=CORNER_COEFFS[k]
    return (pol(*(0,)) if False else
            (sum_polys(pscale(N[j],x[j]) for j in range(4)),
             sum_polys(pscale(N[j],y[j]) for j in range(4))))


def sum_polys(items):
    result=Z
    for p in items:
        result=padd(result,p)
    return result


@dataclass(frozen=True)
class ChartCertificate:
    assignment: tuple
    left: F
    right: F
    denominator: tuple
    numerators: tuple
    sign_checks: int
    left_root: tuple
    right_root: tuple
    left_parameters: tuple
    right_parameters: tuple


def params_at(P,assignment,root,t):
    """Exact t=endpoint edge coordinates for the four square corners."""
    x0,y0,u,v=root
    corners=(pt(x0,y0),pt(x0+u,y0+v),pt(x0+u-v,y0+v+u),
             pt(x0-v,y0+u))
    result=[]
    for k,i in enumerate(assignment):
        a=ppoint_at(P[i],t)
        b=ppoint_at(P[(i+1)%len(P)],t)
        e=sub(b,a)
        result.append(dot(e,sub(corners[k],a))/dot(e,e))
    return tuple(result)


def certify_chart(P,assignment,left,right,side2_limit=S2_LIMIT):
    """Prove one chart is feasible at EVERY time in [left,right].

    Membership includes edge endpoints, critical for chart handoffs.
    The side lower bound is strict.  No tests at a finite time grid.
    """
    left=F(left);right=F(right);side2_limit=F(side2_limit)
    if not (0<=left<right<=1) or side2_limit<=0:
        raise ValueError("invalid chart interval or square cutoff")
    D,N=chart_polynomials(P,assignment)
    try:
        checks=certify_sign(D,left,right,strict=True)
        direction=1
    except ValueError:
        checks=certify_sign(pneg(D),left,right,strict=True)
        direction=-1
    side_sq_num=psub(padd(pmul(N[2],N[2]),pmul(N[3],N[3])),
                     pscale(pmul(D,D),side2_limit))
    checks+=certify_sign(side_sq_num,left,right,strict=True)
    for k,i in enumerate(assignment):
        a=P[i]
        e=pvsub(P[(i+1)%len(P)],a)
        nx,ny=corner_numerators(N,k)
        displacement=(psub(nx,pmul(a[0],D)),
                      psub(ny,pmul(a[1],D)))
        T=pvdot(e,displacement)
        ED=pmul(pvdot(e,e),D)
        # true edge parameter = T/(D*|e|²). Since sign(D) is fixed,
        # these two POLYNOMIAL sign certificates prove 0<=t<=1.
        checks+=certify_sign(pscale(T,direction),left,right)
        checks+=certify_sign(pscale(psub(ED,T),direction),left,right)
        checks+=certify_sign(pvdot(e,e),left,right,strict=True)
    lroot=root_at(D,N,left);rroot=root_at(D,N,right)
    return ChartCertificate(tuple(assignment),left,right,D,N,checks,
                            lroot,rroot,params_at(P,assignment,lroot,left),
                            params_at(P,assignment,rroot,right))


def prove_jordan_family(P):
    """Prove EVERY family member is simple, strictly convex up to the split.

    The original 8-vertex polygon changes only at old vertex 3.
    Non-adjacent original edges remain disjoint by a rational coordinate
    bounding-box separation of > 2delta, delta=1/20. Consecutive turns
    remain strictly positive by Bernstein coefficients. The added moving
    split vertex stays strictly between original vertices 0 and 1.
    """
    validate_jordan_polygon(OCTAGON)
    if len(P)!=9:
        raise ValueError("expected nine family vertices")
    old=[P[0],*P[2:]]
    assert len(old)==8
    delta=F(1,20)
    margin=None
    for i in range(8):
        a=OCTAGON[i];b=OCTAGON[(i+1)%8]
        for j in range(i+1,8):
            if j==i+1 or (i==0 and j==7):
                continue
            c=OCTAGON[j];d=OCTAGON[(j+1)%8]
            def sep(k):
                a1,a2=(a.x,b.x) if k==0 else (a.y,b.y)
                c1,c2=(c.x,d.x) if k==0 else (c.y,d.y)
                return max(min(a1,a2)-max(c1,c2),
                           min(c1,c2)-max(a1,a2))
            gap=max(sep(0),sep(1))
            if gap<=2*delta:
                raise ValueError("moving polygon lacks nonintersection margin")
            margin=gap if margin is None else min(margin,gap)
    checks=0
    for i in range(8):
        p=old[i]
        prev=old[(i-1)%8];nxt=old[(i+1)%8]
        turn=pvcross(pvsub(p,prev),pvsub(nxt,p))
        checks+=certify_sign(turn,F(0),F(1),strict=True)
        checks+=certify_sign(pvdot(pvsub(nxt,p),pvsub(nxt,p)),
                             F(0),F(1),strict=True)
        # All original moving vertices lie in radius-3 ball.
        # Affine endpoints inside convex ball imply entire interpolation
        # and the inserted point lie there.
        for t in (F(0),F(1)):
            v=ppoint_at(p,t)
            if dot(v,v)>=9:
                raise ValueError("radius-3 bounding disk claim fails")
    q=pol(ALPHA-F(1,4),HALF)
    assert 0<peval(q,F(0))<peval(q,F(1))<1
    # Verify supplied split vertex is exactly on the old first edge.
    assert P[1]==pvadd(P[0],(pmul(q,psub(P[2][0],P[0][0])),
                               pmul(q,psub(P[2][1],P[0][1]))))
    return {"minimal_nonadjacent_coordinate_separation":str(margin),
            "max_vertex_motion_Linf":str(delta),
            "turn_and_edge_sign_certificates":checks,
            "all_polygons_jordan":True,
            "all_polygons_inside_radius_3":True}


def certify_cover(P,intervals,side2_limit=S2_LIMIT):
    """Certify a FULL time cover by charts with a continuous common root.

    All intervals are closed. Adjacent roots must coincide exactly at
    seams; a time gap, jump in selected solution, singular chart or a
    sign-indeterminate polynomial fails closed.
    """
    if not intervals:
        raise ValueError("no chart intervals")
    ordered=sorted(intervals,key=lambda x:F(x[0]))
    if F(ordered[0][0])!=0 or F(ordered[-1][1])!=1:
        raise ValueError("cover must reach lambda=0 and lambda=1")
    certificates=[]
    for left,right,assignment in ordered:
        c=certify_chart(P,assignment,left,right,side2_limit)
        if certificates:
            prev=certificates[-1]
            if prev.right!=c.left:
                raise ValueError("gap or overlap in continuation atlas")
            if prev.right_root!=c.left_root:
                raise ValueError("square branch jumps across chart handoff")
        certificates.append(c)
    return tuple(certificates)


def audit():
    P=family_vertices()
    jordan=prove_jordan_family(P)
    low=(1,3,5,7)
    high=(0,3,5,7)
    atlas=certify_cover(P,((F(0),HALF,low),(HALF,F(1),high)))
    # The split is a genuine change of which polygon edge carries a
    # moving square corner, not a false handoff at an interior parameter.
    seam=atlas[0].right_root
    split_vertex=ppoint_at(P[1],HALF)
    assert pt(seam[0],seam[1])==split_vertex
    assert atlas[0].right_parameters[0]==0
    assert atlas[1].left_parameters[0]==1
    assert atlas[0].left_parameters[0]>0
    assert atlas[1].right_parameters[0]<1
    # The geometric square branch is in fact identical algebraically in
    # the two edge charts, not merely equal at the vertex seam.
    for j in range(4):
        assert pmul(atlas[0].numerators[j],atlas[1].denominator)==pmul(
            atlas[1].numerators[j],atlas[0].denominator)
    assert S2_LIMIT>F(27,10)**2-F(1,100) # same strict squared threshold
    return {
        "status":"PASS: exact two-chart continuation across a moving polygon vertex",
        "jordan_family":jordan,
        "chart_intervals":[{"left":str(c.left),"right":str(c.right),
                             "edges":list(c.assignment),
                             "sign_checks":c.sign_checks,
                             "left_edge_parameters":[str(v) for v in c.left_parameters],
                             "right_edge_parameters":[str(v) for v in c.right_parameters]}
                           for c in atlas],
        "handoff_lambda":str(HALF),
        "handoff_corner_on_moving_vertex":True,
        "geometric_square_branch_symbolically_identical_across_charts":True,
        "guaranteed_square_side_squared_strictly_above":str(S2_LIMIT),
        "restricted_family_big_peg_ratio_strictly_above":"9/10",
        "scope":"Restricted affine PL Jordan homotopy with a subdivision-vertex chart event; not general singular-orbit continuation or universal Toeplitz/Big Peg.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
