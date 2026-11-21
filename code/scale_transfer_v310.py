"""Exact refinement-loss bounds and a convex discontinuity test for square parity."""
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
import json

from polygon_orbit_census_v223 import (
    _affine_system, _evaluate, _on_seg, _rank, _solve, _square_residual, add,
    census_polygon, dot, mul, pt, sub, validate_jordan_polygon,
)
from tilted_sweep_v301 import (SAMPLE_TIMES, audit_charts, corner_x, determinant, seed, strict_sign, winning_tags)


def stretch(d):
    return tuple(pt(p.x,(1+d)*p.y) for p in seed(0))


def stretch_census():
    count = 0
    for chart in combinations_with_replacement(range(4),4):
        if len(set(chart))<3:
            continue
        samples = []
        for t in SAMPLE_TIMES:
            p = stretch(t/100)
            a,rhs = _affine_system([(p[i],p[(i+1)%4]) for i in chart],1)
            row = [determinant(a)]
            for j in range(4):
                changed = [list(x) for x in a]
                for line,z in zip(changed,rhs):
                    line[j]=z
                row.append(determinant(changed))
            samples.append(row)
        den = strict_sign([r[0] for r in samples])
        faces = [den*strict_sign(values) for i in range(1,5) for values in
                 ([r[i] for r in samples],[r[0]-r[i] for r in samples])]
        winner = len(set(chart)) == 4
        certified = min(faces) == 1 if winner else -1 in faces
        if not den or not certified:
            raise ArithmeticError("unresolved stretch chart",chart)
        count += 1
    return count


def audit_switch():
    # Clear the positive denominator before checking the degree-four polynomials.
    rows = []
    for t in SAMPLE_TIMES:
        d = t/100
        den = 17 + 12*d + 4*d*d
        xs = corner_x(d)
        a, _, _, e = winning_tags(d)
        vx = xs[1] - xs[0]
        vy = -3*(1-a)/8 - 3*(1-e)/4 + d*vx
        lam = 1+d
        p = stretch(d)
        x = 4/(8+9*lam)
        tags = (2*x, 1-2*x, 2*x, 1-2*x)
        points = _evaluate(list(zip(p,p[1:]+p[:1])), tags)
        v = sub(points[1],points[0])
        if any(_square_residual(points,1)) or dot(v,v)!=(9*lam/(8+9*lam))**2:
            raise ArithmeticError("incorrect vertical-stretch witness")
        rows.append((den*(vx-Q(21,40)), den*(-vy-Q(43,125)),
                     den*(Q(47,100)-(xs[0]+xs[2])/2),
                     Q(3,10)*(8+9*lam)**2-(9*lam)**2))
    if any(strict_sign(col)!=1 for col in zip(*rows)):
        raise ArithmeticError("uncertified threshold")
    if Q(21,40)**2 + Q(43,125)**2 <= Q(39,100):
        raise ArithmeticError("insufficient shear size margin")

    # Each corner is affine in t, each squared side quadratic.
    p = seed(0)
    edges = list(zip(p,p[1:]+p[:1]))
    for t in (Q(5,17),Q(7,17),Q(9,17)):
        tags = (2*t-Q(10,17),2*t-Q(9,17),t-Q(1,17),t)
        points = _evaluate(edges,tags)
        v = sub(points[1],points[0])
        if any(_square_residual(points,1)) or dot(v,v)!=Q(81,289)+Q(9,4)*(t-Q(9,17))**2:
            raise ArithmeticError("incorrect limiting square family")


def parity_intervals(sizes, window):
    lo,hi = window
    cuts = sorted({lo,hi,*(z for z in sizes if lo<z<hi)})
    return [(a,b) for a,b in zip(cuts,cuts[1:]) if sum(z>(a+b)/2 for z in sizes)%2]


def merged(intervals):
    result = []
    for a,b in sorted(intervals):
        if not result or a>result[-1][1]:
            result.append((a,b))
        else:
            result[-1]=(result[-1][0],max(b,result[-1][1]))
    return result


def normal(edge):
    a,b=edge
    n=(b.y-a.y,a.x-b.x)
    if n==(0,0):
        raise ValueError("zero edge")
    return n,n[0]*a.x+n[1]*a.y


def triple_bound(lines, anchor, eps):
    # Images of the side vector at the next, previous and opposite corner.
    transformed = []
    for i,((x,y),h) in enumerate(lines):
        m=((x,y),(y,-x),(x+y,y-x))[i]
        transformed.append((m,h-x*anchor.x-y*anchor.y,(x,y)))
    bounds=[]
    for i,j in combinations(range(3),2):
        mi,bi,ni=transformed[i]
        mj,bj,nj=transformed[j]
        A=(mi,mj)
        center=_solve(A,(bi,bj))
        if center is None:
            continue
        dx=_solve(A,(-ni[0],-nj[0]))
        dy=_solve(A,(-ni[1],-nj[1]))
        lo,hi=Q(0),Q(0)
        for c,x,y in zip(center,dx,dy):
            spread=eps*(abs(x)+abs(y))
            lo+=max(0,abs(c)-spread)**2
            hi+=(abs(c)+spread)**2
        bounds.append((lo,hi))
    if not bounds:
        return None
    return max(a for a,b in bounds),min(b for a,b in bounds)


def flux_cover(edges,anchor,eps,window):
    if eps<0 or not 0<=window[0]<window[1]:
        raise ValueError("invalid patch or window")
    lo,hi=window
    intervals=[(lo,min(hi,8*eps**2))] if lo<=8*eps**2 else []
    for triple in product(map(normal,edges),repeat=3):
        bounds=triple_bound(triple,anchor,eps)
        lo,hi=window if bounds is None else bounds
        lo,hi=max(lo,window[0]),min(hi,window[1])
        if lo<=hi:
            intervals.append((lo,hi))
    return merged(intervals)


def multi_patch_cover(edges, boxes, window):
    """Include squares with one edited corner and all squares with two or more."""
    edges = tuple(edges)
    boxes = tuple(boxes)
    intervals = [interval for anchor, eps in boxes
                 for interval in flux_cover(edges, anchor, eps, window)]
    for (p, e), (q, f) in combinations(boxes, 2):
        # The two corners can be adjacent or opposite. Their coordinate
        # differences are within e+f of the differences of the box centers.
        offsets = (abs(p.x-q.x), abs(p.y-q.y))
        low = sum(max(0, x-e-f)**2 for x in offsets)
        high = sum((x+e+f)**2 for x in offsets)
        for divisor in (1, 2):
            a = max(window[0], Q(low, divisor))
            b = min(window[1], Q(high, divisor))
            if a <= b:
                intervals.append((a, b))
    return merged(intervals)


def displaced_contact_example(triangle, square, eps):
    """Cut out four old contacts; solve the four new middle segments exactly."""
    edges = list(zip(triangle, triangle[1:]+triangle[:1]))
    replacements = [[] for edge in edges]
    middle = []
    for index, corner in enumerate(square.vertices, 1):
        owner, = [i for i,(a,b) in enumerate(edges) if _on_seg(corner,a,b)]
        a,b = edges[owner]
        direction = sub(b,a)
        step = mul(direction, eps/(4*max(abs(direction.x),abs(direction.y))))
        shift = pt(0,index*eps/128)
        points = [sub(corner,mul(step,2)), add(sub(corner,step),shift),
                  add(add(corner,step),shift), add(corner,mul(step,2))]
        replacements[owner].append((dot(sub(corner,a),direction),points))
        middle.append(tuple(points[1:3]))
    polygon = []
    for (a,_), groups in zip(edges,replacements):
        polygon.append(a)
        for _,points in sorted(groups):
            polygon.extend(points)
    polygon = validate_jordan_polygon(polygon)
    boundary = list(zip(polygon,polygon[1:]+polygon[:1]))
    if any(_on_seg(c,a,b) for c in square.vertices for a,b in boundary):
        raise ArithmeticError("an original contact survived")
    tags = _solve(*_affine_system(middle,1))
    if tags is None or not all(0<t<1 for t in tags):
        raise ArithmeticError("new contacts escaped their segments")
    corners = _evaluate(middle,tags)
    v = sub(corners[1],corners[0])
    if any(_square_residual(corners,1)) or dot(v,v)==square.side_squared:
        raise ArithmeticError("replacement square failed its geometric check")
    return dot(v,v)


def integration_example():
    p=(pt(0,0),pt(4,0),pt(1,3))
    edges=list(zip(p,p[1:]+p[:1]))
    squares=census_polygon(p).squares
    for square in squares:
        chart=[edges[int(t)] for t in square.curve_positions]
        if square.has_polygon_vertex or _rank(_affine_system(chart,1)[0])!=4:
            raise ArithmeticError("singular triangle square")
    sizes=sorted(s.side_squared for s in squares)
    window=(Q(1),Q(4))
    initial=parity_intervals(sizes,window)
    eps=Q(1,4096)
    cover=flux_cover(edges,pt(Q(4,7),Q(12,7)),eps,window)
    initial_mass=sum(b-a for a,b in initial)
    loss=sum(b-a for i,j in initial for u,v in cover
             for a,b in [(max(i,u),min(j,v))] if a<b)
    if not 8*eps**2<window[0] or initial_mass-loss<=0:
        raise ArithmeticError("patch budget exhausted")
    lower,upper=Q(29,10),Q(59,20)
    if [s for s in sizes if lower<s<upper]!=[Q(144,49)]:
        raise ArithmeticError("band does not isolate the square at the patch")
    if len(cover)!=1 or not lower<cover[0][0]<=Q(144,49)<=cover[0][1]<upper:
        raise ArithmeticError("possible scale transfer reaches a band boundary")
    axis_square, = [s for s in squares if s.side_squared == Q(144,49)]
    boxes = [(corner,eps) for corner in axis_square.vertices]
    multi_cover = multi_patch_cover(edges,boxes,(lower,upper))
    if len(multi_cover)!=1 or not lower<multi_cover[0][0]<=Q(144,49)<=multi_cover[0][1]<upper:
        raise ArithmeticError("four-patch transfer reaches a band boundary")
    moved = displaced_contact_example(p,axis_square,eps)
    if not multi_cover[0][0]<=moved<=multi_cover[0][1]:
        raise ArithmeticError("replacement square escaped its enclosure")
    return {"four_contact_replacement_squared_side":str(moved),
            "four_patch_centers": [[str(p.x),str(p.y)] for p,e in boxes],
            "four_patch_squared_side_enclosure": list(map(str,multi_cover[0])),
            "isolating_band": [str(lower),str(upper)],"initial_squared_sides":list(map(str,sizes)),
            "possible_flux_intervals":[list(map(str,row)) for row in cover],
            "initial_odd_measure":str(initial_mass),"lost_measure_upper":str(loss),
            "squared_side_floor":str(window[0]+initial_mass-loss)}


def report():
    audit_switch()
    return {"shear_charts":len(audit_charts()),"stretch_charts":stretch_census(),
            "limiting_squared_sides":["117/289","81/289"],
            "uniform_switch_interval":["3/10","39/100"],
            "patch_example":integration_example()}


if __name__=="__main__":
    print(json.dumps(report(),indent=2))
