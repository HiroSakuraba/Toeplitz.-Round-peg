"""v216: boundary-only exact Brouwer degree for four arbitrary polygonal arcs.

All calculations use rational fractions. This is a conditional, finite PL
certificate on FOUR ALREADY CHOSEN arcs with fixed physical lifts, not a
universal square-peg result or an implementation of Floer homology.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product, combinations
from math import isqrt
from typing import Tuple


@dataclass(frozen=True)
class Point:
    x: Q
    y: Q


def pt(x,y):
    return Point(Q(x),Q(y))


def add(a,b):
    return pt(a.x+b.x,a.y+b.y)


def sub(a,b):
    return pt(a.x-b.x,a.y-b.y)


def mul(t,a):
    return pt(t*a.x,t*a.y)


def dot(a,b):
    return a.x*b.x+a.y*b.y


def cross(a,b):
    return a.x*b.y-a.y*b.x


@dataclass(frozen=True)
class Arc:
    vertices: Tuple[Point, ...]

    def __post_init__(self):
        if len(self.vertices)<2:
            raise ValueError('An arc needs at least two vertices')
        if any(a == b for a,b in zip(self.vertices,self.vertices[1:])):
            raise ValueError('Zero-length consecutive polygon edges')
        if any(type(p.x) is not Q or type(p.y) is not Q for p in self.vertices):
            raise ValueError('Use Fraction coordinates')


def arc(*coords):
    return Arc(tuple(pt(*xy) for xy in coords))


def residual(points, orientation):
    p0,p1,p2,p3=points
    e=sub(p1,p0)
    # p2-p1 - sigma J(e), p3-p0 - sigma J(e)
    j=pt(-orientation*e.y,orientation*e.x)
    f=sub(sub(p2,p1),j)
    g=sub(sub(p3,p0),j)
    return (f.x,f.y,g.x,g.y)


def _rank(a):
    m=[list(map(Q,row)) for row in a]
    if not m: return 0
    n=len(m[0]); k=0
    for j in range(n):
        t=next((r for r in range(k,len(m)) if m[r][j]),None)
        if t is None: continue
        m[k],m[t]=m[t],m[k]
        z=m[k][j];m[k]=[v/z for v in m[k]]
        for i in range(k+1,len(m)):
            if m[i][j]:
                z=m[i][j];m[i]=[x-z*y for x,y in zip(m[i],m[k])]
        k+=1
        if k==len(m): break
    return k


def _solve(A,b):
    n=len(b)
    if len(A)!=n or any(len(row)!=n for row in A):
        raise ValueError('Expected square matrix')
    m=[list(map(Q,row))+[Q(t)] for row,t in zip(A,b)]
    for col in range(n):
        k=next((i for i in range(col,n) if m[i][col]),None)
        if k is None:return None
        m[col],m[k]=m[k],m[col]
        p=m[col][col];m[col]=[v/p for v in m[col]]
        for i in range(n):
            if i!=col and m[i][col]:
                z=m[i][col];m[i]=[x-z*y for x,y in zip(m[i],m[col])]
    return tuple(row[-1] for row in m)


def _det(A):
    n=len(A);m=[list(map(Q,r)) for r in A];value=Q(1)
    for col in range(n):
        k=next((i for i in range(col,n) if m[i][col]),None)
        if k is None:return Q(0)
        if k!=col:
            m[col],m[k]=m[k],m[col];value=-value
        p=m[col][col];value*=p
        for i in range(col+1,n):
            z=m[i][col]/p
            for j in range(col+1,n):m[i][j]-=z*m[col][j]
            m[i][col]=Q(0)
    return value


def feasible_box(A,b,n):
    """Exact nonempty affine slice of [0,1]^n, including singular systems."""
    if any(len(row)!=n for row in A) or len(A)!=len(b):
        raise ValueError('Inconsistent dimensions')
    rank=_rank(A);rows=[];rhs=[]
    for row,v in zip(A,b):
        if _rank(rows+[row])>len(rows):
            rows.append(list(row));rhs.append(Q(v))
    for js in combinations(range(n),n-rank):
        for vals in product((Q(0),Q(1)),repeat=n-rank):
            matrix=[list(row) for row in rows];target=list(rhs)
            for j,val in zip(js,vals):
                unit=[Q(0)]*n;unit[j]=Q(1);matrix.append(unit);target.append(val)
            sol=_solve(matrix,target)
            if sol is not None and all(Q(0)<=u<=Q(1) for u in sol) and all(
                sum(x*y for x,y in zip(row,sol))==v for row,v in zip(A,b)):
                return True
    return False


def polygon_segments(a):
    return tuple((u,v) for u,v in zip(a.vertices,a.vertices[1:]))


def _affine_pieces(selected, shifts, orientation):
    starts=[pt(a.x+k,a.y) for (a,_),k in zip(selected,shifts)]
    d=[sub(b,a) for a,b in selected]
    base=residual(starts,orientation)
    cols=[]
    for j in range(4):
        arr=list(starts)
        arr[j]=add(arr[j],d[j])
        nxt=residual(arr,orientation)
        cols.append(tuple(v-w for v,w in zip(nxt,base)))
    return base,tuple(cols)


def _facet_pieces(arcs,shifts,orientation):
    segs=[polygon_segments(a) for a in arcs]
    for fixed in range(4):
        for side in (0,1):
            endpoint=arcs[fixed].vertices[0 if side==0 else -1]
            choices=[((endpoint,endpoint),) if j==fixed else segs[j] for j in range(4)]
            for selected in product(*choices):
                base, cols =_affine_pieces(selected,shifts,orientation)
                free=tuple(j for j in range(4) if j!=fixed)
                yield fixed,side,base,tuple(cols[j] for j in free)


def _base_free_to_matrix(base,cols):
    A=[[cols[j][i] for j in range(len(cols))] for i in range(4)]
    b=[-v for v in base]
    return A,b


def _face_bound(base,cols,max_nodes):
    """Positive exact lower bound on ||F||_inf for a zero-free affine facet."""
    stack=[((Q(0),Q(1)),)*3];nodes=0;minmargin=None
    while stack:
        box=stack.pop();nodes+=1
        if nodes>max_nodes:raise DegreeInconclusive('Boundary subdivision budget exceeded')
        margins=[]
        for row in range(4):
            lo=base[row];hi=base[row]
            for j,(a,b) in enumerate(box):
                x,y=cols[j][row]*a,cols[j][row]*b
                lo+=min(x,y);hi+=max(x,y)
            if lo>0:margins.append(lo)
            if hi<0:margins.append(-hi)
        if margins:
            m=max(margins)
            minmargin=m if minmargin is None else min(minmargin,m)
            continue
        sizes=[(b-a)*max(abs(cols[j][i]) for i in range(4)) for j,(a,b) in enumerate(box)]
        j=max(range(3),key=lambda k:sizes[k]);a,b=box[j]
        if sizes[j]==0 or a==b:raise DegreeInconclusive('Unable to separate facet from zero')
        mid=(a+b)/2
        x=list(box);y=list(box);x[j]=(a,mid);y[j]=(mid,b)
        stack.extend((tuple(x),tuple(y)))
    return minmargin,nodes


def _unit_direction_candidates(count):
    # Moment-curve candidates. When every candidate fails, report inconclusive.
    for n in range(1,count+1):
        yield (Q(1),Q(n),Q(n*n),Q(n*n*n))


class BoundaryRoot(ValueError):
    pass

class DegreeInconclusive(RuntimeError):
    pass


@dataclass(frozen=True)
class DegreeCertificate:
    degree: int
    ray: Tuple[Q,Q,Q,Q]
    intersections: int
    boundary_margin: Q
    opposite_distance_sq: Q
    zero_free_facets: int
    nodes: int
    attempts: int

    @property
    def exists_nondegenerate_square(self):
        return self.degree!=0 and self.opposite_distance_sq>0

    @property
    def perturbation_bound_linf(self):
        # Under ||Δp_i||∞<=eps, each residual coordinate changes <=4eps.
        # Need 4eps < margin; opposite-arc separation squared D².
        # A rational L<=sqrt(D²) gives D>=L, and each endpoint moves <=sqrt2 eps,
        # hence distance remains >0 for eps<L/4.
        L=_rational_sqrt_floor(self.opposite_distance_sq)
        return min(self.boundary_margin/4,L/4)


def _rational_sqrt_floor(q):
    if q<=0:return Q(0)
    # Floor sqrt(q) at 2^-20 resolution, then fall back to 2^-k for tiny q.
    scale=1<<20
    val=isqrt(q.numerator*scale*scale//q.denominator)
    if val: return Q(val,scale)
    k=20
    while Q(1,4**k)>q:k+=1
    return Q(1,2**k)


def segment_distance_squared(a,b,c,d):
    """Minimum squared Euclidean distance between two closed segments (exact)."""
    ab=sub(b,a);cd=sub(d,c)
    def point_to_segment_sq(p,u,v):
        w=sub(v,u);h=dot(w,w)
        t=max(Q(0),min(Q(1),dot(sub(p,u),w)/h))
        diff=sub(p,add(u,mul(t,w)))
        return dot(diff,diff)
    if cross(ab,cd):
        t=cross(sub(c,a),cd)/cross(ab,cd)
        s=cross(sub(c,a),ab)/cross(ab,cd)
        if 0<=t<=1 and 0<=s<=1:return Q(0)
    return min(point_to_segment_sq(p,u,v) for p,u,v in (
        (a,c,d),(b,c,d),(c,a,b),(d,a,b)))


def arc_distance_squared(a,b):
    return min(segment_distance_squared(c,d,e,f) for (c,d),(e,f) in product(polygon_segments(a),polygon_segments(b)))


def boundary_degree(arcs,shifts=(0,0,0,0),orientation=1,
                    max_rays=128,max_boundary_nodes=200000):
    """Boundary-only exact degree. No interior square enumeration is used.

    Raises BoundaryRoot for a source-window boundary zero; DegreeInconclusive
    when the selected rational ray is nontransverse or a resource cap hits.
    """
    if len(arcs)!=4 or any(not isinstance(a,Arc) for a in arcs):
        raise ValueError('Four Arc objects are required')
    if len(shifts)!=4 or shifts[0]!=0 or any(type(k) is not int for k in shifts):
        raise ValueError('Fixed relative physical lifts must be 4 integers, k0=0')
    if orientation not in (1,-1):raise ValueError('Orientation must be +1 or -1')
    lifted_opposite=Arc(tuple(pt(p.x+shifts[2]-shifts[0],p.y) for p in arcs[2].vertices))
    D2=arc_distance_squared(arcs[0],lifted_opposite)
    if D2<=0:raise ValueError('Opposite arcs 0 and 2 need positive physical separation')
    facets=list(_facet_pieces(arcs,shifts,orientation))
    total_nodes=0;margins=[]
    for fixed,side,base,cols in facets:
        A,b=_base_free_to_matrix(base,cols)
        if feasible_box(A,b,3):raise BoundaryRoot('Zero on external source-arc boundary')
    for fixed,side,base,cols in facets:
        m,n=_face_bound(base,cols,max_boundary_nodes-total_nodes)
        total_nodes+=n
        if m is None or m<=0:raise DegreeInconclusive('No positive boundary margin')
        margins.append(m)
    boundary_margin=min(margins)
    for attempt,ray in enumerate(_unit_direction_candidates(max_rays),start=1):
        degree=0;hits=0;generic=True
        for fixed,side,base,cols in facets:
            # Solve F(u)=lambda*ray, lambda>=0.
            M=[[cols[j][i] for j in range(3)]+[-ray[i]] for i in range(4)]
            rhs=[-v for v in base]
            det=_det(M)
            if det:
                sol=_solve(M,rhs)
                u=sol[:3];lam=sol[3]
                if lam<=0 or any(v<0 or v>1 for v in u):continue
                if any(v==0 or v==1 for v in u):
                    generic=False;break
                # Boundary orientation defined by (outward normal, tangent basis).
                face_sign=(-1)**(fixed if side==1 else fixed+1)
                degree+=face_sign*(1 if det>0 else -1)
                hits+=1
            else:
                # A singular affine span may miss the patch. Certify that
                # before declaring the ray unusable.
                B=max(abs(base[i]+sum(cols[j][i]*u[j] for j in range(3)))
                      for u in product((Q(0),Q(1)),repeat=3) for i in range(4))
                L=B/max(abs(t) for t in ray)
                if L<=0:raise BoundaryRoot('Zero residual on entire facet')
                scaled=[[cols[j][i] for j in range(3)]+[-L*ray[i]] for i in range(4)]
                if feasible_box(scaled,rhs,4):
                    generic=False;break
        if generic:
            return DegreeCertificate(degree,ray,hits,boundary_margin,D2,len(facets),total_nodes,attempt)
    raise DegreeInconclusive('No transverse ray found under budget')


def _mid(a,b):
    return pt((a.x+b.x)/2,(a.y+b.y)/2)


def vertex_neighborhood(polygon,i):
    """The two-half-edge marked arc around vertex i of a closed polygon."""
    P=tuple(polygon);n=len(P)
    return Arc((_mid(P[(i-1)%n],P[i]),P[i],_mid(P[i],P[(i+1)%n])))


def candidate_vertex_arcs(polygon):
    """A priori four-vertex-neighborhood choices, no square preselection."""
    P=tuple(polygon)
    if len(P)<4:raise ValueError('At least four polygon vertices')
    for inds in combinations(range(len(P)),4):
        yield inds,tuple(vertex_neighborhood(P,i) for i in inds)


def candidate_edge_arcs(polygon):
    """A priori choices of four whole polygon edges."""
    P=tuple(polygon)
    if len(P)<4:raise ValueError('At least four polygon vertices')
    for inds in combinations(range(len(P)),4):
        yield inds,tuple(Arc((P[i],P[(i+1)%len(P)])) for i in inds)


def _on_segment(a,b,p):
    return cross(sub(b,a),sub(p,a))==0 and (
        min(a.x,b.x)<=p.x<=max(a.x,b.x) and
        min(a.y,b.y)<=p.y<=max(a.y,b.y))


def _intersect_segments(a,b,c,d):
    o1=cross(sub(b,a),sub(c,a));o2=cross(sub(b,a),sub(d,a))
    o3=cross(sub(d,c),sub(a,c));o4=cross(sub(d,c),sub(b,c))
    if o1*o2<0 and o3*o4<0:return True
    return ((o1==0 and _on_segment(a,b,c)) or
            (o2==0 and _on_segment(a,b,d)) or
            (o3==0 and _on_segment(c,d,a)) or
            (o4==0 and _on_segment(c,d,b)))


def validate_simple_polygon(polygon):
    """Reject repeated vertices and nonadjacent edge intersections exactly.

    Collinear adjacent edges are allowed when there is no backtracking;
    adjacent edge backtracking is rejected as an overlapping segment.
    """
    P=tuple(polygon);n=len(P)
    if n<3 or any(type(p.x) is not Q or type(p.y) is not Q for p in P):
        raise ValueError('Need at least three rational polygon vertices')
    if len(set(P))!=n:raise ValueError('Repeated polygon vertex')
    if any(P[i]==P[(i+1)%n] for i in range(n)):
        raise ValueError('Zero-length polygon edge')
    twice_area=sum(cross(P[i],P[(i+1)%n]) for i in range(n))
    if not twice_area:raise ValueError('Zero polygon area')
    for i,j in combinations(range(n),2):
        a,b=P[i],P[(i+1)%n]
        c,d=P[j],P[(j+1)%n]
        if (j-i)%n in (1,n-1):
            # If adjacent and collinear, non-common endpoints must not
            # lie on the opposing segment; this excludes backtracking.
            if (cross(sub(b,a),sub(d,c))==0 and
                ( _on_segment(a,b,d) and d not in (a,b) or
                  _on_segment(c,d,a) and a not in (c,d))):
                raise ValueError('Adjacent edges overlap')
            continue
        if _intersect_segments(a,b,c,d):
            raise ValueError('Self-intersecting polygon')
    return twice_area


def select_certified_four_arcs(polygon,max_candidates=300,max_rays=32):
    """Finite, deterministic polygon source-window scan.

    The four selected subarcs always belong to the same given polygon,
    with physical lift tuple (0,0,0,0). Vertex and edge windows are
    selected without knowing a square. A negative answer is INCONCLUSIVE.
    """
    P=tuple(polygon)
    validate_simple_polygon(P)
    if len(P)<4:raise ValueError('Four polygon edges required')
    attempts=0;errors=[]
    for family,choices in (('vertex',candidate_vertex_arcs(P)),
                           ('edge',candidate_edge_arcs(P))):
        for inds,arcs in choices:
            for sigma in (1,-1):
                if attempts>=max_candidates:
                    return None,attempts,tuple(errors)
                attempts+=1
                try:
                    c=boundary_degree(arcs,orientation=sigma,max_rays=max_rays)
                    if c.degree!=0:return (family,inds,sigma,c),attempts,tuple(errors)
                except (BoundaryRoot,DegreeInconclusive,ValueError) as exc:
                    if len(errors)<10:errors.append((family,inds,sigma,type(exc).__name__))
    return None,attempts,tuple(errors)


if __name__=='__main__':
    p0=pt(Q(1,2),0);e=pt(Q(1,8),Q(1,8));r=pt(-e.y,e.x)
    ps=(p0,add(p0,e),add(add(p0,e),r),add(p0,r))
    slopes=(Q(-2),Q(3),Q(7,2),Q(-1,3));w=Q(1,32)
    arcs=tuple(Arc((pt(p.x-w,p.y-m*w),pt(p.x+w,p.y+m*w))) for p,m in zip(ps,slopes))
    cert=boundary_degree(arcs)
    print('boundary-only signed degree:',cert.degree)
    print('strict position-perturbation bound:',cert.perturbation_bound_linf)
    print('ray intersection count:',cert.intersections)
    polygon=tuple(pt(*pair) for pair in ((0,0),(3,1),(2,3),(-1,2)))
    witness,attempts,_=select_certified_four_arcs(polygon)
    print('a priori polygon scan:',witness[:3] if witness else 'INCONCLUSIVE')
    print('polygon windows tried:',attempts)