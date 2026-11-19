"""v223: exact COMPLETE nondegenerate square-orbit census for a rational Jordan polygon.

Includes repeated-edge vertices and edge endpoints, all curve-order sectors,
all eight geometric D4 labels, and rank-deficient polygon edge products.
A positive-dimensional square family is reported as NONISOLATED, never
silently discarded or interpreted as a finite orbit parity.

For a restricted family of genuine Jordan polygons obtained by a common
positive similarity scale plus translation, gives a complete exact ledger
of square-size cutoff crossings. NOT a general polygon-homotopy cobordism,
not a relative Floer theory, and not an arbitrary-Jordan Toeplitz proof.
Standard library only, exact Fraction arithmetic.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product, combinations
from math import isqrt
from typing import Tuple


@dataclass(frozen=True, order=True)
class Point:
    x: Q
    y: Q


def pt(x, y):
    return Point(Q(x), Q(y))


def sub(a,b):
    return pt(a.x-b.x,a.y-b.y)


def add(a,b):
    return pt(a.x+b.x,a.y+b.y)


def mul(a,s):
    return pt(a.x*Q(s),a.y*Q(s))


def dot(a,b):
    return a.x*b.x+a.y*b.y


def cross(a,b):
    return a.x*b.y-a.y*b.x


def orientation_area(P):
    return sum(cross(P[i],P[(i+1)%len(P)]) for i in range(len(P)))


def _on_seg(p,a,b):
    return cross(sub(p,a),sub(b,a))==0 and (
        min(a.x,b.x)<=p.x<=max(a.x,b.x) and min(a.y,b.y)<=p.y<=max(a.y,b.y))


def _intersect(a,b,c,d):
    c1,c2=cross(sub(b,a),sub(c,a)),cross(sub(b,a),sub(d,a))
    c3,c4=cross(sub(d,c),sub(a,c)),cross(sub(d,c),sub(b,c))
    if c1*c2<0 and c3*c4<0:return True
    return ((c1==0 and _on_seg(c,a,b)) or (c2==0 and _on_seg(d,a,b))
            or (c3==0 and _on_seg(a,c,d)) or (c4==0 and _on_seg(b,c,d)))


def validate_jordan_polygon(vertices):
    """Exact simplicity, including rejection of adjacent backtracking."""
    P=tuple(vertices);n=len(P)
    if n<3 or len(set(P))!=n or any(not isinstance(p,Point) for p in P):
        raise ValueError('need at least three distinct rational polygon vertices')
    if any(type(p.x) is not Q or type(p.y) is not Q for p in P):
        raise ValueError('coordinates must be Fraction')
    if any(P[i]==P[(i+1)%n] for i in range(n)):
        raise ValueError('zero polygon edge')
    if orientation_area(P)==0:raise ValueError('zero-area polygon')
    for i in range(n):
        prev,at,nxt=P[(i-1)%n],P[i],P[(i+1)%n]
        a,b=sub(at,prev),sub(nxt,at)
        if cross(a,b)==0 and dot(a,b)<0:
            raise ValueError('adjacent overlapping/backtracking edges')
    for i,j in combinations(range(n),2):
        if j==i+1 or (i==0 and j==n-1):continue
        if _intersect(P[i],P[(i+1)%n],P[j],P[(j+1)%n]):
            raise ValueError('nonadjacent polygon edges intersect')
    return P


def _rank(rows):
    if not rows:return 0
    m=[list(map(Q,row)) for row in rows]
    rank=0
    for j in range(len(m[0])):
        piv=next((k for k in range(rank,len(m)) if m[k][j]),None)
        if piv is None:continue
        m[rank],m[piv]=m[piv],m[rank]
        a=m[rank][j];m[rank]=[x/a for x in m[rank]]
        for k in range(rank+1,len(m)):
            if m[k][j]:
                fac=m[k][j]
                m[k]=[x-fac*y for x,y in zip(m[k],m[rank])]
        rank+=1
        if rank==len(m):break
    return rank


def _solve(A,b):
    n=len(b)
    if len(A)!=n or any(len(row)!=n for row in A):
        raise ValueError('square system required')
    m=[list(map(Q,row))+[Q(v)] for row,v in zip(A,b)]
    for i in range(n):
        piv=next((j for j in range(i,n) if m[j][i]),None)
        if piv is None:return None
        m[i],m[piv]=m[piv],m[i]
        v=m[i][i];m[i]=[x/v for x in m[i]]
        for j in range(n):
            if j!=i and m[j][i]:
                fac=m[j][i];m[j]=[x-fac*y for x,y in zip(m[j],m[i])]
    return tuple(row[-1] for row in m)


def _box_vertices(A,b):
    """Complete vertices of {Au=b,0<=u<=1}; includes rank deficiency."""
    r=_rank(A);rows=[];targets=[]
    for row,rhs in zip(A,b):
        if _rank(rows+[row])>len(rows):
            rows.append(list(row));targets.append(Q(rhs))
    seen=set()
    for ids in combinations(range(4),4-r):
        for vals in product((Q(0),Q(1)),repeat=4-r):
            M=[list(row) for row in rows];v=list(targets)
            for i,t in zip(ids,vals):
                u=[Q(0)]*4;u[i]=Q(1);M.append(u);v.append(t)
            x=_solve(M,v)
            if x is None or x in seen or any(t<0 or t>1 for t in x):continue
            if all(sum(q*t for q,t in zip(row,x))==rhs for row,rhs in zip(A,b)):
                seen.add(x);yield x


def _square_residual(P,sigma):
    p0,p1,p2,p3=P
    e=sub(p1,p0)
    J=pt(-sigma*e.y,sigma*e.x)
    z=sub(sub(p2,p1),J);w=sub(sub(p3,p0),J)
    return (z.x,z.y,w.x,w.y)


def _affine_system(edges,sigma):
    start=[a for a,b in edges]
    const=_square_residual(start,sigma)
    cols=[]
    for i in range(4):
        s=start.copy();s[i]=edges[i][1]
        nxt=_square_residual(s,sigma)
        cols.append(tuple(a-b for a,b in zip(nxt,const)))
    return [list(row) for row in zip(*cols)],tuple(-z for z in const)


def _evaluate(edges,us):
    return tuple(add(a,mul(sub(b,a),u)) for (a,b),u in zip(edges,us))


def _side_sq(pts):
    v=sub(pts[1],pts[0]);return dot(v,v)


def _vertex_location(i,u,n):
    return (Q(i)+u)%n

@dataclass(frozen=True)
class GeometricSquare:
    vertices: Tuple[Point,Point,Point,Point]
    curve_positions: Tuple[Q,Q,Q,Q]
    cyclic_order_sector: Tuple[int,int,int,int]
    graceful: bool
    side_squared: Q
    labelings_checked: int
    has_polygon_vertex: bool
    uses_repeated_edge: bool


@dataclass(frozen=True)
class PolygonCensus:
    squares: Tuple[GeometricSquare,...]
    cutoff_squared: Q
    checked_edge_charts: int
    singular_charts: int
    isolated_singular_charts: int
    redundant_labelings: int
    @property
    def count(self):return len(self.squares)
    @property
    def parity(self):return self.count%2
    @property
    def graceful_count(self):return sum(s.graceful for s in self.squares)
    @property
    def nongraceful_count(self):return self.count-self.graceful_count


class PositiveDimensionalSquares(ValueError):
    """A continuum of positive-size square orbits prevents finite parity."""
    def __init__(self,edge_tuple,orientation,vertex_count,max_side_sq,sample_square,sample_parameters):
        self.edge_tuple=edge_tuple;self.orientation=orientation
        self.vertex_count=vertex_count;self.max_side_sq=max_side_sq
        self.sample_square=sample_square;self.sample_parameters=sample_parameters
        super().__init__(f'positive-dimensional square solutions on edges {edge_tuple}, '
                         f'orientation {orientation}; cannot assign a finite orbit count')


def canonical_orbit(pts,positions):
    """Unique CCW geometry starting at lexicographically least vertex.

    The curve order sector is measured from THAT designated geometric
    vertex, so all eight labelings give the same one-sector representative.
    """
    pts=tuple(pts);positions=tuple(positions)
    key=tuple(sorted(pts))
    if len(set(key))!=4:raise ValueError('zero or repeated square vertex')
    chosen=None
    for start in range(4):
        if pts[start]!=key[0]:continue
        for direction in (1,-1):
            inds=tuple((start+direction*i)%4 for i in range(4))
            q=tuple(pts[i] for i in inds)
            if _square_residual(q,1)==(0,0,0,0):
                chosen=inds;break
        if chosen is not None:break
    if chosen is None:raise ValueError('not a geometric square')
    coords=tuple(pts[i] for i in chosen)
    locs=tuple(positions[i] for i in chosen)
    return coords,locs


def _can_share_edge(positions, n):
    """Whether any two geometrically distinct vertices lie on one polygon edge.

    Endpoints belong to their two incident edges. This flag is geometric
    and independent of the edge-chart choices made by the enumerator.
    """
    edge_sets=[]
    for t in positions:
        j=int(t)
        if t==j:
            edge_sets.append({(j-1)%n,j%n})
        else:
            edge_sets.append({j})
    return any(edge_sets[i]&edge_sets[j]
               for i,j in combinations(range(4),2))


def census_polygon(vertices,cutoff_squared=Q(0)):
    """Exact exhaustive square-ORBIT census on one rational Jordan polygon.

    Works for endpoint hits, repeated edges, and all cyclic-order sectors.
    If any affine feasible set has positive dimension and a square with
    side^2>cutoff, it RAISES rather than fabricate a finite parity.
    """
    P=validate_jordan_polygon(vertices);n=len(P);cutoff=Q(cutoff_squared)
    if cutoff<0:raise ValueError('nonnegative squared-size cutoff required')
    edges=tuple((P[i],P[(i+1)%n]) for i in range(n))
    found={};checked=singular=isolated_singular=0
    for ids in product(range(n),repeat=4):
        # Three points on the same straight edge cannot form a nonzero square.
        if any(ids.count(i)>=3 for i in set(ids)):continue
        selected=tuple(edges[i] for i in ids)
        for sigma in (1,-1):
            checked+=1
            A,b=_affine_system(selected,sigma)
            rank=_rank(A)
            if rank==4:
                x=_solve(A,b)
                verts=() if (x is None or any(u<0 or u>1 for u in x)) else (x,)
            else:
                singular+=1
                verts=tuple(_box_vertices(A,b))
                if len(verts)==1:isolated_singular+=1
                if len(verts)>1:
                    best=max(verts,key=lambda u:_side_sq(_evaluate(selected,u)))
                    sample=_evaluate(selected,best)
                    maxsq=_side_sq(sample)
                    if maxsq>cutoff:
                        assert _square_residual(sample,sigma)==(0,0,0,0)
                        raise PositiveDimensionalSquares(ids,sigma,len(verts),maxsq,
                                                         sample,best)
                    continue
            for u in verts:
                pts=_evaluate(selected,u);size=_side_sq(pts)
                if size<=cutoff:continue
                if _square_residual(pts,sigma)!=(0,0,0,0):
                    raise AssertionError('invalid algebraic square')
                positions=tuple(_vertex_location(i,t,n) for i,t in zip(ids,u))
                ccw,ccw_t=canonical_orbit(pts,positions)
                order=(0,)+tuple(sorted((1,2,3),key=lambda i:(ccw_t[i]-ccw_t[0])%n))
                graceful=order in ((0,1,2,3),(0,3,2,1))
                key=tuple(sorted(pts))
                info=GeometricSquare(ccw,ccw_t,order,graceful,size,0,
                                     any(t in (Q(i) for i in range(n)) for t in positions),
                                     _can_share_edge(positions,n))
                if key in found:
                    old,count=found[key]
                    if old.vertices!=info.vertices or old.curve_positions!=info.curve_positions:
                        raise AssertionError('same embedded square has different unique curve parameters')
                    found[key]=(old,count+1)
                else:found[key]=(info,1)
    squares=[]
    for info,count in found.values():
        squares.append(GeometricSquare(info.vertices,info.curve_positions,
                  info.cyclic_order_sector,info.graceful,info.side_squared,count,
                  info.has_polygon_vertex,info.uses_repeated_edge))
    squares.sort(key=lambda s:s.vertices)
    return PolygonCensus(tuple(squares),cutoff,checked,singular,
                         isolated_singular,sum(c-1 for _,c in found.values()))


@dataclass(frozen=True)
class CutoffOrbitEvent:
    # For an irrational time, exact_time=None and (lower,upper) is a
    # proven rational root isolating interval. No float substitution.
    exact_time: Q | None
    lower: Q
    upper: Q
    square_key: Tuple[Point,...]
    base_side_squared: Q


@dataclass(frozen=True)
class SimilarityLedger:
    start: Q
    end: Q
    cutoff_squared: Q
    initial_orbits: int
    final_orbits: int
    events: Tuple[CutoffOrbitEvent,...]
    total_unfiltered_orbits: int
    @property
    def balances(self):
        return (self.initial_orbits+self.final_orbits+len(self.events))%2==0
    @property
    def parity_change(self):return (self.initial_orbits+self.final_orbits)%2


def similarity_polygon(P,t,scale0=1,scale_velocity=0,translation0=(0,0),translation_velocity=(0,0)):
    """Single embedded polygon family by common positive scaling and translation."""
    t,s0,sv=Q(t),Q(scale0),Q(scale_velocity)
    z=s0+sv*t
    if z<=0:raise ValueError('similarity must be nonsingular positive')
    q0=pt(*translation0);v=pt(*translation_velocity)
    return tuple(add(add(mul(p,z),q0),mul(v,t)) for p in P)


def similarity_cutoff_ledger(P,start,end,cutoff_squared,scale0=1,scale_velocity=0,
                             translation0=(0,0),translation_velocity=(0,0)):
    """Complete geometric-orbit cutoff ledger for an embedded similarity homotopy.

    All geometric square orbits correspond bijectively at every time;
    uniform positive scale multiplies squared side by scale(t)^2.
    Thus no source/sector or birth-death events are possible in this
    restricted genuine-one-polygon family. Reject nongeneric endpoint hits.
    """
    a,b,k=Q(start),Q(end),Q(cutoff_squared)
    s0,sv=Q(scale0),Q(scale_velocity)
    if a>=b or k<=0:raise ValueError('ordered time and positive cutoff required')
    if min(s0+sv*a,s0+sv*b)<=0:raise ValueError('polygon must stay embedded')
    P=validate_jordan_polygon(P)
    base=census_polygon(P)
    events=[];counts=[]
    for t in (a,b):
        scale=s0+sv*t
        if any(s.side_squared*scale*scale==k for s in base.squares):
            raise ValueError('time slice lies on size cutoff')
        counts.append(sum(s.side_squared*scale*scale>k for s in base.squares))
    for square in base.squares:
        if sv==0:continue
        q=square.side_squared
        va=q*(s0+sv*a)**2-k
        vb=q*(s0+sv*b)**2-k
        if va*vb<0:
            # Unique crossing since scale(t)>0 and scale is strictly
            # monotone. Its exact time is (sqrt(k/q)-s0)/sv.
            rad=k/q
            rootnum=isqrt(rad.numerator)
            rootden=isqrt(rad.denominator)
            exact=None
            if rootnum*rootnum==rad.numerator and rootden*rootden==rad.denominator:
                exact=(Q(rootnum,rootden)-s0)/sv
                if not a<exact<b:raise AssertionError('exact root outside interval')
                lo=hi=exact
            else:
                # Certified rational isolating interval, width <= 2^-40.
                lo,hi=a,b
                while hi-lo>Q(1,2**40):
                    mid=(lo+hi)/2
                    vm=q*(s0+sv*mid)**2-k
                    if vm==0:
                        exact=mid;lo=hi=mid;break
                    if va*vm<0:hi=mid
                    else:lo=mid;va=vm
                if exact is None:
                    assert lo<hi and (q*(s0+sv*lo)**2-k)*(q*(s0+sv*hi)**2-k)<0
            events.append(CutoffOrbitEvent(exact,lo,hi,square.vertices,q))
    events.sort(key=lambda e:(e.lower,e.upper,e.square_key))
    result=SimilarityLedger(a,b,k,counts[0],counts[1],tuple(events),base.count)
    if not result.balances:
        raise AssertionError('size-incidence parity failed')
    return result



def convex_vertex_alignment_polygon(t):
    """A strictly convex single-polygon deformation with a singular vertex event.

    For 0<=t<4/3 the consecutive turning cross products are
    6-3t, (3t^2-13t+16)/2, (t-3)(3t-4)/2, 4-3t,
    all strictly positive. At t=1 one additional square lies exactly
    on a polygon vertex; it does not define a generic mod-2 parity.
    """
    t=Q(t)
    if not Q(0)<=t<Q(4,3):
        raise ValueError('model requires 0<=t<4/3')
    P=(pt(0,0),pt(2,0),pt(3-t,3-Q(3,2)*t),pt(0,2-Q(3,2)*t))
    for i in range(4):
        x=sub(P[(i+1)%4],P[i]);y=sub(P[(i+2)%4],P[(i+1)%4])
        if cross(x,y)<=0:
            raise AssertionError('derived convexity condition failed')
    return validate_jordan_polygon(P)


def vertex_alignment_snapshots():
    """Exact separate-time counts; NOT a complete homotopy event ledger."""
    times=(Q(0),Q(63,64),Q(1),Q(33,32))
    return tuple((t,census_polygon(convex_vertex_alignment_polygon(t)))
                 for t in times)


def report():
    quad=tuple(pt(*xy) for xy in ((0,0),(2,0),(3,3),(0,2)))
    result=census_polygon(quad)
    snaps=vertex_alignment_snapshots()
    return {'polygon_vertices':len(quad), 'geometric_square_orbits':result.count,
            'labeling_duplications':result.redundant_labelings,
            'graceful_orbits':result.graceful_count,
            'nongraceful_orbits':result.nongraceful_count,
            'square_sizes_squared':tuple(map(str,(s.side_squared for s in result.squares))),
            'convex_vertex_alignment_snapshots':tuple((str(t),c.count) for t,c in snaps),
            'scope':'finite rational Jordan polygon; not arbitrary C0 Jordan or full generic homotopy'}

if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
