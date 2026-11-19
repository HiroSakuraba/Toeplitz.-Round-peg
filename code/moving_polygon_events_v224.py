"""v224: exact algebraic event CANDIDATES for linearly moving Jordan polygons.

On one four-edge product chart p_i(t,u_i)=a_i(t)+u_i d_i(t),
F_sigma(t,u)=A(t)u-b(t). Cramer's rule gives u_i=N_i(t)/D(t).
Any change to the isolated in-chart square status, or to its positive
squared-side cutoff, must occur at a root of finitely many exact integer/
rational polynomials. Sturm sequences isolate all real roots, including
irrational and even-multiplicity roots. Singular-in-all-time charts and
persistent endpoint/cutoff identities are reported as unresolved.

This is NOT a global parity theorem, a complete cobordism, or a Toeplitz proof.
Standard library Fraction only, no approximate numerical root oracle.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import permutations, product, combinations
from math import gcd, lcm

from polygon_orbit_census_v223 import (
    Point, pt, add, sub, mul, dot, cross, validate_jordan_polygon,
    _square_residual, _affine_system, _solve, _side_sq, census_polygon,
)


def polynomial(values):
    a = tuple(Q(v) for v in values)
    while len(a)>1 and not a[-1]: a=a[:-1]
    return a or (Q(0),)

ZERO=polynomial((0,))
ONE=polynomial((1,))
T=polynomial((0,1))

def iszero(p):return len(p)==1 and p[0]==0

def plus(a,b):
    return polynomial(((a[i] if i<len(a) else Q(0))+
                       (b[i] if i<len(b) else Q(0)))
                      for i in range(max(len(a),len(b))))

def neg(a):return polynomial(-v for v in a)
def minus(a,b):return plus(a,neg(b))

def times(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):out[i+j]+=u*v
    return polynomial(out)

def derivative(a):
    return polynomial(i*a[i] for i in range(1,len(a)))

def evaluate(a,t):
    t=Q(t);r=Q(0)
    for v in reversed(a):r=r*t+v
    return r

def divide(a,b):
    if iszero(b):raise ZeroDivisionError('zero polynomial divisor')
    rem=list(polynomial(a));quot=[Q(0)]*max(1,len(a)-len(b)+1)
    while len(rem)>=len(b) and any(rem):
        k=len(rem)-len(b);v=rem[-1]/b[-1]
        quot[k]+=v
        for i,c in enumerate(b):rem[k+i]-=v*c
        while len(rem)>1 and not rem[-1]:rem.pop()
    return polynomial(quot),polynomial(rem)

def _primitive_keep_sign(p):
    p=polynomial(p)
    if iszero(p):return p
    denom=lcm(*(v.denominator for v in p))
    v=[int(x*denom) for x in p]
    g=gcd(*v)
    return polynomial(Q(z//g) for z in v)

def poly_gcd(a,b):
    while not iszero(b):
        _,r=divide(a,b)
        a,b=b,r
    return polynomial(x/a[-1] for x in a) if not iszero(a) else ZERO

def squarefree(p):
    p=polynomial(p)
    if len(p)<=1:return p
    g=poly_gcd(p,derivative(p))
    q,r=divide(p,g)
    assert iszero(r)
    return _primitive_keep_sign(q)


def sturm_sequence(p):
    p=squarefree(p)
    if len(p)==1:return ()
    result=[_primitive_keep_sign(p),_primitive_keep_sign(derivative(p))]
    while len(result[-1])>1:
        _,r=divide(result[-2],result[-1])
        if iszero(r):raise AssertionError('non-squarefree input in Sturm chain')
        result.append(_primitive_keep_sign(neg(r)))
    return tuple(result)


def variations_at(sturm,x):
    signs=[]
    for p in sturm:
        v=evaluate(p,x)
        if v:signs.append(1 if v>0 else -1)
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def count_roots_open(sturm,lo,hi):
    """Sturm number of distinct roots in STRICT interval, even when hi root."""
    if not lo<hi:raise ValueError('strict rational bracket required')
    if not sturm:return 0
    return variations_at(sturm,lo)-variations_at(sturm,hi)-int(evaluate(sturm[0],hi)==0)


@dataclass(frozen=True)
class RootBracket:
    # lo==hi denotes a proven exact rational root. Otherwise exactly one
    # real root is contained in the OPEN rational interval (lo,hi).
    lower:Q
    upper:Q
    def __post_init__(self):
        if self.lower>self.upper:raise ValueError('ordered bracket expected')
    @property
    def rational(self):return self.lower==self.upper


def isolate_roots(p,lo,hi,width=Q(1,2**24),max_nodes=100000):
    """All distinct real roots in [lo,hi], certified without floating point.

    Root inclusion at dyadic split points is exact. Irrational roots are
    bracketed into rational open intervals of width <= width. Returns an
    EXACT root when encountered rationally. Multiple roots are isolated
    using the squarefree part (multiplicity is not inferred).
    """
    p=polynomial(p);lo,hi,width=Q(lo),Q(hi),Q(width)
    if not lo<hi or width<=0 or max_nodes<1:raise ValueError('invalid bounds/budget')
    if iszero(p):raise ValueError('identically-zero polynomial has no isolated roots')
    if len(p)==1:return ()
    sq=squarefree(p);st=sturm_sequence(sq)
    output=[]
    if evaluate(sq,lo)==0:output.append(RootBracket(lo,lo))
    if evaluate(sq,hi)==0:output.append(RootBracket(hi,hi))
    todo=[(lo,hi,count_roots_open(st,lo,hi))]
    count=0
    while todo:
        a,b,n=todo.pop();count+=1
        if count>max_nodes:raise RuntimeError('root isolation budget exhausted')
        if n==0:continue
        if n==1 and b-a<=width and evaluate(sq,a)!=0 and evaluate(sq,b)!=0:
            output.append(RootBracket(a,b));continue
        m=(a+b)/2
        v=evaluate(sq,m)
        if v==0:output.append(RootBracket(m,m))
        left=count_roots_open(st,a,m)
        right=count_roots_open(st,m,b)
        if left+right+(v==0)!=n:raise AssertionError('Sturm split identity failed')
        todo.append((m,b,right));todo.append((a,m,left))
    output.sort(key=lambda z:(z.lower,z.upper))
    return tuple(output)


def determinant_poly(M):
    if not M or any(len(row)!=len(M) for row in M):
        raise ValueError('square polynomial matrix required')
    n=len(M);result=ZERO
    for perm in permutations(range(n)):
        term=ONE
        for i,j in enumerate(perm):term=times(term,M[i][j])
        inversions=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        result=plus(result,neg(term) if inversions%2 else term)
    return result


def moving_point(p0,p1):
    return (polynomial((p0.x,p1.x-p0.x)),polynomial((p0.y,p1.y-p0.y)))


def pairplus(a,b):return (plus(a[0],b[0]),plus(a[1],b[1]))
def pairminus(a,b):return (minus(a[0],b[0]),minus(a[1],b[1]))
def pairtimes(a,b):return (times(a,b[0]),times(a,b[1]))
def pairrot(a,orientation):
    return (neg(a[1]),a[0]) if orientation==1 else (a[1],neg(a[0]))
def pairdot(a,b):return plus(times(a[0],b[0]),times(a[1],b[1]))


def residual_poly(pts,orientation):
    p0,p1,p2,p3=pts
    e=pairminus(p1,p0);je=pairrot(e,orientation)
    z=pairminus(pairminus(p2,p1),je)
    w=pairminus(pairminus(p3,p0),je)
    return (*z,*w)


@dataclass(frozen=True)
class ChartPolynomials:
    edges:tuple
    orientation:int
    determinant:tuple
    numerators:tuple
    upper_differences:tuple
    cutoff:tuple
    cutoff_squared:Q
    @property
    def identically_singular(self):return iszero(self.determinant)
    def algebraic_parameters(self,t):
        D=evaluate(self.determinant,t)
        if not D:raise ValueError('singular at given time')
        return tuple(evaluate(N,t)/D for N in self.numerators)
    def eligible(self,t):
        u=self.algebraic_parameters(t)
        return (all(0<x<1 for x in u) and evaluate(self.cutoff,t)>0)
    def candidate_polynomials(self):
        pairs=[('singular',self.determinant)]
        for i,p in enumerate(self.numerators):pairs.append((f'edge-{i}-start',p))
        for i,p in enumerate(self.upper_differences):pairs.append((f'edge-{i}-end',p))
        pairs.append(('size-cutoff',self.cutoff))
        return tuple(pairs)


def chart_polynomials(start_vertices,end_vertices,edge_indices,orientation=1,cutoff_squared=Q(0)):
    """Exact complete event candidates on one moving ordered edge chart.

    Input vertices linearly interpolate in common time t in [0,1].
    Parameters are root ratios N_i(t)/D(t) when D(t)!=0. Interior square
    edge E(t) has polynomial numerator, with side² = |E|² / D².
    """
    P=tuple(start_vertices);R=tuple(end_vertices);n=len(P)
    if len(R)!=n or n<3 or len(edge_indices)!=4 or any(not 0<=i<n for i in edge_indices):
        raise ValueError('matching polygon vertices and four valid edges required')
    if orientation not in (-1,1):raise ValueError('orientation +/-1')
    cutoff=Q(cutoff_squared)
    if cutoff<0:raise ValueError('nonnegative squared-side cutoff')
    verts=tuple(moving_point(a,b) for a,b in zip(P,R))
    starts=tuple(verts[i] for i in edge_indices)
    dirs=tuple(pairminus(verts[(i+1)%n],verts[i]) for i in edge_indices)
    base=residual_poly(starts,orientation)
    columns=[]
    for i in range(4):
        q=list(starts);q[i]=pairplus(q[i],dirs[i])
        r=residual_poly(q,orientation)
        columns.append(tuple(minus(r[j],base[j]) for j in range(4)))
    A=tuple(tuple(columns[j][i] for j in range(4)) for i in range(4))
    b=tuple(neg(v) for v in base)
    D=determinant_poly(A)
    N=[]
    for i in range(4):
        mat=tuple(tuple(b[row] if col==i else A[row][col] for col in range(4)) for row in range(4))
        N.append(determinant_poly(mat))
    U=tuple(minus(D,p) for p in N)
    e_base=pairminus(starts[1],starts[0])
    E=pairplus(pairtimes(D,e_base),pairminus(pairtimes(N[1],dirs[1]),pairtimes(N[0],dirs[0])))
    H=minus(pairdot(E,E),times(polynomial((cutoff,)),times(D,D)))
    return ChartPolynomials(tuple(edge_indices),orientation,D,tuple(N),U,H,cutoff)


@dataclass(frozen=True)
class CandidateEvent:
    chart:tuple
    orientation:int
    kind:str
    polynomial:tuple
    bracket:RootBracket


@dataclass(frozen=True)
class ChartEventCover:
    chart:tuple
    orientation:int
    events:tuple
    unresolved:tuple
    cutoff_squared:Q
    def certified_chambers(self,start,end):
        """Open rational intervals containing no candidate event root."""
        if self.unresolved:
            raise ValueError('persistent identities prevent complete chart cover')
        boundaries=sorted({Q(start),Q(end)}|
                          {e.bracket.lower for e in self.events}|
                          {e.bracket.upper for e in self.events})
        return tuple((a,b) for a,b in zip(boundaries,boundaries[1:])
                     if a<b and not any(e.bracket.lower<b and e.bracket.upper>a
                                       for e in self.events))


def isolate_chart_events(chart,start=0,end=1,width=Q(1,2**24),budget=100000):
    """Every possible chart occupancy/cutoff change lies in these brackets.

    An identically zero event polynomial is marked unresolved; this avoids
    claiming a finite wall list for persistent singular/endpoint roots.
    The result is a *superset* of physical events, not a classification.
    """
    start,end=Q(start),Q(end)
    events=[];unresolved=[]
    for kind,p in chart.candidate_polynomials():
        if iszero(p):
            unresolved.append(kind)
            continue
        for b in isolate_roots(p,start,end,width,budget):
            events.append(CandidateEvent(chart.edges,chart.orientation,kind,p,b))
    return ChartEventCover(chart.edges,chart.orientation,tuple(events),tuple(unresolved),chart.cutoff_squared)


def polygon_event_cover(P0,P1,cutoff_squared,start=0,end=1,
                        width=Q(1,2**16),include_repeated_edges=True,
                        max_charts=100000,require_jordan_isotopy=True):
    """Potential event-time cover, NOT a certified global orbit cobordism.

    Each nonsingular moving four-edge chart contributes all determinant,
    segment-parameter-endpoint and size-cutoff candidate roots. Returns
    explicit unresolved charts when some polynomial is identically zero.
    Complete only on chart strata for which unresolved is empty. A
    candidate may be spurious because the square lies outside its edges.
    By default a separate root-free triple-collinearity certificate proves
    that the whole polygonal homotopy stays embedded; it is sufficient,
    not necessary, and rejects paths with persistent collinear triples.
    """
    P0=validate_jordan_polygon(tuple(P0));P1=validate_jordan_polygon(tuple(P1))
    if require_jordan_isotopy:
        certify_jordan_isotopy(P0,P1,start,end)
    n=len(P0)
    if len(P1)!=n:raise ValueError('same vertex count required')
    results=[];unresolved=[];count=0
    for ids in product(range(n),repeat=4):
        if not include_repeated_edges and len(set(ids))<4:continue
        if any(ids.count(i)>=3 for i in set(ids)):
            continue  # cannot carry nonzero square with three collinear points
        for sigma in (1,-1):
            count+=1
            if count>max_charts:raise RuntimeError('moving-chart budget reached')
            chart=chart_polynomials(P0,P1,ids,sigma,cutoff_squared)
            cover=isolate_chart_events(chart,start,end,width)
            results.append(cover)
            if cover.unresolved:unresolved.append((ids,sigma,cover.unresolved))
    return tuple(results),tuple(unresolved)



@dataclass(frozen=True)
class JordanIsotopyCertificate:
    """Sufficient certificate that every intermediate polygon is Jordan."""
    start:Q
    end:Q
    triples_checked:int
    nonzero_triple_orientation:tuple


def interpolated_polygon(P0,P1,t):
    t=Q(t)
    if len(P0)!=len(P1):raise ValueError('same polygon vertex count required')
    return tuple(add(a,mul(sub(b,a),t)) for a,b in zip(P0,P1))


def three_point_orientation_polynomial(P0,P1,i,j,k):
    a,b,c=(moving_point(P0[h],P1[h]) for h in (i,j,k))
    v=pairminus(b,a);w=pairminus(c,a)
    return minus(times(v[0],w[1]),times(v[1],w[0]))


def certify_jordan_isotopy(P0,P1,start=0,end=1):
    """Sufficient *general-position* isotopy proof, exact rational only.

    If no three distinct polygon vertices become collinear at any time,
    a simple polygon cannot acquire an edge crossing or edge collapse:
    the first failure of simplicity entails a vertex-on-edge, overlapping
    collinear edges, or coinciding vertices, each forces collinearity of
    some triple. A root-free triple orientation test thus certifies that
    the initial Jordan polygon remains embedded across the full interval.

    A triple that stays collinear is REJECTED: it may be harmless, but this
    simple test is not an exhaustive isotopy decision procedure.
    """
    a,b=Q(start),Q(end)
    if not a<b or not 0<=a<b<=1:raise ValueError('require [a,b] subset [0,1]')
    if len(P0)!=len(P1) or len(P0)<3:raise ValueError('matching polygons required')
    validate_jordan_polygon(interpolated_polygon(P0,P1,a))
    results=[]
    for i,j,k in combinations(range(len(P0)),3):
        f=three_point_orientation_polynomial(P0,P1,i,j,k)
        if iszero(f):raise ValueError(f'permanently collinear triple {(i,j,k)}; certificate inconclusive')
        roots=isolate_roots(f,a,b,width=Q(1,2**12))
        if roots:
            raise ValueError(f'triple {(i,j,k)} may become collinear: {roots}')
        v=evaluate(f,a)
        if not v:raise AssertionError('root-free polynomial zero at starting time')
        results.append((i,j,k,1 if v>0 else -1))
    return JordanIsotopyCertificate(a,b,len(results),tuple(results))


@dataclass(frozen=True)
class RationalVertexSquare:
    """A true square ON a polygon vertex at a specified rational time."""
    time:Q
    edge_indices:tuple
    orientation:int
    parameters:tuple
    points:tuple
    squared_side:Q
    endpoint_dimensions:tuple


def rational_vertex_square(P0,P1,chart,t):
    """Verify an actual endpoint square, not just a candidate polynomial root.

    Does not assert that the square's geometric orbit is born or destroyed.
    Candidate roots with parameters outside the four edge charts are rejected.
    """
    t=Q(t)
    if evaluate(chart.determinant,t)==0:
        raise ValueError('chart determinant vanishes: Cramer branch undefined')
    u=chart.algebraic_parameters(t)
    if any(v<0 or v>1 for v in u):
        raise ValueError('algebraic square lies outside its polygon edges')
    ends=tuple(i for i,v in enumerate(u) if v in (0,1))
    if not ends:raise ValueError('no edge-endpoint incidence at supplied time')
    P=interpolated_polygon(P0,P1,t)
    edges=tuple((P[i],P[(i+1)%len(P)]) for i in chart.edges)
    pts=tuple(add(a,mul(sub(b,a),v)) for (a,b),v in zip(edges,u))
    if _square_residual(pts,chart.orientation)!=(0,0,0,0):
        raise AssertionError('Cramer polynomial does not give a square')
    q=_side_sq(pts)
    if q<=chart.cutoff_squared:
        raise ValueError('square below or exactly at the positive size cutoff')
    if len(set(pts))!=4:raise ValueError('collapsed or repeated square point')
    return RationalVertexSquare(t,chart.edges,chart.orientation,u,pts,q,ends)

def report():
    from polygon_orbit_census_v223 import convex_vertex_alignment_polygon
    P0=convex_vertex_alignment_polygon(Q(0))
    P1=convex_vertex_alignment_polygon(Q(9,8))
    chart=chart_polynomials(P0,P1,(0,1,2,3),1,Q(1,100))
    cov=isolate_chart_events(chart)
    return {
        'selected_chart':chart.edges,
        'orientation':chart.orientation,
        'determinant_degree':len(chart.determinant)-1,
        'cutoff_degree':len(chart.cutoff)-1,
        'event_candidates':len(cov.events),
        'unresolved_identities':cov.unresolved,
        'scope':'complete candidate cover for nonsingular moving chart; not global orbit cobordism',
    }

if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
