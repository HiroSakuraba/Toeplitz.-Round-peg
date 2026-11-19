"""v225: exact elimination of persistent-identity flags on moving polygon charts.

Resolves regular Cramer charts with permanently endpoint-locked roots using
half-open ownership of polygon vertices. A square at P_j belongs to its
outgoing edge j at parameter 0, not the incoming edge j-1 at parameter 1.
Also excludes 2+2 configurations on adjacent polygon edges geometrically;
and detects singular chart inconsistency via a nonzero polynomial minor.

Not a complete algebraic moving-polygon orbit cobordism. A persistently
consistent rank-deficient chart is explicitly INCONCLUSIVE unless separately
shown geometrically impossible. Exact Fraction arithmetic throughout.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations,product
from collections import Counter
from functools import lru_cache

from moving_polygon_events_v224 import (
    ZERO, ONE, polynomial, iszero, plus, minus, neg, times, evaluate,
    isolate_roots, determinant_poly, moving_point, pairminus, pairplus,
    residual_poly, chart_polynomials, interpolated_polygon,
    certify_jordan_isotopy,
)
from polygon_orbit_census_v223 import (
    pt, add, sub, mul, _square_residual, _side_sq,
    _affine_system, _box_vertices, _evaluate, _rank,
    validate_jordan_polygon, census_polygon, PositiveDimensionalSquares,
)



@lru_cache(maxsize=8192)
def _cached_chart_polynomials(P0,P1,ids,sigma,cutoff):
    """Memoize exact polynomial charts, not approximate floating-point data."""
    return chart_polynomials(P0,P1,ids,sigma,cutoff)


def _chart(P0,P1,ids,sigma,cutoff):
    return _cached_chart_polynomials(tuple(P0),tuple(P1),tuple(ids),sigma,Q(cutoff))


def adjacent_two_pairs(ids, polygon_size):
    """Two vertices per each of two incident straight edges cannot form a square.

    Two vertices on each of two lines of a square must be an opposite pair
    (diagonals, perpendicular lines with SAME midpoint) or two adjacent
    pairs (parallel supporting lines). Incident polygon edges share an end
    point and each occupies a ray from it. Their two distinct points have
    midpoints strictly inside distinct rays, never a shared midpoint; the
    lines are not distinct parallel lines. If collinear, all four points
    lie on a single line. These exhaust square vertex pairings.
    """
    c=Counter(ids)
    if sorted(c.values())!=[2,2]:return False
    a,b=tuple(c)
    return (a-b)%polygon_size in (1,polygon_size-1)


def affine_system_polynomials(start_vertices,end_vertices,ids,sigma):
    """Exact A(t)u=b(t); A is 4x4, with polynomial entries of degree <=1."""
    P=tuple(start_vertices);R=tuple(end_vertices)
    if len(P)!=len(R) or len(P)<3 or len(ids)!=4:
        raise ValueError('matching polygon endpoints and 4 edges required')
    if any(type(i) is not int or not 0<=i<len(P) for i in ids):
        raise ValueError('bad edge index')
    if sigma not in (-1,1):raise ValueError('orientation +/-1 required')
    verts=tuple(moving_point(a,b) for a,b in zip(P,R))
    starts=tuple(verts[i] for i in ids)
    dirs=tuple(pairminus(verts[(i+1)%len(P)],verts[i]) for i in ids)
    base=residual_poly(starts,sigma)
    cols=[]
    for i in range(4):
        changed=list(starts);changed[i]=pairplus(changed[i],dirs[i])
        f=residual_poly(changed,sigma)
        cols.append(tuple(minus(f[j],base[j]) for j in range(4)))
    A=tuple(tuple(cols[j][i] for j in range(4)) for i in range(4))
    b=tuple(neg(z) for z in base)
    return A,b


def _nonzero_minor(matrix,order):
    if order==0:return ONE
    rows=len(matrix);cols=len(matrix[0])
    for r in combinations(range(rows),order):
        for c in combinations(range(cols),order):
            d=determinant_poly(tuple(tuple(matrix[i][j] for j in c) for i in r))
            if not iszero(d):return d
    return None


def generic_polynomial_rank(matrix):
    """Rank over Q(t), not rank at a particular rational time."""
    for k in range(min(len(matrix),len(matrix[0])),0,-1):
        if _nonzero_minor(matrix,k) is not None:return k
    return 0


@dataclass(frozen=True)
class StratumCertificate:
    edges:tuple
    orientation:int
    status:str
    forced_start:tuple
    forced_end:tuple
    active_polynomials:tuple
    exceptional_roots:tuple
    obstruction:tuple
    rank_A:int
    rank_augmented:int

    @property
    def resolved(self):return self.status!='genuinely-unresolved'

    @property
    def has_regular_branch(self):return self.status in ('regular','vertex-owned')



def singular_inconsistency_certificate(A,b,start=Q(0),end=Q(1),width=Q(1,2**14)):
    """Rank over Q(t) gives a necessary algebraic support set for A(t)u=b(t).

    Let r=rank_Q(t)(A) and s=rank_Q(t)([A|b]). If s>r, a
    nonzero (r+1)-minor M(t) of the augmentation certifies inconsistency
    for EVERY real time with M(t)!=0. This remains true when A's rank
    drops at special times, because rank(A(t))<=r holds identically.

    Roots are exact rationals or isolating intervals. At such roots the
    method makes NO positive feasibility claim; a separate algebraic
    snapshot or more general CAD is required. If s==r, this rank test is
    inconclusive and returns witness=None, not a false existence claim.
    """
    A=tuple(tuple(polynomial(z) for z in row) for row in A)
    b=tuple(polynomial(z) for z in b)
    if len(A)!=len(b) or len(A)!=4 or any(len(row)!=4 for row in A):
        raise ValueError('4x4 affine polynomial system expected')
    a,c=Q(start),Q(end)
    if not a<c or not 0<=a<c<=1:raise ValueError('valid [a,c] subset [0,1] required')
    aug=tuple(tuple(row)+(b[i],) for i,row in enumerate(A))
    ra=generic_polynomial_rank(A)
    rb=generic_polynomial_rank(aug)
    if rb<=ra:return ra,rb,None,()
    witness=_nonzero_minor(aug,ra+1)
    assert witness is not None
    return ra,rb,witness,isolate_roots(witness,a,c,width)


def classify_chart(P0,P1,ids,sigma=1,cutoff_squared=Q(0),
                   start=Q(0),end=Q(1),width=Q(1,2**14)):
    """Prove how a moving chart can affect canonical positive-size orbit count.

    * regular / vertex-owned: root branch exists when D!=0; endpoints at
      u==0 are included; at u==1 the chart is redundant (incoming edge).
    * structurally-empty: 2+2 on incident edges impossible geometrically.
    * cutoff-locked: side^2 identically equal cutoff on regular branch.
    * exceptional-time-only: rank([A|b])>rank(A) generically and nonzero
      obstruction minor proves NO square outside its exact real roots.
    * genuinely-unresolved: identically singular rank-deficient equation
      may have a square continuum; no false generic parity issued.
    """
    P0=tuple(P0);P1=tuple(P1);ids=tuple(ids)
    a,b=Q(start),Q(end)
    if not a<b or not 0<=a<b<=1 or Q(cutoff_squared)<0:
        raise ValueError('invalid time interval or cutoff')
    if adjacent_two_pairs(ids,len(P0)):
        return StratumCertificate(ids,sigma,'structurally-empty',(),(),(),(),ZERO,-1,-1)
    c=_chart(P0,P1,ids,sigma,cutoff_squared)
    if not iszero(c.determinant):
        forced0=tuple(i for i,z in enumerate(c.numerators) if iszero(z))
        forced1=tuple(i for i,z in enumerate(c.upper_differences) if iszero(z))
        if forced1:
            return StratumCertificate(ids,sigma,'incoming-edge-duplicate',forced0,forced1,(),(),ZERO,4,4)
        if iszero(c.cutoff):
            return StratumCertificate(ids,sigma,'cutoff-locked',forced0,forced1,(),(),ZERO,4,4)
        nonzero=tuple((kind,p) for kind,p in c.candidate_polynomials() if not iszero(p))
        events=tuple((kind,r) for kind,p in nonzero for r in isolate_roots(p,a,b,width))
        return StratumCertificate(ids,sigma,'vertex-owned' if forced0 else 'regular',
                                  forced0,forced1,nonzero,events,ZERO,4,4)
    A,v=affine_system_polynomials(P0,P1,ids,sigma)
    r,ra,witness,roots=singular_inconsistency_certificate(A,v,a,b,width)
    if witness is not None:
        # If obstruction has no root in the interval, whole slab is empty.
        return StratumCertificate(ids,sigma,'exceptional-time-only' if roots else 'always-empty',
                                  (),(),(),tuple(('augmented-rank',x) for x in roots),witness,r,ra)
    return StratumCertificate(ids,sigma,'genuinely-unresolved',(),(),(),(),ZERO,r,ra)


def owned_square_at_rational_time(P0,P1,ids,sigma,t,cutoff_squared=Q(0)):
    """A rational-time Cramer witness canonically assigned to outgoing edges.

    Returns None if the chart is singular, root lies off its canonical
    half-open edge windows, or side <= cutoff. Never silently interprets
    singularity as emptiness; use singular_chart_snapshot() in that case.
    """
    c=_chart(P0,P1,ids,sigma,cutoff_squared)
    D=evaluate(c.determinant,t)
    if not D:return None
    us=c.algebraic_parameters(t)
    if not all(0<=u<1 for u in us):return None
    if evaluate(c.cutoff,t)<=0:return None
    P=interpolated_polygon(P0,P1,t)
    n=len(P);edges=tuple((P[j],P[(j+1)%n]) for j in ids)
    points=_evaluate(edges,us)
    if len(set(points))!=4 or _square_residual(points,sigma)!=(0,0,0,0):
        raise AssertionError('Cramer root is not nondegenerate square')
    return points,us,_side_sq(points)


@dataclass(frozen=True)
class SingularSnapshot:
    status:str
    rank:int
    vertices:tuple
    max_side_squared:Q


def singular_chart_snapshot(polygon,ids,sigma=1,cutoff_squared=Q(0)):
    """Exact static rank-deficient feasibility/cutoff classification.

    More than one feasible parameter-polytope vertex and max side²>cutoff
    implies a positive-dimensional family. One vertex may give one
    isolated square at a polygon vertex. No false finite count is assigned.
    """
    P=tuple(polygon);n=len(P)
    edges=tuple((P[i],P[(i+1)%n]) for i in ids)
    A,b=_affine_system(edges,sigma)
    r=_rank(A)
    if r==4:raise ValueError('static chart nonsingular')
    verts=tuple(_box_vertices(A,b))
    mx=max((_side_sq(_evaluate(edges,u)) for u in verts),default=Q(0))
    if mx<=Q(cutoff_squared):status='no-positive-size-square'
    elif len(verts)>1:status='positive-dimensional-family'
    else:status='isolated-boundary-square'
    return SingularSnapshot(status,r,verts,mx)


@dataclass(frozen=True)
class GlobalStrataAudit:
    certificates:tuple
    status_counts:tuple
    raw_candidate_events:int
    unresolved_charts:tuple
    @property
    def completely_stratified(self):return not self.unresolved_charts


def classify_polygon_motion(P0,P1,cutoff_squared=Q(0),
                            start=Q(0),end=Q(1),width=Q(1,2**12),
                            prove_isotopy=True,max_charts=100000):
    """Audit *all* admissible moving charts, with unique vertex ownership.

    The status count is about chart stratification, not a count of square
    orbits/events. Charts with unresolved singular families are reported.
    """
    P0=validate_jordan_polygon(tuple(P0));P1=validate_jordan_polygon(tuple(P1))
    if len(P0)!=len(P1):raise ValueError('same vertex count required')
    if prove_isotopy:certify_jordan_isotopy(P0,P1,start,end)
    n=len(P0)
    result=[];counts=Counter();unresolved=[];total=0
    for ids in product(range(n),repeat=4):
        if any(ids.count(i)>=3 for i in set(ids)):
            continue  # three collinear points cannot be square vertices
        for sigma in (1,-1):
            if len(result)>=max_charts:raise RuntimeError('chart budget exhausted')
            c=classify_chart(P0,P1,ids,sigma,cutoff_squared,start,end,width)
            result.append(c);counts[c.status]+=1;total+=len(c.exceptional_roots)
            if not c.resolved:unresolved.append((ids,sigma))
    return GlobalStrataAudit(tuple(result),tuple(sorted(counts.items())),total,tuple(unresolved))


def generic_snapshot_check(P0,P1,t,cutoff_squared=Q(0)):
    """Compare canonical nonsingular orbit set against v223's complete census.

    This is an audit, NOT a new theorem if singular charts have actual
    positive-size roots. At such times status is INCONCLUSIVE, not false.
    """
    P=interpolated_polygon(P0,P1,Q(t))
    known=census_polygon(P,cutoff_squared)
    target={tuple(sorted(s.vertices)) for s in known.squares}
    found=set();uncertain=[]
    n=len(P)
    for ids in product(range(n),repeat=4):
        if any(ids.count(i)>=3 for i in ids) or adjacent_two_pairs(ids,n):continue
        for sigma in (1,-1):
            c=_chart(P0,P1,ids,sigma,cutoff_squared)
            if evaluate(c.determinant,t)==0:
                snapshot=singular_chart_snapshot(P,ids,sigma,cutoff_squared)
                if snapshot.status!='no-positive-size-square':
                    uncertain.append((ids,sigma,snapshot.status))
                continue
            v=owned_square_at_rational_time(P0,P1,ids,sigma,Q(t),cutoff_squared)
            if v is not None:found.add(tuple(sorted(v[0])))
    return {'known':target,'found':found,'unresolved':tuple(uncertain),'matches':found==target}


def report():
    from polygon_orbit_census_v223 import convex_vertex_alignment_polygon
    P0=convex_vertex_alignment_polygon(Q(0))
    P1=convex_vertex_alignment_polygon(Q(9,8))
    result=classify_polygon_motion(P0,P1,Q(1,100),Q(3,4),Q(1))
    return {
        'total_charts':len(result.certificates),
        'status_counts':dict(result.status_counts),
        'unresolved_charts':len(result.unresolved_charts),
        'raw_event_candidates':result.raw_candidate_events,
        'sample_t_3_4':generic_snapshot_check(P0,P1,Q(3,4),Q(1,100))['matches'],
        'sample_t_11_12':generic_snapshot_check(P0,P1,Q(11,12),Q(1,100))['matches'],
        'scope':'resolves persistent facet flags and incident-edge singular charts; no global parity or Toeplitz proof'
    }

if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
