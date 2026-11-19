"""v228: exact affine feasibility over a *real algebraic event time*.

For an event alpha specified by a squarefree defining polynomial and a
certified isolating bracket, solve A(alpha) u=b(alpha), 0<=u_i<=1, using
finite vertex enumeration over Q(alpha). When A,b come from a moving Jordan
polygon edge chart, test positive squared-side cutoff at every feasible
vertex. Distinguish an empty system, a filtered-out feasible system, an
isolated filtered square, and a positive-dimensional square family.

This resolves event-slice feasibility, NOT one-sided singular Puiseux germs,
branch continuation, or an arbitrary-Jordan square-peg theorem.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product
from collections import Counter

from moving_polygon_events_v224 import (
    ZERO, ONE, polynomial, plus, minus, times, evaluate,
    determinant_poly, interpolated_polygon, moving_point, pairminus,
)
from persistent_strata_v225 import (
    affine_system_polynomials, classify_polygon_motion,
)
from local_orbit_germs_v227 import AlgebraicTime, sign_at_time
from polygon_orbit_census_v223 import validate_jordan_polygon


class AlgebraicFeasibilityInconclusive(RuntimeError):
    """Explicit bounded-resource / unsupported-case failure, not emptiness."""


def _unit_row(j):
    return tuple(ONE if k==j else ZERO for k in range(4))


def _minor(matrix, row_ids, col_ids):
    return determinant_poly(tuple(tuple(matrix[i][j] for j in col_ids)
                                  for i in row_ids))


def rank_and_basis_rows_at(A, alpha):
    """Exact rank over real algebraic number Q(alpha), with independent rows."""
    if len(A)!=4 or any(len(row)!=4 for row in A):
        raise ValueError('expected 4x4 polynomial matrix')
    for k in range(4,0,-1):
        for ri in combinations(range(4),k):
            for ci in combinations(range(4),k):
                if sign_at_time(_minor(A,ri,ci),alpha):
                    return k,ri
    return 0,()


@dataclass(frozen=True)
class AlgebraicQuotient:
    numerator:tuple
    denominator:tuple
    event:AlgebraicTime

    def __post_init__(self):
        object.__setattr__(self,'numerator',polynomial(self.numerator))
        object.__setattr__(self,'denominator',polynomial(self.denominator))
        if sign_at_time(self.denominator,self.event)==0:
            raise ZeroDivisionError('denominator vanishes at algebraic event')

    def sign(self):
        return sign_at_time(self.numerator,self.event)*sign_at_time(self.denominator,self.event)

    def equal(self,other):
        if self.event!=other.event:raise ValueError('incompatible event descriptions')
        return sign_at_time(minus(times(self.numerator,other.denominator),
                                  times(other.numerator,self.denominator)),self.event)==0

    def rational_value(self):
        if not self.event.rational:raise ValueError('event is irrational')
        return evaluate(self.numerator,self.event.lo)/evaluate(self.denominator,self.event.lo)


def _box_feasible_quotients(nums,den,alpha):
    sd=sign_at_time(den,alpha)
    return all(sign_at_time(n,alpha)*sd>=0 and
               sign_at_time(minus(den,n),alpha)*sd>=0 for n in nums)


def _equivalent_vertices(a,b):
    return all(x.equal(y) for x,y in zip(a,b))


def algebraic_box_vertices(A,b,alpha,max_candidates=192):
    """All vertices of {u in [0,1]^4 : A(alpha)u=b(alpha)}.

    Each returned coordinate is a rational-function expression evaluated at
    alpha with a *provably nonzero* denominator. Vertex enumeration is
    complete for a bounded affine slice: at a vertex, 4-r independent box
    constraints can be added to r independent affine equations.
    """
    if len(A)!=4 or any(len(row)!=4 for row in A) or len(b)!=4:
        raise ValueError('4x4 polynomial affine system required')
    A=tuple(tuple(polynomial(z) for z in row) for row in A)
    b=tuple(polynomial(x) for x in b)
    if not isinstance(alpha,AlgebraicTime):raise TypeError('exact algebraic time required')
    rank,ri=rank_and_basis_rows_at(A,alpha)
    out=[];trials=0
    for fixed in combinations(range(4),4-rank):
        for vals in product((0,1),repeat=4-rank):
            trials+=1
            if trials>max_candidates:
                raise AlgebraicFeasibilityInconclusive('active-box vertex budget exhausted')
            M=tuple(A[i] for i in ri)+tuple(_unit_row(j) for j in fixed)
            rhs=tuple(b[i] for i in ri)+tuple(ONE if val else ZERO for val in vals)
            d=determinant_poly(M)
            if sign_at_time(d,alpha)==0:continue
            numerators=tuple(determinant_poly(tuple(tuple(rhs[k] if j==col else M[k][j]
                                        for j in range(4)) for k in range(4)))
                             for col in range(4))
            if not _box_feasible_quotients(numerators,d,alpha):continue
            if any(sign_at_time(minus(
                   sum_polynomials(times(A[i][j],numerators[j]) for j in range(4)),
                   times(b[i],d)),alpha)!=0 for i in range(4)):
                continue
            solution=tuple(AlgebraicQuotient(z,d,alpha) for z in numerators)
            if not any(_equivalent_vertices(solution,old) for old in out):
                out.append(solution)
    return rank,tuple(out)


def sum_polynomials(values):
    result=ZERO
    for p in values:result=plus(result,p)
    return result


def _physical_coordinate_numerators(P0,P1,ids,nums,den):
    n=len(P0)
    verts=tuple(moving_point(x,y) for x,y in zip(P0,P1))
    result=[]
    for i,edge in enumerate(ids):
        a=verts[edge];v=pairminus(verts[(edge+1)%n],a)
        result.append(tuple(plus(times(den,a[k]),times(nums[i],v[k])) for k in range(2)))
    return tuple(result)


def _positive_cutoff_sign(P0,P1,ids,solution,cutoff):
    alpha=solution[0].event
    # Each complete vertex solution has a common Cramer denominator.
    d=solution[0].denominator
    if any(s.denominator!=d for s in solution):
        raise AssertionError('expected common determinant denominator')
    nums=tuple(s.numerator for s in solution)
    coords=_physical_coordinate_numerators(P0,P1,ids,nums,d)
    edge=(minus(coords[1][0],coords[0][0]),
          minus(coords[1][1],coords[0][1]))
    H=minus(plus(times(edge[0],edge[0]),times(edge[1],edge[1])),
            times(polynomial((Q(cutoff),)),times(d,d)))
    return sign_at_time(H,alpha),coords


@dataclass(frozen=True)
class SingularEventSnapshot:
    edges:tuple
    orientation:int
    event:AlgebraicTime
    rank:int
    parameter_vertices:tuple
    status:str
    witness:tuple|None  # four algebraic coordinate pairs, each as AlgebraicQuotient
    positive_vertex_indices:tuple

    @property
    def has_positive_square(self):
        return bool(self.positive_vertex_indices)


def singular_polygon_chart_at_event(P0,P1,ids,sigma,event,cutoff_squared=Q(0),max_candidates=192):
    """Complete exact *instantaneous* decision for ONE moving polygon chart.

    In addition to an exact affine box feasibility certificate, the squared
    side cutoff is checked over the entire affine feasible polytope by
    convexity: a convex quadratic attains its maximum at a vertex.

    All nondegenerate solutions on the chosen physical polygon edges are
    included, even at edge endpoints; orbit deduplication is separate.
    """
    P0=tuple(P0);P1=tuple(P1);ids=tuple(ids)
    if len(P0)!=len(P1) or len(P0)<3 or len(ids)!=4 or any(
        type(i) is not int or not 0<=i<len(P0) for i in ids):
        raise ValueError('matching polygon vertices and valid edge tuple required')
    if sigma not in (-1,1):raise ValueError('orientation +/-1 required')
    if Q(cutoff_squared)<0:raise ValueError('nonnegative cutoff required')
    A,b=affine_system_polynomials(P0,P1,ids,sigma)
    rank,vertices=algebraic_box_vertices(A,b,event,max_candidates)
    positives=[];witness=None
    for j,sol in enumerate(vertices):
        sg,physical_nums=_positive_cutoff_sign(P0,P1,ids,sol,cutoff_squared)
        if sg>0:
            positives.append(j)
            if witness is None:
                common=sol[0].denominator
                witness=tuple(tuple(AlgebraicQuotient(coord,common,event)
                                    for coord in pair) for pair in physical_nums)
    if not vertices:status='infeasible'
    elif not positives:status='no-positive-size-square'
    elif len(vertices)==1:status='isolated-square'
    else:status='positive-dimensional-square-family'
    return SingularEventSnapshot(ids,sigma,event,rank,vertices,status,witness,
                                 tuple(positives))


@dataclass(frozen=True)
class SingularEventAudit:
    event:AlgebraicTime
    cutoff_squared:Q
    checked_singular_charts:int
    statuses:tuple
    witnesses:tuple
    unresolved_chart_count:int
    scope:str


def audit_singular_event_slices(P0,P1,event,cutoff_squared=Q(0),max_charts=600,
                                verify_rational_jordan=True):
    """Inspect *all* charts whose square matrix is singular at event.

    This audits time-slice feasibility, not one-sided branches, and does not
    infer orbit parity. Runtime is bounded by number of examined charts.
    """
    P0=validate_jordan_polygon(tuple(P0));P1=validate_jordan_polygon(tuple(P1))
    if len(P0)!=len(P1):raise ValueError('same polygon vertices required')
    if verify_rational_jordan and event.rational:
        validate_jordan_polygon(interpolated_polygon(P0,P1,event.lo))
    from persistent_strata_v225 import _chart
    n=len(P0);stats=Counter();witnesses=[];checked=0;unresolved=0
    for ids in product(range(n),repeat=4):
        if any(ids.count(i)>=3 for i in set(ids)):continue
        for sigma in (1,-1):
            c=_chart(P0,P1,ids,sigma,Q(cutoff_squared))
            if sign_at_time(c.determinant,event)!=0:continue
            checked+=1
            if checked>max_charts:
                raise AlgebraicFeasibilityInconclusive('singular chart audit budget exhausted')
            try:s=singular_polygon_chart_at_event(P0,P1,ids,sigma,event,cutoff_squared)
            except AlgebraicFeasibilityInconclusive:
                unresolved+=1;continue
            stats[s.status]+=1
            if s.has_positive_square:witnesses.append((ids,sigma,s.status,s.witness))
    return SingularEventAudit(event,Q(cutoff_squared),checked,tuple(sorted(stats.items())),
         tuple(witnesses),unresolved,'event-slice feasibility only, not orbit continuation')


def report():
    from polygon_orbit_census_v223 import pt
    P=tuple(pt(*xy) for xy in ((0,0),(2,0),(2,1),(0,1)))
    event=AlgebraicTime((-Q(1,2),1),Q(1,2),Q(1,2))
    snap=singular_polygon_chart_at_event(P,P,(0,0,2,2),-1,event,Q(1,100))
    return {'rectangle_status':snap.status,'rectangle_rank':snap.rank,
            'feasible_vertices':len(snap.parameter_vertices),
            'scope':'exact event-slice algebraic feasibility; no singular branch Puiseux continuation or Toeplitz proof'}


if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
