"""v229: exact one-sided Cramer germs at a real algebraic rank-drop event.

The standalone core classifies rational-polynomial quotients N(t)/D(t).
It distinguishes removable finite limits, genuine poles, and persistent
singularity, with exact signs on both sides of a rational/irrational event.
The optional polygon bridge uses the v224, v225 and v228 code in the complete
Toeplitz repository. No floating-point step and no asserted global cobordism.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd, lcm


def poly(values):
    p=tuple(Q(x) for x in values)
    if not p:return (Q(0),)
    while len(p)>1 and p[-1]==0:p=p[:-1]
    return p

ZERO=poly((0,))
ONE=poly((1,))

def iszero(p):return poly(p)==ZERO

def plus(a,b):
    a,b=poly(a),poly(b)
    return poly((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0)
                for i in range(max(len(a),len(b))))

def neg(a):return poly(-x for x in a)
def minus(a,b):return plus(a,neg(b))

def times(a,b):
    a,b=poly(a),poly(b)
    result=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):result[i+j]+=x*y
    return poly(result)

def evaluate(p,t):
    s=Q(0)
    for x in reversed(poly(p)):s=s*Q(t)+x
    return s

def diff(p):return poly(i*x for i,x in enumerate(poly(p)) if i) if len(poly(p))>1 else ZERO

def divide(a,b):
    a,b=poly(a),poly(b)
    if iszero(b):raise ZeroDivisionError('zero polynomial')
    q=[Q(0)]*max(1,len(a)-len(b)+1);r=list(a)
    while len(r)>=len(b) and any(r):
        d=len(r)-len(b);v=r[-1]/b[-1];q[d]+=v
        for j,x in enumerate(b):r[d+j]-=v*x
        while len(r)>1 and r[-1]==0:r.pop()
    return poly(q),poly(r)

def polygcd(a,b):
    a,b=poly(a),poly(b)
    while not iszero(b):
        _,r=divide(a,b);a,b=b,r
    return poly(x/a[-1] for x in a) if not iszero(a) else ZERO

def squarefree(p):
    p=poly(p)
    if len(p)<=1:return p
    q,r=divide(p,polygcd(p,diff(p)))
    if not iszero(r):raise AssertionError('squarefree exact division failed')
    return q

def sturm(p):
    p=squarefree(p)
    if len(p)<=1:return ()
    out=[p,diff(p)]
    while len(out[-1])>1:
        _,r=divide(out[-2],out[-1])
        if iszero(r):raise AssertionError('squarefree Sturm remainder zero')
        out.append(neg(r))
    return tuple(out)

def _variations(chain,t):
    signs=[]
    for p in chain:
        v=evaluate(p,t)
        if v:signs.append((v>0)-(v<0))
    return sum(a!=b for a,b in zip(signs,signs[1:]))

def roots_open(chain,a,b):
    if not Q(a)<Q(b):raise ValueError('ordered bracket required')
    if not chain:return 0
    return (_variations(chain,Q(a))-_variations(chain,Q(b))-
            int(evaluate(chain[0],Q(b))==0))

class ContinuationInconclusive(RuntimeError):
    """Exact decision cannot be completed in prescribed scope/budget."""

@dataclass(frozen=True)
class ExactTime:
    factor:tuple
    lo:Q
    hi:Q
    def __post_init__(self):
        f=poly(self.factor);lo,hi=Q(self.lo),Q(self.hi)
        object.__setattr__(self,'factor',f)
        object.__setattr__(self,'lo',lo)
        object.__setattr__(self,'hi',hi)
        if iszero(f) or len(f)==1 or lo>hi:raise ValueError('proper defining polynomial and ordered bracket required')
        if lo==hi:
            if evaluate(f,lo):raise ValueError('exact rational time must be a root')
        elif evaluate(f,lo)==0 or evaluate(f,hi)==0 or roots_open(sturm(f),lo,hi)!=1:
            raise ValueError('open bracket must isolate one interior real root')
    @property
    def rational(self):return self.lo==self.hi


def root_sign(p,alpha,max_steps=160):
    """Sign of p(alpha), where alpha is a certified real algebraic root."""
    p=poly(p)
    if iszero(p):return 0
    if alpha.rational:
        v=evaluate(p,alpha.lo);return (v>0)-(v<0)
    g=polygcd(p,alpha.factor)
    if len(g)>1 and roots_open(sturm(g),alpha.lo,alpha.hi):return 0
    a,b=alpha.lo,alpha.hi;root_chain=sturm(alpha.factor)
    for _ in range(max_steps):
        if not roots_open(sturm(p),a,b):
            mid=(a+b)/2;v=evaluate(p,mid)
            if v==0:
                raise ContinuationInconclusive('sign sample hit an unrelated root')
            return (v>0)-(v<0)
        m=(a+b)/2
        if evaluate(alpha.factor,m)==0:
            v=evaluate(p,m);return (v>0)-(v<0)
        if roots_open(root_chain,a,m)==1:b=m
        elif roots_open(root_chain,m,b)==1:a=m
        else:raise ContinuationInconclusive('could not refine algebraic event isolator')
    raise ContinuationInconclusive('algebraic sign budget exhausted')


def order_at(p,alpha):
    """Vanishing order at alpha; None means polynomial is identically zero."""
    q=poly(p)
    if iszero(q):return None
    m=0
    while root_sign(q,alpha)==0:
        m+=1;q=diff(q)
        if iszero(q):raise AssertionError('nonzero polynomial derivatives all vanished')
    return m


def side_signs(p,alpha):
    """(left,at,right) signs at exact root and infinitesimally close sides."""
    p=poly(p)
    if iszero(p):return (0,0,0)
    v=root_sign(p,alpha)
    if v:return (v,v,v)
    order=order_at(p,alpha)
    q=p
    for _ in range(order):q=diff(q)
    right=root_sign(q,alpha)
    return (right if order%2==0 else -right,0,right)


@dataclass(frozen=True)
class QuotientGerm:
    numerator:tuple
    denominator:tuple
    denominator_order:int|None
    numerator_order:int|None
    classification:str  # persistent-denominator, pole, finite
    limit_numerator:tuple|None
    limit_denominator:tuple|None
    left_sign:int|None
    right_sign:int|None

    def finite_sign(self,alpha):
        if self.classification!='finite':raise ContinuationInconclusive('no finite event limit')
        return root_sign(self.limit_numerator,alpha)*root_sign(self.limit_denominator,alpha)


def quotient_germ(N,D,alpha):
    """Classify N/D by orders and m-th derivatives at an algebraic root.

    If ord(N)<ord(D), N/D has a true pole and cannot stay in any bounded
    source-parameter box approaching alpha. When ord(N)>=ord(D), it has a
    unique removable limit N^(m)(alpha)/D^(m)(alpha) where m=ord(D).
    The zero numerator is treated as having infinite valuation.
    """
    n,d=poly(N),poly(D)
    dn=order_at(d,alpha)
    nn=order_at(n,alpha)
    if dn is None:return QuotientGerm(n,d,None,nn,'persistent-denominator',None,None,None,None)
    if nn is not None and nn<dn:
        return QuotientGerm(n,d,dn,nn,'pole',None,None,None,None)
    a,b=n,d
    for _ in range(dn):a=diff(a);b=diff(b)
    if root_sign(b,alpha)==0:raise AssertionError('not enough derivatives removed')
    ns,ds=side_signs(n,alpha),side_signs(d,alpha)
    return QuotientGerm(n,d,dn,nn,'finite',a,b,ns[0]*ds[0],ns[2]*ds[2])


@dataclass(frozen=True)
class CramerContinuation:
    event:ExactTime
    status:str
    denominator_order:int|None
    parameter_germs:tuple
    left_owned:bool
    right_owned:bool
    left_filtered:bool
    right_filtered:bool
    event_parameters:tuple  # (numer,den) for each finite coordinate
    square_cutoff_at_event:int|None
    event_box_feasible:bool|None


def _ratio_in_closed(num,den,alpha):
    d=root_sign(den,alpha)
    return (root_sign(num,alpha)*d>=0 and
            root_sign(minus(den,num),alpha)*d>=0)


def _side_owned(N,D,alpha,which):
    """Half-open outgoing-edge u in [0,1): N/D>=0 and (D-N)/D>0."""
    sd=side_signs(D,alpha)[which]
    if not sd:return False
    return (side_signs(N,alpha)[which]*sd>=0 and
            side_signs(minus(D,N),alpha)[which]*sd>0)


def cramer_continuation(D,numerators,cutoff_numerator,alpha):
    """Exact local chart decision for a 4-variable Cramer square branch.

    cutoff_numerator is H(t) with q(t)-kappa=H(t)/D(t)^2.
    This is *only* the unique nonsingular branch for t!=alpha. An instant-only
    singular square may also exist; run v228 event-slice feasibility separately.
    """
    if len(numerators)!=4:raise ValueError('four source parameters expected')
    D=poly(D);nums=tuple(poly(x) for x in numerators);H=poly(cutoff_numerator)
    gs=tuple(quotient_germ(n,D,alpha) for n in nums)
    if iszero(D):
        return CramerContinuation(alpha,'persistent-singular',None,gs,False,False,False,False,(),None,None)
    if root_sign(D,alpha):raise ValueError('Cramer matrix is not singular at event')
    if any(g.classification=='pole' for g in gs):
        return CramerContinuation(alpha,'pole-excluded',order_at(D,alpha),gs,
                                  False,False,False,False,(),None,False)
    event_params=tuple((g.limit_numerator,g.limit_denominator) for g in gs)
    inside=all(_ratio_in_closed(n,d,alpha) for n,d in event_params)
    # Side eligibility uses exact germ signs; no rational sample substitutes
    # for an infinitesimal one-sided neighborhood.
    left=all(_side_owned(n,D,alpha,0) for n in nums)
    right=all(_side_owned(n,D,alpha,2) for n in nums)
    size=side_signs(H,alpha)
    size_at=None
    if inside:
        # Finite limits imply q(t) extends continuously through event.
        # Since q-kappa = H/D^2, its event sign is obtained from the first
        # nonzero Laurent coefficient when H vanishes at alpha.
        h=order_at(H,alpha);d=order_at(D,alpha)
        if h is None:size_at=0
        elif h<2*d:
            # A finite box limit implies H/D^2 cannot have a pole; this
            # would contradict the Cramer geometry and input polynomial.
            raise ContinuationInconclusive('cutoff polynomial inconsistent with finite geometric limit')
        elif h>2*d:size_at=0
        else:
            a,b=H,times(D,D)
            for _ in range(2*d):a,b=diff(a),diff(b)
            size_at=root_sign(a,alpha)*root_sign(b,alpha)
    return CramerContinuation(alpha,'removable' if inside else 'outside-event-box',
          order_at(D,alpha),gs,left,right,left and size[0]>0,
          right and size[2]>0,event_params,size_at,inside)


def verify_affine_event_limit(A,b,parameters,alpha):
    """Independent exact test that a finite Cramer limit solves A(alpha)u=b(alpha).

    Each coordinate may have its own rational polynomial denominator. The
    cross-multiplied identity is checked for every original equation; no
    implication of feasibility is inferred merely from a nonempty event box.
    """
    if len(A)!=4 or any(len(row)!=4 for row in A) or len(b)!=4 or len(parameters)!=4:
        raise ValueError('4x4 affine system and four finite parameters expected')
    denominators=[poly(d) for n,d in parameters]
    if any(root_sign(d,alpha)==0 for d in denominators):
        raise ContinuationInconclusive('event coordinate denominator not invertible')
    total=ONE
    for d in denominators:total=times(total,d)
    for i in range(4):
        left=ZERO
        for j,(n,d) in enumerate(parameters):
            others=ONE
            for k,den in enumerate(denominators):
                if k!=j:others=times(others,den)
            left=plus(left,times(times(poly(A[i][j]),poly(n)),others))
        if root_sign(minus(left,times(poly(b[i]),total)),alpha)!=0:
            return False
    return True


@dataclass(frozen=True)
class PersistentSupport:
    status:str
    obstruction_sign:int|None
    event_positive_square:bool|None


def exceptional_instant_support(obstruction,alpha,event_positive_square=None):
    """Classify GENERICALLY INCONSISTENT singular charts by a nonzero minor.

    If a polynomial augmented-rank minor M(t) is not identically zero,
    any square solution must have M(t)=0. A single isolated zero alpha
    supports no nontrivial one-sided branch, even if a square exists exactly
    at alpha. The event-slice existence must be decided separately by v228.
    A zero minor has no exclusion power and remains explicitly inconclusive.
    """
    obstruction=poly(obstruction)
    if iszero(obstruction):
        return PersistentSupport('unknown-persistent-rank',None,event_positive_square)
    sign=root_sign(obstruction,alpha)
    if sign:return PersistentSupport('locally-infeasible',sign,False)
    if event_positive_square is None:
        return PersistentSupport('possible-instant-only',0,None)
    return PersistentSupport('instant-only-positive' if event_positive_square else
                             'instant-only-filtered-out',0,bool(event_positive_square))


@dataclass(frozen=True)
class PolygonRankDropAudit:
    event:ExactTime
    cutoff_squared:Q
    checked_charts:int
    statuses:tuple
    side_filtered_charts:int
    positive_event_singular_charts:int
    complete_slice_certificate:bool


def audit_polygon_rank_drops(P0,P1,event,cutoff_squared,max_charts=600):
    """Audit ALL rank-drop charts of a single polygon at an algebraic event.

    For a persistently singular chart, v228 decides instantaneous feasibility
    and a nonzero augmented-rank minor may forbid punctured-neighborhood
    roots. For temporary rank loss, Cramer valuations prove removable vs
    pole. This is chart bookkeeping, not a geometric square-orbit count.
    """
    from collections import Counter
    from persistent_strata_v225 import _chart, affine_system_polynomials, singular_inconsistency_certificate
    from singular_event_feasibility_v228 import singular_polygon_chart_at_event
    from polygon_orbit_census_v223 import validate_jordan_polygon
    p=validate_jordan_polygon(tuple(P0));q=validate_jordan_polygon(tuple(P1))
    if len(p)!=len(q):raise ValueError('equal vertex count required')
    if max_charts<=0:raise ValueError('positive rank-drop budget required')
    n=len(p);counts=Counter();total=0;side=0;positive=0
    legacy=_to_legacy_time(event)
    for ids in product(range(n),repeat=4):
        if any(ids.count(i)>=3 for i in set(ids)):continue
        for sigma in (-1,1):
            chart=_chart(p,q,ids,sigma,Q(cutoff_squared))
            if root_sign(chart.determinant,event)!=0:continue
            total+=1
            if total>max_charts:raise ContinuationInconclusive('rank-drop chart budget exhausted')
            result=polygon_chart_continuation(p,q,ids,sigma,event,cutoff_squared)
            snapshot=singular_polygon_chart_at_event(p,q,ids,sigma,legacy,cutoff_squared)
            positive+=int(snapshot.has_positive_square)
            side+=int(result.left_filtered)+int(result.right_filtered)
            if result.status=='persistent-singular':
                A,b=affine_system_polynomials(p,q,ids,sigma)
                ra,rb,M,br=singular_inconsistency_certificate(
                    A,b,Q(0),Q(1),width=Q(1,2**10))
                if M is not None:
                    supp=exceptional_instant_support(M,event,snapshot.has_positive_square)
                    counts[supp.status]+=1
                else:
                    counts['persistent-rank-consistent']+=1
            else:
                counts[result.status]+=1
    return PolygonRankDropAudit(event,Q(cutoff_squared),total,
          tuple(sorted(counts.items())),side,positive,True)


def _to_legacy_time(alpha):
    from local_orbit_germs_v227 import AlgebraicTime
    return AlgebraicTime(alpha.factor,alpha.lo,alpha.hi)


def polygon_chart_continuation(P0,P1,edges,orientation,event,cutoff_squared):
    """Bridge to the complete v224/v228 polygonal square implementation.

    The returned chart decision is independently crosschecked at the event
    by v228 algebraic event-slice feasibility whenever it has a finite
    physical limit. Does not claim a full orbit cobordism across all charts.
    """
    from persistent_strata_v225 import _chart, affine_system_polynomials
    from singular_event_feasibility_v228 import (
        singular_polygon_chart_at_event, AlgebraicFeasibilityInconclusive,
    )
    chart=_chart(tuple(P0),tuple(P1),tuple(edges),orientation,Q(cutoff_squared))
    if not iszero(chart.determinant) and not root_sign(chart.determinant,event):
        status=cramer_continuation(chart.determinant,chart.numerators,chart.cutoff,event)
    elif iszero(chart.determinant):
        status=cramer_continuation(chart.determinant,chart.numerators,chart.cutoff,event)
    else:
        raise ValueError('selected chart does not lose rank at event')
    if status.status in ('removable','outside-event-box'):
        legacy=_to_legacy_time(event)
        snapshot=singular_polygon_chart_at_event(P0,P1,edges,orientation,legacy,cutoff_squared)
        if status.event_box_feasible:
            A,b=affine_system_polynomials(P0,P1,edges,orientation)
            if not verify_affine_event_limit(A,b,status.event_parameters,event):
                raise AssertionError('finite Cramer limit fails original rank-drop equations')
            if not snapshot.parameter_vertices:
                raise AssertionError('finite Cramer limit has no event-slice solution')
        if status.square_cutoff_at_event is not None and status.square_cutoff_at_event>0 and not snapshot.has_positive_square:
            raise AssertionError('positive Cramer limit conflicts with v228 square polytope')
    return status


def report():
    a=ExactTime((-Q(1,2),1),Q(1,2),Q(1,2))
    d=poly((-Q(1,2),1))
    nums=(times(d,poly((Q(1,2),))),times(d,poly((Q(1,3),))),
          times(d,poly((Q(2,3),))),times(d,poly((Q(1,4),))))
    h=times(times(d,d),poly((Q(1,8),)))
    r=cramer_continuation(d,nums,h,a)
    return {'example':'four removable Cramer parameters at t=1/2',
            'status':r.status,'event_box_feasible':r.event_box_feasible,
            'left_filtered':r.left_filtered,'right_filtered':r.right_filtered,
            'singular_instant_only_squares':'must be checked with v228 independently',
            'scope':'finite Cramer germs only; no global Toeplitz invariant'}

if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
