"""v227: exact event-germ incidence for isolated polygon square orbits.

The model is a *single* rational Jordan polygon with vertices affine in time.
For an event alpha from v226, test square branches from every nonsingular
four-edge chart using rational-polynomial algebraic signs at alpha and its
one-sided germs. Group descriptions by the actual geometric square at alpha;
then match eligible branch germs to distinct geometric square orbits in the
adjacent event-free chambers, checked independently by v223.

IMPORTANT: The result is fail-closed. Rank drops, positive-dimensional square
families, or a side orbit without a nonsingular event anchor are inconclusive.
No unrestricted Jordan curve invariant or Toeplitz proof is claimed.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations

from algebraic_orbit_events_v226 import (
    EventCoverInconclusive, EventCluster, OrbitEventAudit,
    audit_moving_polygon_orbits, _has_root_in_closed,
)
from moving_polygon_events_v224 import (
    polynomial, ZERO, iszero, derivative, evaluate, poly_gcd,
    sturm_sequence, count_roots_open, isolate_roots,
    plus, minus, times, moving_point, pairminus,
    interpolated_polygon,
)
from persistent_strata_v225 import (
    classify_polygon_motion, _chart, owned_square_at_rational_time,
)
from polygon_orbit_census_v223 import (
    census_polygon, PositiveDimensionalSquares,
)


def _sgn(x):return (x>0)-(x<0)


class GermInconclusive(EventCoverInconclusive):
    """Cannot certify every physical orbit meeting the event."""


@dataclass(frozen=True)
class AlgebraicTime:
    factor: tuple
    lo: Q
    hi: Q
    def __post_init__(self):
        if self.lo>self.hi:raise ValueError('invalid algebraic time bracket')
        if not _has_root_in_closed(self.factor,self.lo,self.hi):
            raise ValueError('bracket does not isolate a root')
        if self.lo<self.hi and count_roots_open(sturm_sequence(self.factor),self.lo,self.hi)!=1:
            raise ValueError('multiple/noninterior roots in isolator')
    @property
    def rational(self):return self.lo==self.hi


@lru_cache(maxsize=65536)
def sign_at_time(p,alpha:AlgebraicTime):
    """Exact sign of rational polynomial p(alpha), using gcd and Sturm.

    Never evaluates at a floating-point estimate of an irrational root.
    """
    p=polynomial(p)
    if iszero(p):return 0
    if alpha.rational:return _sgn(evaluate(p,alpha.lo))
    g=poly_gcd(p,alpha.factor)
    if len(g)>1 and _has_root_in_closed(g,alpha.lo,alpha.hi):return 0
    lo,hi=alpha.lo,alpha.hi
    # Refine the defining root until p is root-free throughout its bracket.
    for _ in range(128):
        if not _has_root_in_closed(p,lo,hi):
            value=evaluate(p,(lo+hi)/2)
            if not value:raise AssertionError('root-free interval has zero midpoint')
            return _sgn(value)
        subbr=isolate_roots(alpha.factor,lo,hi,(hi-lo)/4)
        if len(subbr)!=1:raise GermInconclusive('could not refine algebraic event root')
        nxt=subbr[0]
        if nxt.lower==lo and nxt.upper==hi:
            raise GermInconclusive('root bracket did not shrink')
        lo,hi=nxt.lower,nxt.upper
    raise GermInconclusive('sign separation budget exhausted')


@lru_cache(maxsize=65536)
def germ_signs(p,alpha:AlgebraicTime):
    """(left, event, right) exact signs of one polynomial near alpha."""
    p=polynomial(p)
    if iszero(p):return (0,0,0)
    at=sign_at_time(p,alpha)
    if at:return (at,at,at)
    q=p;m=0
    while sign_at_time(q,alpha)==0:
        q=derivative(q)
        m+=1
        if iszero(q):raise AssertionError('nonzero polynomial has all derivatives zero')
    right=sign_at_time(q,alpha)
    return (right if m%2==0 else -right,0,right)


def _ratio_sign(numer,denom,alpha,side):
    n=germ_signs(numer,alpha)[side]
    d=germ_signs(denom,alpha)[side]
    if d==0:raise GermInconclusive('determinant zero in infinitesimal germ')
    return n*d


def _points_numerators(P0,P1,chart):
    """Four physical points at alpha as rational polynomial numerator/D."""
    n=len(P0)
    vertices=tuple(moving_point(a,b) for a,b in zip(P0,P1))
    points=[]
    for i,edge in enumerate(chart.edges):
        base=vertices[edge]
        direction=pairminus(vertices[(edge+1)%n],base)
        N=chart.numerators[i];D=chart.determinant
        points.append(tuple(plus(times(D,base[j]),times(N,direction[j]))
                            for j in range(2)))
    return tuple(points)


def _equal_rational_functions(a,da,b,db,alpha):
    return sign_at_time(minus(times(a,db),times(b,da)),alpha)==0


def _same_geometric_orbit(root1,root2,alpha):
    """Equality of unlabelled vertex sets over Q(alpha), *exactly*."""
    p,q=root1.point_numerators,root2.point_numerators
    for perm in permutations(range(4)):
        if all(_equal_rational_functions(p[i][k],root1.determinant,
                                         q[perm[i]][k],root2.determinant,alpha)
               for i in range(4) for k in range(2)):
            return True
    return False


def _geometric_key_at_rational(root,alpha):
    if not alpha.rational:raise ValueError('rational time required')
    t=alpha.lo
    d=evaluate(root.determinant,t)
    if not d:raise GermInconclusive('singular event representation')
    return tuple(sorted((evaluate(px,t)/d,evaluate(py,t)/d)
                        for px,py in root.point_numerators))


@dataclass(frozen=True)
class ChartGerm:
    edges:tuple
    orientation:int
    determinant:tuple
    point_numerators:tuple
    side_at_event:str   # "above-cutoff" or "on-cutoff"
    left_key:tuple|None
    right_key:tuple|None
    at_vertex:bool


@dataclass(frozen=True)
class SquareEventGerm:
    event_key:tuple|None
    size_status:str
    chart_presentations:int
    left_keys:tuple
    right_keys:tuple
    classification:str
    vertex_incidence:bool

    @property
    def left_count(self):return len(self.left_keys)
    @property
    def right_count(self):return len(self.right_keys)


@dataclass(frozen=True)
class EventGermLedger:
    algebraic_time:AlgebraicTime
    raw_candidate_entries:int
    groups:tuple
    left_chamber_count:int
    right_chamber_count:int
    event_slice_count:int|None
    def balances(self):
        return (sum(g.left_count for g in self.groups)==self.left_chamber_count and
                sum(g.right_count for g in self.groups)==self.right_chamber_count)
    @property
    def net_change(self):return self.right_chamber_count-self.left_chamber_count


def _classify(left,right,at_cutoff):
    if at_cutoff:
        if left==0 and right==1:return 'size-entry'
        if left==1 and right==0:return 'size-exit'
        if left==right==0:return 'cutoff-tangency-or-isolated'
        return 'multiple-cutoff-germs'
    if left==right==1:return 'continuation'
    if left==0 and right==2:return 'pair-birth'
    if left==2 and right==0:return 'pair-death'
    if left==right==0:return 'isolated-event-root'
    if left==right and left>1:return 'reconnection'
    return 'other-boundary-incidence'


def _chart_germ(P0,P1,cert,alpha,left_time,right_time,cutoff):
    """Return a nonsingular feasible square *at* event and its owned sides."""
    c=_chart(tuple(P0),tuple(P1),cert.edges,cert.orientation,Q(cutoff))
    if iszero(c.determinant) or sign_at_time(c.determinant,alpha)==0:
        return None
    d=sign_at_time(c.determinant,alpha)
    if any(sign_at_time(n,alpha)*d<0 or sign_at_time(u,alpha)*d<0
           for n,u in zip(c.numerators,c.upper_differences)):
        return None
    hs=sign_at_time(c.cutoff,alpha)
    if hs<0:return None
    # For each one-sided germ, outgoing vertex u=0 is OWNED; u=1 is not.
    def get_side(idx,t):
        if germ_signs(c.cutoff,alpha)[idx]<=0:
            # q-k=H/D^2, so its sign is exactly sign(H).
            return None
        if any(_ratio_sign(n,c.determinant,alpha,idx)<0 or
               _ratio_sign(u,c.determinant,alpha,idx)<=0
               for n,u in zip(c.numerators,c.upper_differences)):
            return None
        # Avoid returning nonphysical geometry through a sign-classification bug.
        out=owned_square_at_rational_time(P0,P1,cert.edges,cert.orientation,t,cutoff)
        if out is None:
            raise GermInconclusive('event germ disagrees with full chamber sample')
        pts=out[0]
        return tuple(sorted((p.x,p.y) for p in pts))
    left=get_side(0,left_time)
    right=get_side(2,right_time)
    vertex=any(sign_at_time(x,alpha)==0 for x in c.numerators+c.upper_differences)
    return ChartGerm(cert.edges,cert.orientation,c.determinant,
                     _points_numerators(P0,P1,c),
                     'on-cutoff' if hs==0 else 'above-cutoff',
                     left,right,vertex)


def classify_rational_or_algebraic_event(P0,P1,cutoff,base_audit:OrbitEventAudit,
                                        event_index,strata=None):
    """Match all regular event roots and BOTH nearby geometric-orbit chambers.

    Reuses v226 global candidate cover and its rigorously event-free chambers.
    If a chamber orbit lacks an event-limit anchor, or an event-time square
    lacks nonsingular chart support, fail closed. No numerical root guesses.
    """
    i=int(event_index)
    if not 0<=i<len(base_audit.clusters):raise ValueError('invalid event index')
    ev=base_audit.events[i];c=ev.cluster
    alpha=AlgebraicTime(c.minimal_factor,c.bracket.lower,c.bracket.upper)
    left,right=base_audit.chambers[i],base_audit.chambers[i+1]
    if strata is None:
        strata=classify_polygon_motion(P0,P1,cutoff,base_audit.start,base_audit.end,
                                       Q(1,2**10),True)
    if strata.unresolved_charts:raise GermInconclusive('singular chart family remains')
    if not alpha.rational:
        # A singular event-time chart might carry an instant-only square that
        # cannot be detected from chamber samples. Without exact real-algebraic
        # semialgebraic feasibility, refuse a *complete* event-germ claim.
        for cert in strata.certificates:
            if cert.status=='structurally-empty':continue
            if cert.status in ('always-empty','exceptional-time-only'):
                if cert.obstruction and not iszero(cert.obstruction) and \
                   sign_at_time(cert.obstruction,alpha)==0:
                    raise GermInconclusive('possibly feasible exceptional algebraic-time chart')
                continue
            if cert.status in ('regular','vertex-owned','incoming-edge-duplicate'):
                ch=_chart(tuple(P0),tuple(P1),cert.edges,cert.orientation,Q(cutoff))
                if sign_at_time(ch.determinant,alpha)==0:
                    raise GermInconclusive('singular matrix at irrational event time')
    roots=[]
    for cert in strata.certificates:
        if cert.status not in ('regular','vertex-owned','incoming-edge-duplicate'):
            continue
        germ=_chart_germ(P0,P1,cert,alpha,left.sample,right.sample,cutoff)
        if germ is not None:roots.append(germ)
    groups=[]
    for root in roots:
        match=next((group for group in groups if _same_geometric_orbit(group[0],root,alpha)),None)
        if match is None:groups.append([root])
        else:match.append(root)
    answers=[]
    for g in groups:
        left_keys=tuple(sorted({x.left_key for x in g if x.left_key is not None}))
        right_keys=tuple(sorted({x.right_key for x in g if x.right_key is not None}))
        statuses={x.side_at_event for x in g}
        if len(statuses)!=1:raise GermInconclusive('inconsistent size classification across equivalent charts')
        cutoff_root=statuses=={'on-cutoff'}
        answers.append(SquareEventGerm(_geometric_key_at_rational(g[0],alpha) if alpha.rational else None,
                       statuses.pop(),len(g),left_keys,right_keys,
                       _classify(len(left_keys),len(right_keys),cutoff_root),
                       any(x.at_vertex for x in g)))
    total_left={k for x in answers for k in x.left_keys}
    total_right={k for x in answers for k in x.right_keys}
    if total_left!=set(left.orbit_keys) or total_right!=set(right.orbit_keys):
        raise GermInconclusive('a chamber orbit has no certified nonsingular limit at event')
    if sum(x.left_count for x in answers)!=left.count or sum(x.right_count for x in answers)!=right.count:
        raise GermInconclusive('one chamber orbit assigned to multiple event roots')
    if alpha.rational:
        try:at=census_polygon(interpolated_polygon(P0,P1,alpha.lo),cutoff)
        except PositiveDimensionalSquares as exc:
            raise GermInconclusive('positive-dimensional event-time square family') from exc
        expected={tuple(sorted((p.x,p.y) for p in s.vertices)) for s in at.squares}
        actual={x.event_key for x in answers if x.size_status=='above-cutoff'}
        if actual!=expected:raise GermInconclusive('event polygon contains a square with no regular-germ anchor')
    result=EventGermLedger(alpha,c.raw_entries,tuple(answers),left.count,right.count,ev.at_event_count)
    if not result.balances():raise GermInconclusive('local incidence does not cover chamber orbits')
    return result


def audit_local_event_germs(P0,P1,cutoff,start=0,end=1,width=Q(1,2**10)):
    global_audit=audit_moving_polygon_orbits(P0,P1,cutoff,start,end,width)
    strata=classify_polygon_motion(P0,P1,cutoff,start,end,width,True)
    return global_audit,tuple(classify_rational_or_algebraic_event(
         P0,P1,cutoff,global_audit,i,strata) for i in range(len(global_audit.events)))


def report():
    from polygon_orbit_census_v223 import pt,convex_vertex_alignment_polygon
    P0=convex_vertex_alignment_polygon(Q(0))
    P1=convex_vertex_alignment_polygon(Q(9,8))
    base,events=audit_local_event_germs(P0,P1,Q(1,100),Q(3,4),1,Q(1,2**10))
    return {
      'convex_raw_entries':base.raw_candidates,
      'convex_event_time':str(events[0].algebraic_time.lo),
      'convex_groups':[(g.classification,g.left_count,g.right_count,
                        g.chart_presentations) for g in events[0].groups],
      'net_change':events[0].net_change,
      'scope':'local event germs on single affine-moving polygon; no all-Jordan global invariant'
    }

if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,default=str))
