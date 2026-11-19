"""v226: exact algebraic event-root CLUSTERING and orbit census by chambers.

For a family of rational polygons with vertices affine in time, v224 produces
necessary polynomial event candidates. v225 classifies persistent charts.
This module identifies candidate roots as algebraic numbers by exact gcd +
Sturm root counting, then performs the independent v223 COMPLETE static orbit
census in the open chambers. It does not invent a global invariant or a
cutoff-independent Jordan-curve result.

Finite cases with persistently consistent singular square families, persistent
cutoff strata, or unresolved rank-drop events remain explicit exceptions.
The census at an irrational event instant is NOT inferred from neighboring
rational samples; only chamber orbit counts are certified in that case.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import reduce
from collections import Counter, defaultdict

from moving_polygon_events_v224 import (
    RootBracket, polynomial, iszero, squarefree, poly_gcd,
    sturm_sequence, count_roots_open, isolate_roots, evaluate,
    interpolated_polygon,
)
from persistent_strata_v225 import (
    classify_polygon_motion, GlobalStrataAudit,
)
from polygon_orbit_census_v223 import (
    census_polygon, PositiveDimensionalSquares,
    validate_jordan_polygon,
)


class EventCoverInconclusive(RuntimeError):
    """No assertion of a complete algebraic event cover can be made."""


@dataclass(frozen=True)
class AlgebraicCandidate:
    polynomial: tuple
    bracket: RootBracket
    sources: tuple  # (edge tuple, orientation, kind) triples


@dataclass(frozen=True)
class EventCluster:
    minimal_factor: tuple  # gcd of all defining polynomials; not irreducible in general
    bracket: RootBracket
    sources: tuple
    raw_entries: int

    @property
    def exact_time(self):
        return self.bracket.lower if self.bracket.rational else None


@dataclass(frozen=True)
class OrbitChamber:
    lower: Q
    upper: Q
    sample: Q
    count: int
    orbit_keys: tuple
    square_sizes_squared: tuple
    sectors: tuple


@dataclass(frozen=True)
class OrbitEvent:
    cluster: EventCluster
    left_count: int
    right_count: int
    at_event_count: int | None
    left_parity: int
    right_parity: int
    net_change: int
    has_parity_change: bool


@dataclass(frozen=True)
class OrbitEventAudit:
    start: Q
    end: Q
    cutoff_squared: Q
    charts_total: int
    raw_candidates: int
    clusters: tuple
    chambers: tuple
    events: tuple
    unresolved_chart_count: int

    @property
    def is_complete_chamber_cover(self):
        return self.unresolved_chart_count == 0


def _canonical_polynomial(p):
    p=squarefree(polynomial(p))
    if iszero(p):raise ValueError('cannot isolate identically zero polynomial')
    return polynomial(v/p[-1] for v in p)


def _has_root_in_closed(p, lo, hi):
    """Exact algebraic-existence test on a rational CLOSED interval."""
    if iszero(p):raise ValueError('gcd polynomial cannot vanish identically')
    if lo>hi:return False
    if lo==hi:return evaluate(p,lo)==0
    sq=squarefree(p)
    if len(sq)==1:return False
    return (evaluate(sq,lo)==0 or evaluate(sq,hi)==0 or
            count_roots_open(sturm_sequence(sq),lo,hi)>0)


def _overlap(a,b):
    return max(a.lower,b.lower)<=min(a.upper,b.upper)


def _same_root(p,a,q,b):
    """Two certified one-root intervals represent same algebraic time iff
    their gcd has a root in the intersection. This distinguishes overlapping
    irrational isolators for nearby, *different* roots without float tests.
    """
    if not _overlap(a,b):return False
    g=poly_gcd(p,q)
    return len(g)>1 and _has_root_in_closed(g,max(a.lower,b.lower),
                                            min(a.upper,b.upper))


def _refine(poly, br, target):
    if br.rational:return br
    if target<=0:raise ValueError('positive requested root width required')
    roots=isolate_roots(poly,br.lower,br.upper,width=target)
    if len(roots)!=1:
        raise AssertionError('certified bracket contained other roots')
    return roots[0]


def _joint_isolator(polys, brackets, width):
    g=reduce(poly_gcd, polys)
    if len(g)<=1:raise AssertionError('cluster has no shared algebraic root')
    lo=max(b.lower for b in brackets);hi=min(b.upper for b in brackets)
    if lo>hi or not _has_root_in_closed(g,lo,hi):
        raise AssertionError('cluster intervals do not share defining root')
    if lo==hi:return _canonical_polynomial(g),RootBracket(lo,hi)
    g=_canonical_polynomial(g)
    if len(g)==2:
        t=-g[0]/g[1]
        if lo<=t<=hi:return g,RootBracket(t,t)
    # If an input bracket is exact, its root must be a root of every polynomial.
    for b in brackets:
        if b.rational:
            t=b.lower
            if lo<=t<=hi and all(evaluate(p,t)==0 for p in polys):
                return g,RootBracket(t,t)
    # Clipped interval may contain one root at an endpoint; isolate the
    # interior in the original defining bracket if the endpoint is not a root.
    if evaluate(g,lo)==0:return g,RootBracket(lo,lo)
    if evaluate(g,hi)==0:return g,RootBracket(hi,hi)
    roots=isolate_roots(g,lo,hi,width=width)
    if len(roots)!=1:raise AssertionError('nonunique root in joint isolator')
    return g,roots[0]


def cluster_candidate_roots(atoms, width=Q(1,2**18), max_refinements=64):
    """Merge all exactly equal algebraic candidate roots from many charts.

    Each input AlgebraicCandidate holds a proven *one-root* interval; multiple
    chart labels can share one candidate object. Distinct roots with initially
    overlapping brackets are separated by exact refinement before returning.
    No rational midpoint is substituted for an irrational event time.
    """
    width=Q(width)
    if width<=0:raise ValueError('positive width required')
    atoms=tuple(atoms)
    if not atoms:return ()
    for x in atoms:
        if not isinstance(x,AlgebraicCandidate) or not x.sources:
            raise ValueError('exact candidate atom and nonempty provenance required')
        if not _has_root_in_closed(x.polynomial,x.bracket.lower,x.bracket.upper):
            raise ValueError('input is not a certified root bracket')
        if not x.bracket.rational and count_roots_open(
              sturm_sequence(x.polynomial),x.bracket.lower,x.bracket.upper)!=1:
            raise ValueError('input interval is not a one-root isolator')
    n=len(atoms);parent=list(range(n))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(i,j):
        a,b=find(i),find(j)
        if a!=b:parent[b]=a
    for i in range(n):
        for j in range(i+1,n):
            if _same_root(atoms[i].polynomial,atoms[i].bracket,
                          atoms[j].polynomial,atoms[j].bracket):union(i,j)
    groups=defaultdict(list)
    for i in range(n):groups[find(i)].append(atoms[i])
    clusters=[]
    for group in groups.values():
        pol,br=_joint_isolator(tuple(x.polynomial for x in group),
                               tuple(x.bracket for x in group),width)
        sources=tuple(s for a in group for s in a.sources)
        clusters.append(EventCluster(pol,br,sources,len(sources)))
    # Avoid declaring an order from overlapping intervals for distinct roots.
    # The refining polynomial itself is exact, so a coincident root would
    # have been grouped above. Shrink intervals until disjoint.
    for _ in range(max_refinements):
        clusters.sort(key=lambda c:(c.bracket.lower,c.bracket.upper))
        bad=next((i for i in range(len(clusters)-1)
                  if clusters[i].bracket.upper>=clusters[i+1].bracket.lower),None)
        if bad is None:return tuple(clusters)
        x,y=clusters[bad:bad+2]
        def shrink(c):
            b=c.bracket
            return EventCluster(c.minimal_factor,
                   _refine(c.minimal_factor,b,(b.upper-b.lower)/4),
                   c.sources,c.raw_entries)
        if x.bracket.rational and y.bracket.rational:
            raise AssertionError('distinct clusters have same rational time')
        if not x.bracket.rational:clusters[bad]=shrink(x)
        if not y.bracket.rational:clusters[bad+1]=shrink(y)
    raise EventCoverInconclusive('could not separate distinct algebraic event times in budget')


def candidate_atoms_from_strata(audit:GlobalStrataAudit):
    """Create exact (polynomial, one-root bracket) atoms with provenance.

    An identically cutoff-locked chart may hide singular-time special roots:
    fail closed instead of silently classifying its entire deformation.
    """
    if audit.unresolved_charts:
        raise EventCoverInconclusive(f'{len(audit.unresolved_charts)} persistent singular charts')
    bad=[c for c in audit.certificates if c.status in ('genuinely-unresolved','cutoff-locked')]
    if bad:raise EventCoverInconclusive(f'{len(bad)} special charts cannot be certified')
    duplicates=defaultdict(list)
    total=0
    for c in audit.certificates:
        mp=dict(c.active_polynomials)
        for kind,br in c.exceptional_roots:
            p=c.obstruction if kind=='augmented-rank' else mp.get(kind)
            if p is None:raise AssertionError('candidate polynomial provenance missing')
            p=_canonical_polynomial(p)
            key=(p,br.lower,br.upper)
            duplicates[key].append((c.edges,c.orientation,kind))
            total+=1
    if total!=audit.raw_candidate_events:
        raise AssertionError('not all algebraic candidates retained')
    return tuple(AlgebraicCandidate(k[0],RootBracket(k[1],k[2]),tuple(v))
                 for k,v in duplicates.items())


def _census_chamber(P0,P1,lo,hi,cutoff):
    if not lo<hi:raise AssertionError('zero-length chamber')
    sample=(lo+hi)/2
    P=interpolated_polygon(P0,P1,sample)
    try:result=census_polygon(P,cutoff)
    except PositiveDimensionalSquares as e:
        raise EventCoverInconclusive(f'positive-dimensional square family at t={sample}') from e
    return OrbitChamber(lo,hi,sample,len(result.squares),
             tuple(tuple(sorted((p.x,p.y) for p in s.vertices)) for s in result.squares),
             tuple(s.side_squared for s in result.squares),
             tuple(s.cyclic_order_sector for s in result.squares))


def audit_moving_polygon_orbits(P0,P1,cutoff_squared,start=0,end=1,
                               width=Q(1,2**12),max_charts=100000):
    """Certified orbit counts on EVERY event-free open time chamber.

    Assumes v225 returns complete chart classifications (no persistent
    singular families). v225's Jordan isotopy condition is checked on the
    whole [start,end]. The isolated algebraic roots include *spurious*
    candidates; each chamber count is instead obtained independently from
    the complete v223 geometric-square census, not by summing chart labels.

    At a RATIONAL event time, also census the exact polygon there.
    At IRRATIONAL times, the event slice count is UNKNOWN; do not pretend
    a rational approximation decides it. The chamber counts are valid
    because all potential changes are covered by the exact algebraic roots.
    """
    a,b,k=Q(start),Q(end),Q(cutoff_squared)
    if not 0<=a<b<=1 or k<=0:raise ValueError('time bounds in [0,1] and positive cutoff')
    P0=validate_jordan_polygon(tuple(P0));P1=validate_jordan_polygon(tuple(P1))
    if len(P0)!=len(P1):raise ValueError('same vertex count required')
    strata=classify_polygon_motion(P0,P1,k,a,b,width,True,max_charts)
    atoms=candidate_atoms_from_strata(strata)
    clusters=cluster_candidate_roots(atoms,width)
    # Any event exactly at an endpoint prevents asserting the endpoint
    # count as an ordinary generic count. Chamber counts still well-defined.
    interior=[c for c in clusters if a<c.bracket.lower and c.bracket.upper<b]
    if len(interior)!=len(clusters):
        raise EventCoverInconclusive('event at outer time boundary: shrink time interval')
    chambers=[];last=a
    for c in clusters:
        if last<c.bracket.lower:
            chambers.append(_census_chamber(P0,P1,last,c.bracket.lower,k))
        last=c.bracket.upper
    if last<b:chambers.append(_census_chamber(P0,P1,last,b,k))
    if len(chambers)!=len(clusters)+1:
        raise EventCoverInconclusive('event isolator exhausts a chamber')
    events=[]
    for i,c in enumerate(clusters):
        lc,rc=chambers[i].count,chambers[i+1].count
        at=None
        if c.exact_time is not None:
            try:at=census_polygon(interpolated_polygon(P0,P1,c.exact_time),k).count
            except PositiveDimensionalSquares:at=None
        events.append(OrbitEvent(c,lc,rc,at,lc%2,rc%2,rc-lc,(lc-rc)%2!=0))
    return OrbitEventAudit(a,b,k,len(strata.certificates),strata.raw_candidate_events,
                           clusters,tuple(chambers),tuple(events),len(strata.unresolved_charts))


def report():
    from polygon_orbit_census_v223 import convex_vertex_alignment_polygon
    P0=convex_vertex_alignment_polygon(Q(0))
    P1=convex_vertex_alignment_polygon(Q(9,8))
    r=audit_moving_polygon_orbits(P0,P1,Q(1,100),Q(3,4),Q(1),Q(1,2**11))
    return {
        'raw_candidates':r.raw_candidates,
        'distinct_algebraic_times':len(r.clusters),
        'event_time':str(r.clusters[0].exact_time),
        'event_chamber_orbit_counts':tuple(c.count for c in r.chambers),
        'orbit_count_at_event':r.events[0].at_event_count,
        'net_orbit_change':r.events[0].net_change,
        'parity_changed':r.events[0].has_parity_change,
        'strata_unresolved':r.unresolved_chart_count,
        'scope':'exact event-time clustering and chamber orbit census for finite rational Jordan polygon motions; not universal Toeplitz',
    }


if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
