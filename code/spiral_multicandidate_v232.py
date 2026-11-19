"""v232: bounded multi-candidate discovery and C0 phase-robust degree certificates.

The numerical stage is NOT a theorem. It searches multiple distinct polygon
square charts, then verifies each on genuine smooth spiral arcs using v231's
exact Fraction/Brouwer-degree checker. A failed candidate is not evidence of
absence. The independent phase perturbation checker replays stored v231
rational source endpoints and never uses numerical search or a stored root.

For any continuous real phase perturbation phi with |phi(t)|<=eta on the four
source intervals, |t exp(i(t^-p+phi(t))) - t exp(i t^-p)| <= eta (t<=1).
Thus the v231 chord bound E becomes E+eta, using no derivative bound for phi.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as Q
from math import hypot

from spiral_pruned_v231 import (
    Sample, SpatialEdges, spiral_samples, solve_chart,
    local_arcs_at_newton, certify_on_smooth_spiral,
    STORED_WITNESSES, STORED_STRONG_FLOORS,
    _certified_point, _chord_error, _R_matrix, _inverse_and_root,
)


@dataclass(frozen=True)
class CandidateScan:
    samples: int
    pairs_examined: int
    pair_budget_exhausted: bool
    positive_geometry_charts: int
    distinct_candidates: tuple
    cap: int


def _source_signature(pts):
    """Approximate de-duplication for discovery ONLY; no proof relies on it."""
    return tuple(sorted(round(hypot(*p), 3) for p in pts))


def _midpoint(a, b):
    return ((a.x+b.x)/2, (a.y+b.y)/2)


def scan_candidates(p, inner=.45, step=.12, max_pairs=30000,
                    topk=4, max_candidates=16, min_side=.18,
                    min_radius=.20, min_margin=.003):
    """Rank multiple, spatially pruned square charts within hard budgets.

    Crucially, unlike v231 discover(), this does NOT return the first
    seemingly good chart. It deduplicates numerical square presentations,
    keeps the best-conditioned margin for each, and ranks several distinct
    candidate roots. All geometric tests in here are floating-point,
    diagnostic only.
    """
    if max_pairs<1 or topk<1 or max_candidates<1:
        raise ValueError('positive search budgets required')
    if not (0<min_margin<.5) or min_side<=0 or min_radius<0:
        raise ValueError('invalid numerical discovery filters')
    pts=spiral_samples(p,inner,step)
    index=SpatialEdges(pts)
    n=len(pts)-1
    pairs=0; matches=0; exhausted=False
    distinct={}
    for i in range(n):
        A=_midpoint(pts[i],pts[i+1])
        for j in range(i+2,n):
            if pairs>=max_pairs:
                exhausted=True;break
            pairs+=1
            C=_midpoint(pts[j],pts[j+1])
            vx=(C[0]-A[0])/2;vy=(C[1]-A[1])/2
            center=((A[0]+C[0])/2,(A[1]+C[1])/2)
            B=(center[0]-vy,center[1]+vx)
            D=(center[0]+vy,center[1]-vx)
            bs=index.nearest_edges(B,limit=topk)
            ds=index.nearest_edges(D,limit=topk)
            if not bs or not ds:continue
            for _,k in bs:
                for _,l in ds:
                    for ids,sg in (((i,k,j,l),-1),((i,l,j,k),1)):
                        if len(set(ids))<4:continue
                        sol=solve_chart(pts,ids,sg)
                        if sol is None:continue
                        u,square=sol
                        margin=min(min(x,1-x) for x in u)
                        if margin<min_margin:continue
                        side=hypot(square[1][0]-square[0][0],
                                   square[1][1]-square[0][1])
                        if side<=min_side or min(hypot(*z) for z in square)<=min_radius:
                            continue
                        matches+=1
                        key=_source_signature(square)
                        # Geometry and margin influence search order only.
                        score=margin * min(side,1.)
                        cand={'status':'found','p':p,'inner':inner,'step':step,
                              'samples':n,'pairs':pairs,'ids':ids,'sigma':sg,
                              'parameters':u,'points':square,'side':side,
                              'samples_data':pts,'search_score':score}
                        old=distinct.get(key)
                        if old is None or score>old['search_score']:
                            distinct[key]=cand
                        # Bound memory, not just output length.
                        if len(distinct)>max_candidates*4:
                            keep=sorted(distinct.items(),
                              key=lambda kv:(-kv[1]['search_score'],kv[0]))[:max_candidates*2]
                            distinct=dict(keep)
        if exhausted:break
    ordered=tuple(sorted(distinct.values(),
            key=lambda c:(-c['search_score'],tuple(c['ids'])))[:max_candidates])
    return CandidateScan(n,pairs,exhausted,matches,ordered,max_candidates)


def certify_candidates(candidates, half_angles=(.02,.01,.005,.002),
                       radius_floor=Q(1,5),side_floor=Q(1,8),
                       max_checks=32):
    """Try several independent proposals; never convert failure to nonexistence.

    Exposes an explicit attempt/failure ledger. The only 'certified' return
    is a successful v231 exact Fraction/Brouwer-degree proof.
    """
    if max_checks<1:raise ValueError('max_checks must be positive')
    failures=Counter();attempted=0
    for pos,cand in enumerate(candidates):
        for half in half_angles:
            if attempted>=max_checks:
                return {'certified':False,'status':'budget-inconclusive',
                        'attempted':attempted,'rejections':dict(failures)}
            attempted+=1
            try:
                local=local_arcs_at_newton(cand,half_angle=half)
                proof=certify_on_smooth_spiral(
                    local,radius_floor=Q(radius_floor),side_floor=Q(side_floor))
                if proof.get('certified'):
                    proof.update({'status':'certified','candidate_rank':pos+1,
                                  'half_angle':half,'attempted':attempted,
                                  'rejections':dict(failures)})
                    return proof
                failures[proof.get('reason','uncertified')]+=1
            except (ValueError,OverflowError,ZeroDivisionError,ArithmeticError) as ex:
                failures[type(ex).__name__]+=1
    return {'certified':False,'status':'candidate-exhausted-inconclusive',
            'attempted':attempted,'rejections':dict(failures)}


def adaptive_multicandidate(p, radius_ladder=(.55,.45,.36,.30,.22),
                            step=.12, max_pairs=30000, max_candidates=16,
                            max_checks_per_radius=32,
                            half_angles=(.02,.01,.005,.002)):
    """Bounded numerical discovery + exact certificate; no saved witness lookup.

    Failing to find a certifiable candidate is always INCONCLUSIVE.
    """
    ledger=[]
    for radius in radius_ladder:
        scan=scan_candidates(p,radius,step,max_pairs=max_pairs,
                             max_candidates=max_candidates)
        result=certify_candidates(scan.distinct_candidates,
                                  half_angles=half_angles,
                                  max_checks=max_checks_per_radius)
        ledger.append({'radius':radius,'pairs':scan.pairs_examined,
                       'budget_exhausted':scan.pair_budget_exhausted,
                       'charts':scan.positive_geometry_charts,
                       'distinct':len(scan.distinct_candidates),
                       'attempted':result['attempted'],
                       'rejections':result['rejections']})
        if result['certified']:
            result['search_ladder']=ledger
            result['searched_from_radius']=radius
            return result
    return {'certified':False,'status':'search-inconclusive','p':str(p),
            'search_ladder':ledger}


@dataclass(frozen=True)
class PhasePerturbationCertificate:
    exponent: Q
    phase_amplitude_at_most: Q
    side_strict_floor: Q
    radius_strict_floor: Q
    margin: Q
    four_E_total: Q
    exact_degree_valid: bool


def certify_phase_neighborhood(p, eta=Q(1,10**6)):
    """Replay v231 rational endpoints; certify ANY continuous bounded phase phi.

    On the original four source intervals, |phi|<=eta implies coordinatewise
    perturbation <=t|phi|<=eta, hence 4(E+eta)<m/||A^-1|| retains local degree.
    The v231 side/radius floors are checked anew with E+eta, not merely
    copied from prior float results. No assumptions on phi's derivative.
    """
    p=Q(str(p));eta=Q(eta)
    if p not in STORED_WITNESSES:raise ValueError('no exact rational witness for exponent')
    if eta<0:raise ValueError('phase bound must be nonnegative')
    sg,enc=STORED_WITNESSES[p]
    side_floor,radius_floor=STORED_STRONG_FLOORS[p]
    src=tuple(Sample(Q(x,2**30)**2,0.,0.) for x in enc)
    ends=[(_certified_point(src[j],p),_certified_point(src[j+1],p))
          for j in (0,2,4,6)]
    a=tuple(pair[0][0] for pair in ends)
    b=tuple(pair[1][0] for pair in ends)
    A,rhs=_R_matrix(a,b,sg)
    u,inv_norm=_inverse_and_root(A,rhs)
    m=min(min(t,1-t) for t in u)
    if m<=0:raise ValueError('stored chart root not in box')
    base_E=max(max(e[0][1],e[1][1])+
               _chord_error(src[j].t,src[j+1].t,p)
               for j,e in zip((0,2,4,6),ends))
    total_E=base_E+eta
    margin=m/inv_norm
    square=tuple((a[j][0]+u[j]*(b[j][0]-a[j][0]),
                  a[j][1]+u[j]*(b[j][1]-a[j][1])) for j in range(4))
    side2=(square[1][0]-square[0][0])**2+(square[1][1]-square[0][1])**2
    side_ok=side2>(side_floor+4*total_E)**2
    radius_ok=all(x*x+y*y>(radius_floor+2*total_E)**2 for x,y in square)
    valid=(4*total_E<margin and side_ok and radius_ok)
    return PhasePerturbationCertificate(p,eta,side_floor,radius_floor,
                                        margin,4*total_E,valid)


def report():
    return {'phase_neighborhoods':[
        {'p':str(p),'phase_amplitude':str(c.phase_amplitude_at_most),
         'degree_valid':c.exact_degree_valid,
         'side_floor':str(c.side_strict_floor),
         'radius_floor':str(c.radius_strict_floor),
         'strict_margin':bool(c.four_E_total<c.margin)}
        for p in (Q(2),Q(5,2),Q(3))
        for c in (certify_phase_neighborhood(p),)],
       'scope':'rigorous a posteriori C0 phase neighborhood on three source boxes;'
               ' numerical candidate search remains bounded and incomplete'}

if __name__=='__main__':
    import json
    out=report()
    print(json.dumps(out,indent=2,sort_keys=True))
    if not all(r['degree_valid'] for r in out['phase_neighborhoods']):
        raise SystemExit('FAIL: proposed phase amplitude exceeds exact four-arc margin')
