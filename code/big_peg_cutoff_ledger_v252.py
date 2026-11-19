"""v252: COMPLETE exact high-size square-orbit event ledger, on ONE moving 9-gon.

For the genuinely moving strictly convex family introduced in v250 and
I=[1/2,501/1000], cutoff side² > 94257/10000, rigorously classify ALL
126 possible increasing, distinct-edge square charts, including every
real parameter, NOT just sample times.

A geometric max-edge-length bound excludes repeated-edge corners.
Convexity allows the four boundary-ordered corners to be assigned to
strictly increasing edge indices (126 = C(9,4)). Each quadruple is
excluded by either a global pair-edge diameter bound or an exact
Bernstein polynomial certifying that a corner is off its edge, or the
candidate is below the size cutoff. The only survivors are v251's
two geometric square branches A and B.

At the LEFT endpoint, a square corner lies on a polygon vertex and
has two supporting-edge descriptions. An independent v223 EXACT complete
geometric-orbit census is used to check and deduplicate this endpoint.
The right endpoint is independently audited too.

The squared-side cutoff events of A/B are uniquely isolated by exact
positive/negative derivatives and rational sign brackets. Consequently
the # of geometric squares STRICTLY ABOVE the size cutoff is 1,2,1 as
time increases. Therefore the parity of this size-filtered count is
NOT a homotopy invariant, even in this convex, nondegenerate example.
The change in filtered count is due to crossings of the imposed size
cutoff, NOT births/deaths of actual square orbits.

Exact rational signs throughout; no float proof decisions. This is a
complete ledger ONLY for a FIXED positive cutoff and specific affine
9-gon interval, not a general polygonal homotopy solver or Toeplitz proof.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json

from polygon_orbit_census_v223 import pt, census_polygon, dot
from big_peg_persistence_v248 import square_corner
from big_peg_multichart_v249 import (
    pol, padd, psub, pneg, pmul, pscale, peval, bcoeff,
    certify_sign, chart_polynomials, corner_numerators, ppoint_at,
    pvsub, pvdot, root_at, HALF,
)
from big_peg_corner_v250 import (
    family_vertices, prove_jordan_family,
)
from big_peg_size_atlas_v251 import (
    BRANCH_A, BRANCH_B, WINDOW_START, WINDOW_END, HIGH_S2,
    side_comparison_polynomial, exact_switch_proof,
)

EXPECTED_CLASSIFICATION={
    "geometry":87,
    "off_edge":36,
    "below_cutoff":1,
    "active":2,
}
EXPECTED_ORBITS=(BRANCH_A,BRANCH_B)
B_ENTER_BRACKET=(Q(50053,100000),Q(50054,100000))
A_EXIT_BRACKET=(Q(500551,1000000),Q(50056,100000))


def _strict_sign_by_bernstein(poly,a,b):
    """Prove strictly positive or strictly negative on CLOSED interval."""
    coefficients=bcoeff(poly,a,b)
    if all(x>0 for x in coefficients):
        return 1
    if all(x<0 for x in coefficients):
        return -1
    raise ValueError("determinant sign cannot be certified on entire interval")


def _positive_except_possible_left(poly,a,b):
    """Prove p(t)>0 for t in (a,b], p(a)>=0 by exact Bernstein basis.

    On an interval's interior ALL degree-n Bernstein basis functions
    are strictly positive. Thus all nonnegative coefficients, with at
    least one positive, imply strict positivity there. Requiring the
    last coefficient >0 also guarantees strict positivity at b.
    A zero at the LEFT endpoint is allowed and audited separately.
    """
    coefficients=bcoeff(poly,a,b)
    return (all(x>=0 for x in coefficients)
            and coefficients[-1]>0)


def _edge_box(P,i,a,b):
    """Rational rectangle containing entire moving edge over [a,b]."""
    n=len(P)
    return tuple((
        min(peval(P[j][axis],t)
            for j in (i,(i+1)%n) for t in (a,b)),
        max(peval(P[j][axis],t)
            for j in (i,(i+1)%n) for t in (a,b))
    ) for axis in (0,1))


def _max_pair_distance_squared(box1,box2):
    """Exact UPPER bound on distance of ANY two points in these boxes."""
    dx=max(abs(box1[0][0]-box2[0][1]),
           abs(box1[0][1]-box2[0][0]))
    dy=max(abs(box1[1][0]-box2[1][1]),
           abs(box1[1][1]-box2[1][0]))
    return dx*dx+dy*dy


def cutoff_numerator(D,N,cutoff):
    """D² * (actual side²-cutoff), where D is nonzero on the chart."""
    return psub(padd(pmul(N[2],N[2]),pmul(N[3],N[3])),
                pscale(pmul(D,D),cutoff))


def _edge_parameter_inequalities(P,assignment,D,N,det_sign):
    """Eight exact polynomials for closed-segment feasibility.

    For each square corner k on edge e, signed T>=0 and signed
    (D*|e|²-T)>=0 encode parameter in [0,1]. det_sign is a rigorously
    certified fixed sign of D on the full time interval.
    """
    for k,i in enumerate(assignment):
        a=P[i];e=pvsub(P[(i+1)%len(P)],a)
        x,y=corner_numerators(N,k)
        t=pvdot(e,(psub(x,pmul(a[0],D)),
                   psub(y,pmul(a[1],D))))
        ed=pmul(pvdot(e,e),D)
        yield pscale(t,det_sign)
        yield pscale(psub(ed,t),det_sign)


def _audited_square_vertices(P,assignment,t):
    """Exact CCW square, convenient for independent geometric comparison."""
    D,N=chart_polynomials(P,assignment)
    z=root_at(D,N,t)
    return tuple(square_corner(z,k) for k in range(4))


def exact_endpoint_census(P,cutoff,a,b):
    """Independent v223 all-chart census includes repeated/singular charts.

    At the a endpoint the A square is represented by TWO different
    supporting edge assignments, which this existing geometric census
    deduplicates exactly by square vertices, not by labels.
    """
    samples=[]
    for t,assignment in ((a,BRANCH_A),(b,BRANCH_B)):
        polygon=tuple(ppoint_at(vertex,t) for vertex in P)
        result=census_polygon(polygon,cutoff_squared=cutoff)
        if result.count!=1:
            raise AssertionError("independent v223 endpoint high-square census differs")
        actual=tuple(sorted(result.squares[0].vertices))
        predicted=tuple(sorted(_audited_square_vertices(P,assignment,t)))
        if predicted!=actual:
            raise AssertionError("endpoint geometric square disagrees with orbit census")
        samples.append({
            "time":str(t),
            "geometric_high_square_orbits":result.count,
            "repeated_labelings_deduplicated":result.redundant_labelings,
            "positive_square_on_polygon_vertex":result.squares[0].has_polygon_vertex,
            "exact_square_side_squared":str(result.squares[0].side_squared),
            "checked_edge_charts":result.checked_edge_charts,
        })
    if not samples[0]["positive_square_on_polygon_vertex"]:
        raise AssertionError("left endpoint must exercise vertex chart duplication")
    return tuple(samples)


def exhaustive_cutoff_classification(P,cutoff,a,b):
    """Prove NO OTHER square with side²>cutoff exists for a<t<=b.

    This is complete for any input passing the exact Jordan/edge-length
    hypotheses and classifying all C(n,4) ordered distinct-edge charts.
    When a chart cannot be excluded or certified as one of the listed
    known orbit branches, FAIL CLOSED rather than dropping it.
    """
    if type(a) is not Q or type(b) is not Q or not 0<=a<b<=1:
        raise ValueError("exact rational time bounds required")
    if type(cutoff) is not Q or cutoff<=0:
        raise ValueError("exact positive squared-side cutoff required")
    n=len(P)
    if n<4:
        raise ValueError("need at least four polygon edges")
    checks=0
    # An inscribed square with side²>cutoff has every pair of its four
    # distinct vertices separated by distance>sqrt(cutoff). If two
    # vertices belonged to one edge, that edge would have length²>cutoff.
    # Prove ALL edge lengths² < cutoff for ALL real times.
    for i in range(n):
        e=pvsub(P[(i+1)%n],P[i])
        checks+=certify_sign(
            psub(pol(cutoff),pvdot(e,e)),a,b,strict=True)

    boxes=tuple(_edge_box(P,i,a,b) for i in range(n))
    classes=Counter()
    active=set()
    representative={}
    for ids in combinations(range(n),4):
        # Convexity forces the four CCW square corners to occur in
        # CCW boundary order, so increasing edge quadruples suffice.
        if any(_max_pair_distance_squared(boxes[i],boxes[j])<cutoff
               for i,j in combinations(ids,2)):
            classes["geometry"]+=1
            continue
        D,N=chart_polynomials(P,ids)
        det_sign=_strict_sign_by_bernstein(D,a,b)
        params=tuple(_edge_parameter_inequalities(P,ids,D,N,det_sign))

        # A polynomial strictly negative on (a,b] proves infeasibility;
        # a possible equality at a is handled by exact v223 census.
        witness=next((k for k,q in enumerate(params)
                      if _positive_except_possible_left(pneg(q),a,b)),None)
        if witness is not None:
            classes["off_edge"]+=1
            representative[ids]=("off_edge",witness)
            continue

        G=cutoff_numerator(D,N,cutoff)
        if _positive_except_possible_left(pneg(G),a,b):
            classes["below_cutoff"]+=1
            representative[ids]=("below_cutoff",0)
            continue

        if ids not in EXPECTED_ORBITS:
            raise AssertionError(
                "unaccounted large-square chart; cannot claim exhaustive census: "+
                repr(ids))
        # Independently verify these named branches really lie on all
        # four segments and have strictly positive side throughout W,
        # including a (they may cross the larger cutoff internally).
        for q in params:
            checks+=certify_sign(q,a,b,strict=False)
        active.add(ids)
        classes["active"]+=1
    if sum(classes.values())!=comb(n,4):
        raise AssertionError("not every ordered distinct-edge chart classified")
    if set(active)!=set(EXPECTED_ORBITS):
        raise AssertionError("unexpected surviving geometric square branch")
    return {"classes":dict(classes),
            "total_four_distinct_edge_assignments":comb(n,4),
            "active_charts":tuple(sorted(active)),
            "exact_auxiliary_sign_checks":checks,
            "classification_witnesses":representative}


def poly_derivative(p):
    return pol(*(Q(i)*p[i] for i in range(1,len(p)))) if len(p)>1 else pol(0)


def unique_cutoff_event(P,assignment,cutoff,a,b,bracket,direction):
    """One algebraic size-crossing, isolated in rational open bracket.

    Strict sign of derivative on entire W proves at most one event;
    strict signs on bracket endpoints prove at least one by IVT.
    No floating point or assumed rational value of the root.
    """
    D,N=chart_polynomials(P,assignment)
    _strict_sign_by_bernstein(D,a,b)
    G=cutoff_numerator(D,N,cutoff)
    derivative=pscale(poly_derivative(G),direction)
    monotonicity=certify_sign(derivative,a,b,strict=True)
    lo,hi=bracket
    if not a<lo<hi<b:
        raise ValueError("invalid event isolation bracket")
    if not (direction*peval(G,lo)<0 and direction*peval(G,hi)>0):
        raise AssertionError("wrong cutoff entry/exit sign bracket")
    if not (direction*peval(G,a)<0 and direction*peval(G,b)>0):
        raise AssertionError("no endpoint sign change over whole window")
    return {
        "assignment":tuple(assignment),
        "direction":"enters above cutoff" if direction==1 else "exits above cutoff",
        "bracket":(lo,hi),
        "derivative_sign_checks":monotonicity,
        "cutoff_polynomial":G,
        "unique_algebraic_crossing":True,
    }


def audit():
    P=family_vertices()
    geometry=prove_jordan_family(P)
    if not geometry["strictly_convex_and_jordan_for_all_times"]:
        raise AssertionError("need certified convex Jordan family")
    a,b=WINDOW_START,WINDOW_END
    complete=exhaustive_cutoff_classification(P,HIGH_S2,a,b)
    if complete["classes"]!=EXPECTED_CLASSIFICATION:
        raise AssertionError("unexpected complete chart classifications")
    ends=exact_endpoint_census(P,HIGH_S2,a,b)

    enters=unique_cutoff_event(
        P,BRANCH_B,HIGH_S2,a,b,B_ENTER_BRACKET,direction=1)
    exits=unique_cutoff_event(
        P,BRANCH_A,HIGH_S2,a,b,A_EXIT_BRACKET,direction=-1)
    if enters["bracket"][1]>=exits["bracket"][0]:
        raise AssertionError("size-filtered events not strictly separated")

    # Both possible high-size orbits are distinct throughout W.
    da,na=chart_polynomials(P,BRANCH_A)
    db,nb=chart_polynomials(P,BRANCH_B)
    first_corner_x=psub(pmul(na[0],db),pmul(nb[0],da))
    distinct_sign=_strict_sign_by_bernstein(first_corner_x,a,b)
    if distinct_sign==0:
        raise AssertionError("two roots may be the same geometric orbit")
    # The separately certified v251 competitor-size crossover occurs
    # between these two cutoff events. It changes the identity of
    # which of A/B is larger, NOT the number above the cutoff.
    crossing=exact_switch_proof(P)
    xlo,xhi=crossing["unique_crossover_bracket"]
    if not (enters["bracket"][1]<xlo<xhi<exits["bracket"][0]):
        raise AssertionError("size comparison crossover should be between cutoff events")
    # Exactly one orbit above cutoff outside both roots and exactly two
    # strictly between them. At either root equality does NOT qualify.
    return {
        "status":"PASS: COMPLETE high-size square orbit ledger on convex 9-gon homotopy",
        "interval":[str(a),str(b)],
        "strict_square_side_squared_cutoff":str(HIGH_S2),
        "all_distinct_edge_charts_checked":complete["total_four_distinct_edge_assignments"],
        "classifications":complete["classes"],
        "exact_auxiliary_sign_checks":complete["exact_auxiliary_sign_checks"],
        "only_two_possible_orbit_branches":[list(x) for x in complete["active_charts"]],
        "both_geometric_orbits_distinct_for_all_times":True,
        "independent_complete_v223_endpoint_censuses":ends,
        "event_B_enters_above_cutoff":{
            "exact_time":"unique real algebraic root, not assumed rational",
            "lower":str(enters["bracket"][0]),
            "upper":str(enters["bracket"][1]),
            "certified_derivative_sign_checks":enters["derivative_sign_checks"],
        },
        "event_A_exits_above_cutoff":{
            "exact_time":"unique real algebraic root, not assumed rational",
            "lower":str(exits["bracket"][0]),
            "upper":str(exits["bracket"][1]),
            "certified_derivative_sign_checks":exits["derivative_sign_checks"],
        },
        "geometric_orbit_counts_strictly_above_cutoff_in_time_order":[1,2,1],
        "filtered_orbit_parities_in_time_order":[1,0,1],
        "filter_parity_is_not_invariant":True,
        "separate_same-size_crossing_bracket":[str(xlo),str(xhi)],
        "same_size_crossing_does_not_change_filtered_count":True,
        "scope":"EXHAUSTIVE only for size²>94257/10000 on this affine 9-gon interval; no assertion on smaller square orbits, all-square homotopy parity, generic singular births, universal Big Peg c, or Toeplitz.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
