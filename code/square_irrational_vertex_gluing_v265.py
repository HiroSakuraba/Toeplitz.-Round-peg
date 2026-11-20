"""v265: reusable exact irrational vertex-gluing / pair-annihilation certificate.

Scope: TWO regular, graceful square charts differing in ONE corner's
support on adjacent polygon edges. Certify that their exact Cramer roots
meet at a UNIQUE algebraic vertex time, that both actual positive squares
exist before but not after, that their Jacobian indices cancel, and that
the continuous 4D kink-domain has degree ZERO through the event.

Use exact Fraction polynomial arithmetic, exact Sturm isolation,
Bernstein signed-face and determinant margins. A broader v264 polygon
atlas is necessary before promoting the local event to a complete
whole-polygon singular-time count. No arbitrary-Jordan theorem.
"""
from fractions import Fraction as Q
from dataclasses import dataclass
from itertools import combinations
import json

from polygon_orbit_census_v223 import pt
from moving_polygon_events_v224 import (
    chart_polynomials,poly_gcd,divide,iszero,
    isolate_roots,interpolated_polygon,evaluate,
)
from big_peg_multichart_v249 import (
    pol,psub,padd,pmul,pscale,pneg,peval,bcoeff,certify_sign,
)
from pentagon_complete_ledger_v263 import _affine_polynomial_system
from convex_square_event_atlas_v264 import (
    _valid_endpoints,certify_strict_convex_motion,
    audit_moving_convex_polygon,InconclusiveAtlas,
)
from near_square_transversality_v258 import det4
from square_degree_cube_v259 import square_defect


class UnresolvedVertexEvent(InconclusiveAtlas):
    """Unsupported event form; no event count or local index assigned."""


def _strict_sign(p):
    if iszero(p):
        return 0
    for sign,v in ((1,p),(-1,pneg(p))):
        try:
            certify_sign(v,Q(0),Q(1),strict=True,max_depth=14)
            return sign
        except ValueError:
            pass
    if not isolate_roots(p,Q(0),Q(1),width=Q(1,2**20)):
        z=peval(p,Q(0))
        if z==0:
            raise AssertionError("root-free polynomial vanishes at t=0")
        return 1 if z>0 else -1
    return None


def _positive_min(p,max_depth=14):
    """Provable strictly positive global lower bound via Bernstein cones."""
    if iszero(p) or _strict_sign(p)!=1:
        raise UnresolvedVertexEvent("positive polynomial floor unavailable")
    todo=[(Q(0),Q(1),0)]
    lo=None
    while todo:
        a,b,depth=todo.pop()
        cs=bcoeff(p,a,b)
        if all(x>0 for x in cs):
            candidate=min(cs)
            lo=candidate if lo is None else min(lo,candidate)
        elif depth<max_depth:
            mid=(a+b)/2
            todo.extend(((a,mid,depth+1),(mid,b,depth+1)))
        else:
            raise UnresolvedVertexEvent("positive Bernstein floor budget exhausted")
    if lo is None or lo<=0:
        raise AssertionError("invalid strict positivity floor")
    return lo


def _positive_max(p):
    """Uniform positive polynomial upper bound from full interval Bernstein."""
    if _strict_sign(p)!=1:
        raise UnresolvedVertexEvent("not uniformly positive")
    return max(bcoeff(p,Q(0),Q(1)))


def _exact_gcd(p,q):
    g=poly_gcd(p,q)
    if len(g)<=1:
        raise UnresolvedVertexEvent("contact factors have no algebraic common root")
    return g


def _exact_quotient(p,q):
    f,r=divide(p,q)
    if not iszero(r):
        raise UnresolvedVertexEvent("expected algebraic root divisibility failed")
    return f


def _signed_contact(chart,j,kind,det_sign):
    if kind=="start":
        p=chart.numerators[j]
    elif kind=="end":
        p=chart.upper_differences[j]
    else:
        raise ValueError("contact kind must be start/end")
    return pscale(p,det_sign)


def _parameter_bound_polynomials(chart,j,lo,hi,sign):
    D,N=chart.determinant,chart.numerators[j]
    return (
        pscale(psub(N,pscale(D,lo)),sign),
        pscale(psub(pscale(D,hi),N),sign),
    )


def _matrix_entry_bound(A):
    return max(abs(peval(z,t)) for row in A for z in row
               for t in (Q(0),Q(1)))


@dataclass(frozen=True)
class VertexPairCertificate:
    left_edges:tuple
    right_edges:tuple
    corner_index:int
    shared_vertex_index:int
    exact_event_polynomial:tuple
    event_bracket:tuple
    left_orientation:int
    right_orientation:int
    external_defect_floor:Q
    source_root_face_gap:Q
    positive_square_size_cutoff:Q
    exact_divisibilities_checked:int
    exact_strict_face_signs_checked:int


def certify_vertex_pair(P0,P1,left,right,corner,cutoff=Q(0)):
    """Prove a genuine pair *death* at a unique algebraic polygon vertex.

    Input endpoint P0/P1 images are CCW STRICTLY convex at every t.
    Both 4-edge charts are increasing, identical except their corner-th
    edge: left supports edge j-1; right supports edge j. The local
    second-parameter chart is glued at left u=1, right u=0, using
    h=u_left-1 on the left and h=u_right on the right.
    The other three parameters use the same edge coordinates.

    The kink source domain is:
       other corners: (0,1), hinge h: (-1,1).
    Extended affine roots of BOTH sides are strictly inside its outer
    box for ALL t. Therefore its Brouwer degree is constant, and since
    it is zero before the event it stays ZERO at the irrational event.
    """
    A,B=_valid_endpoints(P0,P1)
    certify_strict_convex_motion(A,B)
    n=len(A)
    L,R=tuple(left),tuple(right)
    if (len(L)!=4 or len(R)!=4
            or type(corner) is not int or not 0<=corner<4
            or any(type(i) is not int or not 0<=i<n for i in L+R)
            or tuple(sorted(L))!=L or tuple(sorted(R))!=R
            or any(L[i]!=R[i] for i in range(4) if i!=corner)):
        raise ValueError("two aligned ordered edge charts differing in one corner required")
    j=R[corner]
    if L[corner]!=(j-1)%n or j==0:
        raise ValueError("must cross one non-wrap adjacent polygon vertex")
    if len(set(L))!=4 or len(set(R))!=4:
        raise UnresolvedVertexEvent("repeated-edge interaction requires extra order analysis")
    if type(cutoff) is not Q or cutoff<0:
        raise ValueError("nonnegative exact squared-side cutoff required")
    cL=chart_polynomials(A,B,L,orientation=1,cutoff_squared=cutoff)
    cR=chart_polynomials(A,B,R,orientation=1,cutoff_squared=cutoff)
    sL,sR=_strict_sign(cL.determinant),_strict_sign(cR.determinant)
    if sL not in (-1,1) or sR!=-sL:
        raise UnresolvedVertexEvent("nonsingular opposite-index charts not certified")

    # Shared algebraic contact time: u_L[corner]=1 and u_R[corner]=0.
    HL=_signed_contact(cL,corner,"end",sL)
    HR=_signed_contact(cR,corner,"start",sR)
    q=_exact_gcd(HL,HR)
    roots=isolate_roots(q,Q(0),Q(1),width=Q(1,2**24))
    if len(roots)!=1 or roots[0].rational:
        raise UnresolvedVertexEvent("not one genuine irrational contact event")
    br=roots[0]
    if not (Q(0)<br.lower<br.upper<Q(1)):
        raise UnresolvedVertexEvent("event touches time boundary")
    # Both charts are physically feasible before and infeasible after.
    # Sign of the critical polynomial is controlled by event q and
    # a strictly signed quotient over the ENTIRE real time interval.
    quotL=_exact_quotient(HL,q)
    quotR=_exact_quotient(HR,q)
    if _strict_sign(quotL) not in (-1,1) or _strict_sign(quotR) not in (-1,1):
        raise UnresolvedVertexEvent("unresolved event contact multiplicity")
    if not (peval(HL,Q(0))>0 and peval(HR,Q(0))>0
            and peval(HL,Q(1))<0 and peval(HR,Q(1))<0):
        raise UnresolvedVertexEvent("event not a simultaneous two-square death")
    if not (peval(q,Q(0))*peval(q,Q(1))<0):
        raise UnresolvedVertexEvent("contact is tangential, not a simple sign-changing event")
    for p in (cL.cutoff,cR.cutoff):
        if _strict_sign(p)!=1:
            raise UnresolvedVertexEvent("extended square side cutoff not positive")
    # Exact equality of all 3 *other* unwrapped corner parameters at q.
    checks=0
    for i in range(4):
        if i==corner:
            continue
        p=psub(pmul(cL.numerators[i],cR.determinant),
               pmul(cR.numerators[i],cL.determinant))
        _exact_quotient(p,q)
        checks+=1
    # The hinge corner is exactly the same polygon vertex since the
    # two incidence coordinates are 1 (left) and 0 (right) at q.
    _exact_quotient(cL.upper_differences[corner],q)
    _exact_quotient(cR.numerators[corner],q)
    checks+=2

    # Compute time-uniform enclosure of both extended roots in the
    # RECTANGULAR outer kink box, including infeasible extensions.
    # h=uL-1 on left, h=uR on right; both h in (-1,1).
    mins=[];strict_face_checks=0
    for chart,sign,is_left in ((cL,sL,True),(cR,sR,False)):
        d_abs=pscale(chart.determinant,sign)
        det_min=_positive_min(d_abs)
        det_max=_positive_max(d_abs)
        A_poly,_=_affine_polynomial_system(A,B,chart.edges)
        matrix_bound=_matrix_entry_bound(A_poly)
        if matrix_bound<=0:
            raise AssertionError("zero affine Jacobian magnitude")
        # Conservative adjugate inverse infinity norm:
        # each 3x3 cofactor <= 6 B^3 and 4 terms per row.
        inv_bound=Q(24)*matrix_bound**3/det_min
        root_gaps=[]
        for i in range(4):
            lower,upper=((-Q(0),Q(2)) if (i==corner and is_left)
                         else ((-Q(1),Q(1)) if i==corner
                               else (Q(0),Q(1))))
            for p in _parameter_bound_polynomials(
                    chart,i,lower,upper,sign):
                root_gaps.append(_positive_min(p)/det_max)
                strict_face_checks+=1
        gap=min(root_gaps)
        if gap<=0:
            raise AssertionError("lost strict root-to-exterior gap")
        mins.append((gap,gap/inv_bound))
    source_gap=min(x[0] for x in mins)
    residual_floor=min(x[1] for x in mins)
    if not source_gap>0 or not residual_floor>0:
        raise AssertionError("missing external boundary proof")
    return VertexPairCertificate(
        L,R,corner,j,tuple(q),(br.lower,br.upper),
        sL,sR,residual_floor,source_gap,cutoff,checks,strict_face_checks)


def complete_v264_vertex_event(P0,P1,left,right,corner,width=Q(1,2**18)):
    """Combine a COMPLETE v264 chart cover with a certified unique pair death.

    If the v264 event cover has any other candidate roots, this function
    REFUSES to infer an event count. It is deliberately not an arbitrary
    multiple-event resolver.
    """
    A,B=_valid_endpoints(P0,P1)
    atlas=audit_moving_convex_polygon(A,B,width=width)
    if atlas["all_time_exact_collinear_refinement_vertices_removed"]:
        raise UnresolvedVertexEvent(
            "pass compressed polygon endpoints to the pair verifier")
    floor=atlas["uniform_positive_square_size_certificate"]
    cutoff=Q(floor["safe_side_squared_cutoff"])
    cert=certify_vertex_pair(A,B,left,right,corner,cutoff)
    events=atlas["possible_event_clusters"]
    if len(events)!=1:
        raise UnresolvedVertexEvent("whole polygon does not have one event cluster")
    e=events[0]
    factors=set(tuple(z) for z in e["candidate_factors"])
    expected={( ",".join(map(str,left)),f"corner-{corner}-end"),
              ( ",".join(map(str,right)),f"corner-{corner}-start")}
    if factors!=expected:
        raise UnresolvedVertexEvent(
            "other algebraic chart events may coexist in this bracket")
    lo,hi=cert.event_bracket
    if not (Q(e["lower"])<hi and lo<Q(e["upper"])):
        raise UnresolvedVertexEvent("pair contact root outside global event cluster")
    chamber=atlas["rational_event_free_chambers"]
    if len(chamber)!=2:
        raise UnresolvedVertexEvent("not exactly two event-free chambers")
    before=chamber[0]["geometric_square_orbits"]
    after=chamber[1]["geometric_square_orbits"]
    if before<2 or after!=before-2:
        raise UnresolvedVertexEvent("global census disagrees with isolated 2-to-0 pair death")
    if (not chamber[0]["all_sample_squares_regular_in_open_positive_sector"]
        or not chamber[1]["all_sample_squares_regular_in_open_positive_sector"]):
        raise UnresolvedVertexEvent("unresolved regular index at a rational chamber")
    # The shared pair is ONE positive square at the singular instant;
    # every remaining root lies on a stable chart, with no other walls.
    event_count=before-1
    out={
        "status":"PASS: complete three-slice whole-polygon positive-square census from exact generic irrational vertex gluing",
        "event_defining_polynomial":[str(x) for x in cert.exact_event_polynomial],
        "event_isolating_bracket_normalized":[str(x) for x in cert.event_bracket],
        "genuine_vertex_index":cert.shared_vertex_index,
        "pair_signed_regular_indices":[cert.left_orientation,cert.right_orientation],
        "glued_pair_relative_degree":0,
        "uniform_outer_boundary_square_defect_floor":str(cert.external_defect_floor),
        "uniform_extended_root_distance_from_outer_faces":str(cert.source_root_face_gap),
        "event_before_geometric_squares":before,
        "event_exact_geometric_squares":event_count,
        "event_after_geometric_squares":after,
        "all_real_time_pre_event_count":before,
        "all_real_time_post_event_count":after,
        "relative_unlabeled_parity_before":before%2,
        "relative_unlabeled_parity_after":after%2,
        "complete_whole_polygon_all_time_ledger":True,
        "local_exterior_C0_defect_margin_positive":True,
        "all_real_time_uniform_positive_square_side_floor":
            floor["every_positive_square_side_squared_lower_bound"],
        "assumptions":"Exact affine convex polygon with one complete event cluster; two distinct regular square charts differing at one genuine vertex; two actual roots before, one at event, none after; exact Cramer identities and zero-free outer kink-domain. Nonconvex, multi-event, determinant-singular and tangency cases remain inconclusive.",
        "unrestricted_toeplitz_proved":False,
    }
    return out


def demo():
    from big_peg_pair_death_v254 import (
        family_vertices,ROOT_LEFT,ROOT_RIGHT,
        COLLIDING_LEFT,COLLIDING_RIGHT,
    )
    from big_peg_multichart_v249 import ppoint_at
    P=family_vertices()
    a=tuple(ppoint_at(v,ROOT_LEFT) for v in P)
    b=tuple(ppoint_at(v,ROOT_RIGHT) for v in P)
    return complete_v264_vertex_event(a,b,COLLIDING_LEFT,COLLIDING_RIGHT,1)


if __name__=="__main__":
    print(json.dumps(demo(),sort_keys=True,indent=2))
