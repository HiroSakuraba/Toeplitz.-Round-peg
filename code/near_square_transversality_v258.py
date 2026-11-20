"""v258: conditional quantitative lift of a near-square to an EXACT square.

Given an ideal Euclidean square z0 and four ACTUAL polygon-boundary points
q_k close to its four corners (e.g. v257), assume:
  - q_k strictly interior to four distinct selected polygon edges;
  - the corresponding 4x4 incidence matrix M is invertible;
  - a rational bound delta on all |q_k-C_k(z0)| is small relative to the
    M-inverse condition number and four edge-contact margins.

Then the UNIQUE supporting-line square z=M^-1 b lies strictly on all
four closed edges, and has a rigorous positive squared-side lower bound.
This is a finite-dimensional implicit correction, NOT a theorem that
approximate squares alone imply exact ones. No convexity is required.

A separate v248 Brouwer certificate proves the reference square also
survives ALL sufficiently small continuous rough normal-graph changes of
these four selected boundary arcs, without smoothness assumptions.

Exactly calculated 10-gon example: centrally symmetric nonconvex reference,
perturb 9 of its 10 coordinates to a genuinely asymmetric polygon, obtain
an approximate square via v257's quarter-turn and antipodal pairing, then
lift it to an EXACT square, with rational quantitative certificates.

All math decisions are exact Fraction operations. Toeplitz still open.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import permutations
import json

from polygon_orbit_census_v223 import (
    pt,add,sub,mul,dot,cross,_on_seg,validate_jordan_polygon,
)
from big_peg_persistence_v248 import (
    CORNER_COEFFS, invert_and_solve, square_corner, certify_persistence,
)
from central_symmetry_square_v256 import (
    J,minus,disk_certificate,certify_near_square,
    is_centrally_symmetric,
)

BASE_TEN=tuple(pt(*p) for p in (
    (4,-1),(4,1),(Q(3,2),Q(3,2)),(2,5),(-2,3),
    (-4,1),(-4,-1),(-Q(3,2),-Q(3,2)),(-2,-5),(2,-3),
))
BASE_ASSIGNMENT=(0,3,5,8)
R2=Q(1)
DELTA=Q(1,30)
GRAPH_EPSILON=Q(1,100)


def asymmetric_ten():
    """Rational normal shifts of four supporting lines plus a notch change.

    All points are explicit rational translations of BASE_TEN.
    Half of the perturbations used during a preliminary search are
    retained to make the v257 global antipodal defect less than 1/30.
    """
    P=list(BASE_TEN)
    shifts={
        0:(Q(1,160),Q(0)),
        1:(Q(1,160),Q(0)),
        2:(Q(1,200),-Q(1,240)),
        3:(Q(1,140),-Q(1,70)),
        4:(Q(1,140),-Q(1,70)),
        5:(-Q(1,180),Q(0)),
        6:(-Q(1,180),Q(0)),
        8:(Q(1,190),-Q(1,95)),
        9:(Q(1,190),-Q(1,95)),
    }
    for i,(x,y) in shifts.items():
        P[i]=add(P[i],pt(x,y))
    return tuple(P)


def det4(M):
    """Exact determinant with an independent 24-permutation expansion."""
    if len(M)!=4 or any(len(row)!=4 for row in M):
        raise ValueError("expected four-by-four matrix")
    d=Q(0)
    for order in permutations(range(4)):
        inversions=sum(order[i]>order[j] for i in range(4) for j in range(i+1,4))
        term=Q(-1 if inversions%2 else 1)
        for i,j in enumerate(order):
            term*=M[i][j]
        d+=term
    return d


def _point_parameter(p,a,b):
    e=sub(b,a)
    ee=dot(e,e)
    if ee<=0:
        raise ValueError("zero-length boundary edge")
    return dot(e,sub(p,a))/ee


def _square_parameter(corners):
    """Invert four-corner square linear parametrization exactly."""
    if len(corners)!=4:
        raise ValueError("four ideal square corners required")
    p=corners[0]
    s=sub(corners[1],p)
    z=(p.x,p.y,s.x,s.y)
    if dot(s,s)<=0 or any(square_corner(z,k)!=corners[k] for k in range(4)):
        raise ValueError("ideal points are not a positive-size oriented square")
    return z


def _candidate_assignments(P,actual):
    """Unique strictly interior edge containing each ordered boundary point."""
    ids=[]
    for q in actual:
        match=[]
        for i in range(len(P)):
            a,b=P[i],P[(i+1)%len(P)]
            t=_point_parameter(q,a,b)
            if 0<t<1 and _on_seg(q,a,b):
                match.append(i)
        if len(match)!=1:
            raise ValueError("near-square contact not interior to a unique edge")
        ids.append(match[0])
    return tuple(ids)


@dataclass(frozen=True)
class TransverseLiftCertificate:
    assignment:tuple
    determinant:Q
    ideal_square_parameters:tuple
    exact_square_parameters:tuple
    ideal_side_squared:Q
    exact_side_squared:Q
    guaranteed_side_squared_lower:Q
    corner_error_bound:Q
    maximum_observed_corner_error_squared:Q
    condition_response_bound:Q
    guaranteed_parameter_displacement:Q
    actual_parameter_displacement:Q
    approximate_edge_parameters:tuple
    exact_edge_parameters:tuple
    maximum_relative_edge_parameter_bounds:tuple
    verified_interior_contact_margins:tuple
    exact_square_corners:tuple


def certify_polygon_square_lift(P,ideal,actual,delta,
                                assignment=None,*,disk_r2=None):
    """Finite, exact transversality sufficient condition.

    The input q_k are ACTUAL boundary points on four distinct edge
    interiors; ideal[k] are corners of ONE positive ideal square.
    Each |q_k-ideal_k|<=delta. The incidence matrix M of their
    selected supporting lines must be invertible.

    Let L=max_j sum_i |M^-1[j,i]|*(|E_ix|+|E_iy|), R=delta*L.
    Cramer gives |z_exact-z_ideal|_infty<=R, because
    |n_i.(ideal_i-q_i)|<=delta*(|E_ix|+|E_iy|).

    For each edge contact with parameter t_i:
      |t_exact_i-t_i|
         <= delta*||E_i||_1/||E_i||_2²
            + R*sum_j |E_i.C_k_basis_j|/||E_i||_2².
    If the RHS is < min(t_i,1-t_i), every corrected corner remains
    strictly inside its assigned edge. The solution is an EXACT
    boundary-inscribed square, regardless of convexity.

    Fails closed if geometry, incidence, contact margin or positive
    squared-side guarantee cannot be certified.
    """
    P=validate_jordan_polygon(tuple(P))
    ideal=tuple(ideal);actual=tuple(actual)
    z0=_square_parameter(ideal)
    if len(actual)!=4:
        raise ValueError("four ordered boundary contacts required")
    if type(delta) is not Q or delta<0:
        raise ValueError("nonnegative exact Fraction displacement required")
    if disk_r2 is not None:
        disk_certificate(P,disk_r2)
        if dot(ideal[0],ideal[0])<disk_r2 or (
            ideal[2]!=minus(ideal[0]) or ideal[3]!=minus(ideal[1])):
            raise ValueError("ideal centered square must clear certified disk")

    error2=tuple(dot(sub(a,b),sub(a,b)) for a,b in zip(actual,ideal))
    if max(error2)>delta*delta:
        raise ValueError("approximate contacts exceed claimed displacement")
    chosen=tuple(assignment) if assignment is not None else _candidate_assignments(P,actual)
    if len(chosen)!=4 or len(set(chosen))!=4 or any(
            type(i) is not int or not 0<=i<len(P) for i in chosen):
        raise ValueError("four distinct valid edge indices required")

    M=[];b=[];norms=[];params=[];sensitivities=[];param_base=[]
    for k,i in enumerate(chosen):
        a=P[i]
        e=sub(P[(i+1)%len(P)],a)
        n=J(e)
        cx,cy=CORNER_COEFFS[k]
        M.append(tuple(n.x*cx[j]+n.y*cy[j] for j in range(4)))
        b.append(dot(n,a))
        q=actual[k]
        t=_point_parameter(q,a,P[(i+1)%len(P)])
        if not (_on_seg(q,a,P[(i+1)%len(P)]) and Q(0)<t<Q(1)):
            raise ValueError("approximate point is not strictly interior to assigned edge")
        params.append(t)
        e2=dot(e,e)
        norms.append(abs(e.x)+abs(e.y))
        sensitivities.append(sum(abs(e.x*cx[j]+e.y*cy[j])
                                 for j in range(4))/e2)
        param_base.append(norms[-1]/e2)

    determinant=det4(M)
    if determinant==0:
        raise ValueError("nontransverse near-square: singular incidence matrix")
    inverse,z_exact=invert_and_solve(M,b)
    response=max(sum(abs(inverse[j][i])*norms[i] for i in range(4))
                 for j in range(4))
    R=delta*response
    actual_R=max(abs(z_exact[j]-z0[j]) for j in range(4))
    if actual_R>R:
        raise AssertionError("derived inverse residual bound failed")
    bounds=tuple(param_base[k]*delta+sensitivities[k]*R for k in range(4))
    margins=tuple(min(t,Q(1)-t) for t in params)
    if any(not bounds[k]<margins[k] for k in range(4)):
        raise ValueError("transversality conditioning cannot protect edge margins")
    exact_corners=tuple(square_corner(z_exact,k) for k in range(4))
    exact_params=[]
    for k,i in enumerate(chosen):
        a=P[i];end=P[(i+1)%len(P)]
        t=_point_parameter(exact_corners[k],a,end)
        if not (Q(0)<t<Q(1) and _on_seg(exact_corners[k],a,end)):
            raise AssertionError("corrected square corner is not strictly on its edge")
        if abs(t-params[k])>bounds[k]:
            raise AssertionError("corrected edge parameter exceeds certified bound")
        exact_params.append(t)

    side2=z_exact[2]**2+z_exact[3]**2
    lower=sum(max(abs(z0[j])-R,Q(0))**2 for j in (2,3))
    if lower<=0 or side2<lower:
        raise ValueError("no certified positive square-size lower bound")
    return TransverseLiftCertificate(
        chosen,determinant,z0,z_exact,z0[2]**2+z0[3]**2,
        side2,lower,delta,max(error2),response,R,actual_R,
        tuple(params),tuple(exact_params),bounds,margins,exact_corners)


def certify_v257_near_square_lift(P,r2,delta):
    """Apply transversality correction to v257's explicit near-square."""
    P=validate_jordan_polygon(tuple(P))
    near=certify_near_square(P,r2)
    if max(dot(sub(a,b),sub(a,b)) for a,b in zip(
            near.exact_boundary_points,near.ideal_centered_square))>delta*delta:
        raise ValueError("v257 witness exceeds proposed exact-lift radius")
    cert=certify_polygon_square_lift(
        P,near.ideal_centered_square,near.exact_boundary_points,
        delta,disk_r2=r2)
    if cert.guaranteed_side_squared_lower<=0:
        raise AssertionError("no positive exact square bound")
    return near,cert


def audit():
    base=validate_jordan_polygon(BASE_TEN)
    if is_centrally_symmetric(base) is False:
        raise AssertionError("reference lacks central symmetry")
    if not any(cross(sub(base[i],base[i-1]),
                     sub(base[(i+1)%len(base)],base[i]))<0
               for i in range(len(base))):
        raise AssertionError("reference must be NONCONVEX")
    disk_certificate(base,R2)
    rough=certify_persistence(base,BASE_ASSIGNMENT,GRAPH_EPSILON)
    if rough.minimum_square_side_squared<=0:
        raise AssertionError("rough normal graph class has no protected side")

    P=asymmetric_ten()
    disk_certificate(P,R2)
    if is_centrally_symmetric(P):
        raise AssertionError("target is not genuinely asymmetric")
    near,cert=certify_v257_near_square_lift(P,R2,DELTA)
    if cert.guaranteed_side_squared_lower<Q(28):
        raise AssertionError("quantitative exact-square size floor unexpectedly weak")
    if cert.determinant==0 or not all(Q(0)<t<Q(1) for t in cert.exact_edge_parameters):
        raise AssertionError("exact square transversality/corners failed")
    return {
        "status":"PASS: v257 near-square lifts to an EXACT positive inscribed square under quantitative edge-chart transversality",
        "reference_vertices":len(base),
        "reference_centrally_symmetric":True,
        "reference_nonconvex":True,
        "reference_normal_graph_rough_epsilon":str(GRAPH_EPSILON),
        "reference_v248_rough_minimum_side_squared":str(rough.minimum_square_side_squared),
        "target_vertices":len(P),
        "target_centrally_symmetric":False,
        "target_certified_disk_radius_squared":str(R2),
        "v257_actual_approximate_max_error_squared":str(near.max_corner_error_squared),
        "v258_proposed_delta":str(DELTA),
        "transverse_four_edge_assignment":list(cert.assignment),
        "four_edge_incidence_determinant":str(cert.determinant),
        "inverse_condition_response_L":str(cert.condition_response_bound),
        "maximum_guaranteed_square_parameter_correction":str(cert.guaranteed_parameter_displacement),
        "actual_exact_correction":str(cert.actual_parameter_displacement),
        "ideal_square_side_squared":str(cert.ideal_side_squared),
        "verified_exact_square_side_squared":str(cert.exact_side_squared),
        "certified_exact_square_side_squared_lower":str(cert.guaranteed_side_squared_lower),
        "exact_square_points":[[str(p.x),str(p.y)] for p in cert.exact_square_corners],
        "guarantee":"For this exact NONCONVEX asymmetric polygon, four v257 actual near-square boundary points have a quantitatively conditioned unique exact four-edge square continuation.",
        "limitation":"Conditional on four distinct interior contacts, invertible square-incidence matrix and conditioning margin. Small symmetry Hausdorff defect ALONE is insufficient for this certificate. No universal Big Peg bound or unrestricted Toeplitz proof.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
