"""v260: EXACT relative-degree gluing across v254's positive-square pair death.

Consider the four-arc square-defect F(s,t) of v259 on a RECTANGULAR
4D parameter domain, where its second boundary parameter s_1 crosses
the genuine convex pentagon vertex V_2(t). The map is continuous and
piecewise affine: negative s_1 uses edge 1, positive s_1 uses edge 2.

The v254 event has TWO genuine positive squares before its unique
algebraic collision, ONE shared vertex square AT the event, and NONE
in this local four-arc box afterward. The two pre-event square roots
have OPPOSITE Jacobian determinant signs (-1,+1), so their total
local Brouwer degree is ZERO.

The essential new v260 fact is uniform OUTER BOUNDARY NONVANISHING:
for EVERY REAL time in v254's exact isolating interval, both affine
Cramer solutions (including their formally infeasible extensions)
are proved to stay strictly inside the outer parameter box, with a
positive RATIONAL gap. This proves F has no zeros on its 3D exterior
faces, even at the nonsmooth irrational vertex collision. Exact
determinant and entry bounds give an explicit uniform positive lower
bound on ||F||_infty at the exterior boundary, robust to SMALL arbitrary
continuous coordinatewise arc displacements.

The standard degree additivity law permits summing the two side-cell
indices only while F is nonzero on their SHARED interface. At the
collision instant the interface contains a zero, so the separate
indices are undefined, but their common outer-union degree stays zero
by homotopy invariance.

This proves an actual local gluing/cancellation mechanism, NOT a
nonzero protected global square index. Degree ZERO cannot force square
existence; the example shows exactly why. No all-size orbit ledger on
all real times, universal Big Peg, or Toeplitz proof is claimed.

Every certificate decision is an exact Fraction/Bernstein sign, or
an exact algebraic polynomial divisibility check.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
import json

from polygon_orbit_census_v223 import (
    pt, add, sub, mul, dot, validate_jordan_polygon, census_polygon,
)
from big_peg_pair_death_v254 import (
    family_vertices, certify_jordan_family, certify_collision,
    PERSISTENT, COLLIDING_LEFT, COLLIDING_RIGHT,
    QPOLY, ROOT_LEFT, ROOT_RIGHT, EVENT,
)
from big_peg_multichart_v249 import (
    pol, padd, psub, pmul, pneg, pscale, peval,
    ppoint_at, certify_sign,
)
from moving_polygon_events_v224 import (
    chart_polynomials as edge_chart_polynomials,
    divide, iszero,
)
from near_square_transversality_v258 import det4
from square_degree_cube_v259 import square_defect, column_for_arc

START=ROOT_LEFT
END=ROOT_RIGHT
BOX=((Q(9,10),Q(93,100)),
     (Q(-1,1000),Q(1,1000)),
     (Q(43,100),Q(46,100)),
     (Q(73,100),Q(78,100)))
INTERIOR_GAP=Q(1,10000)
ENTRY_BOUND=Q(2)
DETERMINANT_FLOOR=Q(1,10)
UNIFORM_INVERSE_NORM_BOUND=Q(24)*ENTRY_BOUND**3/DETERMINANT_FLOOR
BOUNDARY_RESIDUAL_FLOOR=INTERIOR_GAP/UNIFORM_INVERSE_NORM_BOUND
C0_EPSILON=Q(1,10**9)
# Every component of the square-defect map is a signed sum of no
# more than FOUR separate input x/y coordinates.
PERTURBATION_RESIDUAL_RESPONSE=Q(4)
# Exact all-time separation of the first two square corners:
SIDE_X_FLOOR=Q(74,100)
SIDE_Y_FLOOR=Q(79,100)
CHARTS=(COLLIDING_LEFT,COLLIDING_RIGHT)


def _snapshot(P,t):
    return tuple(ppoint_at(v,t) for v in P)


def _square_points(P,t,s):
    """The ACTUAL four continuous boundary arcs, with a genuine kink.

    s=(u0,h,u2,u3), h<=0 on edge 1 measured back from V2, h>=0
    on edge 2 measured forward from V2. Edge-interior parameters
    u0,u2,u3 belong to edges 0,3,4 respectively.
    """
    if len(s)!=4:
        raise ValueError("four boundary parameters required")
    p=_snapshot(P,t)
    p0=add(p[0],mul(sub(p[1],p[0]),s[0]))
    p1=(add(p[2],mul(sub(p[2],p[1]),s[1])) if s[1]<=0
        else add(p[2],mul(sub(p[3],p[2]),s[1])))
    p2=add(p[3],mul(sub(p[4],p[3]),s[2]))
    p3=add(p[4],mul(sub(p[0],p[4]),s[3]))
    return (p0,p1,p2,p3)


def _map(P,t,s):
    return square_defect(_square_points(P,t,s))


def _matrix(P,t,assignment):
    p=_snapshot(P,t)
    columns=tuple(column_for_arc(k,sub(p[(j+1)%len(p)],p[j]))
                  for k,j in enumerate(assignment))
    return tuple(tuple(columns[k][i] for k in range(4))
                 for i in range(4))


def _signed_sign(p,a=Q(0),b=Q(1)):
    try:
        certify_sign(p,a,b,strict=True)
        return 1
    except ValueError:
        certify_sign(pneg(p),a,b,strict=True)
        return -1


def _local_numerators(chart,assignment):
    """Cramer root in the SHARED 4D coordinates, including h=u1-1 on e1."""
    numerators=list(chart.numerators)
    if assignment==COLLIDING_LEFT:
        numerators[1]=psub(numerators[1],chart.determinant)
    elif assignment!=COLLIDING_RIGHT:
        raise ValueError("only the two certified adjacent square charts")
    return tuple(numerators)


def _chart_at(P,assignment):
    if tuple(P)!=family_vertices() or assignment not in CHARTS:
        raise ValueError("v260 is scoped to its proved affine pentagon")
    start,end=_snapshot(P,START),_snapshot(P,END)
    return edge_chart_polynomials(
        start,end,assignment,orientation=1,cutoff_squared=Q(1))


def _normalized_event_polynomial():
    """Q(START+(END-START)*tau); unique root in normalized [0,1]."""
    q0,q1,q2=QPOLY
    d=END-START
    return pol(q0+q1*START+q2*START**2,
               (q1+2*q2*START)*d,
               q2*d*d)


@dataclass(frozen=True)
class CellDegreeCertificate:
    assignment:tuple
    orientation_sign:int
    determinant_floor:Q
    entry_bound:Q
    gap_from_outer_faces:Q
    parameter_coordinate_sign_checks:int
    event_factor:tuple
    event_quotient_sign:int
    signed_roots_at_start:tuple
    signed_roots_at_end:tuple


def certify_cell(P,assignment,box=BOX,gap=INTERIOR_GAP):
    """Uniform exact root enclosure, determinant orientation and event germ.

    This proves the EXTENDED affine branch root lies at least gap INSIDE
    every outer coordinate bound for all real t in [START,END], even
    when the root's kink parameter is of the wrong sign and hence not
    an actual square of the glued curve.

    Any genuine square-defect root in a halfcell must be this root.
    Consequently there can be NO roots on the OUTER box boundary.
    """
    if tuple(P)!=family_vertices() or assignment not in CHARTS:
        raise ValueError("uncertified pentagon family or chart")
    if tuple(tuple(x) for x in box)!=BOX or gap!=INTERIOR_GAP:
        raise ValueError("only the audited v260 box and gap are certified")
    chart=_chart_at(P,assignment)
    D=chart.determinant
    N=_local_numerators(chart,assignment)
    sign=_signed_sign(D)
    checks=certify_sign(psub(pscale(D,sign),pol(DETERMINANT_FLOOR)),
                         Q(0),Q(1),strict=True)

    # Each matrix entry is AFFINE in normalized time: the pentagon has
    # one moving vertex, and each boundary edge depends affinely on t.
    # Endpoint absolute-value checks therefore bound EVERY real time.
    for t,normalized in ((START,Q(0)),(END,Q(1))):
        M=_matrix(P,t,assignment)
        if any(abs(x)>ENTRY_BOUND for row in M for x in row):
            raise AssertionError("uniform matrix entry bound fails")
        if det4(M)!=peval(D,normalized):
            raise AssertionError("independent determinant formula disagrees")

    # The rational Cramer root coordinate N_j/D must stay in the
    # strictly smaller rectangle (lo+gap,hi-gap). Multiplying by
    # sign(D) converts each claim to strict polynomial positivity.
    for j,(lo,hi) in enumerate(BOX):
        if not (type(lo) is Q and type(hi) is Q
                and hi-lo>2*INTERIOR_GAP):
            raise AssertionError("invalid rational outer source box")
        lower=pscale(psub(N[j],pscale(D,lo+gap)),sign)
        upper=pscale(psub(pscale(D,hi-gap),N[j]),sign)
        checks+=certify_sign(lower,Q(0),Q(1),strict=True)
        checks+=certify_sign(upper,Q(0),Q(1),strict=True)

    # Both critical kink coordinates factor by EXACTLY the same
    # algebraic event polynomial Q(t). The quotient has a fixed sign,
    # so their h(t) change sign ONCE, with no hidden interface root.
    normalized_q=_normalized_event_polynomial()
    factor,remainder=divide(N[1],normalized_q)
    if not iszero(remainder):
        raise AssertionError("kink root missing the exact event factor")
    factor_sign=_signed_sign(factor)
    if factor_sign<=0:
        raise AssertionError("unexpected Cramer kink root factor sign")
    if not (peval(normalized_q,0)>0 and peval(normalized_q,1)<0):
        raise AssertionError("wrong event orientation")

    # Use the same determinant's sign to derive h signs at both
    # RATIONAL endpoints. Left branch h<0 before and h>0 after;
    # right branch h>0 before and h<0 after.
    endpoints=tuple(tuple(peval(n,t)/peval(D,t) for n in N)
                    for t in (Q(0),Q(1)))
    if assignment==COLLIDING_LEFT:
        if not (endpoints[0][1]<0<endpoints[1][1] and sign==-1):
            raise AssertionError("negative branch did not die at vertex")
    else:
        if not (endpoints[0][1]>0>endpoints[1][1] and sign==1):
            raise AssertionError("positive branch did not die at vertex")
    return CellDegreeCertificate(
        assignment,sign,DETERMINANT_FLOOR,ENTRY_BOUND,
        INTERIOR_GAP,checks,factor,factor_sign,*endpoints)


def certify_interface_gluing(P):
    """Prove all FOUR parameter coordinates coincide at the kink root."""
    if tuple(P)!=family_vertices():
        raise ValueError("unproved geometry")
    a=_chart_at(P,COLLIDING_LEFT)
    b=_chart_at(P,COLLIDING_RIGHT)
    A=_local_numerators(a,COLLIDING_LEFT)
    B=_local_numerators(b,COLLIDING_RIGHT)
    q=_normalized_event_polynomial()
    quotients=[]
    for k in range(4):
        numerator=psub(pmul(A[k],b.determinant),
                       pmul(B[k],a.determinant))
        quotient,remainder=divide(numerator,q)
        if not iszero(remainder):
            raise AssertionError("two glued roots differ at algebraic event")
        quotients.append(quotient)
    if not any(quotients):
        raise AssertionError("branches are identically the same root")
    if not (peval(A[1],0)/peval(a.determinant,0)<0
            and peval(B[1],0)/peval(b.determinant,0)>0):
        raise AssertionError("two actual roots must exist before event")
    # The collision is a unique inner-interface zero. The v254
    # certificate additionally establishes positivity and agreement of
    # all square's GEOMETRIC coordinates, not just source parameters.
    collision=certify_collision(P)
    if not collision["root_identity_divisible_by_event_in_four_components"]:
        raise AssertionError("geometric square roots do not glue")
    return tuple(quotients)


def certify_sides_and_robust_boundary(P,epsilon=C0_EPSILON):
    """Uniform external-boundary nonvanishing and square-side noncollapse.

    On either closed halfbox, F=M(t)(s-r(t)). The certified root r(t)
    is at least GAP away from every OUTER face. Every 3x3 cofactor of
    a 4x4 matrix with entries <=ENTRY_BOUND has magnitude <=6*B^3.
    With |det M|>DET_FLOOR, each inverse row sum is <24*B^3/d.
    Therefore every exterior point satisfies

      ||F||inf >= GAP / inverse_bound >0.

    Arbitrary continuous coordinatewise displacements |h_k.x|,|h_k.y|
    <=epsilon change any square-defect component by at most 4epsilon.
    If 4epsilon < residual_floor, the degree of F on the entire glued
    outer domain stays unchanged (ZERO). This does not guarantee a root.

    Finally, a multi-AFFINE endpoint rectangle check bounds all possible
    first-to-second side vectors away from zero, including at the kink.
    """
    if tuple(P)!=family_vertices():
        raise ValueError("unproved geometry")
    if type(epsilon) is not Q or epsilon<0:
        raise ValueError("nonnegative exact Fraction epsilon required")
    if PERTURBATION_RESIDUAL_RESPONSE*epsilon>=BOUNDARY_RESIDUAL_FLOOR:
        raise ValueError("continuous perturbation may create an outer boundary zero")
    checked=0
    # For each s1 branch, both first square-corner locations are
    # multi-affine in time t and their own edge parameters. Therefore
    # their coordinate differences reach extrema at these vertices.
    for t,u,h in product((START,END),BOX[0],
                         (BOX[1][0],Q(0),BOX[1][1])):
        p0,p1,*_= _square_points(P,t,(u,h,Q(43,100),Q(73,100)))
        side=sub(p1,p0)
        if not (side.x < -SIDE_X_FLOOR and side.y < -SIDE_Y_FLOOR):
            raise AssertionError("nonzero square-side floor fails")
        checked+=1
    side_floor=(SIDE_X_FLOOR-2*epsilon)**2+(SIDE_Y_FLOOR-2*epsilon)**2
    if side_floor<=Q(1):
        raise AssertionError("noncollapsed side bound too weak")
    return {
        "checked_multiaffine_side_vertices":checked,
        "rigorous_unperturbed_boundary_defect_floor":BOUNDARY_RESIDUAL_FLOOR,
        "max_continuous_perturbed_defect":PERTURBATION_RESIDUAL_RESPONSE*epsilon,
        "safe_continuous_coordinate_epsilon":epsilon,
        "certified_side_squared_floor":side_floor,
    }


def _verify_independent_snapshots(P):
    """Independent v223 exhaustive unfiltered geometric orbit counts."""
    samples=[]
    for t,expected in ((START,3),(END,1)):
        points=_snapshot(P,t)
        validate_jordan_polygon(points)
        census=census_polygon(points)
        if census.count!=expected or any(s.side_squared<=Q(1) for s in census.squares):
            raise AssertionError("independent complete geometric orbit count differs")
        samples.append((t,census.count))
    return samples


@dataclass(frozen=True)
class RelativeDegreeCertificate:
    interval:tuple
    box:tuple
    left_orientation:int
    right_orientation:int
    signed_degree_before_event:int
    signed_degree_after_event:int
    glued_degree_through_event:int
    strict_boundary_nonvanishing_floor:Q
    uniform_side_squared_lower:Q
    exact_algebraic_event:tuple
    interface_root_is_inner_at_event:bool
    all_time_root_box_sign_checks:int
    robustness_epsilon:Q
    independent_global_snapshot_counts:tuple


def audit():
    P=family_vertices()
    jordan=certify_jordan_family(P)
    if not jordan["convex_Jordan_for_all_real_times"]:
        raise AssertionError("not a Jordan isotopy")
    collision=certify_collision(P)
    if collision["determinants_nonzero_at_event"] is not True:
        raise AssertionError("not a true nonsingular-chart vertex gluing")
    before=certify_cell(P,COLLIDING_LEFT)
    after=certify_cell(P,COLLIDING_RIGHT)
    if (before.orientation_sign,after.orientation_sign)!=(-1,+1):
        raise AssertionError("pair death must have opposite signed indices")
    equality=certify_interface_gluing(P)
    if len(equality)!=4:
        raise AssertionError("incomplete four-coordinate gluing proof")
    boundary=certify_sides_and_robust_boundary(P)
    snapshots=_verify_independent_snapshots(P)
    if [n for t,n in snapshots]!=[3,1]:
        raise AssertionError("expected independent 3-to-1 event")

    # For t < alpha: one regular zero in each of the two halfboxes,
    # with derivative signs -1,+1 => relative degree zero.
    # For t > alpha: neither halfbox root has the correct sign => no
    # roots in this local box => relative degree zero.
    # For t = alpha: the two solutions coincide at the INNER seam.
    # Their individual halfbox degrees are UNDEFINED because their
    # boundaries contain the root, while the root-free outer boundary
    # preserves degree zero by homotopy invariance.
    certificate=RelativeDegreeCertificate(
        (START,END),BOX,-1,+1,0,0,0,
        boundary["rigorous_unperturbed_boundary_defect_floor"],
        boundary["certified_side_squared_floor"],
        tuple(QPOLY),True,
        before.parameter_coordinate_sign_checks+
        after.parameter_coordinate_sign_checks,
        C0_EPSILON,tuple(snapshots))
    return {
        "status":"PASS: exact relative Brouwer-degree cancellation through a genuine positive-size square pair death",
        "time_interval":[str(START),str(END)],
        "algebraic_collision_polynomial":[str(x) for x in QPOLY],
        "algebraic_collision_strictly_inside_time_bracket":True,
        "outer_parameter_box":[[str(z) for z in bounds] for bounds in BOX],
        "strict_root_to_outer_face_margin":str(INTERIOR_GAP),
        "left_cell_square_orbit_index":before.orientation_sign,
        "right_cell_square_orbit_index":after.orientation_sign,
        "signed_pre_collision_sum":certificate.signed_degree_before_event,
        "signed_post_collision_sum":certificate.signed_degree_after_event,
        "glued_degree_at_nonsmooth_algebraic_event":certificate.glued_degree_through_event,
        "inner_interface_undefined_cell_degrees_at_event":True,
        "certified_uniform_outer_boundary_residual_lower":str(certificate.strict_boundary_nonvanishing_floor),
        "uniform_continuous_perturbation_epsilon":str(certificate.robustness_epsilon),
        "uniform_perturbed_square_side_squared_lower_if_a_root_exists":str(certificate.uniform_side_squared_lower),
        "exact_polynomial_parameter_sign_checks":certificate.all_time_root_box_sign_checks,
        "independent_all_square_census_regular_times":[
            {"time":str(t),"total_geometric_square_orbits":n}
            for t,n in snapshots],
        "positive_square_collision_but_degree_zero":True,
        "scope":"Proves degree ADDITIVITY/CANCELLATION on one explicit, full-time, nonsmooth 4D local event box. Degree zero does NOT force a square, and robustness under small C0 perturbation preserves only the zero degree, not square existence. No general event atlas, universal Big Peg c, or unrestricted Toeplitz proof.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
