"""v261: an oriented two-region atlas covering v254's three tracked squares.

The disappearing pair forms a SINGLE glued region with total degree zero
(v260); the persistent third square is isolated in a DISJOINT rational
four-box and has Brouwer degree +1. Their union has degree +1 on the
fixed, labeled, counterclockwise boundary-parameter lift.

This is NOT a global degree on unlabeled square configurations. The
cyclic shift of four labeled corners reverses the 4D source orientation
but preserves the orientation of the four square-defect equations.
Hence the integer index changes sign under this relabeling; only its
mod-2 reduction is immediately invariant under that Z4 identification.

Every time-uniform existence claim below follows from exact rational
Bernstein signs and the pre-existing v260 pair-death certificate. The
full moving-pentagon square census at every real time is NOT claimed.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
import json

from polygon_orbit_census_v223 import pt,add,sub,mul
from big_peg_pair_death_v254 import family_vertices,PERSISTENT
from big_peg_multichart_v249 import pol,psub,pscale,peval,certify_sign
from moving_polygon_events_v224 import chart_polynomials as edge_chart_polynomials
from near_square_transversality_v258 import det4
from square_degree_cube_v259 import square_defect
from square_relative_degree_v260 import (
    START,END,BOX as PAIR_BOX,C0_EPSILON,
    BOUNDARY_RESIDUAL_FLOOR as PAIR_RESIDUAL_FLOOR,
    PERTURBATION_RESIDUAL_RESPONSE,
    _snapshot,_matrix,audit as pair_audit,
)

# Coordinates u=(u0,u1,u2,u3) run in the positive directions along
# polygon edges (0,1,2,4), respectively. The associated unwrapped
# global ordered boundary parameters are (u0,1+u1,2+u2,4+u3).
PERSISTENT_BOX=(
    (Q(14,100),Q(17,100)),
    (Q(30,100),Q(34,100)),
    (Q(78,100),Q(84,100)),
    (Q(27,100),Q(31,100)),
)
PERSISTENT_GAP=Q(1,1000)
PERSISTENT_DET_FLOOR=Q(1,5)
PERSISTENT_ENTRY_BOUND=Q(2)
PERSISTENT_INVERSE_BOUND=(
    Q(24)*PERSISTENT_ENTRY_BOUND**3/PERSISTENT_DET_FLOOR
)
PERSISTENT_RESIDUAL_FLOOR=PERSISTENT_GAP/PERSISTENT_INVERSE_BOUND
SIDE_X_FLOOR=Q(1)
SIDE_Y_FLOOR=Q(9,100)

# F(p0,p1,p2,p3)=(p2-p1-J(p1-p0), p3-p0-J(p1-p0)).
# F after cyclically relabeling the corners is R*F before relabeling.
CYCLIC_RESIDUAL_TRANSFORM=(
    (Q(-1),Q(1),Q(1),Q(0)),
    (Q(-1),Q(-1),Q(0),Q(1)),
    (Q(0),Q(1),Q(0),Q(0)),
    (Q(-1),Q(0),Q(0),Q(0)),
)
CYCLIC_SOURCE_PERMUTATION=(
    (Q(0),Q(1),Q(0),Q(0)),
    (Q(0),Q(0),Q(1),Q(0)),
    (Q(0),Q(0),Q(0),Q(1)),
    (Q(1),Q(0),Q(0),Q(0)),
)


def persistent_points(P,t,s):
    """Four actual polygon boundary points, in the fixed CCW order."""
    if tuple(P)!=family_vertices() or len(s)!=4:
        raise ValueError("v261 only certifies the v254 pentagon and four parameters")
    polygon=_snapshot(P,t)
    return tuple(
        add(polygon[e],mul(sub(polygon[(e+1)%5],polygon[e]),u))
        for e,u in zip(PERSISTENT,s)
    )


def _unwrapped_q_persistent(s):
    return (s[0],Q(1)+s[1],Q(2)+s[2],Q(4)+s[3])


def _unwrapped_q_pair(s):
    return (s[0],Q(2)+s[1],Q(3)+s[2],Q(4)+s[3])


def global_boundary_point(P,t,q):
    """Continuous CCW polygon boundary point for 0<=q<5."""
    if tuple(P)!=family_vertices() or not (Q(0)<=q<Q(5)):
        raise ValueError("invalid pentagon or unwrapped parameter")
    vertices=_snapshot(P,t)
    j=int(q)
    return add(vertices[j],mul(sub(vertices[(j+1)%5],vertices[j]),q-j))


def _apply(matrix,vector):
    return tuple(sum(a*x for a,x in zip(row,vector)) for row in matrix)


def certify_cyclic_relabeling():
    """Check a universal R4 linear identity on the eight input basis vectors.

    Linearity gives F(cyclic(p))=R*F(p) for ALL four input points.
    The relabeling q -> (q1,q2,q3,q0+5) has derivative P, with
    det(P)=-1, while det(R)=+1. Therefore the signed LOCAL degree
    changes sign on the cyclic quotient. This prevents interpreting
    the atlas's +1 as a canonical integer index of an unlabeled square.
    """
    R=CYCLIC_RESIDUAL_TRANSFORM
    P=CYCLIC_SOURCE_PERMUTATION
    if det4(R)!=Q(1) or det4(P)!=Q(-1):
        raise AssertionError("wrong cyclic domain/target orientation")
    for k in range(8):
        p=tuple(pt(Q(int(k==2*i)),Q(int(k==2*i+1)))
                for i in range(4))
        defect=square_defect(p)
        cycled=square_defect(p[1:]+p[:1])
        if cycled!=_apply(R,defect):
            raise AssertionError("cyclic square-defect identity fails")
    return {"source_jacobian_sign":-1,
            "defect_transform_determinant":1,
            "local_degree_relabeling_factor":-1,
            "number_of_exact_basis_checks":8}


@dataclass(frozen=True)
class PersistentCertificate:
    orientation_sign:int
    parameter_sign_checks:int
    exact_endpoint_roots:tuple
    residual_floor:Q
    side_squared_floor:Q
    inverse_bound:Q


def certify_persistent_box(P,box=PERSISTENT_BOX,gap=PERSISTENT_GAP,
                           epsilon=C0_EPSILON):
    """Certify the persistent square has +1 degree for EVERY real time.

    The source is one affine 4D cell. Compute its unique Cramer root
    with the v224 polynomial chart engine; prove each coordinate
    strictly lies between the fixed box faces with margin gap by
    exact rational Bernstein polynomial positivity on [0,1].
    The root never touches the box boundary, so the +1 persists
    through the nonsmooth pair-death event elsewhere in configuration
    space. The C0 budget is a conservative direct cofactor estimate.
    """
    if tuple(P)!=family_vertices():
        raise ValueError("only the fully certified v254 family")
    if tuple(tuple(x) for x in box)!=PERSISTENT_BOX or gap!=PERSISTENT_GAP:
        raise ValueError("uncertified persistent square isolating region")
    if type(epsilon) is not Q or epsilon<Q(0):
        raise ValueError("nonnegative exact Fraction epsilon required")
    if PERTURBATION_RESIDUAL_RESPONSE*epsilon>=PERSISTENT_RESIDUAL_FLOOR:
        raise ValueError("C0 perturbation could reach the outer boundary")
    if (SIDE_X_FLOOR-2*epsilon)<=0 or (SIDE_Y_FLOOR-2*epsilon)<=0:
        raise ValueError("side-component lower bound no longer positive")

    start,end=_snapshot(P,START),_snapshot(P,END)
    chart=edge_chart_polynomials(
        start,end,PERSISTENT,orientation=1,cutoff_squared=Q(1))
    D=chart.determinant
    N=chart.numerators
    checks=certify_sign(
        psub(D,pol(PERSISTENT_DET_FLOOR)),Q(0),Q(1),strict=True)
    for t,normalized in ((START,Q(0)),(END,Q(1))):
        M=_matrix(P,t,PERSISTENT)
        if any(abs(x)>PERSISTENT_ENTRY_BOUND for row in M for x in row):
            raise AssertionError("uniform affine matrix entry bound fails")
        if det4(M)!=peval(D,normalized):
            raise AssertionError("independent determinant disagrees")
    for j,(lo,hi) in enumerate(PERSISTENT_BOX):
        lower=psub(N[j],pscale(D,lo+gap))
        upper=psub(pscale(D,hi-gap),N[j])
        checks+=certify_sign(lower,Q(0),Q(1),strict=True)
        checks+=certify_sign(upper,Q(0),Q(1),strict=True)

    roots=tuple(tuple(peval(n,tau)/peval(D,tau) for n in N)
                for tau in (Q(0),Q(1)))
    for t,root in zip((START,END),roots):
        if not all(lo+gap<u<hi-gap
                   for u,(lo,hi) in zip(root,PERSISTENT_BOX)):
            raise AssertionError("endpoint root escaped rational box")
        if square_defect(persistent_points(P,t,root))!=(Q(0),)*4:
            raise AssertionError("exact root is not a geometric square")

    # The first-to-second square-side vector is multi-affine in
    # (t,u0,u1). Its coordinate extrema on the whole 3-box occur at
    # the 2^3 rational vertices; u2,u3 do not enter this vector.
    for t,u0,u1 in product(
            (START,END),PERSISTENT_BOX[0],PERSISTENT_BOX[1]):
        vertices=_snapshot(P,t)
        first=add(vertices[0],mul(sub(vertices[1],vertices[0]),u0))
        second=add(vertices[1],mul(sub(vertices[2],vertices[1]),u1))
        side=sub(second,first)
        if not (side.x < -SIDE_X_FLOOR and side.y > SIDE_Y_FLOOR):
            raise AssertionError("positive-size uniform box floor fails")
    side_floor=(SIDE_X_FLOOR-2*epsilon)**2+(
        SIDE_Y_FLOOR-2*epsilon)**2
    if side_floor<=Q(1):
        raise AssertionError("noncollapse bound is insufficient")

    return PersistentCertificate(
        1,checks,roots,PERSISTENT_RESIDUAL_FLOOR,side_floor,
        PERSISTENT_INVERSE_BOUND,
    )


def audit():
    P=family_vertices()
    persistent=certify_persistent_box(P)
    pair=pair_audit()
    if pair["glued_degree_at_nonsmooth_algebraic_event"]!=0:
        raise AssertionError("v260 pair region degree unexpectedly nonzero")
    if pair["left_cell_square_orbit_index"]!=-1 or (
            pair["right_cell_square_orbit_index"]!=1):
        raise AssertionError("v260 orientations no longer match the atlas")
    if persistent.orientation_sign!=1:
        raise AssertionError("persistent index not +1")
    if not (PERSISTENT_BOX[0][1]<PAIR_BOX[0][0]):
        raise AssertionError("two regions overlap in first global coordinate")
    if not (Q(2)+PERSISTENT_BOX[2][1] <
            Q(3)+PAIR_BOX[2][0]):
        raise AssertionError("ordered global configuration regions overlap")

    # On each chart q_i=edge-index+u_i, with dq/du=+identity.
    # Pair's h has q1=2+h on BOTH sides of the real polygon kink.
    # Thus orientations and degrees are directly comparable in the
    # same labeled, unwrapped ordered-configuration space.
    for t in (START,END):
        for u in persistent.exact_endpoint_roots:
            q=_unwrapped_q_persistent(u)
            if not (Q(0)<q[0]<q[1]<q[2]<q[3]<Q(5)):
                raise AssertionError("persistent chart not cyclically ordered")
            direct=tuple(global_boundary_point(P,t,w) for w in q)
            if direct!=persistent_points(P,t,u):
                raise AssertionError("inconsistent global chart")
        for h in (Q(-1,1000),Q(0),Q(1,1000)):
            s=(Q(91,100),h,Q(44,100),Q(75,100))
            q=_unwrapped_q_pair(s)
            if not (Q(0)<q[0]<q[1]<q[2]<q[3]<Q(5)):
                raise AssertionError("pair chart not cyclically ordered")
    relabel=certify_cyclic_relabeling()
    total=persistent.orientation_sign+pair[
        "glued_degree_at_nonsmooth_algebraic_event"]
    if total!=1:
        raise AssertionError("disjoint region index should be +1")
    residual=min(persistent.residual_floor,PAIR_RESIDUAL_FLOOR)
    if not (PERTURBATION_RESIDUAL_RESPONSE*C0_EPSILON<residual):
        raise AssertionError("no robust union boundary")
    return {
        "status":"PASS: nonzero +1 relative degree in two separated, oriented four-parameter regions",
        "interval":[str(START),str(END)],
        "persistent_edges":list(PERSISTENT),
        "persistent_box":[[str(x) for x in v] for v in PERSISTENT_BOX],
        "colliding_pair_box":[[str(x) for x in v] for v in PAIR_BOX],
        "persistent_index":persistent.orientation_sign,
        "colliding_left_index":pair["left_cell_square_orbit_index"],
        "colliding_right_index":pair["right_cell_square_orbit_index"],
        "pair_glued_index":pair["glued_degree_at_nonsmooth_algebraic_event"],
        "oriented_disconnected_union_degree":total,
        "all_real_time_persistent_root_sign_checks":
            persistent.parameter_sign_checks,
        "persistent_uniform_boundary_defect_floor":
            str(persistent.residual_floor),
        "combined_uniform_outer_defect_floor":str(residual),
        "certified_C0_coordinatewise_epsilon":str(C0_EPSILON),
        "persistent_noncollapse_side_squared_floor":
            str(persistent.side_squared_floor),
        "surviving_square_endpoint_parameters":[
            [str(v) for v in root] for root in persistent.exact_endpoint_roots],
        "independent_whole_polygon_census_at_two_rational_slices":
            pair["independent_all_square_census_regular_times"],
        "cyclic_unlabeled_quotient_audit":relabel,
        "canonical_unlabeled_integer_degree_proved":False,
        "complete_all_real_time_global_square_census_proved":False,
        "scope":"Two selected disjoint labeled regions, each with an all-real-time external-boundary certificate. The surviving square has +1 local degree; the disappearing pair has 0. The +1 is NOT a canonical global integer degree on unlabeled square-orbit space because cyclic relabeling reverses sign. No arbitrary-Jordan noncollapse or unrestricted Toeplitz theorem.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
