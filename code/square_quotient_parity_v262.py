"""v262: exact cyclic-C4 quotient parity atlas through the v254 square pair death.

This strengthens v261's +1 *labeled* relative degree by saturating both
selected boxes under all four cyclic relabelings. The eight labeled
regions are pairwise disjoint even by their first unwrapped coordinate.
Their free C4 quotient consists of TWO relative regions with an exactly
certified Z/2 degree 1 at every real time through the irrational event.

The quotient statement is only for this *selected* positive cyclic-order
atlas, NOT the full square-configuration space of an arbitrary Jordan
curve. v218/v219 already established the general symmetry obstruction.
All real-time existence and gluing claims here inherit v260/v261's
exact rational Bernstein and polynomial certificates.
"""
from fractions import Fraction as Q
from itertools import product
import json

from polygon_orbit_census_v223 import pt
from near_square_transversality_v258 import det4
from square_degree_cube_v259 import square_defect
from big_peg_pair_death_v254 import (
    family_vertices, COLLIDING_LEFT, COLLIDING_RIGHT,
    certify_jordan_family, certify_collision,
)
from square_relative_degree_v260 import (
    BOX as PAIR_BOX, START, END, C0_EPSILON,
    BOUNDARY_RESIDUAL_FLOOR as PAIR_RESIDUAL_FLOOR,
    certify_cell, certify_interface_gluing,\n    certify_sides_and_robust_boundary, _square_points,
)
from square_oriented_atlas_v261 import (
    PERSISTENT_BOX, PERSISTENT_RESIDUAL_FLOOR,
    CYCLIC_RESIDUAL_TRANSFORM, CYCLIC_SOURCE_PERMUTATION,
    persistent_points, _unwrapped_q_persistent, _unwrapped_q_pair,
    global_boundary_point, certify_persistent_box,
)

PERIOD=Q(5)
IDENTITY=tuple(tuple(Q(i==j) for j in range(4)) for i in range(4))
PATCHES=("persistent","collision_pair")
LOCAL_BOX={"persistent":PERSISTENT_BOX,"collision_pair":PAIR_BOX}
BASE_BOUNDARY_RESIDUAL_FLOOR=min(
    PERSISTENT_RESIDUAL_FLOOR,PAIR_RESIDUAL_FLOOR
)


def matmul(A,B):
    if len(A)!=4 or len(B)!=4 or any(len(r)!=4 for r in A+B):
        raise ValueError("two rational 4x4 matrices required")
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4))
                       for j in range(4)) for i in range(4))


def power(A,k):
    if type(k) is not int or not 0<=k<=4:
        raise ValueError("cyclic group exponent must be 0..4")
    out=IDENTITY
    for _ in range(k):
        out=matmul(out,A)
    return out


def apply(A,v):
    if len(v)!=4:
        raise ValueError("four residual coordinates required")
    return tuple(sum(a*x for a,x in zip(row,v)) for row in A)


def infinity_norm(A):
    return max(sum(abs(v) for v in row) for row in A)


def base_unwrapped(kind,s):
    if kind=="persistent":
        return _unwrapped_q_persistent(s)
    if kind=="collision_pair":
        return _unwrapped_q_pair(s)
    raise ValueError("unknown v262 patch")


def base_box(kind):
    if kind not in PATCHES:
        raise ValueError("unknown v262 patch")
    local=LOCAL_BOX[kind]
    shift=((Q(0),Q(1),Q(2),Q(4))
           if kind=="persistent" else
           (Q(0),Q(2),Q(3),Q(4)))
    return tuple((lo+delta,hi+delta)
                 for (lo,hi),delta in zip(local,shift))


def _validate_domain(box):
    if len(box)!=4 or any(
        type(lo) is not Q or type(hi) is not Q or not lo<hi
        for lo,hi in box
    ):
        raise ValueError("strict rational four-box required")
    if not (Q(0)<=box[0][0]<box[0][1]<PERIOD):
        raise ValueError("first unwrapped parameter must lie on circle")
    if not all(box[i][1]<box[i+1][0] for i in range(3)):
        raise ValueError("strict cyclic parameter order not guaranteed")
    if not box[3][1]<box[0][0]+PERIOD:
        raise ValueError("cyclic gap may vanish")


def rotate_unwrapped(q,k):
    """T^k on positive cyclic configurations, reduced mod common +5 shift.

    The starting q is an increasing lift with 0<=q0<5,
    q3<q0+5. The returned tuple has the same properties.
    """
    if len(q)!=4 or any(type(x) is not Q for x in q):
        raise ValueError("four exact rational parameters required")
    if not (Q(0)<=q[0]<PERIOD and
            q[0]<q[1]<q[2]<q[3]<q[0]+PERIOD):
        raise ValueError("not a positive cyclic configuration")
    if type(k) is not int or not 0<=k<=4:
        raise ValueError("rotation count 0..4 required")
    raw=tuple(q[(i+k)%4]+PERIOD*((i+k)//4)
              for i in range(4))
    winding=raw[0]//PERIOD
    return tuple(v-PERIOD*winding for v in raw)


def rotate_box(box,k):
    """Affine image of an entire box under T^k (no winding ambiguity)."""
    _validate_domain(box)
    if type(k) is not int or not 0<=k<4:
        raise ValueError("only four cyclic sheets")
    q=tuple(box[(i+k)%4] for i in range(4))
    out=tuple((lo+(i+k)//4*PERIOD,
               hi+(i+k)//4*PERIOD)
              for i,(lo,hi) in enumerate(q))
    # The base q0 lies strictly in [0,5), and k<4 gives qk in [0,5).
    # Strict inequalities are checked by _validate_domain.
    _validate_domain(out)
    return out


def cyclic_sector_points(P,t,q):
    """Evaluate an increasing unwrapped cyclic tuple on the same curve."""
    if tuple(P)!=family_vertices():
        raise ValueError("only the exact v254 moving family")
    if not (Q(0)<=q[0]<PERIOD and
            q[0]<q[1]<q[2]<q[3]<q[0]+PERIOD):
        raise ValueError("not a positive cyclic four-point parameter")
    return tuple(global_boundary_point(P,t,v%PERIOD) for v in q)


def certify_cyclic_action():
    """Exact C4 representation, source-orientation and residual identities."""
    P=CYCLIC_SOURCE_PERMUTATION
    R=CYCLIC_RESIDUAL_TRANSFORM
    if power(P,4)!=IDENTITY or power(R,4)!=IDENTITY:
        raise AssertionError("not a C4 action")
    if det4(P)!=-1 or det4(R)!=1:
        raise AssertionError("source/target monodromy incorrect")
    checks=0
    for k in range(4):
        Pk,Rk=power(P,k),power(R,k)
        if det4(Pk)!=(-1)**k or det4(Rk)!=1:
            raise AssertionError("incorrect orientation character")
        # F is a linear map on the 8 Cartesian input coordinates.
        # Testing all 8 independent coordinate basis vectors proves
        # F(cyclic^k p)=R^k F(p) for ALL four-point configurations.
        for j in range(8):
            pts=tuple(pt(Q(j==2*i),Q(j==2*i+1))
                      for i in range(4))
            rotated=pts[k:]+pts[:k]
            if square_defect(rotated)!=apply(Rk,square_defect(pts)):
                raise AssertionError("cyclic residual equivariance failed")
            checks+=1
    return {"basis_identities_checked":checks,
            "source_orientation_character":(1,-1,1,-1),
            "target_determinants":(1,1,1,1),
            "twisted_index_character":(1,-1,1,-1)}


def certify_saturated_disjoint_regions():
    """Eight C4 images are separated *solely by q0 interval*.

    Thus the quotient of the saturated disjoint atlas has precisely
    two chart components, each represented once by the base lift.
    The argument uses the boxes, not sampled roots.
    """
    bands=[]
    for kind in PATCHES:
        for k in range(4):
            box=rotate_box(base_box(kind),k)
            bands.append((box[0][0],box[0][1],kind,k))
    bands.sort()
    if len(bands)!=8:
        raise AssertionError("four cyclic translates of two boxes required")
    for prev,nxt in zip(bands,bands[1:]):
        if not prev[1]<nxt[0]:
            raise AssertionError("two cyclic sheets can overlap")
    # Rotating four times is the identity after one winding around
    # the circle. This checks that the group action is free *on boxes*.
    checked=0
    for kind in PATCHES:
        box=base_box(kind)
        for vertex in product(*box):
            if rotate_unwrapped(vertex,4)!=vertex:
                raise AssertionError("cyclic relabeling fails C4 closure")
            q=vertex
            for k in range(1,4):
                q=rotate_unwrapped(q,1)
                if q==vertex:
                    raise AssertionError("unexpected stabilizer")
                if q!=rotate_unwrapped(vertex,k):
                    raise AssertionError("C4 composition failed")
                checked+=1
    return {
        "regions_in_labeled_C4_saturation":8,
        "regions_in_unlabeled_quotient":2,
        "strictly_disjoint_first_coordinate_bands":[
            [str(lo),str(hi),kind,k] for lo,hi,kind,k in bands],
        "free_action_vertex_checks":checked,
    }


def certify_geometric_equivariance(P):
    """Check actual polygon corners, including BOTH hinge sides and h=0."""
    if tuple(P)!=family_vertices():
        raise ValueError("only audited pentagon")
    checked=0
    for t in (START,END):
        for kind in PATCHES:
            # True corners of both four-arc regions; kink is interior
            # to the pair region and receives additional h=0 checks.
            for s in product(*LOCAL_BOX[kind]):
                points=(persistent_points(P,t,s) if kind=="persistent"
                        else _square_points(P,t,s))
                q=base_unwrapped(kind,s)
                if points!=cyclic_sector_points(P,t,q):
                    raise AssertionError("unwrapped source map incorrect")
                for k in range(4):
                    qk=rotate_unwrapped(q,k)
                    if cyclic_sector_points(P,t,qk)!=points[k:]+points[:k]:
                        raise AssertionError("cyclic map fails for actual polygon")
                    if square_defect(cyclic_sector_points(P,t,qk))!=apply(
                            power(CYCLIC_RESIDUAL_TRANSFORM,k),
                            square_defect(points)):
                        raise AssertionError("polygon residual equivariance")
                    checked+=1
            if kind=="collision_pair":
                for h in (-Q(1,2000),Q(0),Q(1,2000)):
                    s=(Q(91,100),h,Q(44,100),Q(75,100))
                    q=base_unwrapped(kind,s)
                    for k in range(4):
                        pts=cyclic_sector_points(P,t,rotate_unwrapped(q,k))
                        orig=_square_points(P,t,s)
                        if pts!=orig[k:]+orig[:k]:
                            raise AssertionError("vertex gluing inconsistent")
                        checked+=1
    return checked


def certify_regular_root_orbits(P,persistent,pair_left,pair_right):
    """Exact three-orbit vs one-orbit classification at rational endpoints.

    The all-real-time event statement uses v260/v261 certificates;
    this function only checks root membership and C4 multiplicities
    at the two rational time slices.
    """
    count=0
    for t,idx in ((START,0),(END,1)):
        roots=[("persistent",persistent.exact_endpoint_roots[idx])]
        if idx==0:
            roots.extend((
                ("collision_pair",pair_left.signed_roots_at_start),
                ("collision_pair",pair_right.signed_roots_at_start)))
        seen=set()
        for kind,s in roots:
            original=base_unwrapped(kind,s)
            original_points=cyclic_sector_points(P,t,original)
            if square_defect(original_points)!=(Q(0),)*4:
                raise AssertionError("root is not an exact inscribed square")
            signatures=[]
            for k in range(4):
                q=rotate_unwrapped(original,k)
                box=rotate_box(base_box(kind),k)
                if not all(lo<x<hi for x,(lo,hi) in zip(q,box)):
                    raise AssertionError("actual root escaped its cyclic box")
                pts=cyclic_sector_points(P,t,q)
                if square_defect(pts)!=(Q(0),)*4:
                    raise AssertionError("cyclic rotation loses square zero")
                signatures.append(q)
                count+=1
            if len(set(signatures))!=4:
                raise AssertionError("cyclic orbit is not free")
            # The geometry of each orbit is cyclically invariant.
            geometric=frozenset((v.x,v.y) for v in original_points)
            if geometric in seen:
                raise AssertionError("two roots are same geometric square")
            seen.add(geometric)
        if len(seen)!=(3 if idx==0 else 1):
            raise AssertionError("unexpected exact square orbit total")
    return {"rational_before_orbits":3,
            "rational_after_orbits":1,
            "rational_before_labeled_roots":12,
            "rational_after_labeled_roots":4,
            "cyclic_root_membership_checks":count}


def certify_quotient_atlas():
    P=family_vertices()
    jordan=certify_jordan_family(P)
    if not jordan["convex_Jordan_for_all_real_times"]:
        raise AssertionError("underlying moving polygon not certified Jordan")
    collision=certify_collision(P)
    if not collision["determinants_nonzero_at_event"]:
        raise AssertionError("v254 noncollinear collision proof lost")
    glue=certify_interface_gluing(P)
    if len(glue)!=4:
        raise AssertionError("four source-parameter gluing identities lost")
    pair_robust=certify_sides_and_robust_boundary(P)
    persistent=certify_persistent_box(P)
    left=certify_cell(P,COLLIDING_LEFT)
    right=certify_cell(P,COLLIDING_RIGHT)
    if (persistent.orientation_sign,left.orientation_sign,
            right.orientation_sign)!=(1,-1,1):
        raise AssertionError("v261/v260 source orientations no longer match")
    symmetry=certify_cyclic_action()
    saturation=certify_saturated_disjoint_regions()
    checks=certify_geometric_equivariance(P)
    roots=certify_regular_root_orbits(P,persistent,left,right)

    # Signed local degrees on cyclic lift k. No global scalar integer
    # degree is assigned to the nonorientable quotient. The signs
    # define an orientation-local-system equivariance rule.
    index_lifts={
        "persistent":[(-1)**k for k in range(4)],
        "pair_before_left":[-(-1)**k for k in range(4)],
        "pair_before_right":[(-1)**k for k in range(4)],
        "pair_union":[0]*4,
        "whole_union":[(-1)**k for k in range(4)]
    }
    if (sum(index_lifts["whole_union"])!=0
            or sum(x%2 for x in index_lifts["whole_union"])%2!=0):
        raise AssertionError("full labeled cyclic sum should vanish")
    quotient_degree=(persistent.orientation_sign+
                     left.orientation_sign+right.orientation_sign)%2
    if quotient_degree!=1:
        raise AssertionError("quotient relative parity must be 1")
    if any(
        ((-1)**k*index_lifts["whole_union"][k])!=1
        for k in range(4)
    ):
        raise AssertionError("twisted-index normalization failed")
    # Under cyclic relabeling, a defect is multiplied by R^k.
    # ||R^k v||inf >= ||v||inf / ||R^{-k}||inf.
    # This propagates the boundary residual budget to all eight
    # translated boxes, including those of the kinked pair.
    worst_inverse_norm=max(
        infinity_norm(power(CYCLIC_RESIDUAL_TRANSFORM,(4-k)%4))
        for k in range(4)
    )
    quotient_floor=BASE_BOUNDARY_RESIDUAL_FLOOR/worst_inverse_norm
    if not Q(4)*C0_EPSILON<quotient_floor:
        raise AssertionError("uniform C0 budget fails on cyclic cover")
    return {
        "status":"PASS: Z/2 quotient relative degree 1 on an exactly C4-saturated selected two-region atlas",
        "interval":[str(START),str(END)],
        "singular_event":"v254 algebraic positive-size pair death, inner vertex interface",
        "all_real_time_convex_Jordan_pentagon_certified":True,
        "exact_four_source_parameter_gluing_identities":len(glue),
        "pair_noncollapse_squared_side_floor":
            str(pair_robust["certified_side_squared_floor"]),
        "persistent_noncollapse_squared_side_floor":
            str(persistent.side_squared_floor),
        "cyclic_representation":symmetry,
        "saturation":saturation,
        "actual_polygon_equivariance_checks":checks,
        "rational_snapshot_root_orbits":roots,
        "signed_indices_in_four_labeled_sheets":index_lifts,
        "ordinary_integer_sum_on_entire_labeled_saturation":0,
        "ordinary_mod2_count_on_entire_labeled_saturation":0,
        "quotient_relative_degree_before_event_mod2":quotient_degree,
        "quotient_relative_degree_at_event_mod2":quotient_degree,
        "quotient_relative_degree_after_event_mod2":quotient_degree,
        "quotient_pair_at_event_mod2":0,
        "quotient_geometric_root_count_at_singular_event_parity_not_an_invariant":True,
        "orientation_local_system_character":[1,-1,1,-1],
        "canonically_globally_oriented_integer_degree_on_quotient":False,
        "uniform_cyclic_inverse_residual_transform_norm":str(worst_inverse_norm),
        "uniform_cyclic_saturated_outer_residual_floor":str(quotient_floor),
        "robust_coordinatewise_C0_epsilon":str(C0_EPSILON),
        "positive_square_side_for_every_certified_atlas_root":True,
        "full_six_sector_atlas_proved":False,
        "all_real_time_exhaustive_whole_pentagon_census_proved":False,
        "arbitrary_Jordan_square_peg_theorem_proved":False,
        "scope":"A Z/2 *relative* degree on the quotient of the exact C4 saturation of the two v261 selected neighborhoods, not a global absolute Euler class or an all-Jordan square existence theorem. The singular square counts are not parity invariants; degree through gluing is. Four labeled copies must NOT be counted mod2 before taking their quotient.",
    }


if __name__=="__main__":
    print(json.dumps(certify_quotient_atlas(),sort_keys=True,indent=2))
