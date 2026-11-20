"""v256-v257: exact quarter-turn boundary pairs and near-symmetry square certificates.

v256 theorem (topological, not dependent on polygons):
  Every Jordan curve whose bounded domain contains the origin meets its
  own 90-degree rotation about the origin. The proof uses Jordan nesting:
  two disjoint Jordan boundaries enclosing 0 would be strictly nested;
  rotating the inclusion FOUR times yields an impossible strict cycle.
  Hence it has two boundary points v and -Jv with |v|=|-Jv|>0.
  If the curve is centrally symmetric, v,Jv,-v,-Jv all belong to the
  boundary and constitute a CENTERED square. A disk centered at the
  symmetry point of squared radius r2 gives square side²>=2*r2.
  No convexity, differentiability, or star-shapedness is required.

v257 quantitative approximate version (not an exact Toeplitz theorem):
  For any Jordan curve with a centered disk of squared radius r2, let
  delta=d_H(boundary,-boundary) be the Hausdorff symmetry defect.
  The same v and -Jv are exact boundary points; two other boundary
  points are within delta of Jv and -v. Four boundary points are
  delta-close to a centered square of side >=sqrt(2)*r.
  An even-vertex polygon has a computable upper defect bound delta2:
  max_i |P_i+P_{i+n/2}|². The corresponding affine edge pairing
  gives EXPLICIT rational approximate corners, exact deviation checks,
  and distinctness whenever delta2<r2/4.

The theorem is a short Jordan-domain argument; computations below
verify examples/construct witnesses with exact Fraction arithmetic.
Never claim that approximate corners form an exact square.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
import json

from polygon_orbit_census_v223 import (
    pt, Point, add, sub, mul, dot, cross, _on_seg, validate_jordan_polygon
)


def J(p):
    return pt(-p.y,p.x)


def minus(p):
    return pt(-p.x,-p.y)


def point_location(P,p):
    """Exact ray parity: +1 inside, 0 on boundary, -1 outside."""
    n=len(P)
    if any(_on_seg(p,P[i],P[(i+1)%n]) for i in range(n)):
        return 0
    crossings=0
    for i in range(n):
        a,b=P[i],P[(i+1)%n]
        if (a.y>p.y)==(b.y>p.y):
            continue
        x=a.x+(p.y-a.y)*(b.x-a.x)/(b.y-a.y)
        if x>p.x:
            crossings+=1
    return 1 if crossings%2 else -1


def distance_squared_point_segment(p,a,b):
    e=sub(b,a)
    length2=dot(e,e)
    if length2==0:
        raise ValueError("zero segment")
    t=dot(sub(p,a),e)/length2
    t=max(Q(0),min(Q(1),t))
    q=add(a,mul(e,t))
    return dot(sub(p,q),sub(p,q))


def disk_certificate(P,r2):
    """Prove CLOSED disk centered at 0, squared radius r2, lies in domain.

    If 0 is strictly inside a Jordan polygon and all closed edges stay
    at Euclidean distance >=sqrt(r2) from 0, the entire closed disk is
    in the bounded closed domain. A path from 0 to a point outside must
    cross an edge within that disk, a contradiction.
    """
    P=validate_jordan_polygon(tuple(P))
    if type(r2) is not Q or r2<=0:
        raise ValueError("positive exact squared radius required")
    origin=pt(0,0)
    if point_location(P,origin)!=1:
        raise ValueError("symmetry center is not strictly interior")
    distances=tuple(
        distance_squared_point_segment(origin,P[i],P[(i+1)%len(P)])
        for i in range(len(P))
    )
    if min(distances)<r2:
        raise ValueError("centered disk is not contained")
    return min(distances)


@dataclass(frozen=True)
class QuarterTurnIntersection:
    point: Point
    rotated_source: Point  # w = -J(point), an ACTUAL point on P
    edge_index: int
    other_edge_index: int
    first_parameter: Q
    other_parameter: Q


def _parameter(p,a,b):
    edge=sub(b,a)
    l2=dot(edge,edge)
    if not l2:
        raise ValueError("zero edge")
    return dot(sub(p,a),edge)/l2


def _segment_intersections(a,b,c,d):
    """Exact intersections, with two segment parameters; overlap endpoints.

    Returns all isolated intersections and endpoints for collinear
    overlaps. This suffices to find a boundary point common to both
    polygonal Jordan curves, without numerical root finding.
    """
    r=sub(b,a);s=sub(d,c)
    determinant=cross(r,s)
    delta=sub(c,a)
    if determinant:
        t=cross(delta,s)/determinant
        u=cross(delta,r)/determinant
        if 0<=t<=1 and 0<=u<=1:
            return ((add(a,mul(r,t)),t,u),)
        return ()
    if cross(delta,r)!=0:
        return ()
    candidates={}
    for p in (a,b,c,d):
        if _on_seg(p,a,b) and _on_seg(p,c,d):
            candidates[p]=(_parameter(p,a,b),_parameter(p,c,d))
    return tuple((p,t,u) for p,(t,u) in candidates.items())


def quarter_turn_pairs(P):
    """All discrete exact boundary intersection witnesses of P and J(P).

    If edges overlap, return their finite endpoints as witnesses, not
    a fictitious finite enumeration of an intersection continuum.
    """
    P=validate_jordan_polygon(tuple(P))
    # Precondition center interior, needed for intersection existence.
    if point_location(P,pt(0,0))!=1:
        raise ValueError("rotation center must lie in polygon interior")
    rotated=tuple(J(p) for p in P)
    found={}
    n=len(P)
    for i in range(n):
        a,b=P[i],P[(i+1)%n]
        for j in range(n):
            c,d=rotated[j],rotated[(j+1)%n]
            for v,t,u in _segment_intersections(a,b,c,d):
                w=minus(J(v))  # v=J(w); thus w=J^{-1}(v)
                if not (_on_seg(v,a,b)
                        and _on_seg(w,P[j],P[(j+1)%n])):
                    raise AssertionError("quarter-turn intersection invalid")
                found.setdefault(v,QuarterTurnIntersection(v,w,i,j,t,u))
    if not found:
        raise AssertionError("Jordan rotation-intersection theorem contradicted")
    return tuple(sorted(found.values(),key=lambda z:(z.point.x,z.point.y)))


def _on_boundary(P,p):
    return any(_on_seg(p,P[i],P[(i+1)%len(P)])
               for i in range(len(P)))


def _antipodal_shift(P):
    n=len(P)
    if n%2 or n<4:
        raise ValueError("requires even vertex count >=4")
    return n//2


def is_centrally_symmetric(P):
    half=_antipodal_shift(P)
    return all(P[i+half]==minus(P[i]) for i in range(half))


@dataclass(frozen=True)
class CenteredSquareCertificate:
    vertices:tuple
    side_squared:Q
    disk_radius_squared:Q
    quarter_turn_intersection_count:int
    polygon_centrally_symmetric:bool


def certify_centered_square(P,r2):
    """Exact witness for v256, sharp lower bound for centered disk."""
    P=validate_jordan_polygon(tuple(P))
    disk_certificate(P,r2)
    if not is_centrally_symmetric(P):
        raise ValueError("central symmetry prerequisite fails")
    pair=max(quarter_turn_pairs(P),
             key=lambda z:dot(z.point,z.point))
    v=pair.point
    corners=(v,J(v),minus(v),minus(J(v)))
    if not all(_on_boundary(P,p) for p in corners):
        raise AssertionError("symmetry does not close square")
    h2=dot(v,v)
    side2=2*h2
    if side2<2*r2 or not h2:
        raise AssertionError("sharp centered disk bound fails")
    return CenteredSquareCertificate(
        corners,side2,r2,len(quarter_turn_pairs(P)),True)


@dataclass(frozen=True)
class ApproximateSquareCertificate:
    exact_boundary_points:tuple   # all four actual points on P
    ideal_centered_square:tuple    # an exact geometric square (not necessarily on P)
    paired_vertex_defect_squared:Q
    max_corner_error_squared:Q
    ideal_side_squared:Q
    min_centered_disk_radius_squared:Q
    quarter_turn_intersection_count:int
    exact_square_claimed:bool


def _interpolate_edge(P,i,t):
    n=len(P)
    return add(P[i],mul(sub(P[(i+1)%n],P[i]),t))


def certify_near_square(P,r2):
    """v257 rational polygon certificate using half-cycle symmetry defect.

    For any even-vertex Jordan polygon, index-opposite edge pairing
    gives d_H(P,-P) <= delta, where delta² <= max_i|P_i+P_{i+m}|².
    C0-close approximate square corners are found explicitly.
    """
    P=validate_jordan_polygon(tuple(P))
    disk_certificate(P,r2)
    m=_antipodal_shift(P)
    n=len(P)
    delta2=max(dot(add(P[i],P[(i+m)%n]),add(P[i],P[(i+m)%n]))
               for i in range(n))
    if delta2>=r2/4:
        raise ValueError("symmetry defect too large for certified distinct corners")
    pairs=quarter_turn_pairs(P)
    witness=max(pairs,key=lambda z:dot(z.point,z.point))
    v=witness.point
    w=witness.rotated_source
    # v on edge i; w = -J(v) on edge j. Their opposite partners
    # are the affinely matched points on i+m and j+m.
    near_minus_v=_interpolate_edge(P,(witness.edge_index+m)%n,
                                   witness.first_parameter)
    near_minus_w=_interpolate_edge(P,(witness.other_edge_index+m)%n,
                                   witness.other_parameter)
    ideal=(v,J(v),minus(v),w)
    actual=(v,near_minus_w,near_minus_v,w)
    if w!=minus(J(v)):
        raise AssertionError("orthogonal radius pair is not exact")
    if not all(_on_boundary(P,p) for p in actual):
        raise AssertionError("approximate corner not on actual curve")
    errors=tuple(dot(sub(a,b),sub(a,b)) for a,b in zip(actual,ideal))
    if max(errors)>delta2:
        raise AssertionError("antipodal affine edge correspondence bound failed")
    if len(set(actual))!=4:
        raise AssertionError("defect bound failed to ensure four distinct corners")
    return ApproximateSquareCertificate(
        actual,ideal,delta2,max(errors),2*dot(v,v),r2,len(pairs),False)


# Nonconvex, centrally symmetric polygon with TWO reentrant rectangular
# notches. It encloses the closed disk of radius 2 centered at (0,0),
# but the ray 0 -> (4,4) passes through the upper notch: NOT star-shaped
# about its center. Its even-vertex opposite index pairing is exact.
NOTCHED_12=tuple(pt(*xy) for xy in (
    (-4,-4),(-3,-4),(-3,-2),(-1,-2),(-1,-4),(4,-4),
    (4,4),(3,4),(3,2),(1,2),(1,4),(-4,4),
))


def audit():
    P=NOTCHED_12
    if point_location(P,pt(Q(5,2),Q(5,2)))!=-1:
        raise AssertionError("non-star-shaped adversarial example failed")
    central=certify_centered_square(P,Q(4))
    perturbed=list(P)
    perturbed[3]=add(perturbed[3],pt(Q(1,20),Q(1,40)))
    perturbed=tuple(perturbed)
    if is_centrally_symmetric(perturbed):
        raise AssertionError("perturbation did not break central symmetry")
    near=certify_near_square(perturbed,Q(4))
    if near.max_corner_error_squared>near.paired_vertex_defect_squared:
        raise AssertionError("near-square error bound failed")
    return {
        "status":"PASS: v256 exact square on non-star-shaped symmetric Jordan polygon; v257 rigorously bounded near-square on asymmetric polygon",
        "central_example_vertices":len(P),
        "central_example_nonconvex_and_not_star_shaped_about_centre":True,
        "central_example_certified_disk_radius_squared":"4",
        "central_example_centered_square_side_squared":str(central.side_squared),
        "central_example_square_coordinates":[[str(p.x),str(p.y)]
                                               for p in central.vertices],
        "near_example_symmetry_broken":True,
        "near_example_max_antipodal_vertex_mismatch_squared":str(near.paired_vertex_defect_squared),
        "near_example_max_actual_corner_error_squared":str(near.max_corner_error_squared),
        "near_example_ideal_square_side_squared":str(near.ideal_side_squared),
        "near_example_all_four_points_distinct_on_actual_boundary":True,
        "near_example_exact_square_claimed":False,
        "topological_scope":"For EVERY Jordan curve enclosing the rotation center, gamma meets J gamma; with central symmetry this yields an exact centered square, side²>=2r² for a disk centered there.",
        "approximate_scope":"Any Jordan curve with symmetry Hausdorff defect delta and a centered disk has four distinct boundary points delta-close to square corners if delta<r/2. This does NOT establish an exact inscribed square.",
        "unrestricted_Toeplitz":"OPEN",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
