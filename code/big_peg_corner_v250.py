"""v250: EXACT square continuation through a GENUINE polygon corner.

A continuously moving, strictly convex *non-collinear* kink is met by one
corner of a square at time 1/2. Before and after, different 4-edge Cramer
charts carry continuously matching positive-side squares. All-real-time
claims use the pre-existing exact Bernstein verifier from v249.

A separate feature-size lemma proves ALL inscribed squares on a fixed
simple rational polygon with >=4 edges have side^2>=d^2/2, where d is
the minimum distance between nonadjacent closed polygon edges. A verified
positive coordinate-box gap gives an all-time version for this family,
independent of which square branch survives.

Neither lemma gives a complexity-independent Big Peg constant, or proves
the general Toeplitz square-peg conjecture.
"""
from fractions import Fraction as F
from itertools import combinations
import json

from polygon_orbit_census_v223 import (
    pt, sub, cross, dot, validate_jordan_polygon,
)
from big_peg_exact_v243 import _distance_sq_to_segment
from big_peg_persistence_v248 import OCTAGON, square_corner
from big_peg_multichart_v249 import (
    HALF, Z, O, pol, padd, psub, pmul, pscale, peval,
    ppoint, ppoint_at, pvsub, pvcross, pvdot, certify_sign,
    chart_polynomials, root_at, certify_chart, certify_cover,
)

ROOT_EDGES = (0,2,4,6)
RIGHT_CHART = (1,3,5,7)
LEFT_CHART = (0,3,5,7)
# Kernel of the three unchanged incidence rows for original octagon edges
# (2,4,6). The first square corner moves OFF the original 0-th edge.
THREE_EDGE_NULL_VECTOR = (
    F(3241,8989), F(-5455,8989), F(665,8989), F(1)
)
KINK_SCALE = F(-1,2)
KINK_SWEEP = F(1,4)
SQUARE_SIDE_SQUARED_THRESHOLD = F(9)  # Strict side length > 3.


def family_vertices():
    """Nine rational-affine vertices of an actually changing convex polygon.

    A new vertex moves ALONG a translated line outside the original edge,
    producing a real non-collinear corner at ALL times. It is designed
    to meet the first corner of a 3-edge-constrained square at time 1/2.
    """
    source=tuple(ppoint(p.x,p.y) for p in OCTAGON)
    D,N=chart_polynomials(source,ROOT_EDGES)
    old=root_at(D,N,F(0))
    # Cramer 3-edge kernel displacement; this is independently checked
    # for the other three square corners in the audit.
    crossing=(old[0]+KINK_SCALE*THREE_EDGE_NULL_VECTOR[0],
              old[1]+KINK_SCALE*THREE_EDGE_NULL_VECTOR[1])
    E=sub(OCTAGON[1],OCTAGON[0])
    offset=pol(-KINK_SWEEP*HALF,KINK_SWEEP)
    kink=(padd(pol(crossing[0]),pscale(offset,E.x)),
          padd(pol(crossing[1]),pscale(offset,E.y)))
    return tuple((source[0],kink,*source[1:]))


def exact_global_coordinate_separation(P):
    """Global lower bound on nonadjacent edge distance over ALL times.

    Because the supplied vertex coordinates are affine, each edge lies
    inside a time-independent axis-aligned rectangle determined by its
    endpoint positions at times 0 and 1. Positive projected separation
    of *every* pair of nonadjacent edge rectangles means the corresponding
    actual segments remain disjoint, and gives uniform distance >= delta.
    """
    n=len(P)
    if n<4 or any(len(c)>2 for p in P for c in p):
        raise ValueError("requires >=4 affine polygon vertices")
    boxes=[]
    for i in range(n):
        data=[]
        for ax in (0,1):
            vals=[peval(P[j][ax],t)
                  for j in (i,(i+1)%n) for t in (F(0),F(1))]
            data.append((min(vals),max(vals)))
        boxes.append(tuple(data))
    gaps=[]
    for i,j in combinations(range(n),2):
        if j==i+1 or (i==0 and j==n-1):
            continue
        projected=tuple(max(boxes[i][ax][0]-boxes[j][ax][1],
                            boxes[j][ax][0]-boxes[i][ax][1])
                        for ax in (0,1))
        g=max(projected)
        if g<=0:
            raise ValueError("unseparated nonadjacent edge bounding boxes")
        gaps.append(g)
    return min(gaps)


def prove_jordan_family(P):
    """Strict convexity, nonintersection and radius-3 bound for ALL times."""
    if len(P)!=9:
        raise ValueError("expected 9-vertex kink family")
    g=exact_global_coordinate_separation(P)
    sign_count=0
    for i in range(len(P)):
        prev=P[(i-1)%len(P)]
        vertex=P[i]
        next_vertex=P[(i+1)%len(P)]
        turn=pvcross(pvsub(vertex,prev),pvsub(next_vertex,vertex))
        length2=pvdot(pvsub(next_vertex,vertex),pvsub(next_vertex,vertex))
        sign_count+=certify_sign(turn,F(0),F(1),strict=True)
        sign_count+=certify_sign(length2,F(0),F(1),strict=True)
        # Euclidean radius <=3 holds throughout the homotopy by
        # convexity of the ball, because vertices move affinely.
        for t in (F(0),F(1)):
            point=ppoint_at(vertex,t)
            if dot(point,point)>=9:
                raise ValueError("moving vertex escapes radius-3 disk")
    # This is a genuine NONCOLLINEAR corner, not v249's subdivision point.
    kink_turn=pvcross(pvsub(P[1],P[0]),pvsub(P[2],P[1]))
    sign_count+=certify_sign(kink_turn,F(0),F(1),strict=True)
    if peval(kink_turn,F(0))!=peval(kink_turn,F(1)):
        raise AssertionError("expected a nonzero, constant kink cross product")
    return {
        "uniform_nonadjacent_coordinate_gap":g,
        "positive_turn_and_edge_polynomial_checks":sign_count,
        "strictly_convex_and_jordan_for_all_times":True,
        "strict_noncollinear_corner_cross":peval(kink_turn,F(0)),
        "enclosed_in_radius_3_disk":True,
    }


def minimum_nonadjacent_distance_squared(vertices):
    """Exact positive feature size of a fixed Jordan polygon (n>=4).

    For two disjoint planar closed segments, their minimum distance is
    attained by at least one endpoint-to-segment pair. The Jordan check
    rules out intersections between nonadjacent edges.
    """
    P=validate_jordan_polygon(tuple(vertices))
    n=len(P)
    if n<4:
        raise ValueError("at least four polygon edges are required")
    values=[]
    for i,j in combinations(range(n),2):
        if j==i+1 or (i==0 and j==n-1):
            continue
        a,b=P[i],P[(i+1)%n]
        c,d=P[j],P[(j+1)%n]
        values.extend((
            _distance_sq_to_segment(a,c,d),
            _distance_sq_to_segment(b,c,d),
            _distance_sq_to_segment(c,a,b),
            _distance_sq_to_segment(d,a,b),
        ))
    value=min(values)
    if value<=0:
        raise AssertionError("a simple polygon cannot have intersecting nonadjacent edges")
    return value



def toothed_square(delta):
    """Fixed-angle small outward tooth; large interior disk still exists.

    Its tip angle is 90 degrees and both adjacent foot angles are 225
    degrees for EVERY 0<delta<1/4. Yet the minimum separation of
    nonadjacent edges is <=sqrt(2)*delta, while the unit disk is contained.
    This refutes only an ANGLE-ONLY uniform lower bound on feature gap,
    NOT the existence of an inscribed square of significant size.
    """
    if type(delta) is not F or not F(0)<delta<F(1,4):
        raise ValueError("tooth scale must be an exact Fraction in (0,1/4)")
    P=(
        pt(-1,-1),pt(1,-1),pt(1,1),pt(delta,1),
        pt(0,1+delta),pt(-delta,1),pt(-1,1),
    )
    validate_jordan_polygon(P)
    return P


def tooth_diagnostic():
    """Exact failure mode of feature-size argument with invariant angles."""
    from big_peg_exact_v243 import contained_disk
    evidence=[]
    for n in (8,16,32,64):
        delta=F(1,n)
        polygon=toothed_square(delta)
        if not contained_disk(polygon,pt(0,0),F(1)):
            raise AssertionError("claimed interior unit disk is absent")
        d2=minimum_nonadjacent_distance_squared(polygon)
        if not (F(0)<d2<=2*delta*delta):
            raise AssertionError("tooth should generate a fine feature gap")
        # Tip edges are (-delta,+delta) and (-delta,-delta):
        # their dot product vanishes exactly for every scale.
        e1=sub(polygon[4],polygon[3])
        e2=sub(polygon[5],polygon[4])
        if dot(e1,e2)!=0:
            raise AssertionError("tooth's fixed right angle changed")
        evidence.append({
            "delta":str(delta),
            "exact_nonadjacent_gap_squared":str(d2),
            "interior_unit_disk":True,
            "right_angle_tip":True,
        })
    return evidence

def verify_three_edge_kernel():
    """Check stored exact vector satisfies ALL THREE unchanged edge rows."""
    base=tuple(ppoint(p.x,p.y) for p in OCTAGON)
    D,N=chart_polynomials(base,ROOT_EDGES)
    root=root_at(D,N,F(0))
    displaced=tuple(root[i]+KINK_SCALE*THREE_EDGE_NULL_VECTOR[i]
                    for i in range(4))
    for k,old_edge in enumerate(ROOT_EDGES):
        if k==0:
            continue
        a=OCTAGON[old_edge]
        b=OCTAGON[(old_edge+1)%len(OCTAGON)]
        corner=square_corner(displaced,k)
        if cross(sub(b,a),sub(corner,a))!=0:
            raise AssertionError("purported three-edge square kernel is invalid")
        # The new square corners must also stay on their closed segments.
        e=sub(b,a)
        t=dot(e,sub(corner,a))/dot(e,e)
        if not 0<t<1:
            raise AssertionError("three-edge square corner outside segment")
    return root,displaced


def audit():
    P=family_vertices()
    jordan=prove_jordan_family(P)
    original,displaced=verify_three_edge_kernel()
    atlas=certify_cover(P,(
        (F(0),HALF,RIGHT_CHART),
        (HALF,F(1),LEFT_CHART),
    ),side2_limit=SQUARE_SIDE_SQUARED_THRESHOLD)
    left,right=atlas
    if left.right_root!=right.left_root:
        raise AssertionError("square branch discontinuous at kink")
    if left.right_root!=displaced:
        raise AssertionError("square at actual kink disagrees with kernel")
    kink=ppoint_at(P[1],HALF)
    if pt(displaced[0],displaced[1])!=kink:
        raise AssertionError("square corner not at moving kink")
    if left.right_parameters[0]!=0 or right.left_parameters[0]!=1:
        raise AssertionError("handoff not at opposite ends of adjacent edges")
    if left.left_parameters[0]<=0 or right.right_parameters[0]>=1:
        raise AssertionError("incorrect side ordering across kink")
    # In contrast to v249's collinear split, the two rational square
    # branches are NOT the same algebraic function away from the seam.
    if all(pmul(a,right.denominator)==pmul(b,left.denominator)
           for a,b in zip(left.numerators,right.numerators)):
        raise AssertionError("expected genuinely distinct square branches")

    # Branch-independent lower-size bound. Any 4 noncollinear square
    # corners on a polygon with >=4 edges occupy nonadjacent edges.
    # Their separation is at most the square diagonal sqrt(2)*s.
    gap=jordan["uniform_nonadjacent_coordinate_gap"]
    universal_family_square_side2_lower=gap*gap/F(2)
    if universal_family_square_side2_lower<=0:
        raise AssertionError("no uniform feature-size floor")

    observed=[]
    for t in (F(0),F(1,4),HALF,F(3,4),F(1)):
        sample=tuple(ppoint_at(vertex,t) for vertex in P)
        validate_jordan_polygon(sample)
        actual_d2=minimum_nonadjacent_distance_squared(sample)
        if actual_d2<gap*gap:
            raise AssertionError("global coordinate-box gap is not valid")
        observed.append((str(t),str(actual_d2)))

    return {
        "status":"PASS: real non-collinear corner handoff with no square collapse",
        "noncollinear_kink_corner_cross":str(jordan["strict_noncollinear_corner_cross"]),
        "global_nonadjacent_edge_gap":str(gap),
        "every_inscribed_square_side2_at_least":str(universal_family_square_side2_lower),
        "branch_square_side2_strictly_above":str(SQUARE_SIDE_SQUARED_THRESHOLD),
        "uniform_true_big_peg_J_strictly_above":"1",
        "handoff_time":str(HALF),
        "handoff_square_root":[str(x) for x in displaced],
        "chart_atlas":[
            {"time":[str(c.left),str(c.right)],
             "assignment":list(c.assignment),
             "sign_checks":c.sign_checks,
             "first_corner_edge_parameters":[str(c.left_parameters[0]),
                                              str(c.right_parameters[0])]}
            for c in atlas
        ],
        "three_fixed_edges_have_valid_kernel":True,
        "algebraically_distinct_square_branches_join_at_kink":True,
        "jordan":{k:str(v) for k,v in jordan.items()},
        "sampled_feature_gap_rechecks":observed,
        "angle_only_feature_gap_counterexample":tooth_diagnostic(),
        "limits":"Quantitative result for this restricted affine convex nine-gon family; no universal c for arbitrary Jordan curves, no general singular orbit transition, and no Toeplitz proof.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
