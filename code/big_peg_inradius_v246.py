"""v246: EXACT global inradius intervals for the original and v245 Big Peg decagons.

A rational quadtree covers the entire axis-aligned polygon bounding rectangle.
Every leaf proves either (A) the box cannot meet the polygon's interior or
(B) each point in the box is within target radius of the polygon boundary.

Proof: distance to a closed set is 1-Lipschitz. The half diagonal of a
rational box is <= its rational L1 half-width + half-height, called H.
If d(midpoint,boundary) + H <= R then no point of the box can be a centre
of an inscribed disk of radius > R. If midpoint is exterior and
d(midpoint,boundary) > H, the box cannot meet the interior.
Otherwise subdivide. All strict/non-strict predicates use Fraction.

This gives a *global upper bound*, not an exact algebraic formula for the
inradius. An explicit rational centre and exact segment distances give a
lower bound. Together these prove the actual v245 ratio is strictly smaller.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
import json

from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from big_peg_exact_v243 import (
    _point_inside, _distance_sq_to_segment, maximum_square, read_polygon,
)

ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = ROOT / "research" / "v244-pochuev-original-ten-gon.json"
IMPROVED = ROOT / "research" / "v245-improved-ten-gon.json"

# Rational disks near the three-contact medial-axis maximizers.
# Their containment is PROVED by the existing point-in-polygon and exact
# distances-to-segments predicates, never inferred from optimization.
ORIGINAL_CENTER = pt(Q(45637, 100000000), Q(10485, 100000000))
IMPROVED_CENTER = pt(Q(302622629, 10**12), Q(2143102, 10**13))

# Global upper bounds on the radius (not on squared radius).
# Both are < 8e-9 above the corresponding independently checked disks.
ORIGINAL_RADIUS_UPPER = Q(100100053, 100000000)
IMPROVED_RADIUS_UPPER = Q(100109166, 100000000)

ORIGINAL_SIDE2 = Q(
    56381275521625791352241234309,
    28343589374868261752462760000,
)
IMPROVED_SIDE2 = Q(
    2263154862866037715637393248853633089,
    1137974336077532830058267809922880000,
)

# Exact radius-squared lower bounds for the particular rational centres.
ORIGINAL_RADIUS2_LOWER = Q(
    246132835013175286089001, 245641050000000000000000
)
IMPROVED_RADIUS2_LOWER = Q(
    68382681809394993823195766186089,
    68233625000000000000000000000000
)


@dataclass(frozen=True)
class QuadtreeProof:
    radius_upper: Q
    visited_boxes: int
    covered_boxes: int
    exterior_boxes: int
    split_boxes: int
    maximum_depth: int
    bounding_box: tuple

    @property
    def complete(self):
        return self.visited_boxes == (self.covered_boxes+self.exterior_boxes+
                                      self.split_boxes)


def min_boundary_distance_squared(P, center):
    """Minimum exact squared distance from centre to the CLOSED boundary."""
    return min(
        _distance_sq_to_segment(center, P[i], P[(i+1)%len(P)])
        for i in range(len(P))
    )


def certify_global_radius_upper(vertices, radius_upper, *,
                                max_depth=48, max_boxes=200000):
    """Exhaustively prove the inradius is <= the positive rational R.

    Raises, never returns 'success', if the rational subdivision fails
    within its prescribed limits. No floating-point computation is used.
    The interior test is exact and the entire polygon is within the root
    rectangle, so no possible disk centre is omitted.
    """
    P = validate_jordan_polygon(tuple(vertices))
    if type(radius_upper) is not Q or radius_upper<=0:
        raise ValueError("supply an exact positive Fraction radius")
    if type(max_depth) is not int or max_depth<0 or type(max_boxes) is not int or max_boxes<1:
        raise ValueError("max_depth and max_boxes must be positive integers")
    xmin=min(p.x for p in P); xmax=max(p.x for p in P)
    ymin=min(p.y for p in P); ymax=max(p.y for p in P)
    R=radius_upper
    stack=[(xmin,xmax,ymin,ymax,0)]
    visited=covered=exterior=split=deepest=0

    while stack:
        if visited>=max_boxes:
            raise RuntimeError("global inradius upper certificate exceeded box limit")
        xl,xh,yl,yh,depth=stack.pop()
        visited+=1
        deepest=max(deepest,depth)
        mx=(xl+xh)/2; my=(yl+yh)/2
        centre=pt(mx,my)
        # Every box point is at most sqrt(hx^2+hy^2) <= hx+hy from midpoint.
        H=(xh-xl+yh-yl)/2
        d2=min_boundary_distance_squared(P,centre)
        # Terminal A: the ENTIRE box is within distance R of boundary.
        # This proves NO interior disk of radius greater than R can be here.
        if H<=R and d2<=(R-H)**2:
            covered+=1
            continue
        # Terminal B: the midpoint is exterior AND the entire box stays
        # strictly farther from boundary than its radius H. Hence no
        # interior point can be reached continuously inside this box.
        if not _point_inside(P,centre) and d2>H**2:
            exterior+=1
            continue
        if depth>=max_depth:
            raise RuntimeError("global inradius proof failed closed at depth "+
                               str(depth))
        split+=1
        xmid=mx; ymid=my
        stack.extend((
            (xl,xmid,yl,ymid,depth+1),
            (xmid,xh,yl,ymid,depth+1),
            (xl,xmid,ymid,yh,depth+1),
            (xmid,xh,ymid,yh,depth+1),
        ))

    certificate=QuadtreeProof(R,visited,covered,exterior,split,deepest,
                              (xmin,xmax,ymin,ymax))
    assert certificate.complete
    assert visited==4*split+1
    return certificate


def _decimal_sqrt(q,places=18):
    """Only for display. Exact comparisons always use rational squares."""
    with localcontext() as ctx:
        ctx.prec=places+16
        return str((Decimal(q.numerator)/Decimal(q.denominator)).sqrt())


def audit():
    """Recompute two GLOBAL radius bounds and strictly separated J intervals."""
    original,_,_=read_polygon(ORIGINAL)
    improved,_,_=read_polygon(IMPROVED)
    proofs=[]
    results=[]
    for P,centre,R,expected_r2,expected_s2 in (
        (original,ORIGINAL_CENTER,ORIGINAL_RADIUS_UPPER,
         ORIGINAL_RADIUS2_LOWER,ORIGINAL_SIDE2),
        (improved,IMPROVED_CENTER,IMPROVED_RADIUS_UPPER,
         IMPROVED_RADIUS2_LOWER,IMPROVED_SIDE2),
    ):
        validate_jordan_polygon(P)
        if not _point_inside(P,centre):
            raise AssertionError("claimed near-inradius rational centre is exterior")
        actual_r2=min_boundary_distance_squared(P,centre)
        assert actual_r2 == expected_r2
        assert actual_r2<R**2
        # Same exact complete square certificate as v243 and v244/v245,
        # so a mistaken 'largest square' assumption cannot enter this audit.
        mx=maximum_square(P)
        assert mx.side_squared==expected_s2
        proof=certify_global_radius_upper(P,R)
        proofs.append(proof)
        results.append({
            "inradius_squared_lower":actual_r2,
            "inradius_radius_upper":R,
            "J_squared_lower":expected_s2/R**2,
            "J_squared_upper":expected_s2/actual_r2,
            "center":centre,
            "square_side_squared":mx.side_squared,
            "proof":proof,
        })
    original_result,improved_result=results
    # A strict EXACT ordering of the TRUE RATIOS, not merely two upper
    # bounds from arbitrarily chosen disks.
    assert improved_result["J_squared_upper"] < original_result["J_squared_lower"]
    assert original_result["J_squared_lower"] <= original_result["J_squared_upper"]
    assert improved_result["J_squared_lower"] <= improved_result["J_squared_upper"]

    def serialize(x):
        return {
            "rational_disk_center":[str(x["center"].x),str(x["center"].y)],
            "inradius_lower_squared_exact":str(x["inradius_squared_lower"]),
            "inradius_upper_exact":str(x["inradius_radius_upper"]),
            "J_squared_lower_exact":str(x["J_squared_lower"]),
            "J_squared_upper_exact":str(x["J_squared_upper"]),
            "inradius_lower_decimal":_decimal_sqrt(x["inradius_squared_lower"]),
            "inradius_upper_decimal":str(x["inradius_radius_upper"].numerator/
                                         x["inradius_radius_upper"].denominator),
            "ratio_lower_decimal":_decimal_sqrt(x["J_squared_lower"]),
            "ratio_upper_decimal":_decimal_sqrt(x["J_squared_upper"]),
            "quadtree_boxes_checked":x["proof"].visited_boxes,
            "quadtree_exterior_leaves":x["proof"].exterior_boxes,
            "quadtree_coverage_leaves":x["proof"].covered_boxes,
            "quadtree_maximum_depth":x["proof"].maximum_depth,
        }
    return {
        "theorem":"TRUE v245 Big Peg ratio strictly smaller than Pochuev's original true ratio",
        "status":"PASS: global rational-covering inradius certificates",
        "original":serialize(original_result),
        "improved":serialize(improved_result),
        "strict_j_squared_gap_lower_exact":str(
            original_result["J_squared_lower"] -
            improved_result["J_squared_upper"]),
        "independence":"Full polygon bounding-box partition; no assumptions about active medial-axis edges.",
        "scope":"Rigorous inradius enclosures, not exact symbolic radical incenter formulas or a universal Big Peg constant.",
    }


if __name__ == "__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
