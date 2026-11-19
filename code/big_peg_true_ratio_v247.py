"""v247: exact true Big Peg ratio improvement for a rational 10-gon.

All theorem decisions use fractions.Fraction and the existing independently
audited v243 exhaustive square maximum + v246 global radius covering engine.
The local numerical search only supplied candidate rational vertex coordinates.

The exact strict comparison establishes:
    J(v247) < J(v245)
where J(P) = largest inscribed square side / true polygon inradius.
It does not assert optimality among polygons or resolve the general problem.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json

from polygon_orbit_census_v223 import pt, validate_jordan_polygon, census_polygon
from big_peg_exact_v243 import read_polygon, contained_disk, maximum_square
from big_peg_inradius_v246 import (
    certify_global_radius_upper, min_boundary_distance_squared,
)

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "research" / "v245-improved-ten-gon.json"
NEW_PATH = ROOT / "research" / "v247-refined-ten-gon.json"

GRID_DENOMINATOR = 1200000
NEW_DISK_CENTER = pt(Q(114334619, 100000000000),
                     Q(-7507522, 100000000000))
NEW_DISK_RADIUS2 = Q(
    984883927824703700502571706090888809,
    982265423769090000000000000000000000)
NEW_GLOBAL_RADIUS_UPPER = Q(100133200316, 100000000000)
NEW_MAX_SQUARE_SIDE2 = Q(
    407200306394224630512330297939647665150234737979752892405301,
    204751870350371073185004616046232753607402236547022400000000)

OLD_MAX_SQUARE_SIDE2 = Q(
    2263154862866037715637393248853633089,
    1137974336077532830058267809922880000)
OLD_GLOBAL_RADIUS_UPPER = Q(100109166, 100000000)
EXPECTED_NEW_COVER_BOXES = 1257
EXPECTED_OLD_COVER_BOXES = 917


def _sqrt_display(q, prec=30):
    """Decimal approximation for presentation only, never for decisions."""
    with localcontext() as ctx:
        ctx.prec=prec
        return str((Decimal(q.numerator)/Decimal(q.denominator)).sqrt())


def audit(crosscheck_v223=False):
    """Fail closed unless every exact polygonal certificate is recomputed."""
    base,_,_ = read_polygon(BASE_PATH)
    new,center,r2 = read_polygon(NEW_PATH)
    raw = json.loads(NEW_PATH.read_text(encoding="utf-8"))
    provenance = raw["provenance"]
    if provenance["derived_from"] != "research/v245-improved-ten-gon.json":
        raise AssertionError("changed source provenance")
    if provenance["coordinate_denominator"] != GRID_DENOMINATOR:
        raise AssertionError("changed rational grid denominator")
    if len(base)!=len(new) or len(new)!=10:
        raise AssertionError("expected two simple ten-gons")
    deltas=provenance["vertex_deltas_from_v245_in_grid_units"]
    if len(deltas)!=10:
        raise AssertionError("wrong number of grid displacements")
    for oldp,newp,pair in zip(base,new,deltas):
        if not (len(pair)==2 and all(type(v) is int for v in pair)):
            raise AssertionError("grid displacements must be exact integers")
        dx,dy=pair
        if newp != pt(oldp.x+Q(dx,GRID_DENOMINATOR),
                      oldp.y+Q(dy,GRID_DENOMINATOR)):
            raise AssertionError("candidate differs from recorded grid moves")
        if any((coord*GRID_DENOMINATOR).denominator != 1 for coord in
               (newp.x,newp.y)):
            raise AssertionError("off-grid candidate coordinate")

    validate_jordan_polygon(base)
    validate_jordan_polygon(new)
    if center != NEW_DISK_CENTER or r2 != NEW_DISK_RADIUS2:
        raise AssertionError("altered exact disk certificate")
    if not contained_disk(new,center,r2):
        raise AssertionError("new rational disk is not contained")
    clearance = min_boundary_distance_squared(new,center)
    if clearance != NEW_DISK_RADIUS2:
        raise AssertionError("claimed squared disk clearance is incorrect")
    if r2<=0 or r2 >= NEW_GLOBAL_RADIUS_UPPER**2:
        raise AssertionError("bad inradius interval")

    # The square maximum is checked over all 2*n^4 ordered edge/orientation
    # charts, including positive-dimensional feasible parameter sets.
    new_max = maximum_square(new)
    old_max = maximum_square(base)
    if new_max.charts != 20000 or old_max.charts != 20000:
        raise AssertionError("incomplete exact square chart enumeration")
    if new_max.side_squared != NEW_MAX_SQUARE_SIDE2:
        raise AssertionError("candidate maximum square differs from claim")
    if old_max.side_squared != OLD_MAX_SQUARE_SIDE2:
        raise AssertionError("baseline maximum square differs from claim")
    if not new_max.positive or new_max.witness is None:
        raise AssertionError("no attained nondegenerate maximum square")

    # Whole-polygon, exact rational rectangle covers prove upper bounds
    # on the GLOBAL inradius. No choice of active medial-axis edges assumed.
    new_cover = certify_global_radius_upper(new,NEW_GLOBAL_RADIUS_UPPER)
    old_cover = certify_global_radius_upper(base,OLD_GLOBAL_RADIUS_UPPER)
    if not (new_cover.complete and old_cover.complete):
        raise AssertionError("incomplete global inradius cover")
    if new_cover.visited_boxes != EXPECTED_NEW_COVER_BOXES:
        raise AssertionError("unexpected new polygon cover statistics")
    if old_cover.visited_boxes != EXPECTED_OLD_COVER_BOXES:
        raise AssertionError("unexpected v245 cover statistics")

    # NEW J^2 <= s_new^2/r_new_lower^2,
    # OLD J^2 >= s_old^2/r_old_global_upper^2.
    # Their strictly disjoint rational bounds prove a true improvement.
    new_j2_lower = NEW_MAX_SQUARE_SIDE2 / NEW_GLOBAL_RADIUS_UPPER**2
    new_j2_upper = NEW_MAX_SQUARE_SIDE2 / r2
    old_j2_lower = OLD_MAX_SQUARE_SIDE2 / OLD_GLOBAL_RADIUS_UPPER**2
    if not (new_j2_lower <= new_j2_upper < old_j2_lower):
        raise AssertionError("claimed strict TRUE ratio improvement not proved")

    orbit_count=None
    if crosscheck_v223:
        # This 4-parameter exact geometry census predates v247 and is
        # independent of the convex 2-parameter maximum-only implementation.
        census = census_polygon(new)
        if census.count != 1 or census.squares[0].side_squared != NEW_MAX_SQUARE_SIDE2:
            raise AssertionError("independent v223 geometric-square audit failed")
        orbit_count=census.count

    invariants=raw["claimed_exact_invariants"]
    for key,truth in (
        ("maximum_square_side_squared",NEW_MAX_SQUARE_SIDE2),
        ("minimum_distance_squared_at_disk_center",clearance),
        ("global_radius_cover_boxes",EXPECTED_NEW_COVER_BOXES),
        ("previous_v245_global_radius_upper",OLD_GLOBAL_RADIUS_UPPER),
        ("previous_v245_maximum_square_side_squared",OLD_MAX_SQUARE_SIDE2),
    ):
        if Q(invariants[key]) != Q(truth):
            raise AssertionError("discrepant JSON claimed invariant: "+key)

    return {
        "status":"PASS: exact true Big Peg ratio improved (v247 < v245)",
        "vertices":10,
        "grid_denominator":GRID_DENOMINATOR,
        "new_max_square_squared":str(NEW_MAX_SQUARE_SIDE2),
        "new_contained_disk_center":[str(center.x),str(center.y)],
        "new_contained_disk_squared_radius":str(r2),
        "new_global_inradius_upper":str(NEW_GLOBAL_RADIUS_UPPER),
        "new_global_cover_boxes":new_cover.visited_boxes,
        "new_global_cover_max_depth":new_cover.maximum_depth,
        "new_global_cover_exterior_leaves":new_cover.exterior_boxes,
        "new_global_cover_boundary_leaves":new_cover.covered_boxes,
        "previous_global_cover_boxes":old_cover.visited_boxes,
        "new_j_squared_lower":str(new_j2_lower),
        "new_j_squared_upper":str(new_j2_upper),
        "previous_j_squared_lower":str(old_j2_lower),
        "new_true_J_interval_display":[_sqrt_display(new_j2_lower),
                                        _sqrt_display(new_j2_upper)],
        "previous_true_J_lower_display":_sqrt_display(old_j2_lower),
        "strict_squared_ratio_separation_lower":str(old_j2_lower-new_j2_upper),
        "new_geometric_square_orbits":orbit_count,
        "new_maximizing_edge_assignment":list(new_max.edge_assignment),
        "new_maximizing_square_vertices":[[str(p.x),str(p.y)] for p in new_max.witness],
        "scope":"Strict improvement over v245's true J, not ten-gon optimality or a general positive-c Big Peg result; Toeplitz remains open.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosscheck-v223",action="store_true",
                        help="verify unique geometric square using independent old census")
    args=parser.parse_args()
    print(json.dumps(audit(args.crosscheck_v223),indent=2,sort_keys=True))


if __name__=="__main__":
    main()
