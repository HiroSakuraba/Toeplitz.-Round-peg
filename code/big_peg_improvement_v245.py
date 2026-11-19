"""v245: certify a NEW rational 10-gon with a sharper Big Peg upper bound.

The search was numerical; every mathematical acceptance decision is made
with exact Fraction arithmetic by the pre-existing v243 polynomial edge-chart
engine and exact disk containment. Not a claim of the global optimal ratio.

The improvement certified here is over the published v244 *explicit*
source-origin disk upper bound, not necessarily over the unknown exact
inradius-optimized ratio of the original polygon.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path

from polygon_orbit_census_v223 import pt, census_polygon, validate_jordan_polygon
from big_peg_exact_v243 import (
    read_polygon, contained_disk, maximum_square, _distance_sq_to_segment,
)

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "research" / "v244-pochuev-original-ten-gon.json"
CANDIDATE = ROOT / "research" / "v245-improved-ten-gon.json"

EXPECTED_BOUNDARY_D2 = Q(1352506475516588161, 1349561502900000000)
EXPECTED_MAX_S2 = Q(
    2263154862866037715637393248853633089,
    1137974336077532830058267809922880000,
)
ORIGINAL_L2 = Q(
    56381275521625791352241234309,
    28343589374868261752462760000,
)
ORIGINAL_ORIGIN_DISK_D2 = Q(2403854890969, 2400326560000)
EXPECTED_RATIO_SQ = EXPECTED_MAX_S2 / EXPECTED_BOUNDARY_D2
BASELINE_CERTIFIED_RATIO_SQ = ORIGINAL_L2 / ORIGINAL_ORIGIN_DISK_D2


def verify(crosscheck_v223=False):
    """Completely replay the exact candidate; never trust the JSON claims."""
    source,_,_ = read_polygon(BASE)
    candidate, center, r2 = read_polygon(CANDIDATE)
    raw = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    provenance = raw["provenance"]
    assert provenance["original_author"] == "Grisha Pochuev"
    assert provenance["baseline"] == "research/v244-pochuev-original-ten-gon.json"
    assert provenance["original_commit"] == "b9e0c95cbab888405692c9efa7e592f2a0807ae4"
    assert len(source)==len(candidate)==10

    deltas = provenance["integer_deltas_from_original_times_10"]
    assert len(deltas)==10
    for original, new, pair in zip(source, candidate, deltas):
        dx, dy = pair
        assert type(dx) is int and type(dy) is int
        assert new == pt(original.x+Q(dx,12000),original.y+Q(dy,12000))
        assert (new.x*12000).denominator == (new.y*12000).denominator == 1

    validate_jordan_polygon(candidate)
    assert center == pt(Q(3,10000),Q(0))
    assert r2 == EXPECTED_BOUNDARY_D2
    assert contained_disk(candidate,center,r2)
    clearance = min(_distance_sq_to_segment(
        center,candidate[i],candidate[(i+1)%10]) for i in range(10))
    assert clearance == EXPECTED_BOUNDARY_D2
    assert clearance>1

    maximum = maximum_square(candidate)
    assert maximum.charts == 20000
    assert maximum.side_squared == EXPECTED_MAX_S2
    assert maximum.witness is not None
    assert maximum.side_squared < 2*r2
    actual_verified_ratio_sq = maximum.side_squared / r2
    assert actual_verified_ratio_sq == EXPECTED_RATIO_SQ
    assert actual_verified_ratio_sq < BASELINE_CERTIFIED_RATIO_SQ

    claimed = raw["claimed_exact_invariants"]
    assert Q(claimed["minimum_boundary_distance_squared"]) == clearance
    assert Q(claimed["maximum_square_side_squared"]) == maximum.side_squared
    assert Q(claimed["candidate_squared_ratio_bound"]) == actual_verified_ratio_sq
    assert Q(claimed["original_source_origin_disk_squared_ratio_bound"]) == BASELINE_CERTIFIED_RATIO_SQ

    geometric_orbits = None
    if crosscheck_v223:
        # Independent old four-parameter census, includes degenerate charts.
        result = census_polygon(candidate)
        assert result.count == 1
        assert result.squares[0].side_squared == maximum.side_squared
        geometric_orbits = result.count

    return {
        "status":"PASS: independent exact rational improved candidate",
        "vertex_count":len(candidate),
        "edge_chart_cases":maximum.charts,
        "feasible_charts":maximum.feasible_charts,
        "singular_feasible_charts":maximum.singular_feasible_charts,
        "exact_disk_center":[str(center.x),str(center.y)],
        "certified_disk_radius_squared":str(r2),
        "exact_maximum_square_side_squared":str(maximum.side_squared),
        "candidate_squared_ratio_upper_bound":str(actual_verified_ratio_sq),
        "original_source_origin_disk_squared_ratio_upper_bound":str(BASELINE_CERTIFIED_RATIO_SQ),
        "strict_upper_bound_improved":True,
        "sharp_sqrt2_big_peg_counterexample":True,
        "maximizing_square_edge_assignment":list(maximum.edge_assignment),
        "maximizing_square_vertices":[[str(p.x),str(p.y)] for p in maximum.witness],
        "v223_geometric_orbit_count":geometric_orbits,
        "scope":"Improved certified upper bound relative to source origin disk; does not prove a globally minimal ratio or resolve positive-c Big Peg / Toeplitz.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosscheck-v223",action="store_true")
    args=parser.parse_args()
    print(json.dumps(verify(args.crosscheck_v223),indent=2,sort_keys=True))


if __name__=="__main__":
    main()
