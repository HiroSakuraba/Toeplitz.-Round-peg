"""v244: independent exact verification of Pochuev's original sharp Big Peg 10-gon.

Uses the project's v243 two-parameter convex-polytope square census, which
is mathematically independent of Pochuev's three published 4x4 solvers.
All decisions use fractions.Fraction. No numerical search or float comparisons.
Optional --crosscheck-v223 reruns an independent four-parameter orbit census.

Source and pinned version:
https://github.com/Grisha-Pochuev/sharp-big-peg-counterexample
commit b9e0c95cbab888405692c9efa7e592f2a0807ae4
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json
import sys

from polygon_orbit_census_v223 import pt, census_polygon, validate_jordan_polygon
from big_peg_exact_v243 import (
    read_polygon, maximum_square, contained_disk, _distance_sq_to_segment,
)

ROOT = Path(__file__).resolve().parents[1]
WITNESS = ROOT / "research" / "v244-pochuev-original-ten-gon.json"
SOURCE_COMMIT = "b9e0c95cbab888405692c9efa7e592f2a0807ae4"
SOURCE_BLOB = "0f8b70e9aab72ca2c62e2e9713049e1558ba013f"
EXPECTED_VERTICES_GRID = (
    (4114,-558), (2022,1568), (-1234,1006), (-868,2516),
    (-1089,2814), (-3769,1169), (-3459,81), (141,-1348),
    (571,-1449), (1105,-485),
)
EXPECTED_MIN_BOUNDARY_D2 = Q(2403854890969, 2400326560000)
EXPECTED_MAX_SIDE2 = Q(
    56381275521625791352241234309,
    28343589374868261752462760000
)
EXPECTED_GAP = Q(
    305903228110732152684285691,
    28343589374868261752462760000
)


def audit(crosscheck_v223=False):
    """A complete exact verification of a known polygon, with optional crosscheck.

    The two-parameter v243 verifier is separate from the upstream author's
    center/offset and vertex/side 4x4 linear-system checkers.
    """
    data = json.loads(WITNESS.read_text(encoding="utf-8"))
    provenance = data["provenance"]
    assert provenance["pinned_commit"] == SOURCE_COMMIT
    assert provenance["source_blob_sha"] == SOURCE_BLOB
    P, center, radius2 = read_polygon(WITNESS)
    assert len(P) == 10
    assert tuple((int(p.x*1200), int(p.y*1200)) for p in P) == EXPECTED_VERTICES_GRID
    assert all((p.x*1200).denominator == 1 and (p.y*1200).denominator == 1 for p in P)
    assert center == pt(0,0) and radius2 == Q(1)

    validate_jordan_polygon(P)
    assert contained_disk(P, center, radius2)
    boundary_d2 = min(
        _distance_sq_to_segment(center,P[i],P[(i+1)%len(P)])
        for i in range(len(P))
    )
    assert boundary_d2 == EXPECTED_MIN_BOUNDARY_D2
    assert boundary_d2 > radius2

    result = maximum_square(P)
    assert result.charts == 2*10**4
    assert result.positive
    assert result.side_squared == EXPECTED_MAX_SIDE2
    assert 2*radius2 - result.side_squared == EXPECTED_GAP
    assert result.side_squared < 2*radius2
    assert result.witness is not None

    claimed = data["claimed_exact_invariants"]
    assert Q(claimed["minimum_boundary_distance_squared"]) == boundary_d2
    assert Q(claimed["maximum_square_side_squared"]) == result.side_squared
    assert Q(claimed["gap_to_two"]) == EXPECTED_GAP

    orbit_count = None
    if crosscheck_v223:
        # Separate pre-existing exact four-source-variable census. This has
        # a different implementation and treats singular edge charts itself.
        census = census_polygon(P)
        assert census.count == 1, ("unexpected geometric square count", census.count)
        assert census.squares[0].side_squared == result.side_squared
        orbit_count = census.count

    return {
        "status":"PASS: independent exact Pochuev Big Peg certificate",
        "source_repo":provenance["repository"],
        "source_commit":provenance["pinned_commit"],
        "vertex_count":len(P),
        "exact_edge_charts":result.charts,
        "feasible_edge_charts":result.feasible_charts,
        "singular_feasible_charts":result.singular_feasible_charts,
        "unit_disk_certified":True,
        "minimum_boundary_distance_squared":str(boundary_d2),
        "maximum_inscribed_square_side_squared":str(result.side_squared),
        "two_minus_maximum_side_squared":str(EXPECTED_GAP),
        "sharp_big_peg_refuted_for_polygon":True,
        "maximizing_square_edge_assignment":result.edge_assignment,
        "maximizing_square_vertices":[[str(v.x),str(v.y)] for v in result.witness],
        "independent_v223_geometric_orbit_count":orbit_count,
        "limits":"Refutes sharp sqrt(2) for polygons; does not refute an unspecified universal positive Big Peg constant or Toeplitz.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosscheck-v223",action="store_true",
                        help="also run independent, more expensive geometric-square orbit census")
    args=parser.parse_args()
    result=audit(crosscheck_v223=args.crosscheck_v223)
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
