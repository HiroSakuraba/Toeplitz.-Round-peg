"""Independently enumerate actual large squares after breaking polygon degeneracies."""
import argparse
from fractions import Fraction as Q
import json

from deep_dents_v308 import certify, example, pt
from polygon_orbit_census_v223 import _affine_system, _rank, census_polygon


def census_example(name):
    polygon, lanes, low, high, cap = example(name)
    eta = Q(1,1024)
    points = tuple(pt(p.x+Q((i*i+3*i+1)%19-9,32768),
                      p.y+Q((2*i*i+i+4)%23-11,32768)) for i,p in enumerate(polygon))
    if any((a.x-b.x)**2+(a.y-b.y)**2>eta**2 for a,b in zip(polygon,points)):
        raise ArithmeticError("sample exceeds its registered perturbation budget")
    bands = [(a-eta,b+eta) for a,b in lanes]
    record = certify(points,bands,low-eta,high+eta,cap+2*eta)
    # q=1 lies strictly in the certified gap. The v223 cutoff uses squared SIDE length.
    if not Q(record["forbidden_half_diagonal_squared_lower_open"])<1<Q(record["forbidden_half_diagonal_upper_open"]):
        raise ArithmeticError("sample cutoff left the gap")
    result = census_polygon(points,cutoff_squared=Q(2))
    for square in result.squares:
        ids = [int(t) for t in square.curve_positions]
        edges = [(points[i],points[(i+1)%len(points)]) for i in ids]
        matrix, _ = _affine_system(edges,1)
        if square.has_polygon_vertex or _rank(matrix)!=4 or square.labelings_checked!=8:
            raise ArithmeticError("expected eight regular interior labelings per square")
    if result.parity!=record["parity"]:
        raise ArithmeticError("independent orbit count disagrees with cell parity")
    return {"fixture":name,"actual_squares_above_cutoff":result.count,
            "predicted_parity":record["parity"],"edge_charts":result.checked_edge_charts,
            "nongraceful":result.nongraceful_count}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixture",choices=("opposite_dents","two_runs","three_runs"))
    args = parser.parse_args()
    print(json.dumps(census_example(args.fixture),indent=2))
