"""Four-cut parity, consecutive contacts, and an all-refinement notch certificate."""
from fractions import Fraction as Q
from itertools import product
import json

from accumulating_notches_v220 import polygon
from polygon_orbit_census_v223 import _affine_system, _rank, census_polygon, pt


def four_cut(cycles, marked, period):
    """Count mixed squares for each cut pair, including source-order wraparound."""
    i0, i1, j0, j1 = marked
    if period <= 0 or any(not 0 <= t < period for t in marked):
        raise ValueError("marks outside source circle")
    offsets = [(t-i0) % period for t in marked]
    if not 0 < offsets[1] < offsets[2] < offsets[3]:
        raise ValueError("need disjoint windows in cyclic order")

    counts = [[0, 0], [0, 0]]
    adjacency = 0
    for cycle in cycles:
        if (len(cycle) != 4 or len(set(cycle)) != 4 or set(cycle) & set(marked)
                or any(not 0 <= t < period for t in cycle)):
            raise ValueError("invalid cycle or cut at a contact")
        tags = sorted((t-i0) % period for t in cycle)
        for a, b in product(range(2), repeat=2):
            inside = sum(offsets[a] < t < offsets[2+b] for t in tags)
            counts[a][b] += int(0 < inside < 4)
        for u, v in zip(tags, tags[1:]+tags[:1]):
            a, b = sorted((u, v))
            adjacency += int(a < offsets[1] and offsets[2] < b < offsets[3])

    matrix = [[z % 2 for z in row] for row in counts]
    bit = sum(map(sum, matrix)) % 2
    if bit != adjacency % 2:
        raise ArithmeticError("four-cut identity failed")
    return {"counts": counts, "matrix": matrix, "consecutive_pairs": adjacency, "bit": bit}


def detecting_windows(squares, period):
    """Recover an odd contact pair and isolating rational parameter intervals."""
    squares = tuple(squares)
    if not squares:
        raise ValueError("need at least one geometric square")
    points = {}
    odd = set()
    observed = set()
    for square in squares:
        tags = sorted(square.curve_positions)
        points.update(zip(square.curve_positions, square.vertices))
        for a,b in zip(tags, tags[1:]+tags[:1]):
            edge = tuple(sorted((a,b)))
            odd.symmetric_difference_update([edge])
            observed.add(edge)
    if not odd:
        raise ArithmeticError("finite square collection has no exposed contact pair")

    def separation_squared(edge):
        a,b = (points[t] for t in edge)
        return (a.x-b.x)**2 + (a.y-b.y)**2

    left,right = max(sorted(odd), key=separation_squared)
    if separation_squared((left,right)) != max(map(separation_squared, observed)):
        raise ArithmeticError("longest contact length lost under cancellation")
    tags = sorted(points)
    min_gap = min((b-a) % period for a,b in zip(tags,tags[1:]+tags[:1]))
    margin = Q(min_gap)/3
    marks = tuple((t+d) % period for t in (left,right) for d in (-margin,margin))
    witness = four_cut((s.curve_positions for s in squares), marks, period)
    if witness["bit"] != 1:
        raise ArithmeticError("isolating intervals failed")
    return {"marks":list(map(str,marks)), "contact_pair":[str(left),str(right)],
            "contact_distance_squared":str(separation_squared((left,right))), "bit":1}


def notch_certificate(eta=Q(1,512)):
    """Bounds in the coordinate sup norm. All tests are strict wall exclusions."""
    eta = Q(eta)
    if eta < 0:
        raise ValueError("negative displacement")
    # Margins for the four faces after eliminating the two horizontal contacts.
    checks = {
        "first_window_endpoints": (Q(1,32), 8),
        "bottom_window_endpoints": (Q(1,8), 4),
        "upper_window_endpoints": (Q(1,8), 8),
        "last_window_endpoints": (Q(1,16), 8),
        # A diagonal completion falls below the entire perturbed curve.
        "diagonal_lower_vertex": (Q(1,16), 3),
        "clockwise_lower_vertex": (Q(5,16), 4),
        # Keep the upper right corner off either side wall and above the bottom.
        "upper_right_side_clearance": (Q(1,4), 4),
        "upper_right_bottom_clearance": (Q(3,8), 4),
        # A top contact of the last corner has reference x <= 1/2 + 12 eta.
        "last_corner_top_extent": (Q(1,8), 12),
    }
    if any(margin <= factor*eta for margin,factor in checks.values()):
        raise ValueError("arc or completion wall reached")
    return {"eta":str(eta), "side_floor":str(1-2*eta),
            "cut_bit":1, "reference_crossing":["2/5","1"],
            "remaining_margins":{key:str(gap-k*eta) for key,(gap,k) in checks.items()}}


def notched_census(n, seed=1):
    if n < 1:
        raise ValueError("the first notch contains a marked cut")
    raw = polygon(n)
    # Deterministic generic perturbation. The registered norm is <= 15/2^17.
    vertices = tuple(pt(p.x + Q((i*i+7*i+seed)%31-15, 2**17),
                        p.y + Q((3*i*i+5*i+7*seed)%31-15, 2**17))
                     for i,p in enumerate(raw))
    period = len(raw)
    downward = period-7
    bottom = period-4
    def bottom_tag(x):
        return bottom + (x-Q(9,16))/(2-Q(9,16))
    marks = (period-Q(1,8), downward+Q(1,2),
             bottom_tag(Q(7,8)), bottom_tag(Q(13,8)))
    squares = census_polygon(vertices).squares
    edges = list(zip(vertices, vertices[1:]+vertices[:1]))
    for square in squares:
        chart = [edges[int(t)] for t in square.curve_positions]
        if square.has_polygon_vertex or _rank(_affine_system(chart,1)[0]) != 4:
            raise ArithmeticError("nontransverse census witness")
    result = four_cut((s.curve_positions for s in squares), marks, period)
    result.update({"notches":n,"total_squares":len(squares),"seed":seed})
    if result["bit"] != 1:
        raise ArithmeticError("census disagrees with the uniform degree")
    return result


def report():
    triangle = (pt(0,0),pt(4,0),pt(1,3))
    cuts = tuple(t+d for t in (Q(1,7),Q(4,7)) for d in (-Q(1,100),Q(1,100)))
    triangular_squares = census_polygon(triangle).squares
    return {"triangle_automatic_windows":detecting_windows(triangular_squares,3),
            "triangle":four_cut((s.curve_positions for s in triangular_squares),
                                cuts,3),
            "uniform_notch_certificate":notch_certificate(),
            "independent_notch_censuses":[notched_census(1,k) for k in (1,4)]
                                         + [notched_census(2)]}


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
