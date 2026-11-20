"""Exact rectangle atlas for rational Jordan polygons with nonparallel edges.

The chart parameter is tan(theta/2)/(1+tan(theta/2)); CCW labels are
kept until gluing. Unresolved positive-dimensional fibers fail closed.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path

from sage.all import AA, QQ, PolynomialRing, RealIntervalField, matrix, vector

from moving_polygon_events_v224 import count_roots_open, sturm_sequence
from polygon_orbit_census_v223 import census_polygon, pt, validate_jordan_polygon


R = PolynomialRing(QQ, "x")
X = R.gen()
K = R.fraction_field()
P = ((0, 0), (4, 0), (1, 1), (0, 3))
ARTIFACT = Path(__file__).resolve().parents[1] / "results" / "v274" / "concave-quad-atlas.json"


class AtlasInconclusive(ValueError):
    pass


@dataclass(frozen=True)
class Chart:
    edges: tuple
    system: object
    rhs: object
    det: object
    parameters: tuple
    events: tuple
    singular: tuple


def _cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def _roots(p):
    if p.is_zero() or p.degree() <= 0:
        return ()
    return tuple(t for t in p.roots(AA, multiplicities=False) if 0 < t < 1)


def _residual(points):
    a, b, c, d = points
    midpoint = a+c-b-d
    u, v = a-c, b-d
    cos_num = 1-2*X
    sin_num = 2*X*(1-X)
    denominator = (1-X)**2+X**2
    rotated = vector(R, (cos_num*u[0]-sin_num*u[1], sin_num*u[0]+cos_num*u[1]))
    return vector(R, tuple(midpoint)+tuple(denominator*v-rotated))


def chart(points, tangents, ids):
    origins = [vector(R, points[i]) for i in ids]
    constant = _residual(origins)
    columns = []
    for j, i in enumerate(ids):
        displaced = origins.copy()
        displaced[j] = displaced[j] + vector(R, tangents[i])
        columns.append(_residual(displaced)-constant)
    system = matrix(R, columns).transpose()
    determinant = system.det()
    if not determinant:
        raise AtlasInconclusive(f"persistently singular chart {ids}")
    rhs = -constant
    parameters = tuple(system.change_ring(K).solve_right(rhs.change_ring(K)))
    singular = _roots(determinant)
    events = set(singular)
    for value in parameters:
        events.update(_roots(value.numerator()))
        events.update(_roots((value-1).numerator()))
    return Chart(ids, system, rhs, determinant, parameters, tuple(sorted(events)), singular)


def _inner_sample(a, b):
    bits = 32
    while bits <= 4096:
        field = RealIntervalField(bits)
        lo = field(a).upper().exact_rational()
        hi = field(b).lower().exact_rational()
        if lo < hi:
            return (lo+hi)/2
        bits *= 2
    raise AtlasInconclusive("could not separate distinct angle events")


def _at(f, x):
    return f.numerator()(x)/f.denominator()(x)


def _evaluate(points, tangents, ids, t):
    return tuple(tuple(a+x*b for a, b in zip(points[i], tangents[i])) for i, x in zip(ids, t))


def _radius_sq(corners):
    a, _, c, _ = corners
    return sum((x-y)**2 for x, y in zip(a, c))/4


def _positions(ids, params, n):
    return tuple(AA(0) if i+t == n else AA(i+t) for i, t in zip(ids, params))


def _polytope_vertices(A, b):
    """All vertices of an affine slice of [0,1]^4 over the real algebraics."""
    independent = tuple(A.transpose().pivots())
    base = [tuple(A.row(i)) for i in independent]
    rhs = [b[i] for i in independent]
    dimension = 4-len(base)
    found = set()
    for walls in combinations(range(4), dimension):
        system = matrix(AA, base+[tuple(AA(i == j) for i in range(4)) for j in walls])
        if not system.det():
            continue
        inverse = system.inverse()
        for fixed in product((0, 1), repeat=dimension):
            u = inverse*vector(AA, rhs+list(fixed))
            if all(0 <= v <= 1 for v in u) and A*u == b:
                found.add(tuple(u))
    return tuple(sorted(found))


def _matrix_at(A, b, x):
    return matrix(AA, [[p(x) for p in row] for row in A.rows()]), vector(AA, [p(x) for p in b])


def _candidates(n):
    for ids in product(range(n), repeat=4):
        # At least three edges are needed when none of their lines is parallel.
        if len(set(ids)) >= 3 and max(ids.count(i) for i in ids) <= 2:
            yield ids


def rotate_node(key):
    angle, positions = key
    return (1-angle, positions[1:]+positions[:1])


def _kind(angle, points):
    if not _radius_sq(points):
        return "total_collision"
    if angle in (0, 1):
        return "two_pair_flat"
    return "square" if angle == QQ(1)/2 else "genuine"


def build_atlas(vertices=P):
    raw = validate_jordan_polygon(tuple(pt(*p) for p in vertices))
    points = tuple((QQ(p.x), QQ(p.y)) for p in raw)
    n = len(points)
    tangents = tuple(tuple(b-a for a, b in zip(points[i], points[(i+1)%n])) for i in range(n))
    if any(not _cross(a, b) for a, b in combinations(tangents, 2)):
        raise AtlasInconclusive("parallel edge lines are outside this atlas")
    charts = tuple(chart(points, tangents, ids) for ids in _candidates(n))
    radii = {c.edges: _radius_sq(_evaluate(points, tangents, c.edges, c.parameters)) for c in charts}
    events = ({AA(0), AA(QQ(1)/2), AA(1)} | {x for c in charts for x in c.events}
              | {x for radius in radii.values() for x in _roots(radius.numerator())})
    angles = sorted(events | {1-x for x in events})
    chambers = [(a, b, _inner_sample(a, b)) for a, b in zip(angles, angles[1:])]
    nodes = {}
    edges = []
    singular_counts = {"empty": 0, "collision_only": 0, "isolated_positive": 0}

    def node(c, angle, parameters):
        position = _positions(c.edges, parameters, n)
        corners = _evaluate(points, tangents, c.edges, parameters)
        key = (angle, position)
        prior = nodes.get(key)
        if prior is None:
            nodes[key] = (corners, _kind(angle, corners))
        elif prior[0] != corners:
            raise AssertionError("same curve parameters give different geometry")
        return key

    for c in charts:
        ids = c.edges
        # Discard identically collapsed charts, using affine diagonal numerators.
        radius = radii[ids]
        if radius:
            critical = None
            for a, b, t in chambers:
                params = tuple(_at(f, t) for f in c.parameters)
                if not all(0 <= u < 1 for u in params) or _at(radius, t) <= 0:
                    continue
                ends = tuple(tuple(_at(f, s) for f in c.parameters) for s in (a, b))
                if any(not 0 <= u <= 1 for us in ends for u in us):
                    raise AssertionError("feasible chamber has infeasible closure")
                begin, end = (node(c, s, us) for s, us in zip((a, b), ends))
                if critical is None:
                    critical = _roots(radius.derivative().numerator())
                floor = min(_at(radius, s) for s in (a, b, *(x for x in critical if a < x < b)))
                pos = _positions(ids, params, n)
                sector = (0,)+tuple(sorted(range(1, 4),
                                           key=lambda i: pos[i]-pos[0]+(n if pos[i] < pos[0] else 0)))
                edges.append((begin, end, ids, floor, sector))
        # A singular angle can contain extra configurations missed by Cramer branches.
        for angle in c.singular:
            A, rhs = _matrix_at(c.system, c.rhs, angle)
            feasible = _polytope_vertices(A, rhs)
            positive = [u for u in feasible if _radius_sq(_evaluate(points, tangents, ids, u))]
            if not feasible:
                singular_counts["empty"] += 1
            elif not positive:
                singular_counts["collision_only"] += 1
            elif len(feasible) == 1:
                singular_counts["isolated_positive"] += 1
                node(c, angle, positive[0])
            else:
                raise AtlasInconclusive(f"positive-dimensional singular fiber {ids}, angle {angle}")
        # Isolated feasible regular events can occur when box walls only touch.
        for angle in c.events if radius else ():
            if angle in c.singular:
                continue
            params = tuple(_at(f, angle) for f in c.parameters)
            if all(0 <= u < 1 for u in params) and _at(radius, angle) > 0:
                node(c, angle, params)

    # Rotating the labels complements x. Check edges and nodes before quotienting.
    for key in nodes:
        if rotate_node(key) not in nodes:
            raise AtlasInconclusive("missing cyclic rotation of a chart node")
    edge_keys = {(a, b): ids for a, b, ids, _, _ in edges}
    if len(edge_keys) != len(edges):
        raise AtlasInconclusive("multiple branches share graph endpoints")
    if any((rotate_node(b), rotate_node(a)) not in edge_keys for a, b in edge_keys):
        raise AtlasInconclusive("missing cyclic rotation of a chart edge")
    by_ids = {c.edges: c.parameters for c in charts}
    for ids, shifted in {(ids, edge_keys[rotate_node(b), rotate_node(a)])
                         for (a, b), ids in edge_keys.items()}:
        if shifted != ids[1:]+ids[:1]:
            raise AtlasInconclusive("rotated branch has the wrong source ownership")
        for f, g in zip(by_ids[ids][1:]+by_ids[ids][:1], by_ids[shifted]):
            if f.numerator()*g.denominator()(1-X) != f.denominator()*g.numerator()(1-X):
                raise AtlasInconclusive("cyclic transition changed the interior geometry")
    squares = {tuple(sorted(pts)) for pts, kind in nodes.values() if kind == "square"}
    census = census_polygon(raw)
    independent = {tuple(sorted((AA(QQ(p.x)), AA(QQ(p.y))) for p in s.vertices))
                   for s in census.squares}
    if squares != independent:
        raise AtlasInconclusive("independent square census disagrees")
    return {"points": points, "charts": charts, "angles": angles, "nodes": nodes,
            "edges": edges, "singular_counts": singular_counts, "squares": squares}


def _tangents(points):
    return tuple(tuple(b-a for a, b in zip(points[i], points[(i+1)%len(points)]))
                 for i in range(len(points)))


def atlas_slice(data, x):
    """Geometric rectangle orbits at one nonflat angle, including isolated nodes."""
    x = AA(QQ(x)) if isinstance(x, (int,)) else AA(x)
    if not 0 < x < 1:
        raise ValueError("slice requires a genuine angle")
    result = {tuple(sorted(pts)) for (angle, _), (pts, kind) in data["nodes"].items()
              if angle == x and kind != "total_collision"}
    if x not in data["angles"]:
        by_ids = {c.edges: c for c in data["charts"]}
        points = data["points"]
        tangents = _tangents(points)
        for a, b, ids, _, _ in data["edges"]:
            if a[0] < x < b[0]:
                corners = _evaluate(points, tangents, ids, tuple(_at(f, x) for f in by_ids[ids].parameters))
                result.add(tuple(sorted(corners)))
    return result


def independent_slice(points, x):
    """Direct adjacent-side ratio census, without chart polynomials or gluing."""
    x = AA(x)
    if not 0 < x < 1:
        raise ValueError("slice requires a genuine angle")
    ratio = (1-x)/x
    tangents = _tangents(points)

    def residual(pts):
        a, b, c, d = pts
        side = b-a
        rotated = vector(AA, (-ratio*side[1], ratio*side[0]))
        return vector(AA, tuple(c-b-rotated)+tuple(d-a-rotated))

    results = set()
    for ids in product(range(len(points)), repeat=4):
        if max(ids.count(i) for i in ids) > 2:
            continue
        origins = [vector(AA, points[i]) for i in ids]
        constant = residual(origins)
        cols = []
        for j, i in enumerate(ids):
            displaced = origins.copy()
            displaced[j] = displaced[j]+vector(AA, tangents[i])
            cols.append(residual(displaced)-constant)
        system = matrix(AA, cols).transpose()
        feasible = _polytope_vertices(system, -constant)
        positive = [pts for u in feasible
                    if _radius_sq(pts := _evaluate(points, tangents, ids, u))]
        if positive and len(feasible) != 1:
            raise AtlasInconclusive("independent slice is positive-dimensional")
        results.update(tuple(sorted(pts)) for pts in positive)
    return results


@lru_cache(maxsize=None)
def algebraic_data(value):
    a = AA(value)
    p = tuple(QQ(x) for x in a.minpoly().list())
    if len(p) == 2:
        return str(-p[0]/p[1])
    sturm = sturm_sequence(tuple(Fraction(str(x)) for x in p))
    for bits in (64, 128, 256, 512, 1024):
        scale = 2**bits
        cell = (a*scale).floor()
        lo, hi = QQ(cell)/scale, QQ(cell+1)/scale
        if count_roots_open(sturm, Fraction(str(lo)), Fraction(str(hi))) == 1:
            return {"minpoly": list(map(str, p)), "isolator": [str(lo), str(hi)]}
    raise AtlasInconclusive("could not serialize isolated real algebraic value")


def artifact(data):
    keys = sorted(data["nodes"])
    index = {key: i for i, key in enumerate(keys)}
    return {
        "polygon": [[str(x) for x in p] for p in data["points"]],
        "ordered_charts": len(data["charts"]),
        "angle_events": [algebraic_data(x) for x in data["angles"]],
        "singular_fibers": data["singular_counts"],
        "rotation": [index[rotate_node(key)] for key in keys],
        "nodes": [{"angle": algebraic_data(x), "parameters": [algebraic_data(t) for t in pos],
                   "corners": [[algebraic_data(v) for v in p] for p in data["nodes"][x, pos][0]],
                   "dot_sign": (0 if data["nodes"][x, pos][1] == "total_collision"
                                else int(x < QQ(1)/2)-int(x > QQ(1)/2)),
                   "kind": data["nodes"][x, pos][1]} for x, pos in keys],
        "branches": [{"ends": [index[a], index[b]], "source_edges": list(ids),
                      "radius_squared_infimum": algebraic_data(r), "source_sector": list(s)}
                     for a, b, ids, r, s in data["edges"]],
        "squares": [[[algebraic_data(x) for x in p] for p in square]
                    for square in sorted(data["squares"])],
    }


def artifact_text(result):
    """Keep generated records on individual lines for review."""
    fields = []
    for key, value in result.items():
        if isinstance(value, list):
            rows = ",\n".join("    "+json.dumps(row, separators=(",", ":")) for row in value)
            encoded = "[\n"+rows+"\n  ]"
        else:
            encoded = json.dumps(value, separators=(",", ":"))
        fields.append(f"  {json.dumps(key)}: {encoded}")
    return "{\n"+",\n".join(fields)+"\n}\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = build_atlas()
    print(f"Glued {len(data['nodes'])} labeled nodes and {len(data['edges'])} branches.",
          flush=True)
    result = artifact(data)
    if args.write:
        ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
        ARTIFACT.write_text(artifact_text(result), encoding="utf-8")
    elif result != json.loads(ARTIFACT.read_text(encoding="utf-8")):
        raise AtlasInconclusive("stored atlas differs from exact replay")
    print(f"{result['ordered_charts']} charts, {len(result['nodes'])} labeled nodes, "
          f"{len(result['branches'])} labeled branches, {len(result['squares'])} geometric squares.")


if __name__ == "__main__":
    main()
