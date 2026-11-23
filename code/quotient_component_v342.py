"""Cyclic component parity via an exact, symmetry-preserving polygon deformation."""
from dataclasses import dataclass
from fractions import Fraction as Q

from polygon_component_v341 import PolygonComponents, segment_floor
from polygon_orbit_census_v223 import _affine_system, pt


def rotations(ids):
    return tuple(ids[j:]+ids[:j] for j in range(4))


def component_stabilizer(cells):
    keys = frozenset(cells)
    return tuple(j for j in range(4)
                 if frozenset(ids[j:]+ids[:j] for ids in keys) == keys)


@dataclass(frozen=True)
class OrbitReplay:
    stabilizer: tuple
    marked_charts: int
    chart_orbits: int
    regular_charts: tuple
    floor_squared: Q

    @property
    def parity(self):
        return len(self.regular_charts) % 2


def orbit_replay(vertices, seed, velocities):
    """Compute the mod-two local index with P(t)=P+t*velocities at t=0+.

    The input direction must regularize every surviving orbit. Reject ambiguous
    chart hits, since counting shared faces or families as distinct gives false parity.
    """
    import sympy as sp
    from terrain_component_degree_v337 import feasible_vertices, sign_near_zero

    engine = PolygonComponents(vertices)
    if len(velocities) != len(engine.edges):
        raise ValueError("expected one velocity per polygon vertex")
    velocity = [pt(*v) for v in velocities]
    edges_v = list(zip(velocity,velocity[1:]+velocity[:1]))
    cells = engine.trace(tuple(pt(*p) for p in seed))
    floor = min(segment_floor(c.squares) for c in cells.values())
    # One tuple per chart orbit in the saturation, also when rotations
    # take the chosen labeled component to a disjoint component.
    representatives = sorted({min(rotations(ids)) for ids in cells})
    t = sp.Symbol("t", positive=True)
    hits = []
    for ids in representatives:
        A,b = _affine_system([engine.edges[i] for i in ids],1)
        dA,db = _affine_system([edges_v[i] for i in ids],1)
        matrix = sp.Matrix(A)+t*sp.Matrix(dA)
        rhs = sp.Matrix(b)+t*sp.Matrix(db)
        feasible = feasible_vertices(matrix,rhs,t)
        if not feasible:
            continue
        if len(feasible) != 1 or sp.cancel(matrix.det()) == 0:
            raise ValueError(f"deformation leaves a nonregular chart {ids}")
        u, = feasible
        if any(min(sign_near_zero(x,t),sign_near_zero(1-x,t)) <= 0 for x in u):
            raise ValueError(f"deformation leaves a boundary zero on {ids}")
        if len(set(rotations(ids))) != 4:
            raise ValueError("periodic chart tuple at a regular square")
        hits.append(ids)
    return OrbitReplay(component_stabilizer(cells),len(cells),len(representatives),
                       tuple(hits),floor)


def recorded_stabilizers(fixtures, records):
    """Independently trace the recorded v341 seeds without repeating its search."""
    from collections import Counter

    sizes = Counter()
    for kind, n, _, _, _, i, boundary in fixtures:
        row = records[kind,n]
        seed = tuple(pt(*p.split(":")) for p in row["marked_square"].split())
        cells = PolygonComponents(boundary).trace(seed)
        floor = min(segment_floor(c.squares) for c in cells.values())
        if floor != Q(row["floor2"]) or 2*floor < i*i or int(row["degree"]) % 2 != 1:
            raise ValueError(f"unexpected v341 certificate on {kind} {n}")
        sizes[kind,len(component_stabilizer(cells))] += 1
    return sizes


if __name__ == "__main__":
    import csv
    from pathlib import Path
    from terrain_replay_v341 import fixtures

    polygon = ((0,0),(2,0),(2,2),(0,2))
    print("square:", orbit_replay(polygon,polygon,[(y,x) for x,y in polygon]))
    data = Path(__file__).resolve().parents[1]/"results/v341/component-certificates.csv"
    with data.open() as f:
        records = {(r["kind"],int(r["seed"])):r for r in csv.DictReader(f)}
    print("v341 component stabilizers:", dict(recorded_stabilizers(fixtures(),records)))
