"""Relative sign line of perpendicularity on a verified compact graph."""

from collections import deque
from fractions import Fraction
import json
from pathlib import Path

from polygon_orbit_census_v223 import _affine_system, _evaluate, _side_sq, _solve, pt


ARTIFACT = Path(__file__).resolve().parents[1] / "results" / "v274" / "concave-quad-atlas.json"


def relative_obstructions(vertices, edges, boundary):
    """Solve δs = twist with the prescribed boundary values; return witnesses."""
    adjacent = {v: [] for v in vertices}
    seen_edges = set()
    for a, b, twist in edges:
        if twist not in (0, 1):
            raise ValueError("transitions must be binary")
        key = frozenset((a, b))
        if key in seen_edges:
            raise ValueError("parallel edges need distinct identifiers")
        seen_edges.add(key)
        adjacent[a].append((b, twist))
        adjacent[b].append((a, twist))
    if any(v not in adjacent or bit not in (0, 1) for v, bit in boundary.items()):
        raise ValueError("invalid relative boundary")
    potential = {}
    parent = {}

    def path(a, b):
        chains = []
        for v in (a, b):
            chain = []
            while v is not None:
                chain.append(v)
                v = parent[v]
            chains.append(chain)
        first, second = chains
        join = next(v for v in first if v in second)
        return first[:first.index(join)+1] + list(reversed(second[:second.index(join)]))

    components = []
    for root in sorted(adjacent):
        if root in potential:
            continue
        potential[root] = 0
        parent[root] = None
        todo = deque([root])
        component = []
        witness = None
        while todo:
            a = todo.popleft()
            component.append(a)
            for b, twist in adjacent[a]:
                want = potential[a]^twist
                if b not in potential:
                    potential[b], parent[b] = want, a
                    todo.append(b)
                elif potential[b] != want and witness is None:
                    witness = {"kind": "odd_monodromy", "path": path(a, b)+[a]}
        fixed = [v for v in sorted(component) if v in boundary]
        if witness is None and fixed:
            first = fixed[0]
            opposite = next((v for v in fixed[1:]
                             if potential[v]^boundary[v] != potential[first]^boundary[first]), None)
            if opposite is not None:
                witness = {"kind": "relative_boundary", "path": path(first, opposite)}
        components.append({"vertices": sorted(component), "boundary": fixed, "witness": witness})
    return components


def _lower(value):
    return Fraction(value if isinstance(value, str) else value["isolator"][0])


def quotient_index(atlas):
    """Cyclic quotient of positive-radius chart edges, with flat boundary frames."""
    nodes, rot = atlas["nodes"], atlas["rotation"]
    if sorted(rot) != list(range(len(nodes))):
        raise ValueError("rotation is not a permutation")
    positive = {i for i, node in enumerate(nodes) if node["kind"] != "total_collision"}

    def orbit(i):
        result = [i]
        for _ in range(3):
            result.append(rot[result[-1]])
        return result

    by_vertex = {}
    for i in positive:
        lifts = orbit(i)
        if (len(set(lifts)) != 4 or rot[lifts[-1]] != i
                or set(lifts)-positive or nodes[i]["dot_sign"] != -nodes[rot[i]]["dot_sign"]
                or any(nodes[j]["kind"] != nodes[i]["kind"] for j in lifts)):
            raise ValueError("positive node lacks a free compatible rotation orbit")
        first = min(lifts)
        by_vertex[i] = (first, lifts.index(first)%2)
    vertices = {item[0] for item in by_vertex.values()}
    boundary = {i: int(nodes[i]["dot_sign"] < 0) for i in vertices
                if nodes[i]["kind"] == "two_pair_flat"}
    if any(nodes[i]["dot_sign"] == 0 for i in boundary):
        raise ValueError("flat boundary sign vanishes")
    lifted = {}
    for branch in atlas["branches"]:
        a, b = branch["ends"]
        key = tuple(sorted((a, b)))
        if key in lifted:
            raise ValueError("ambiguous labeled branch")
        if a in positive and b in positive:
            lifted[key] = _lower(branch["radius_squared_infimum"])
            if lifted[key] <= 0:
                raise ValueError("positive branch lacks a certified radius floor")
    seen = set()
    edges = []
    floors = {}
    degree = dict.fromkeys(vertices, 0)
    for a, b in sorted(lifted):
        if (a, b) in seen:
            continue
        lifts = {tuple(sorted((i, j))) for i, j in zip(orbit(a), orbit(b))}
        if len(lifts) != 4 or lifts-set(lifted):
            raise ValueError("positive branch lacks its cyclic orbit")
        seen.update(lifts)
        u, p = by_vertex[a]
        v, q = by_vertex[b]
        edge = tuple(sorted((u, v)))
        if edge in floors:
            raise ValueError("parallel quotient edges need distinct identifiers")
        floors[edge] = min(lifted[e] for e in lifts)
        edges.append((u, v, p^q))
        degree[u] += 1
        degree[v] += 1
    result = relative_obstructions(vertices, edges, boundary)
    for c in result:
        if witness := c["witness"]:
            word = witness["path"]
            witness["squared_radius_floor"] = str(min(floors[tuple(sorted(e))]
                                                      for e in zip(word, word[1:])))
            witness["square_nodes"] = sorted({i for i in word if nodes[i]["kind"] == "square"})
            if not witness["square_nodes"]:
                raise ValueError("certified obstruction lacked a square on its glued witness")
    return {"universality_status": "restricted_theorem", "components": result,
            "quotient_vertices": len(vertices), "quotient_edges": len(edges),
            "relative_boundary": sorted(boundary), "dropped_labeled_branches": len(atlas["branches"])-len(lifted),
            "singular_vertices": [[i, d] for i, d in sorted(degree.items())
                                  if d != (1 if i in boundary else 2)]}


def local_protection(atlas, info, half_width=Fraction(1, 100)):
    """A uniform C0 degree box based at a square on a nonzero class witness."""
    if not 0 < half_width < Fraction(1, 2):
        raise ValueError("half-width must lie in (0,1/2)")
    points = [pt(x, y) for x, y in atlas["polygon"]]
    for component in info["components"]:
        witness = component["witness"]
        if witness is None:
            continue
        for i in witness["square_nodes"]:
            try:
                positions = tuple(Fraction(t) for t in atlas["nodes"][i]["parameters"])
            except (TypeError, ValueError):
                continue
            ids = tuple(t.numerator//t.denominator for t in positions)
            root = tuple(t-j for t, j in zip(positions, ids))
            if not all(half_width < u < 1-half_width for u in root):
                continue
            edges = [(points[j], points[(j+1)%len(points)]) for j in ids]
            system, rhs = _affine_system(edges, 1)
            if _solve(system, rhs) != root:
                continue
            inverse_cols = [_solve(system, tuple(Fraction(i == j) for i in range(4)))
                            for j in range(4)]
            inverse_norm = max(sum(map(abs, row)) for row in zip(*inverse_cols))
            eps = min(half_width/(8*inverse_norm), Fraction(1, 1000))
            side_loss = half_width*sum(abs(a.x-b.x)+abs(a.y-b.y) for a, b in edges[:2])+4*eps
            if _side_sq(_evaluate(edges, root)) <= (1+side_loss)**2:
                continue
            return {"node": i, "edge_ids": list(ids), "edge_parameters": list(map(str, root)),
                    "half_width": str(half_width), "inverse_inf_norm": str(inverse_norm),
                    "coordinate_error": str(eps), "square_side_strictly_greater_than": "1"}
    raise ValueError("no protected rational interior square on a nonzero component")


if __name__ == "__main__":
    atlas = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    info = quotient_index(atlas)
    info["local_protection"] = local_protection(atlas, info)
    print(json.dumps(info, indent=2))
