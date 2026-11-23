"""Trace a labeled polygonal square component and certify its local degree over Q."""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, count, product
from math import gcd

from polygon_orbit_census_v223 import (
    _affine_system, _box_vertices, _evaluate, _on_seg, _rank, _side_sq, _solve,
    _square_residual, add, dot, mul, pt, sub, validate_jordan_polygon,
)


class UnsupportedComponent(ValueError):
    """The trace reaches collapse or needs a higher-dimensional size minimum."""


def segment_floor(squares):
    """Minimum squared edge length on an affine point or segment of squares."""
    if len(squares) not in (1, 2):
        raise UnsupportedComponent("component replay supports only point and segment cells")
    first, last = (sub(s[1], s[0]) for s in (squares[0], squares[-1]))
    delta = sub(last, first)
    if dot(delta, delta) == 0:
        return dot(first, first)
    t = max(0, min(1, -dot(first, delta)/dot(delta, delta)))
    e = add(first, mul(delta, t))
    return dot(e, e)


def determinant(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    return sum((-1)**j * c * determinant([row[:j]+row[j+1:] for row in matrix[1:]])
               for j,c in enumerate(matrix[0]) if c)


@dataclass
class Cell:
    matrix: list
    parameters: tuple
    squares: tuple


class PolygonComponents:
    def __init__(self, vertices):
        polygon = validate_jordan_polygon([pt(*p) for p in vertices])
        self.edges = tuple(zip(polygon, polygon[1:] + polygon[:1]))
        self.cells = {}

    def incident(self, p):
        return tuple(i for i,(a,b) in enumerate(self.edges) if _on_seg(p,a,b))

    def cell(self, ids):
        if ids not in self.cells:
            edges = [self.edges[i] for i in ids]
            matrix, rhs = _affine_system(edges, 1)
            parameters = tuple(_box_vertices(matrix, rhs))
            self.cells[ids] = Cell(matrix, parameters,
                                   tuple(_evaluate(edges,u) for u in parameters))
        return self.cells[ids]

    def trace(self, seed):
        if any(_square_residual(seed,1)) or _side_sq(seed) == 0:
            raise ValueError("expected a positive counterclockwise square")
        pending, visited, cells = [seed], set(), {}
        while pending:
            square = pending.pop()
            if square in visited:
                continue
            visited.add(square)
            incident = [self.incident(p) for p in square]
            if not all(incident):
                raise ValueError("square vertex off polygon")
            for ids in product(*incident):
                if ids not in cells:
                    cell = self.cell(ids)
                    if not cell.parameters:
                        raise ArithmeticError("incident reference chart has no feasible vertex")
                    if any(_side_sq(s) == 0 for s in cell.squares):
                        raise UnsupportedComponent("component meets the collapsed diagonal")
                    cells[ids] = cell
                    pending.extend(cell.squares)
        return cells

    def lattice_seeds(self):
        """Find search seeds on integral polygons. Positive proofs use trace(), not exhaustion."""
        points = set()
        for a,b in self.edges:
            if any(x.denominator != 1 for x in (a.x,a.y,b.x,b.y)):
                raise ValueError("lattice search requires integral vertices")
            delta = sub(b,a)
            n = gcd(int(delta.x), int(delta.y))
            points.update(add(a,mul(delta,Q(j,n))) for j in range(n))
        found = set()
        for a,b in combinations(sorted(points),2):
            e = sub(b,a)
            for turn in (1,-1):
                v = pt(-turn*e.y,turn*e.x)
                c,d = add(b,v),add(a,v)
                if c in points and d in points:
                    square = (a,b,c,d) if turn == 1 else (a,d,c,b)
                    found.add(min(square[j:]+square[:j] for j in range(4)))
        return sorted(found, key=lambda s: (-_side_sq(s),s))


def degree(cells):
    """Exact signed count for a sufficiently small generic constant residual shift."""
    for k in count(2):
        target = [1,k,k*k,k**3]
        signed_hits = []
        for ids, cell in sorted(cells.items()):
            matrix = cell.matrix
            velocity = _solve(matrix, target)
            if velocity is None:
                augmented = [row+[v] for row,v in zip(matrix,target)]
                if _rank(augmented) == _rank(matrix):
                    break  # ray lies in this singular image; choose another direction
                continue
            u, = cell.parameters
            if any(q in (0,1) and v == 0 for q,v in zip(u,velocity)):
                break
            if all((q != 0 or v > 0) and (q != 1 or v < 0) for q,v in zip(u,velocity)):
                signed_hits.append((ids, 1 if determinant(matrix)>0 else -1))
        else:
            return k, signed_hits
