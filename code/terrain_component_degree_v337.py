"""Exact edge-chart and infinitesimal shear certificates for the v337 fixtures."""
from fractions import Fraction
from itertools import combinations, product

from polygon_orbit_census_v223 import (
    _affine_system, _box_vertices, _evaluate, _side_sq, _solve, pt,
    validate_jordan_polygon,
)

TERRAIN = ((0, 0), (2, 0), (2, 7), (1, 7), (1, 2), (0, 2))
GLUED = ((0, 0), (6, 0), (6, -1), (8, -1), (8, 3),
         (7, 3), (7, 7), (3, 7), (3, 8), (0, 8))


def chart_cells(vertices):
    """Complete positive-size affine cells, up to rotation of CCW labels.

    Each cell is returned in the vertex representation of its convex polytope.
    Endpoint hits may occur in more than one chart.
    """
    polygon = validate_jordan_polygon([pt(*p) for p in vertices])
    edges = tuple(zip(polygon, polygon[1:] + polygon[:1]))
    for ids in product(range(len(edges)), repeat=4):
        if any(ids.count(i) > 2 for i in ids):
            continue
        if any(ids[k:] + ids[:k] < ids for k in (1, 2, 3)):
            continue
        selected = [edges[i] for i in ids]
        matrix, rhs = _affine_system(selected, 1)
        solution = _solve(matrix, rhs)
        if solution is None:
            parameters = _box_vertices(matrix, rhs)
        elif all(0 <= u <= 1 for u in solution):
            parameters = (solution,)
        else:
            continue
        corners = tuple(_evaluate(selected, u) for u in parameters)
        if any(_side_sq(square) > 0 for square in corners):
            yield ids, corners


def large_cells(vertices, cutoff_squared):
    # Include equality so a separate component touching the cutoff is detected.
    return {ids: corners for ids, corners in chart_cells(vertices)
            if max(map(_side_sq, corners)) >= cutoff_squared}


def sign_near_zero(expr, t):
    """Exact sign in the ordered field of rational germs at t=0+."""
    import sympy as sp

    numerator, denominator = sp.fraction(sp.cancel(expr))

    def lowest_sign(poly):
        if poly == 0:
            return 0
        return sp.sign(sp.Poly(poly, t).terms()[-1][1])

    return lowest_sign(numerator) * lowest_sign(denominator)


def shear_system(vertices, ids, direction, t):
    import sympy as sp

    polygon = [sp.Matrix([x + direction*t*y, y]) for x, y in vertices]
    u = sp.symbols('u:4')
    points = [polygon[i] + v * (polygon[(i+1) % len(polygon)] - polygon[i])
              for i, v in zip(ids, u)]
    p, q, r, s = points
    dx, dy = q - p
    quarter_turn = sp.Matrix([-dy, dx])
    residual = (r - q - quarter_turn).col_join(s - p - quarter_turn)
    matrix, rhs = sp.linear_eq_to_matrix(residual, u)
    return matrix, rhs, (q-p).dot(q-p), u


def feasible_vertices(matrix, rhs, t):
    """Vertices of a rational-function affine slice of the closed unit box."""
    import sympy as sp

    independent = matrix.T.rref()[1]
    rows = list(independent)
    rank = len(rows)
    vertices = set()
    for axes in combinations(range(4), 4-rank):
        basis = matrix[rows, :].col_join(sp.eye(4)[list(axes), :])
        if sp.cancel(basis.det()) == 0:
            continue
        inverse = basis.inv()
        for bounds in product((0, 1), repeat=len(axes)):
            target = rhs[rows, :].col_join(sp.Matrix(len(bounds), 1, bounds))
            solution = (inverse * target).applyfunc(sp.cancel)
            if any(sp.cancel(e) != 0 for e in matrix*solution - rhs):
                continue
            if all(sign_near_zero(u, t) >= 0 and sign_near_zero(1-u, t) >= 0
                   for u in solution):
                vertices.add(tuple(solution))
    return vertices


def shear_witnesses(vertices, cells, direction):
    """Check every limiting candidate chart for all sufficiently small t>0."""
    import sympy as sp

    t = sp.Symbol('t', positive=True)
    found = {}
    for ids in cells:
        matrix, rhs, size, variables = shear_system(vertices, ids, direction, t)
        feasible = feasible_vertices(matrix, rhs, t)
        if not feasible:
            continue
        if len(feasible) != 1 or sp.cancel(matrix.det()) == 0:
            raise ValueError(f'nonregular sheared solution on {ids}')
        coordinates, = feasible
        if any(sign_near_zero(u, t) <= 0 or sign_near_zero(1-u, t) <= 0
               for u in coordinates):
            raise ValueError(f'sheared solution on a chart boundary: {ids}')
        found[ids] = (sp.factor(matrix.det()), sp.cancel(size.subs(dict(zip(variables, coordinates)))))
    return found


if __name__ == '__main__':
    for name, polygon, cutoff, direction in (
        ('(2,7)', TERRAIN, Fraction(3, 2), 1), ('glued A', GLUED, 24, -1)
    ):
        cells = large_cells(polygon, cutoff)
        print(name, 'candidate charts', len(cells))
        for ids, (det, size) in shear_witnesses(polygon, cells, direction).items():
            print('  regular chart', ids, 'determinant', det, 'squared side', size)
