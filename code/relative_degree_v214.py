"""v214: exact relative Brouwer degree for four marked periodic PL arcs.

The domain is four *chosen real* arc windows and a *fixed relative deck
translation tuple*. It intentionally excludes connector points and any other
periodic lifts. The reported integer is a finite-dimensional Brouwer degree,
NOT a Floer invariant and NOT an unconditional Toeplitz result.

Python 3.9+, standard library only; all certificate arithmetic is rational.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from typing import Optional, Sequence, Tuple

from periodic_bridge_v212 import PeriodicGraph, Point, rotate90, squared_length
from periodic_square_cert_v213 import (
    _affine_square_system, _extreme_solutions, _interval_gap,
    _solve_square, segments_of, window_height_range,
)


def _determinant(matrix):
    """Exact determinant by elimination; no floating point."""
    a = [list(map(F, row)) for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("square matrix required")
    value = F(1)
    for j in range(n):
        k = next((i for i in range(j, n) if a[i][j] != 0), None)
        if k is None:
            return F(0)
        if k != j:
            a[j], a[k] = a[k], a[j]
            value = -value
        pivot = a[j][j]
        value *= pivot
        for i in range(j+1, n):
            t = a[i][j] / pivot
            for col in range(j+1, n):
                a[i][col] -= t * a[j][col]
            a[i][j] = 0
    return value


def _distance_to_unit_interval(u):
    return max(F(0), -u, u - 1)


def _inverse_rows(A):
    n = len(A)
    cols = []
    for j in range(n):
        e = [F(0)] * n
        e[j] = F(1)
        x = _solve_square(A, e)
        if x is None:
            return None
        cols.append(x)
    return tuple(tuple(cols[j][i] for j in range(n)) for i in range(n))


def _cell_face_residual_bound(A, b):
    """Certified >0 lower bound for ||Au-b||∞ on ALL box faces.

    If A invertible, u*=A^-1 b, then each row j of A^-1 gives
    |u_j-u*_j| <= ||row_j||_1 ||Au-b||∞. Restricting u to a fixed face
    u_k=0 or 1 yields the exact elementary lower bound below. It is
    positive precisely when u* lies off every face of [0,1]^4.
    Note: positive lower bound need not be sharp.
    """
    root = _solve_square(A, b)
    rows = _inverse_rows(A)
    if root is None or rows is None:
        return None
    norms = tuple(sum(abs(v) for v in row) for row in rows)
    assert all(n > 0 for n in norms)
    distances = tuple(_distance_to_unit_interval(v) / norms[j]
                      for j, v in enumerate(root))
    per_face = []
    for j in range(4):
        for side in (F(0), F(1)):
            candidates = [abs(side - root[j]) / norms[j]]
            candidates.extend(distances[k] for k in range(4) if k != j)
            per_face.append(max(candidates))
    return min(per_face)


@dataclass(frozen=True)
class DegreeWitness:
    segment_indices: Tuple[int, int, int, int]
    parameters: Tuple[F, F, F, F]
    points: Tuple[Point, Point, Point, Point]
    orientation: int
    shifts: Tuple[int, int, int, int]
    jacobian_det: F
    side_squared: F


@dataclass(frozen=True)
class RelativeDegreeCertificate:
    degree: int
    roots: Tuple[DegreeWitness, ...]
    orientation: int
    relative_lifts: Tuple[int, int, int, int]
    opposite_height_gap: F
    boundary_residual_lower_bound: Optional[F]
    singular_empty_cells: int
    regular_cells: int

    @property
    def has_nonzero_degree(self):
        return self.degree != 0

    @property
    def vertical_perturbation_bound(self):
        """Strict epsilon threshold for all continuous vertical perturbations.

        None means the boundary bound was NOT certified, not zero degree.
        """
        if self.boundary_residual_lower_bound is None:
            return None
        return min(self.boundary_residual_lower_bound / 2,
                   self.opposite_height_gap / 2)


class InadmissibleDegreeDomain(ValueError):
    """Raised when a zero hits a cell face or the size guard fails."""


def relative_degree(graphs: Sequence[PeriodicGraph],
                    windows: Sequence[Tuple[F, F]],
                    relative_lifts=(0,0,0,0), orientation=1) -> RelativeDegreeCertificate:
    """Compute the signed relative degree when the finite certificate applies.

    The square residual is a continuous map on the four-dimensional
    product of the chosen source arcs. Each PL chart has affine residual
    A u-b. If all zeros lie strictly in regular chart interiors, the
    Brouwer degree is sum sign(det A). We *reject* any zero at an internal
    knot or source-window endpoint, and any non-isolated affine family;
    the mathematical degree may exist there but this evaluator does not
    certify it. The input also requires fixed physical relative deck lifts.

    A required positive vertical separation between opposite arcs rules
    out degenerate fourfold-point zeros at every deformation stage.
    """
    if len(graphs) != 4 or len(windows) != 4:
        raise ValueError("exactly four graphs/windows required")
    if (len(relative_lifts) != 4 or relative_lifts[0] != 0 or
            any(type(k) is not int for k in relative_lifts)):
        raise ValueError("relative lifts must be 4 integers with k0=0")
    if orientation not in (1, -1):
        raise ValueError("orientation must be +1 or -1")
    ws = tuple(tuple(map(F, w)) for w in windows)
    if any(not (0 <= lo < hi <= 1) for lo, hi in ws):
        raise ValueError("strict positive-length windows in [0,1] required")
    height0 = window_height_range(graphs[0], ws[0])
    height2 = window_height_range(graphs[2], ws[2])
    g = _interval_gap(height0, height2)
    if g <= 0:
        raise InadmissibleDegreeDomain(
            "no positive opposite-arc height gap; degeneracy not excluded")
    arcs = tuple(segments_of(graph, window) for graph, window in zip(graphs, ws))
    if any(not segs for segs in arcs):
        raise ValueError("empty source arc")
    roots = []
    margin = None
    singular_empty = regular = 0
    for indices in product(*(range(len(parts)) for parts in arcs)):
        segs = tuple(arcs[i][k] for i, k in enumerate(indices))
        if any(seg.start.x == seg.end.x for seg in segs):
            raise InadmissibleDegreeDomain("zero-length cell")
        A,b = _affine_square_system(segs, relative_lifts, orientation)
        det = _determinant(A)
        if det == 0:
            if any(True for _ in _extreme_solutions(A,b)):
                raise InadmissibleDegreeDomain(
                    "non-isolated/singular square on a PL cell")
            singular_empty += 1
            margin = None  # Cannot certify uniform all-cell margin with inverse bound.
            continue
        regular += 1
        u = _solve_square(A,b)
        assert u is not None
        if all(0 <= t <= 1 for t in u):
            if any(t == 0 or t == 1 for t in u):
                raise InadmissibleDegreeDomain(
                    "square hits a segment boundary or a source-window endpoint")
            pts = tuple(s.point(t,k) for s,t,k in zip(segs,u,relative_lifts))
            side = squared_length(pts[1]-pts[0])
            if side <= 0:
                raise InadmissibleDegreeDomain("degenerate square despite guard")
            roots.append(DegreeWitness(tuple(indices), tuple(u), pts, orientation,
                                       tuple(relative_lifts), det, side))
        # This formula additionally excludes zeros at all *internal* cell
        # interfaces, not merely the external boundary of the full domain.
        lower = _cell_face_residual_bound(A,b)
        if lower is None or lower <= 0:
            raise InadmissibleDegreeDomain("unresolved zero on a PL cell face")
        # Once a singular empty cell occurs, margin stays uncertified.
        if singular_empty == 0:
            margin = lower if margin is None else min(margin,lower)
    if singular_empty:
        margin = None
    degree = sum(1 if r.jacobian_det > 0 else -1 for r in roots)
    return RelativeDegreeCertificate(degree, tuple(roots), orientation,
                                     tuple(relative_lifts), g, margin,
                                     singular_empty, regular)


if __name__ == '__main__':
    # A standalone deterministic example, not imported from the test suite.
    origin = Point(F(1, 2), F(0))
    edge = Point(F(1, 8), F(1, 8))
    re = rotate90(edge)
    planted = (origin, origin + edge, origin + edge + re, origin + re)
    slopes = (F(-2), F(3), F(7, 2), F(-1, 3))
    graphs, windows = [], []
    hw = F(1, 32)
    for p, slope in zip(planted, slopes):
        lo, hi = p.x - hw, p.x + hw
        yl, yr = p.y - slope * hw, p.y + slope * hw
        graphs.append(PeriodicGraph((Point(F(0), yl), Point(lo, yl),
                                     Point(hi, yr), Point(F(1), yl))))
        windows.append((lo, hi))
    c = relative_degree(graphs, windows)
    print('v214 signed relative degree:',c.degree)
    print('isolated labeled squares:',len(c.roots))
    print('fixed deck tuple:',c.relative_lifts)
    print('strict allowed vertical perturbation:',c.vertical_perturbation_bound)
