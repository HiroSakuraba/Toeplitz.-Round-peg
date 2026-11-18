"""Exact decision procedure for labeled squares on four periodic PL graphs.

Period = 1 horizontally. No floating point, no sampling, no unbounded
search over lifts. This solves the *finite PL model*, not Toeplitz for
arbitrary Jordan curves and not any Floer theorem.

The admissible domain for each vertex is either a full graph or a closed
subarc [lo, hi] of its periodic parameter. Both square orientations are
considered. If the result is absent, it certifies that no nondegenerate
labeled square exists with vertices in the requested four subarcs.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from typing import Optional, Sequence, Tuple

from periodic_bridge_v212 import Point, PeriodicGraph, rotate90, squared_length


@dataclass(frozen=True)
class Segment:
    start: Point
    end: Point

    def __post_init__(self):
        if self.start.x > self.end.x:
            raise ValueError("graph segment must not go backwards in x")

    def point(self, u: F, shift: int = 0) -> Point:
        return Point(self.start.x + (self.end.x - self.start.x) * u + shift,
                     self.start.y + (self.end.y - self.start.y) * u)

    def ranges(self, shift: int = 0):
        return ((self.start.x + shift, self.end.x + shift),
                (min(self.start.y, self.end.y), max(self.start.y, self.end.y)))


def segments_of(g: PeriodicGraph, window: Optional[Tuple[F, F]] = None) -> Tuple[Segment, ...]:
    """All nonempty segment pieces within the supplied closed parameter arc."""
    lo, hi = (F(0), F(1)) if window is None else tuple(map(F, window))
    if not (F(0) <= lo <= hi <= F(1)):
        raise ValueError("window must satisfy 0 <= lo <= hi <= 1")
    result = []
    for a, b in zip(g.vertices, g.vertices[1:]):
        left, right = max(lo, a.x), min(hi, b.x)
        if left < right or (lo == hi and left == right == lo and a.x <= lo <= b.x and not result):
            result.append(Segment(Point(left, g.height(left)),
                                  Point(right, g.height(right))))
    return tuple(result)


def _rank(rows: Sequence[Sequence[F]]) -> int:
    if not rows:
        return 0
    m = [list(map(F, row)) for row in rows]
    n = len(m[0]); pivot = 0
    for col in range(n):
        q = next((i for i in range(pivot, len(m)) if m[i][col]), None)
        if q is None:
            continue
        m[pivot], m[q] = m[q], m[pivot]
        factor = m[pivot][col]
        m[pivot] = [x / factor for x in m[pivot]]
        for i in range(pivot + 1, len(m)):
            if m[i][col]:
                t = m[i][col]
                m[i] = [a - t*b for a, b in zip(m[i], m[pivot])]
        pivot += 1
        if pivot == len(m):
            break
    return pivot


def _solve_square(A: Sequence[Sequence[F]], b: Sequence[F]) -> Optional[Tuple[F, ...]]:
    """Solve one exactly determined rational square linear system."""
    n = len(b)
    m = [list(map(F, row)) + [F(rhs)] for row, rhs in zip(A, b)]
    for j in range(n):
        i = next((k for k in range(j, n) if m[k][j]), None)
        if i is None:
            return None
        m[i], m[j] = m[j], m[i]
        t = m[j][j]
        m[j] = [z / t for z in m[j]]
        for k in range(n):
            if k != j and m[k][j]:
                t = m[k][j]
                m[k] = [a - t*z for a, z in zip(m[k], m[j])]
    return tuple(row[-1] for row in m)


def _extreme_solutions(A: Sequence[Sequence[F]], b: Sequence[F]):
    """All candidate vertices of {A u=b, 0<=u<=1} for four variables.

    A nonempty affine slice of a box has a vertex. A vertex has at least
    (4-rank(A)) active coordinate bounds. We enumerate all such bounds,
    with no numerical optimization dependency or generic-rank assumption.
    """
    rank = _rank(A)
    independent = []
    for row, rhs in zip(A, b):
        if _rank([r for r, _ in independent] + [row]) > len(independent):
            independent.append((list(row), F(rhs)))
    if len(independent) != rank:
        raise AssertionError("rank selection inconsistent")
    seen = set()
    for indexes in combinations(range(4), 4 - rank):
        for vals in product((F(0), F(1)), repeat=4-rank):
            rows = [r for r, _ in independent]
            target = [rhs for _, rhs in independent]
            for j, val in zip(indexes, vals):
                row = [F(0)]*4
                row[j] = F(1)
                rows.append(row)
                target.append(val)
            sol = _solve_square(rows, target)
            if (sol is None or sol in seen or any(u < 0 or u > 1 for u in sol)):
                continue
            if all(sum(a*u for a,u in zip(row,sol)) == rhs for row,rhs in zip(A,b)):
                seen.add(sol)
                yield sol


def _square_residual(points: Sequence[Point], orientation: int) -> Tuple[F, F, F, F]:
    p0, p1, p2, p3 = points
    e = p1 - p0
    r = rotate90(e)
    r = Point(orientation*r.x, orientation*r.y)
    f = p2 - p1 - r
    g = p3 - p0 - r
    return f.x, f.y, g.x, g.y




def _affine_square_system(segs: Sequence[Segment], shifts, orientation):
    if len(segs) != 4 or len(shifts) != 4 or orientation not in (-1,1):
        raise ValueError("four segments and shifts, orientation +/-1 required")
    starts = [s.point(F(0), shift) for s,shift in zip(segs,shifts)]
    base = _square_residual(starts, orientation)
    cols = []
    for j in range(4):
        pts = starts.copy()
        pts[j] = segs[j].point(F(1),shifts[j])
        z = _square_residual(pts,orientation)
        cols.append(tuple(a-b for a,b in zip(z,base)))
    A = [list(row) for row in zip(*cols)]
    b = tuple(-z for z in base)
    return A,b

@dataclass(frozen=True)
class SquareWitness:
    points: Tuple[Point, ...]
    orientation: int
    shifts: Tuple[int, ...]
    parameters: Tuple[F, ...]
    side_squared: F


def solve_four_segments(segs: Sequence[Segment], shifts=(0,0,0,0),
                        orientation=1) -> Optional[SquareWitness]:
    """Complete exact feasible-square solver for one segment/lift choice.

    Includes rank-deficient families and endpoints. A nondegenerate solution
    exists iff at least one feasible polytope vertex has nonzero edge: the
    edge is affine in u, and every feasible point is a convex combination
    of vertices. Hence testing only exact vertices is complete.
    """
    if len(segs) != 4 or len(shifts) != 4 or orientation not in (-1,1):
        raise ValueError("four segments and shifts, orientation +/-1 required")
    A, b = _affine_square_system(segs,shifts,orientation)
    for u in _extreme_solutions(A,b):
        points = tuple(s.point(t,k) for s,t,k in zip(segs,u,shifts))
        sq = squared_length(points[1]-points[0])
        if sq > 0 and _square_residual(points,orientation) == (0,0,0,0):
            return SquareWitness(points,orientation,tuple(shifts),u,sq)
    return None




@dataclass(frozen=True)
class TransverseSquareCertificate:
    """Sufficient open-set existence certificate, NOT a global Toeplitz proof."""
    witness: SquareWitness
    parameter_margin: F
    inverse_inf_norm: F
    opposite_height_gap: F
    epsilon_bound: F


def certify_transverse_square(segs: Sequence[Segment], shifts=(0,0,0,0),
                               orientation=1) -> TransverseSquareCertificate:
    """Certifies persistence under arbitrary continuous vertical perturbations.

    For each segment p_i(u), allow p_i(u)+(0,h_i(u)) with |h_i(u)|<=eps.
    If eps < epsilon_bound, Brouwer applied to u* - A^{-1}E(u) gives
    a square with parameters remaining in (0,1). Vertical separation
    of opposite segment ranges ensures this square is nondegenerate.

    Requires exact full-rank linearized square system, interior witness,
    and disjoint y-ranges of opposite labeled arcs. Strictly sufficient.
    """
    A,b = _affine_square_system(segs,shifts,orientation)
    sol = _solve_square(A,b)
    if sol is None or any(not (0 < z < 1) for z in sol):
        raise ValueError("no unique transverse square with interior parameters")
    points = tuple(s.point(u,k) for s,u,k in zip(segs,sol,shifts))
    side2=squared_length(points[1]-points[0])
    if side2 == 0 or any(sum(a*u for a,u in zip(row,sol)) != rhs for row,rhs in zip(A,b)):
        raise ValueError("not a nondegenerate square")
    gap=_interval_gap(segs[0].ranges()[1],segs[2].ranges()[1])
    if gap == 0:
        raise ValueError("opposite segments lack required vertical separation")
    # Matrix columns of A^{-1}, assembled using exact Gaussian elimination.
    cols=[]
    for j in range(4):
        unit=[F(0)]*4; unit[j]=F(1)
        col=_solve_square(A,unit)
        if col is None:
            raise ValueError("singular constraint matrix")
        cols.append(col)
    inverse_inf=max(sum(abs(cols[j][i]) for j in range(4)) for i in range(4))
    margin=min(min(u,1-u) for u in sol)
    limit=min(margin/(2*inverse_inf),gap/2)
    return TransverseSquareCertificate(
        SquareWitness(points,orientation,tuple(shifts),sol,side2),
        margin,inverse_inf,gap,limit)

def _add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _sub(a,b):
    return a[0]-b[1], a[1]-b[0]


def _neg(a):
    return -a[1], -a[0]


def _overlap(a,b):
    return max(a[0],b[0]) <= min(a[1],b[1])


def _intervals_allow_square(segs, shifts, orientation):
    """Necessary interval identities; pruning only, never claims existence."""
    boxes = [s.ranges(k) for s,k in zip(segs,shifts)]
    x = [z[0] for z in boxes]
    y = [z[1] for z in boxes]
    if not _overlap(_sub(x[0],x[2]),
                    _sub(y[1],y[3]) if orientation == 1 else _neg(_sub(y[1],y[3]))):
        return False
    if not _overlap(_sub(y[0],y[2]),
                    _sub(x[3],x[1]) if orientation == 1 else _neg(_sub(x[3],x[1]))):
        return False
    return (_overlap(_add(x[0],x[2]), _add(x[1],x[3])) and
            _overlap(_add(y[0],y[2]), _add(y[1],y[3])))




def _interval_gap(a: Tuple[F, F], b: Tuple[F, F]) -> F:
    """Minimum distance between two ordinary real closed intervals."""
    if a[1] < b[0]:
        return b[0] - a[1]
    if b[1] < a[0]:
        return a[0] - b[1]
    return F(0)


def window_horizontal_gap(a: Tuple[F, F], b: Tuple[F, F]) -> F:
    """Exact circular distance between two closed intervals in R/Z."""
    a = tuple(map(F,a)); b = tuple(map(F,b))
    if not (F(0) <= a[0] <= a[1] <= 1 and F(0) <= b[0] <= b[1] <= 1):
        raise ValueError("unit-period closed windows required")
    return min(_interval_gap(a,(b[0]+k,b[1]+k)) for k in (-1,0,1))


def window_height_range(graph: PeriodicGraph, window: Tuple[F, F]) -> Tuple[F, F]:
    """Exact range of linear graph height on a closed window."""
    lo,hi = tuple(map(F,window))
    if not F(0) <= lo <= hi <= F(1):
        raise ValueError("window outside full period")
    samples = [graph.height(lo),graph.height(hi)]
    samples += [p.y for p in graph.vertices if lo < p.x < hi]
    return min(samples),max(samples)


@dataclass(frozen=True)
class ExclusionCertificate:
    horizontal_separation: F
    maximum_opposite_height_difference: F

    @property
    def margin(self):
        return self.horizontal_separation - self.maximum_opposite_height_difference

    @property
    def proves_absence(self):
        return self.margin > 0


def source_window_exclusion(graphs: Sequence[PeriodicGraph],
                            windows: Sequence[Tuple[F, F]]) -> ExclusionCertificate:
    """Sufficient analytic certificate excluding labeled squares of both orientations.

    Every square obeys |x0-x2|=|y1-y3| in planar lifts.
    LHS is at least the circular gap between W0,W2, while RHS is at most
    the farthest y-difference between graph1|W1 and graph3|W3.
    The strict positive margin proves *absence* without enumeration.
    """
    if len(graphs) != 4 or len(windows) != 4:
        raise ValueError("four graphs and four windows required")
    dx = window_horizontal_gap(windows[0],windows[2])
    a,b = window_height_range(graphs[1],windows[1]),window_height_range(graphs[3],windows[3])
    dy = max(abs(a[0]-b[1]),abs(a[1]-b[0]))
    return ExclusionCertificate(dx,dy)

@dataclass(frozen=True)
class SearchResult:
    witness: Optional[SquareWitness]
    solved_cases: int
    rejected_cases: int
    lift_bound: int

    @property
    def no_square(self):
        return self.witness is None


def find_labeled_square(graphs: Sequence[PeriodicGraph],
                        windows: Optional[Sequence[Tuple[F, F]]] = None,
                        orientations=(1,-1),
                        relative_lifts: Optional[Tuple[int, int, int, int]] = None) -> SearchResult:
    """Complete search over all period lifts and four PL segment choices.

    A square has x-extent equal to its y-extent. If all graph heights lie
    between Ymin and Ymax, then its |x_i-x_0| <= H=Ymax-Ymin. Each x_i
    has a representative in [0,1], hence after shifting x0 into [0,1],
    every integer lift k_i obeys |k_i| <= ceil(H)+1. Finite completeness.

    If relative_lifts is set, the search is confined to that specified
    deck-translation tuple (normalized k0=0). This is necessary when the
    original source arcs occupy *fixed planar lifts*, rather than any
    translated copy of their cylindrical projections.

    Returns first exact nondegenerate square or an exhaustive negative answer.
    """
    if len(graphs) != 4:
        raise ValueError("four graphs required")
    orientations = tuple(orientations)
    if not orientations or any(o not in (1, -1) for o in orientations):
        raise ValueError("orientations must be nonempty and +/-1")
    if windows is not None and len(windows) != 4:
        raise ValueError("four closed windows required")
    if relative_lifts is not None:
        if (len(relative_lifts) != 4 or relative_lifts[0] != 0 or
                any(type(k) is not int for k in relative_lifts)):
            raise ValueError("relative_lifts must be four integers, first zero")
    arcs = [segments_of(g, None if windows is None else windows[i])
            for i,g in enumerate(graphs)]
    heights = [p.y for g in graphs for p in g.vertices]
    H = max(heights)-min(heights)
    K = -(-H.numerator // H.denominator) + 1
    if H == 0 or any(not a for a in arcs):
        return SearchResult(None, 0, 0, K)
    if relative_lifts is not None and any(abs(k)>K for k in relative_lifts):
        return SearchResult(None, 0, 0, K)
    solved = rejected = 0
    # The first vertex lift is normalized to k0=0 by integer translation.
    for segs in product(*arcs):
        choices_k2 = range(-K,K+1) if relative_lifts is None else (relative_lifts[2],)
        for k2 in choices_k2:
            # This x0/x2 projection must match a y1/y3 difference.
            rx0, ry0 = segs[0].ranges()
            rx2, ry2 = segs[2].ranges(k2)
            ry1, ry3 = segs[1].ranges()[1], segs[3].ranges()[1]
            difference = _sub(rx0,rx2)
            ydifference = _sub(ry1,ry3)
            if not (_overlap(difference,ydifference) or
                    _overlap(difference,_neg(ydifference))):
                rejected += ((2*K+1)**2 if relative_lifts is None else 1)*len(orientations)
                continue
            # A fresh iterator is mandatory for each k2. Do not reuse a
            # consumed itertools.product: that would miss valid lifts.
            choices_k13 = (product(range(-K,K+1),repeat=2)
                           if relative_lifts is None else
                           ((relative_lifts[1],relative_lifts[3]),))
            for k1,k3 in choices_k13:
                shifts = (0,k1,k2,k3)
                for orientation in orientations:
                    if not _intervals_allow_square(segs,shifts,orientation):
                        rejected += 1
                        continue
                    solved += 1
                    witness = solve_four_segments(segs,shifts,orientation)
                    if witness is not None:
                        return SearchResult(witness,solved,rejected,K)
    return SearchResult(None,solved,rejected,K)


if __name__ == "__main__":
    from periodic_bridge_v212 import four_strands, SOURCE_WINDOWS
    fs = four_strands()
    orig = find_labeled_square(fs, SOURCE_WINDOWS)
    full = find_labeled_square(fs)
    certificate = source_window_exclusion(fs, SOURCE_WINDOWS)
    print("v212 original windows:","NO LABELED SQUARE" if orig.no_square else "square found")
    print("exact exclusion margin:", certificate.margin)
    print("completed periodic graphs:", "square side²=" + str(full.witness.side_squared)
          if full.witness else "no square")

    levels = tuple(PeriodicGraph((Point(F(0), F(h)), Point(F(1), F(h))))
                   for h in (1, 0, -1, 0))
    window = ((F(7,16), F(9,16)),) * 4
    periodic = find_labeled_square(levels, window)
    physical = find_labeled_square(levels, window, relative_lifts=(0,0,0,0))
    print("marked quotient windows:",periodic.witness.shifts if periodic.witness else "none")
    print("same physical lifts:","NO SQUARE" if physical.no_square else "square found")