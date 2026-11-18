"""v217: exact signed boundary-flux ledger for subdivided square configurations.

Four ORIENTED rational polygonal arcs, one physical lift tuple, and one
square orientation. Each child box inherits the same residual. An internal
face is counted twice with opposite orientations, even when it contains a
true square. A child's ray flux is NOT its degree if its boundary has zeros.
No Floer construction or general Toeplitz proof is claimed.

Python standard library only (3.9+). Fractions are exact.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product
from typing import Tuple


@dataclass(frozen=True)
class Point:
    x: Q
    y: Q


def pt(x, y):
    return Point(Q(x), Q(y))


@dataclass(frozen=True)
class Arc:
    vertices: Tuple[Point, ...]

    def __post_init__(self):
        if len(self.vertices) < 2:
            raise ValueError("An arc needs two vertices")
        if any(a == b for a, b in zip(self.vertices, self.vertices[1:])):
            raise ValueError("Consecutive polygon vertices must differ")
        if any(type(v.x) is not Q or type(v.y) is not Q for v in self.vertices):
            raise ValueError("Coordinates must be Fraction")


def arc(*pairs):
    return Arc(tuple(pt(x, y) for x, y in pairs))


def _at(a, t):
    t = Q(t)
    if not Q(0) <= t <= Q(1):
        raise ValueError("Arc parameter outside [0,1]")
    n = len(a.vertices) - 1
    j = min(int(t * n), n - 1)
    u = t * n - j
    p, q = a.vertices[j:j+2]
    return pt(p.x + u*(q.x-p.x), p.y + u*(q.y-p.y))


def subarc(a, lo, hi):
    """Orientation-preserving restriction in equal-segment parameterization."""
    lo, hi = Q(lo), Q(hi)
    if not Q(0) <= lo < hi <= Q(1):
        raise ValueError("Require 0<=lo<hi<=1")
    n = len(a.vertices) - 1
    parts = [_at(a, lo)]
    parts.extend(a.vertices[i] for i in range(1, n) if lo < Q(i, n) < hi)
    parts.append(_at(a, hi))
    return Arc(tuple(parts))


def partition_arcs(arcs, cuts):
    """Return four tuples of subarcs and their exact parameter intervals."""
    if len(arcs) != 4 or len(cuts) != 4:
        raise ValueError("Four arcs and four cut schedules required")
    pieces, domains = [], []
    for a, inner in zip(arcs, cuts):
        if not isinstance(a, Arc):
            raise ValueError("Expected Arc")
        c = (Q(0),) + tuple(Q(t) for t in inner) + (Q(1),)
        if any(x >= y for x, y in zip(c, c[1:])):
            raise ValueError("Cuts must be strictly ordered and internal")
        domains.append(tuple(zip(c, c[1:])))
        pieces.append(tuple(subarc(a, x, y) for x, y in zip(c, c[1:])))
    return tuple(pieces), tuple(domains)


def _residual(points, orientation):
    p0, p1, p2, p3 = points
    ex, ey = p1.x-p0.x, p1.y-p0.y
    jx, jy = -orientation*ey, orientation*ex
    return (p2.x-p1.x-jx, p2.y-p1.y-jy,
            p3.x-p0.x-jx, p3.y-p0.y-jy)


def _rank(rows):
    if not rows:
        return 0
    m = [list(map(Q, row)) for row in rows]
    n, rank = len(m[0]), 0
    for j in range(n):
        k = next((i for i in range(rank, len(m)) if m[i][j]), None)
        if k is None:
            continue
        m[rank], m[k] = m[k], m[rank]
        t = m[rank][j]
        m[rank] = [v/t for v in m[rank]]
        for i in range(rank+1, len(m)):
            if m[i][j]:
                t = m[i][j]
                m[i] = [x-t*y for x, y in zip(m[i], m[rank])]
        rank += 1
        if rank == len(m):
            break
    return rank


def _solve(a, b):
    n = len(b)
    if len(a) != n or any(len(row) != n for row in a):
        raise ValueError("Square system required")
    m = [list(map(Q, row)) + [Q(target)] for row, target in zip(a, b)]
    for j in range(n):
        k = next((i for i in range(j, n) if m[i][j]), None)
        if k is None:
            return None
        m[j], m[k] = m[k], m[j]
        t = m[j][j]
        m[j] = [v/t for v in m[j]]
        for i in range(n):
            if i != j and m[i][j]:
                t = m[i][j]
                m[i] = [x-t*y for x, y in zip(m[i], m[j])]
    return tuple(row[-1] for row in m)


def _det(a):
    m = [list(map(Q, row)) for row in a]
    n, sign = len(m), Q(1)
    for j in range(n):
        k = next((i for i in range(j, n) if m[i][j]), None)
        if k is None:
            return Q(0)
        if k != j:
            m[j], m[k] = m[k], m[j]
            sign = -sign
        t = m[j][j]
        sign *= t
        for i in range(j+1, n):
            c = m[i][j]/t
            for h in range(j+1, n):
                m[i][h] -= c*m[j][h]
            m[i][j] = Q(0)
    return sign


def _box_vertices(a, b, n):
    """Enumerate rational vertices of {Ax=b, 0<=x<=1}."""
    if len(a) != len(b) or any(len(row) != n for row in a):
        raise ValueError("Bad rectangular system")
    rank = _rank(a)
    rows, targets = [], []
    for row, target in zip(a, b):
        if _rank(rows + [row]) > len(rows):
            rows.append(list(row))
            targets.append(Q(target))
    for indices in combinations(range(n), n-rank):
        for values in product((Q(0), Q(1)), repeat=n-rank):
            square = [list(row) for row in rows]
            rhs = list(targets)
            for j, value in zip(indices, values):
                row = [Q(0)]*n
                row[j] = Q(1)
                square.append(row)
                rhs.append(value)
            sol = _solve(square, rhs)
            if sol is not None and all(0 <= v <= 1 for v in sol) and all(
                sum(x*y for x, y in zip(row, sol)) == target
                for row, target in zip(a, b)):
                yield sol


def _feasible_box(a, b, n):
    return next(_box_vertices(a, b, n), None) is not None


def _facets(arcs, shifts, orientation):
    segments = [tuple(zip(a.vertices, a.vertices[1:])) for a in arcs]
    for fixed in range(4):
        for side in (0, 1):
            endpoint = arcs[fixed].vertices[0 if side == 0 else -1]
            choices = [((endpoint, endpoint),) if j == fixed else segments[j]
                       for j in range(4)]
            for segments4 in product(*choices):
                starts = [pt(p.x+k, p.y) for (p, _), k in zip(segments4, shifts)]
                base = _residual(starts, orientation)
                columns = []
                for j in range(4):
                    modified = list(starts)
                    p, q = segments4[j]
                    modified[j] = pt(starts[j].x+q.x-p.x,
                                     starts[j].y+q.y-p.y)
                    v = _residual(modified, orientation)
                    columns.append(tuple(v[i]-base[i] for i in range(4)))
                free = tuple(columns[j] for j in range(4) if j != fixed)
                yield fixed, side, base, free


def _sub(a, b):
    return pt(a.x-b.x, a.y-b.y)


def _dot(a, b):
    return a.x*b.x + a.y*b.y


def _cross(a, b):
    return a.x*b.y - a.y*b.x


def _segment_distance_sq(a, b, c, d):
    ab, cd = _sub(b,a), _sub(d,c)
    crossval = _cross(ab,cd)
    if crossval:
        t = _cross(_sub(c,a),cd)/crossval
        u = _cross(_sub(c,a),ab)/crossval
        if 0 <= t <= 1 and 0 <= u <= 1:
            return Q(0)
    def point_dist(p, v, w):
        vw = _sub(w,v)
        t = max(Q(0),min(Q(1),_dot(_sub(p,v),vw)/_dot(vw,vw)))
        diff = pt(p.x-v.x-t*vw.x,p.y-v.y-t*vw.y)
        return _dot(diff,diff)
    return min(point_dist(a,c,d),point_dist(b,c,d),
               point_dist(c,a,b),point_dist(d,a,b))


def opposite_gap_squared(arcs, shifts=(0,0,0,0)):
    """Exact Euclidean separation of physical opposing source arcs 0 and 2."""
    a,b=arcs[0],arcs[2]
    delta=shifts[2]-shifts[0]
    return min(_segment_distance_sq(p,q,pt(r.x+delta,r.y),pt(t.x+delta,t.y))
               for p,q in zip(a.vertices,a.vertices[1:])
               for r,t in zip(b.vertices,b.vertices[1:]))


class BoundaryRoot(ValueError):
    pass


class NonGenericRay(RuntimeError):
    pass


class Inconclusive(RuntimeError):
    pass


def boundary_has_square(arcs, shifts=(0,0,0,0), orientation=1):
    for _, _, base, cols in _facets(arcs, shifts, orientation):
        matrix = [[cols[j][i] for j in range(3)] for i in range(4)]
        if _feasible_box(matrix, [-v for v in base], 3):
            return True
    return False


@dataclass(frozen=True)
class Flux:
    total: int
    face_counts: Tuple[Tuple[int, int, int], ...]
    positive_ray_hits: int

    def face(self, i, side):
        return next(count for a,b,count in self.face_counts if a==i and b==side)


def ray_flux(arcs, ray, shifts=(0,0,0,0), orientation=1):
    """Oriented ray intersections with F(∂product(arcs)).

    This is a FORMAL boundary-chain flux. It equals Brouwer degree only
    when the particular box has no zero on its EXTERNAL boundary.
    """
    if len(arcs) != 4 or any(not isinstance(a, Arc) for a in arcs):
        raise ValueError("Four Arc objects required")
    if len(shifts) != 4 or shifts[0] != 0 or any(type(k) is not int for k in shifts):
        raise ValueError("Four physical lift integers, normalized at first=0")
    if orientation not in (-1, 1):
        raise ValueError("Orientation must be +/-1")
    if len(ray) != 4 or not any(ray):
        raise ValueError("Nonzero four-dimensional ray required")
    ray = tuple(map(Q, ray))
    counts = {(i, s): 0 for i in range(4) for s in (0, 1)}
    hits = 0
    for i, side, base, cols in _facets(arcs, shifts, orientation):
        matrix = [[cols[j][k] for j in range(3)] + [-ray[k]] for k in range(4)]
        rhs = [-v for v in base]
        det = _det(matrix)
        if det:
            solution = _solve(matrix, rhs)
            assert solution is not None
            parameters, lam = solution[:3], solution[3]
            if lam < 0 or any(v < 0 or v > 1 for v in parameters):
                continue
            if lam == 0:
                # Origin is an internal child-wall root: the positive ray
                # deliberately excludes lambda=0. Flux is not child degree.
                continue
            if any(v in (0, 1) for v in parameters):
                raise NonGenericRay("Ray meets a positive-distance facet seam")
            face_sign = (-1)**(i if side == 1 else i+1)
            counts[(i, side)] += face_sign * (1 if det > 0 else -1)
            hits += 1
        else:
            span = max(abs(base[k] + sum(cols[j][k]*u[j] for j in range(3)))
                       for u in product((Q(0), Q(1)), repeat=3)
                       for k in range(4))
            limit = span / max(abs(t) for t in ray)
            if limit <= 0:
                # This facet maps only to 0. It contributes no positive hits.
                continue
            scaled = [[cols[j][k] for j in range(3)] + [-limit*ray[k]]
                      for k in range(4)]
            if any(sol[3] > 0 for sol in _box_vertices(scaled, rhs, 4)):
                raise NonGenericRay("Positive ray lies in a singular facet image")
    return Flux(sum(counts.values()),
                tuple((i, s, counts[(i,s)]) for i in range(4) for s in (0,1)),
                hits)


@dataclass(frozen=True)
class Ledger:
    ray: Tuple[Q, Q, Q, Q]
    parent_degree: int
    opposite_distance_squared: Q
    sum_child_flux: int
    children: Tuple[Tuple[Tuple[int, ...], int, bool], ...]
    internal_wall_pairs: int
    boundary_root_children: int
    generic_attempts: int


def multiscale_ledger(arcs, cuts=((),(),(),()), shifts=(0,0,0,0),
                      orientation=1, max_rays=64):
    """Certified sum of subdivided boundary-chain fluxes.

    The parent's boundary MUST exclude square zeros. Children may have
    square zeros on internal cuts, so their fluxes are NOT called degrees
    unless their own boundary is zero-free.
    """
    gap2 = opposite_gap_squared(arcs, shifts)
    if gap2 <= 0:
        raise ValueError("Opposite physical arcs intersect: square may collapse")
    if boundary_has_square(arcs, shifts, orientation):
        raise BoundaryRoot("The parent has a square on its external boundary")
    pieces, domains = partition_arcs(arcs, cuts)
    indexed = tuple(product(*(range(len(x)) for x in pieces)))
    children = tuple((index, tuple(pieces[j][k] for j,k in enumerate(index)))
                     for index in indexed)
    for n in range(1, max_rays+1):
        ray = (Q(1), Q(n), Q(n*n), Q(n*n*n))
        try:
            parent = ray_flux(arcs, ray, shifts, orientation)
            child_flows = tuple((idx, ray_flux(child, ray, shifts, orientation))
                                for idx, child in children)
        except NonGenericRay:
            continue
        break
    else:
        raise Inconclusive("No common generic ray within the search budget")

    interface = {}
    for idx, flow in child_flows:
        for i in range(4):
            for side in (0,1):
                neighbor = idx[i] - 1 if side == 0 else idx[i] + 1
                if neighbor < 0 or neighbor >= len(pieces[i]):
                    continue
                cut_at = domains[i][idx[i]][side]
                key = (i, cut_at, tuple(idx[j] for j in range(4) if j != i))
                interface.setdefault(key, []).append(flow.face(i,side))
    for key, values in interface.items():
        if len(values) != 2 or sum(values) != 0:
            raise AssertionError(("Internal boundary cancellation failed", key, values))
    subtotal = sum(flow.total for _, flow in child_flows)
    if subtotal != parent.total:
        raise AssertionError(("Telescoping degree failed", subtotal, parent.total))
    certified = tuple((idx, flow.total,
                       not boundary_has_square(child, shifts, orientation))
                      for (idx, child), (_,flow) in zip(children, child_flows))
    return Ledger(ray, parent.total, gap2, subtotal, certified, len(interface),
                  sum(not ok for _,_,ok in certified), n)


def planted_example():
    """An isolated transverse square with rational source arcs."""
    p0 = pt(Q(1,2), Q(0))
    edge = pt(Q(1,8), Q(1,8))
    j = pt(-edge.y, edge.x)
    p1 = pt(p0.x+edge.x, p0.y+edge.y)
    p2 = pt(p1.x+j.x, p1.y+j.y)
    p3 = pt(p0.x+j.x, p0.y+j.y)
    slopes = (Q(-2), Q(3), Q(7,2), Q(-1,3))
    hw = Q(1,32)
    return tuple(arc((p.x-hw,p.y-m*hw),(p.x+hw,p.y+m*hw))
                 for p,m in zip((p0,p1,p2,p3),slopes))


def cancellation_example():
    """Two transverse square zeros with opposite signs on one bent source arc."""
    parts = list(planted_example())
    lo, mid, hi = Q(15,32), Q(1,2), Q(17,32)
    bump = Q(1,100)
    target = lambda x: Q(1,4)-2*(x-Q(1,2))
    parts[2] = arc((lo,target(lo)-bump),
                   (mid,target(mid)+bump),
                   (hi,target(hi)-bump))
    return tuple(parts)


def audit():
    """Exact representative ledger, reproduced in tests."""
    good = multiscale_ledger(planted_example(), cuts=((Q(1,3),),(),(),()))
    wall = multiscale_ledger(planted_example(), cuts=((Q(1,2),),(),(),()))
    cancel = multiscale_ledger(cancellation_example(), cuts=((),(),(Q(1,2),),()))
    return {
        "nonzero": {"degree":good.parent_degree,"child_flux":good.sum_child_flux,
                    "children":good.children,"internal_pairs":good.internal_wall_pairs},
        "root_on_internal_cut": {"degree":wall.parent_degree,
                    "child_flux":wall.sum_child_flux,"wall_children":wall.boundary_root_children},
        "opposite_sign_pair": {"degree":cancel.parent_degree,
                    "child_flux":cancel.sum_child_flux,"children":cancel.children},
    }


if __name__ == '__main__':
    import json
    print(json.dumps(audit(), indent=2, default=str))