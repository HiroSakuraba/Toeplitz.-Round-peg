"""v215: boundary-certified Brouwer degree on four PL source arcs.

The degree is for the four-dimensional square residual restricted to four
physically marked planar graph arcs and a fixed integer lift tuple.
Exact Fraction arithmetic; no sampling, no Floer-homology calculation.
A nonzero degree is sufficient for a real square ONLY under the
separation and boundary hypotheses. Does not prove the Toeplitz conjecture.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations, product
from typing import Sequence, Tuple


@dataclass(frozen=True)
class Point:
    x: Q
    y: Q


def point(x, y):
    return Point(Q(x), Q(y))


@dataclass(frozen=True)
class Arc:
    vertices: Tuple[Point, ...]

    def __post_init__(self):
        if len(self.vertices) < 2:
            raise ValueError("an arc needs at least two vertices")
        if any(a.x >= b.x for a, b in zip(self.vertices, self.vertices[1:])):
            raise ValueError("x coordinates must increase strictly")
        if any(not isinstance(p.x, Q) or not isinstance(p.y, Q)
               for p in self.vertices):
            raise ValueError("use Fraction coordinates")


def arc(*coords):
    return Arc(tuple(point(x, y) for x, y in coords))


def arcs_from_periodic_windows(graphs, windows):
    """Convert v212 PeriodicGraph windows into v215 source-arc coordinates.

    This is a convenience adapter. It retains the physical lift selection
    as a separate explicit argument to source_lift_degree().
    """
    if len(graphs)!=4 or len(windows)!=4:
        raise ValueError("four graphs and four windows required")
    arcs=[]
    for g,(lo,hi) in zip(graphs,windows):
        lo,hi=Q(lo),Q(hi)
        if not Q(0)<=lo<hi<=Q(1):
            raise ValueError("windows must lie inside a single period")
        pts=[point(lo,g.height(lo))]
        pts.extend(point(v.x,v.y) for v in g.vertices if lo<v.x<hi)
        pts.append(point(hi,g.height(hi)))
        arcs.append(Arc(tuple(pts)))
    return tuple(arcs)


@dataclass(frozen=True)
class Segment:
    a: Point
    b: Point

    def at(self, u: Q, shift: int):
        return Point(self.a.x + (self.b.x - self.a.x)*u + shift,
                     self.a.y + (self.b.y - self.a.y)*u)


def _residual(p, orientation):
    a, b, c, d = p
    vx, vy = b.x-a.x, b.y-a.y
    s = orientation
    # c-b-s R90(b-a), d-a-s R90(b-a)
    return (c.x-b.x+s*vy, c.y-b.y-s*vx,
            d.x-a.x+s*vy, d.y-a.y-s*vx)


def _cell(segs, shifts, orientation):
    starts = [s.at(Q(0), k) for s, k in zip(segs, shifts)]
    constant = _residual(starts, orientation)
    cols = []
    for i in range(4):
        pts = list(starts)
        pts[i] = segs[i].at(Q(1), shifts[i])
        cols.append(tuple(v-c for v,c in zip(_residual(pts,orientation),constant)))
    a = tuple(tuple(col[j] for col in cols) for j in range(4))
    b = tuple(-v for v in constant)
    return a,b


def _rank(a):
    if not a:
        return 0
    m = [list(map(Q, row)) for row in a]
    n = len(m[0]); r = 0
    for col in range(n):
        piv = next((i for i in range(r, len(m)) if m[i][col]), None)
        if piv is None:
            continue
        m[r],m[piv] = m[piv],m[r]
        v=m[r][col]
        m[r]=[x/v for x in m[r]]
        for i in range(r+1,len(m)):
            if m[i][col]:
                v=m[i][col]
                m[i]=[x-v*y for x,y in zip(m[i],m[r])]
        r+=1
        if r==len(m):
            break
    return r


def _solve(a, b):
    n=len(b)
    mat=[list(map(Q,row))+[Q(y)] for row,y in zip(a,b)]
    if len(mat)!=n or any(len(row)!=n+1 for row in mat):
        raise ValueError("square system expected")
    for j in range(n):
        piv=next((i for i in range(j,n) if mat[i][j]), None)
        if piv is None:
            return None
        mat[j],mat[piv]=mat[piv],mat[j]
        v=mat[j][j]
        mat[j]=[x/v for x in mat[j]]
        for i in range(n):
            if i!=j and mat[i][j]:
                v=mat[i][j]
                mat[i]=[x-v*y for x,y in zip(mat[i],mat[j])]
    return tuple(row[-1] for row in mat)


def _det(a):
    n=len(a)
    m=[list(map(Q,row)) for row in a]
    sign=Q(1)
    for j in range(n):
        piv=next((i for i in range(j,n) if m[i][j]),None)
        if piv is None:
            return Q(0)
        if piv!=j:
            m[j],m[piv]=m[piv],m[j]
            sign=-sign
        v=m[j][j]
        sign*=v
        for i in range(j+1,n):
            t=m[i][j]/v
            for k in range(j+1,n):
                m[i][k]-=t*m[j][k]
            m[i][j]=Q(0)
    return sign


def _feasible_box(a,b):
    """Exact feasibility of A*u=b and 0<=u<=1, including singular families."""
    rank=_rank(a)
    rows=[]; rhs=[]
    for row,target in zip(a,b):
        if _rank(rows+[row])>len(rows):
            rows.append(list(row))
            rhs.append(target)
    for inds in combinations(range(4),4-rank):
        for vals in product((Q(0),Q(1)),repeat=4-rank):
            constraints=[list(row) for row in rows]
            targets=list(rhs)
            for k,v in zip(inds,vals):
                unit=[Q(0)]*4;unit[k]=Q(1)
                constraints.append(unit);targets.append(v)
            sol=_solve(constraints,targets)
            if sol is not None and all(0<=u<=1 for u in sol) and all(
                sum(x*y for x,y in zip(row,sol))==v for row,v in zip(a,b)):
                return True
    return False


def _segments(arcs):
    return tuple(tuple(Segment(a,b) for a,b in zip(arc.vertices,arc.vertices[1:]))
                 for arc in arcs)


def _boundary_cells(segs,arcs):
    for i in range(4):
        for side in (0,1):
            p=arcs[i].vertices[0 if side==0 else -1]
            fixed=Segment(p,p)
            collections=[(fixed,) if j==i else segs[j] for j in range(4)]
            for choice in product(*collections):
                yield i,choice


def _residual_inf_lower_on_boundary_cell(a,b,fixed,max_nodes):
    """Adaptive exact dyadic interval certificate on 3D affine facet.

    No boundary root must have been verified first. For an affine residual,
    each coordinate's interval enclosure is exact on a parameter box.
    If one coordinate excludes zero, it certifies a positive l_infty
    residual bound. Subdivision terminates on a compact zero-free facet.
    """
    free=tuple(j for j in range(4) if j!=fixed)
    stack=[tuple((Q(0),Q(1)) for _ in free)]
    nodes=0
    bound=None
    while stack:
        box=stack.pop()
        nodes+=1
        if nodes>max_nodes:
            raise CertificateInconclusive("dyadic boundary certificate budget exceeded")
        estimates=[]
        for row,target in zip(a,b):
            lo=-target; hi=-target
            for col,(u,v) in zip(free,box):
                lo+=min(row[col]*u,row[col]*v)
                hi+=max(row[col]*u,row[col]*v)
            if lo>0:
                estimates.append(lo)
            if hi<0:
                estimates.append(-hi)
        if estimates:
            value=max(estimates)
            bound=value if bound is None else min(bound,value)
            continue
        sizes=[(v-u)*max(abs(row[j]) for row in a) for j,(u,v) in zip(free,box)]
        j=max(range(3),key=lambda k:sizes[k])
        u,v=box[j]
        if sizes[j]==0 or u==v:
            raise CertificateInconclusive("no separation found on a boundary cell")
        mid=(u+v)/2
        lower=list(box);upper=list(box)
        lower[j]=(u,mid);upper[j]=(mid,v)
        stack.append(tuple(lower));stack.append(tuple(upper))
    return bound,nodes


class BoundarySquare(ValueError):
    """The square-residual map has a zero on a marked source boundary."""


class CertificateInconclusive(RuntimeError):
    """An exact certificate could not be finished under its resource budget."""


@dataclass(frozen=True)
class RegularPreimage:
    indices: Tuple[int,int,int,int]
    parameters: Tuple[Q,Q,Q,Q]
    determinant: Q


@dataclass(frozen=True)
class DegreeCertificate:
    degree: int
    preimages_of_regular_target: Tuple[RegularPreimage,...]
    regular_target: Tuple[Q,Q,Q,Q]
    boundary_margin: Q
    opposite_vertical_gap: Q
    candidate_attempts: int
    checked_boundary_cells: int
    subdivision_nodes: int

    @property
    def vertical_perturbation_bound(self):
        return min(self.boundary_margin/2,self.opposite_vertical_gap/2)


def source_lift_degree(arcs: Sequence[Arc], shifts=(0,0,0,0),
                       orientation=1, max_subdivision_nodes=200000,
                       max_target_attempts=128):
    """Exact degree for a four-arc PL residual, including knot roots.

    Returns a conservative positive margin and signed count at a small
    *regular target*. The target preimages are NOT exact square witnesses;
    a nonzero degree ensures at least one true square F=0 exists.

    Raises BoundarySquare if a zero touches any external source endpoint;
    raises CertificateInconclusive if resource caps prevent certification.
    It deliberately requires a positive y-gap between opposite arcs 0,2.
    """
    if len(arcs)!=4 or any(not isinstance(a,Arc) for a in arcs):
        raise ValueError("four Arc objects required")
    if (len(shifts)!=4 or shifts[0]!=0 or
            any(type(k) is not int for k in shifts)):
        raise ValueError("relative lifts must be four integers starting with zero")
    if orientation not in (-1,1):
        raise ValueError("orientation must be +/-1")
    ymin0=min(v.y for v in arcs[0].vertices);ymax0=max(v.y for v in arcs[0].vertices)
    ymin2=min(v.y for v in arcs[2].vertices);ymax2=max(v.y for v in arcs[2].vertices)
    gap=max(Q(0),ymin2-ymax0,ymin0-ymax2)
    if gap<=0:
        raise ValueError("opposite arc vertical ranges must be separated")
    segs=_segments(arcs)
    # 1. Detect *any* boundary root, even those lying on internal PL knots.
    #    The facet equation has three free parameters but we solve it
    #    exactly as a rank-deficient four-variable box system.
    margins=[]; checked=nodes=0
    boundary=list(_boundary_cells(segs,arcs))
    for i,choice in boundary:
        a,b=_cell(choice,shifts,orientation)
        if _feasible_box(a,b):
            raise BoundarySquare("square touches a marked source endpoint")
    for i,choice in boundary:
        a,b=_cell(choice,shifts,orientation)
        lower,count=_residual_inf_lower_on_boundary_cell(
            a,b,i,max_subdivision_nodes-nodes)
        checked+=1;nodes+=count
        if lower is None or lower<=0:
            raise CertificateInconclusive("no positive boundary margin")
        margins.append(lower)
    m=min(margins)
    # 2. A moment-curve rational target avoids all internal chart faces
    #    and all singular images after finitely many attempts. Because
    #    ||q||_inf < m, deg(F,Q,0) == deg(F,Q,q).
    cells=[]
    for indices in product(*(range(len(parts)) for parts in segs)):
        a,b=_cell([segs[j][v] for j,v in enumerate(indices)],shifts,orientation)
        cells.append((tuple(indices),a,b,_det(a)))
    for n in range(2,max_target_attempts+2):
        delta=Q(1,2**n)
        if delta>=m/2:
            continue
        q=tuple(delta**(j+1) for j in range(4))
        preimages=[]
        regular=True
        for indices,a,b,det in cells:
            target=tuple(y+z for y,z in zip(b,q))
            if det:
                sol=_solve(a,target)
                assert sol is not None
                if any(u in (0,1) for u in sol):
                    regular=False;break
                if all(0<u<1 for u in sol):
                    preimages.append(RegularPreimage(indices,sol,det))
            else:
                augmented=[tuple(row)+(target[j],) for j,row in enumerate(a)]
                if _rank(augmented)==_rank(a):
                    regular=False;break
        if regular:
            deg=sum(1 if z.determinant>0 else -1 for z in preimages)
            return DegreeCertificate(deg,tuple(preimages),q,m,gap,n-1,
                                     checked,nodes)
    raise CertificateInconclusive("could not find regular target in budget")


if __name__=="__main__":
    # Four short transverse lines with a planted square.
    a=point(Q(1,2),0);e=point(Q(1,8),Q(1,8))
    r=point(-e.y,e.x)
    p=(a,point(a.x+e.x,a.y+e.y),
       point(a.x+e.x+r.x,a.y+e.y+r.y),point(a.x+r.x,a.y+r.y))
    slopes=(Q(-2),Q(3),Q(7,2),Q(-1,3))
    w=Q(1,32)
    arcs=tuple(arc((v.x-w,v.y-s*w),(v.x+w,v.y+s*w)) for v,s in zip(p,slopes))
    c=source_lift_degree(arcs)
    print('degree:',c.degree)
    print('boundary margin:',c.boundary_margin)
    print('opposite vertical gap:',c.opposite_vertical_gap)
    print('regular target:',c.regular_target)
    print('target preimages:',len(c.preimages_of_regular_target))
    print('strict vertical perturbation:',c.vertical_perturbation_bound)
