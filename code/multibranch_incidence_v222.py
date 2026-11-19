"""v222 exact multi-chart filtered square-incidence audit on moving polygonal arcs.

Each marked arc is a rigidly translated rational polygonal chain. Product
cells have a full-rank constant affine square Jacobian. Their one-dimensional
zero branches glue across internal PL knots. The ledger pairs shared-knot
endpoints *by actual time and ordered physical square*, so internal knot
crossings and two-branch PL folds contribute zero exposed boundary parity.

This is a conditional finite PL result. It does not claim that generic smooth
births/deaths are represented by this restricted model, that the four arcs
come from one Jordan curve, or that all geometric square orbits are counted.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product, combinations
from collections import defaultdict
from typing import Tuple

from cutoff_incidence_v221 import (
    Point, MovingSegment, LinearSquareBranch, NonGenericIncidence,
    pt, add, sub, scale, dot, residual, linear_square_branch, solve,
    _crossings_of_quadratic, polynomial_value,
)

@dataclass(frozen=True)
class MovingArc:
    vertices: Tuple[Point, ...]
    velocity: Point

    def __post_init__(self):
        if len(self.vertices)<2 or any(not isinstance(v,Point) for v in self.vertices):
            raise ValueError('each arc needs two or more exact vertices')
        if not isinstance(self.velocity,Point):
            raise ValueError('velocity must be an exact Point')
        if any(type(q.x) is not Q or type(q.y) is not Q
               for q in self.vertices+(self.velocity,)):
            raise ValueError('rational Fraction coordinates required')
        if any(a==b for a,b in zip(self.vertices,self.vertices[1:])):
            raise ValueError('zero-length segment in moving arc')


def moving_arc(coords, velocity=(0,0)):
    return MovingArc(tuple(pt(*xy) for xy in coords),pt(*velocity))


def moving_segments(arc, deck_shift=0):
    return tuple(MovingSegment(pt(p.x+deck_shift,p.y),pt(q.x+deck_shift,q.y),arc.velocity)
                 for p,q in zip(arc.vertices,arc.vertices[1:]))


def _physical_points_key(points):
    return tuple((p.x,p.y) for p in points)


def _branch_admissible(branch,t,cutoff):
    return all(0<u<1 for u in branch.parameters(t)) and branch.side_squared(t)>cutoff


def _active_closed_time(branch,a,b):
    """Intersection of four [0,1] affine-parameter strips and [a,b]."""
    left,right=a,b
    for p,v in zip(branch.params0,branch.params1):
        if v==0:
            if p in (0,1):
                # Coincidence along an entire cell boundary would need
                # special gluing and is not a generic isolated branch.
                raise NonGenericIncidence('persistent internal/source-face root')
            if not 0<p<1:return None
        else:
            x,y=sorted((-p/v,(1-p)/v))
            left=max(left,x);right=min(right,y)
            if left>right:return None
    return left,right


def _rank(rows):
    if not rows:return 0
    a=[list(map(Q,row)) for row in rows]
    r=0
    for j in range(len(a[0])):
        k=next((i for i in range(r,len(a)) if a[i][j]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]
        p=a[r][j];a[r]=[x/p for x in a[r]]
        for i in range(r+1,len(a)):
            if a[i][j]:
                z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a):break
    return r


def _box_feasible(A,b,n):
    """Exact rational feasibility for A x=b inside unit n-cube."""
    if len(A)!=len(b) or any(len(row)!=n for row in A):
        raise ValueError('invalid affine system')
    rank=_rank(A)
    rows=[];rhs=[]
    for row,target in zip(A,b):
        if _rank(rows+[row])>len(rows):
            rows.append(list(row));rhs.append(Q(target))
    for js in combinations(range(n),n-rank):
        for values in product((Q(0),Q(1)),repeat=n-rank):
            mat=[list(row) for row in rows]
            out=list(rhs)
            for j,v in zip(js,values):
                unit=[Q(0)]*n;unit[j]=Q(1)
                mat.append(unit);out.append(v)
            sol=solve(mat,out)
            if (sol is not None and all(0<=u<=1 for u in sol) and
                all(sum(c*u for c,u in zip(row,sol))==val
                    for row,val in zip(A,b))):
                return True
    return False


def _singular_chart_can_hold_square(segs,orientation,a,b):
    """Exact check of F(t,u)=0 for u in [0,1]^4 and t in [a,b].

    Use fifth normalized time coordinate. This can reject a singular chart
    only when it contains *no square at all*; feasible singular charts
    remain explicitly inconclusive.
    """
    starts=[s.at(a,0) for s in segs]
    const=residual(starts,orientation)
    cols=[]
    for j in range(4):
        moved=list(starts);moved[j]=segs[j].at(a,1)
        nxt=residual(moved,orientation)
        cols.append(tuple(x-y for x,y in zip(nxt,const)))
    finish=residual([s.at(b,0) for s in segs],orientation)
    cols.append(tuple(x-y for x,y in zip(finish,const)))
    A=tuple(tuple(col[i] for col in cols) for i in range(4))
    return _box_feasible(A,tuple(-x for x in const),5)


@dataclass(frozen=True)
class ChartBranch:
    indices: Tuple[int,int,int,int]
    branch: LinearSquareBranch
    active_start: Q
    active_end: Q
    cutoff_crossings: int


@dataclass(frozen=True)
class IncidenceEvent:
    time: Q
    kind: str
    square: Tuple[Point,...]
    chart_indices: Tuple[Tuple[int,int,int,int],...]


@dataclass(frozen=True)
class MultiBranchLedger:
    start: Q
    end: Q
    cutoff_squared: Q
    start_count: int
    end_count: int
    source_events: Tuple[IncidenceEvent,...]
    internal_events: Tuple[IncidenceEvent,...]
    cutoff_crossings: int
    chart_branches: Tuple[ChartBranch,...]
    orientation: int
    shifts: Tuple[int,int,int,int]

    @property
    def balances(self):
        return (self.start_count+self.end_count+self.cutoff_crossings+
                len(self.source_events))%2==0

    @property
    def parity(self):
        return (self.start_count%2,self.end_count%2)


def multibranch_ledger(arcs, start, end, cutoff_squared, orientation=1,
                       shifts=(0,0,0,0)):
    """Complete finite generic ledger over all four-arc segment products.

    Assumptions enforced: every product-cell square matrix invertible;
    one-dimensional feasible branch strata; each internal endpoint has
    exactly one paired adjacent chart with same square/time; no simultaneous
    cutoff/source events; time slices generic. Each crossing of the positive
    size cutoff is simple, with exact (possibly irrational) root counting.
    A pair meeting at a PL knot may be born/die, or a root may just cross the
    knot. Both have two chart half-branches and cancel internally.
    """
    a,b,k=Q(start),Q(end),Q(cutoff_squared)
    if not a<b or k<=0:raise ValueError('ordered time interval and positive cutoff required')
    if len(arcs)!=4 or any(not isinstance(x,MovingArc) for x in arcs):
        raise ValueError('four moving polygonal arcs required')
    if orientation not in (-1,1):raise ValueError('orientation +/-1 required')
    if (len(shifts)!=4 or shifts[0]!=0 or
        any(type(t) is not int for t in shifts)):
        raise ValueError('four integer physical lifts with first zero required')
    lists=tuple(moving_segments(arc,shift) for arc,shift in zip(arcs,shifts))
    segments_per_arc=tuple(len(x) for x in lists)
    live=[]
    raw=defaultdict(list)
    ends=[[],[]]
    cutoffs=0
    for indices in product(*(range(n) for n in segments_per_arc)):
        segs=tuple(lists[i][j] for i,j in enumerate(indices))
        try:
            branch=linear_square_branch(segs,orientation)
        except NonGenericIncidence as exc:
            if 'singular residual Jacobian' not in str(exc):raise
            if _singular_chart_can_hold_square(segs,orientation,a,b):
                raise NonGenericIncidence('feasible singular chart: nonisolated squares') from exc
            continue  # exact infeasibility certificate, no dropped roots
        interval=_active_closed_time(branch,a,b)
        if interval is None:continue
        left,right=interval
        if left==right:
            if branch.side_squared(left)>=k:
                raise NonGenericIncidence('isolated chart-corner root, not 1D generic')
            continue
        poly=branch.cutoff_polynomial(k)
        for t,side in ((left,0),(right,1)):
            if polynomial_value(poly,t)==0:
                raise NonGenericIncidence('cutoff coincides with a chart/time endpoint')
            if t in (a,b):
                if branch.side_squared(t)>k:
                    if any(u in (0,1) for u in branch.parameters(t)):
                        raise NonGenericIncidence('endpoint time slice hits a source/knot wall')
                    if _branch_admissible(branch,t,k):
                        global_u=tuple(indices[i]+u for i,u in enumerate(branch.parameters(t)))
                        ends[int(t==b)].append((indices,global_u))
                continue
            if branch.side_squared(t)<=k:continue
            u=branch.parameters(t)
            boundary_dims=tuple(i for i,v in enumerate(u) if v in (0,1))
            if len(boundary_dims)!=1:
                raise NonGenericIncidence('simultaneous multiple chart-face crossing')
            j=boundary_dims[0]
            if not all(0<x<1 for i,x in enumerate(u) if i!=j):
                raise NonGenericIncidence('chart corner involving another arc')
            chart_end=(u[j]==0 and indices[j]==0) or (
                u[j]==1 and indices[j]==segments_per_arc[j]-1)
            kind='source' if chart_end else 'internal'
            global_u=tuple(indices[i]+v for i,v in enumerate(u))
            key=(t,global_u)
            raw[key].append((kind,indices,j,int(u[j]),_physical_points_key(branch.square(t))))
        if polynomial_value(poly,left)==0 or polynomial_value(poly,right)==0:
            raise NonGenericIncidence('cutoff at active interval endpoint')
        n=_crossings_of_quadratic(poly,left,right)
        cutoffs+=n
        live.append(ChartBranch(indices,branch,left,right,n))
    for e in ends:
        # A root on different cells at a TIME boundary is non-generic.
        if len({key for _,key in e})!=len(e):
            raise NonGenericIncidence('time slice square duplicated across charts')
    source=[];internal=[]
    for (t,global_u),events in sorted(raw.items()):
        # Equal geometric points do not imply equal parameter configurations
        # on self-crossing arcs. Internal-knot gluing uses *global arc
        # parameters*, not just square vertices in the ambient plane.
        coords=events[0][4]
        if any(e[4]!=coords for e in events):
            raise AssertionError('same marked parameters yield different planar points')
        square=tuple(pt(*x) for x in coords)
        kinds={e[0] for e in events}
        if kinds=={'source'} and len(events)==1:
            source.append(IncidenceEvent(t,'source',square,(events[0][1],)))
        elif kinds=={'internal'} and len(events)==2:
            e0,e1=events
            i,j=e0[2],e1[2]
            if i!=j or {e0[3],e1[3]}!={0,1}:
                raise NonGenericIncidence('unpaired internal arc-knot endpoints')
            x,y=e0[1],e1[1]
            if (abs(x[i]-y[i])!=1 or any(x[z]!=y[z] for z in range(4) if z!=i)):
                raise NonGenericIncidence('internal endpoints not on adjacent charts')
            internal.append(IncidenceEvent(t,'internal',square,(x,y)))
        else:
            raise NonGenericIncidence('unmatched or multiple square boundary incidence')
    result=MultiBranchLedger(a,b,k,len(ends[0]),len(ends[1]),
                             tuple(source),tuple(internal),cutoffs,
                             tuple(live),orientation,tuple(shifts))
    if not result.balances:
        raise AssertionError('incidence identity failed: a chart endpoint was missed')
    return result


def canonical_square_orbit(points):
    """Geometric vertex set key, independent of D4 labeling.

    Rejects non-squares and collapsed configurations. An orbit counts once,
    not four cyclic labels (or eight labels including reversal).
    """
    if len(points)!=4 or len(set(points))!=4:
        raise ValueError('four distinct vertices required')
    if (residual(points,1)!=(0,0,0,0) and
        residual(points,-1)!=(0,0,0,0)):
        raise ValueError('ordered vertices do not form a square')
    return tuple(sorted((p.x,p.y) for p in points))


def dihedral_labels(points):
    """The eight geometric labels of one noncollapsed oriented square."""
    canonical_square_orbit(points)
    p=tuple(points)
    reverse=(p[0],p[3],p[2],p[1])
    return tuple(tuple(seq[i:]+seq[:i]) for seq in (p,reverse) for i in range(4))


def polygon_example_two_roots(delta=Q(1,100)):
    """Four independent marked arcs: one tent with TWO transverse squares.

    Translate the tent downward by 2*delta over t in [0,1]. At t=1/2
    its two roots coalesce at the internal knot, and then disappear.
    They contribute an even internal half-branch count and no exposed event.
    """
    d=Q(delta)
    if d<=0:raise ValueError('positive tent height required')
    h=Q(1,32)
    points=(pt(Q(1,2),0),pt(Q(5,8),Q(1,8)),
            pt(Q(1,2),Q(1,4)),pt(Q(3,8),Q(1,8)))
    slopes=(-Q(2),Q(3),Q(7,2),-Q(1,3))
    result=[]
    for i,(p,m) in enumerate(zip(points,slopes)):
        lo,hi=p.x-h,p.x+h
        if i==2:
            target=lambda x:Q(1,4)-2*(x-Q(1,2))
            coords=((lo,target(lo)-d),(p.x,target(p.x)+d),(hi,target(hi)-d))
            vel=(Q(0),-2*d)
        else:
            coords=((lo,p.y-m*h),(hi,p.y+m*h))
            vel=(Q(0),Q(0))
        result.append(moving_arc(coords,vel))
    return tuple(result)



def fold_model_root_count(t):
    """Abstract smooth fold F_t=(x²-t,y,z,w), not a square-arc model.

    In the full five-dimensional (t,x,y,z,w) space its zero locus is a
    smooth parabola, t=x². For negative t no roots; positive t has two
    regular roots with opposite Jacobian signs. This is a separate
    local topological model of pair birth/death, not a new claim about
    the four-arc translating-polygon specialization.
    """
    t=Q(t)
    return 0 if t<0 else (1 if t==0 else 2)

def report():
    arcs=polygon_example_two_roots()
    d=multibranch_ledger(arcs,0,1,Q(1,1024))
    example=dihedral_labels(d.chart_branches[0].branch.square(Q(0))) if d.chart_branches else ()
    return {'two_root_model':{'start':d.start_count,'end':d.end_count,
                              'source_events':len(d.source_events),
                              'internal_pairs':len(d.internal_events),
                              'size_crossings':d.cutoff_crossings,'balances':d.balances},
            'dihedral_orbit_labels':len(example),
            'scope':'translating PL arc products, not full Jordan curve or all-orbits cobordism'}

if __name__=='__main__':
    import json
    print(json.dumps(report(),sort_keys=True,indent=2))
