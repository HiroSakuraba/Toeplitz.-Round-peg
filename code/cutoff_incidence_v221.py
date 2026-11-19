"""v221: exact mod-2 size-cutoff incidence in one *transverse linear chart*.

This module treats four independently marked segments undergoing translations
linear in a homotopy parameter. Their square equations have ONE affine root
branch when the 4x4 residual Jacobian is invertible. It enumerates exits
through source-segment ends and crossings of a *positive* squared-side cutoff.
It does not count every square on a polygon or Jordan curve, construct a
relative Floer invariant, or establish a global square-parity theorem.

All event counts are exact rational decisions (quadratic roots can be
irrational; their count is obtained from signs and the exact vertex).
"""

from dataclasses import dataclass
from fractions import Fraction as Q
from typing import Tuple


@dataclass(frozen=True)
class Point:
    x: Q
    y: Q


def pt(x, y):
    return Point(Q(x), Q(y))


def add(p, q):
    return pt(p.x + q.x, p.y + q.y)


def sub(p, q):
    return pt(p.x - q.x, p.y - q.y)


def scale(t, p):
    return pt(Q(t)*p.x, Q(t)*p.y)


def dot(p, q):
    return p.x*q.x+p.y*q.y


def rot(p, orientation):
    return pt(-orientation*p.y, orientation*p.x)


def residual(points, orientation=1):
    if len(points) != 4 or orientation not in (1, -1):
        raise ValueError("four points and an orientation +/-1 required")
    p0,p1,p2,p3=points
    r=rot(sub(p1,p0), orientation)
    a=sub(sub(p2,p1),r)
    b=sub(sub(p3,p0),r)
    return (a.x,a.y,b.x,b.y)


@dataclass(frozen=True)
class MovingSegment:
    start: Point
    end: Point
    velocity: Point

    def at(self, t, u):
        return add(add(self.start,scale(t,self.velocity)),
                   scale(u,sub(self.end,self.start)))


def segment_through(point, direction, half_width, velocity=None):
    """Segment centered at point, with fixed direction during translation."""
    v=velocity if velocity is not None else pt(0,0)
    h=Q(half_width)
    if h<=0 or dot(direction,direction)==0:
        raise ValueError("nontrivial segment and positive half-width required")
    return MovingSegment(sub(point,scale(h,direction)),
                         add(point,scale(h,direction)), v)


def solve(A, b):
    """Exact 4x4 rational Gaussian elimination; None for singular A."""
    n=len(b)
    if len(A)!=n or any(len(row)!=n for row in A):
        raise ValueError("square system required")
    m=[list(map(Q,row))+[Q(y)] for row,y in zip(A,b)]
    for j in range(n):
        pivot=next((k for k in range(j,n) if m[k][j]),None)
        if pivot is None: return None
        m[j],m[pivot]=m[pivot],m[j]
        fac=m[j][j]
        m[j]=[v/fac for v in m[j]]
        for k in range(n):
            if k!=j and m[k][j]:
                fac=m[k][j]
                m[k]=[v-fac*w for v,w in zip(m[k],m[j])]
    return tuple(r[-1] for r in m)


class NonGenericIncidence(ValueError):
    """Branch, cutoff, or boundary crossing is not sufficiently transverse."""


@dataclass(frozen=True)
class LinearSquareBranch:
    segments: Tuple[MovingSegment, ...]
    orientation: int
    params0: Tuple[Q, ...]
    params1: Tuple[Q, ...]
    edge0: Point
    edge1: Point

    def parameters(self, t):
        return tuple(a+Q(t)*b for a,b in zip(self.params0,self.params1))

    def square(self, t):
        return tuple(s.at(t,u) for s,u in zip(self.segments,self.parameters(t)))

    def side_squared(self, t):
        e=add(self.edge0,scale(t,self.edge1))
        return dot(e,e)

    def cutoff_polynomial(self, cutoff_squared):
        k=Q(cutoff_squared)
        return (dot(self.edge1,self.edge1),
                2*dot(self.edge0,self.edge1),
                dot(self.edge0,self.edge0)-k)


def linear_square_branch(segments, orientation=1):
    """Solve the unique affine square root branch over all real t.

    A(t) is constant because each segment moves by a uniform translation.
    Parameter and square-edge coordinates are affine in time. The returned
    solution need not lie within the four segment interiors for all times.
    """
    if len(segments)!=4 or any(not isinstance(s,MovingSegment) for s in segments):
        raise ValueError("four moving segments required")
    if orientation not in (1,-1):raise ValueError("orientation must be +/-1")
    starts=[s.start for s in segments]
    base=residual(starts,orientation)
    columns=[]
    for j in range(4):
        q=starts.copy();q[j]=segments[j].end
        r=residual(q,orientation)
        columns.append(tuple(a-b for a,b in zip(r,base)))
    A=tuple(tuple(col[i] for col in columns) for i in range(4))
    b0=tuple(-v for v in base)
    vel=residual([s.velocity for s in segments],orientation)
    b1=tuple(-v for v in vel)
    u0=solve(A,b0);u1=solve(A,b1)
    if u0 is None or u1 is None:
        raise NonGenericIncidence("singular residual Jacobian: no unique branch")
    branch=LinearSquareBranch(tuple(segments),orientation,u0,u1,pt(0,0),pt(0,0))
    e0=sub(branch.square(0)[1],branch.square(0)[0])
    e1=sub(sub(branch.square(1)[1],branch.square(1)[0]),e0)
    result=LinearSquareBranch(tuple(segments),orientation,u0,u1,e0,e1)
    for t in (Q(0),Q(1,3),Q(1)):
        if residual(result.square(t),orientation)!=(0,0,0,0):
            raise AssertionError("linear branch algebra inconsistent")
    return result


def polynomial_value(poly,t):
    A,B,C=poly;t=Q(t)
    return (A*t+B)*t+C


def _crossings_of_quadratic(poly,lo,hi):
    """Count *simple* roots in (lo,hi) over Q, without extracting radicals.

    Endpoints must not be roots. A tangency at the interior critical time
    is non-generic and rejected even though it does not change parity.
    """
    lo,hi=Q(lo),Q(hi)
    if not lo<hi:raise ValueError("ordered interval required")
    A,B,C=poly
    if A==0 and B==0 and C==0:
        raise NonGenericIncidence("entire branch lies on cutoff")
    va,vb=polynomial_value(poly,lo),polynomial_value(poly,hi)
    if va==0 or vb==0:
        raise NonGenericIncidence("cutoff root coincides with a time/source boundary")
    if A==0:
        if B==0:
            if C==0:raise NonGenericIncidence("entire branch lies on cutoff")
            return 0
        return int(lo< -C/B <hi)
    vertex=-B/(2*A)
    if lo<vertex<hi:
        vv=polynomial_value(poly,vertex)
        if vv==0:
            raise NonGenericIncidence("tangential cutoff contact, not a regular crossing")
        return int(va*vv<0)+int(vv*vb<0)
    return int(va*vb<0)


@dataclass(frozen=True)
class IncidenceLedger:
    start: Q
    end: Q
    cutoff_squared: Q
    start_count: int
    end_count: int
    source_wall_times: Tuple[Q,...]
    cutoff_crossings: int
    active_cutoff_intervals: Tuple[Tuple[Q,Q,int],...]

    @property
    def source_events(self):
        return len(self.source_wall_times)

    @property
    def event_parity(self):
        return (self.source_events+self.cutoff_crossings)%2

    @property
    def slice_parity_change(self):
        return (self.start_count+self.end_count)%2

    @property
    def balances(self):
        return self.event_parity==self.slice_parity_change


def _eligible(branch,t,cutoff):
    return (all(0<u<1 for u in branch.parameters(t)) and
            branch.side_squared(t)>cutoff)


def incidence_ledger(branch,start,end,cutoff_squared):
    """Certified mod-2 boundary accounting for one regular square branch.

    Under the explicit genericity conditions, roots on the cutoff and on
    the four source-parameter walls are exactly the boundary of the 1D
    eligible branch locus. The two slice counts differ by their number mod2.
    This is NOT an all-curve mod-2 orbit count.
    """
    a,b,k=Q(start),Q(end),Q(cutoff_squared)
    if not a<b or k<=0:raise ValueError("time interval and positive cutoff needed")
    poly=branch.cutoff_polynomial(k)
    candidates={}
    for i,(u,v) in enumerate(zip(branch.params0,branch.params1)):
        if not v:
            if u in (0,1):
                raise NonGenericIncidence("branch lives on a source-segment boundary")
            continue
        for side in (0,1):
            t=(Q(side)-u)/v
            if a<t<b:candidates.setdefault(t,set()).add(i)
    for t in (a,b):
        if polynomial_value(poly,t)==0 or any(u in (0,1) for u in branch.parameters(t)):
            raise NonGenericIncidence("time slice meets a filter boundary")
    wall_times=[]
    for t,indices in sorted(candidates.items()):
        side_value=polynomial_value(poly,t)
        if side_value==0:
            raise NonGenericIncidence("simultaneous cutoff/source-wall event")
        if side_value<0:
            continue  # outside size filter at this source boundary
        # Only a wall with ALL other coordinates interior contributes.
        if len(indices)>1:
            if all(0<=u<=1 for u in branch.parameters(t)):
                raise NonGenericIncidence("multiple simultaneous active source walls")
            continue
        i=next(iter(indices))
        if all(0<value<1 for j,value in enumerate(branch.parameters(t)) if j!=i):
            wall_times.append(t)
    all_cuts=(a,)+tuple(sorted(candidates))+(b,)
    crossing_total=0;active_intervals=[]
    for lo,hi in zip(all_cuts,all_cuts[1:]):
        middle=(lo+hi)/2
        if all(0<u<1 for u in branch.parameters(middle)):
            n=_crossings_of_quadratic(poly,lo,hi)
            crossing_total+=n
            active_intervals.append((lo,hi,n))
    led=IncidenceLedger(a,b,k,int(_eligible(branch,a,k)),
                        int(_eligible(branch,b,k)),tuple(wall_times),
                        crossing_total,tuple(active_intervals))
    if not led.balances:
        raise AssertionError("mod-2 incidence law failed; inspect genericity hypothesis")
    return led


def planted_square_family(edge0,edge1,half_width=Q(1,1024)):
    """Four transverse segments centered on the vertices of a moving square.

    All four chosen source arcs are *independent*, not extracted from one
    Jordan curve. This intentionally tests only the local filtered ledger.
    """
    base=pt(Q(1,2),0)
    def corners(e):
        r=rot(e,1)
        return (base,add(base,e),add(add(base,e),r),add(base,r))
    p0,p1=corners(edge0),corners(add(edge0,edge1))
    dirs=(pt(1,-2),pt(1,3),pt(1,Q(7,2)),pt(1,-Q(1,3)))
    segs=tuple(segment_through(p,d,half_width,sub(q,p))
               for p,q,d in zip(p0,p1,dirs))
    return linear_square_branch(segs)


def shrinking_notch_tagged_count(cutoff_squared):
    """Count only the explicitly designated notch-corner square per index.

    This is NOT an exhaustive count of square orbits on the Jordan curve.
    Each fixed positive threshold admits finitely many named notches;
    passing one notch side² flips the tagged parity.
    """
    k=Q(cutoff_squared)
    if k<=0:raise ValueError("positive cutoff required; zero counts infinitely many")
    n=1
    while Q(1,16**(n+1))>=k:
        n+=1
    return n-1


def report():
    branch=planted_square_family(pt(0,0),pt(Q(1,8),Q(1,8)))
    cut=incidence_ledger(branch,Q(1,4),Q(3,4),Q(1,128))
    motion=planted_square_family(pt(Q(1,8),Q(1,8)),pt(0,0))
    moved=list(motion.segments)
    moved[1]=MovingSegment(moved[1].start,moved[1].end,pt(0,1))
    source=incidence_ledger(linear_square_branch(moved),Q(0),Q(1,5),Q(1,1024))
    return {
        "scaling_cutoff_start_end":(cut.start_count,cut.end_count),
        "scaling_cutoff_crossings":cut.cutoff_crossings,
        "scaling_source_wall_events":cut.source_events,
        "scaling_balances":cut.balances,
        "source_exit_start_end":(source.start_count,source.end_count),
        "source_exit_times":tuple(map(str,source.source_wall_times)),
        "source_exit_cutoff_crossings":source.cutoff_crossings,
        "source_balances":source.balances,
        "named_notches_cutoff_at_n5":shrinking_notch_tagged_count(Q(1,16**6)),
        "scope":"one transverse branch and named notch squares only; no all-Jordan parity invariant",
    }


if __name__=="__main__":
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
