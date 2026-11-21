"""Exact size-wall census and perturbation clearance for the half-turn degree."""
from fractions import Fraction as Q
from itertools import product
import json

from polygon_orbit_census_v223 import (
    _affine_system, _box_vertices, _evaluate, _side_sq, _solve,
    add, mul, pt, sub, validate_jordan_polygon,
)
from rough_arc_winding_v312 import whole_curve_winding
from boundary_or_winding_v313 import dyadic_arc, scan_level, separation_squared, window_alternative


def profile():
    return [pt(-1,0),pt(Q(-1,8),0),pt(Q(-3,8),Q(1,2)),
            pt(Q(7,16),Q(17,32)),pt(Q(1,8),0),pt(1,0)]


def ear_census():
    points = profile()
    edges = list(zip(points,points[1:]))
    squares = {}
    for ids in product(range(len(edges)),repeat=4):
        if ids[0]!=min(ids) or max(ids.count(i) for i in ids)>2:
            continue
        cell = [edges[i] for i in ids]
        a,b = _affine_system(cell,1)
        tags = _solve(a,b)
        if tags is None:
            # Convexity of squared norm: zero on all vertices implies
            # zero throughout the feasible tag polytope.
            if any(_side_sq(_evaluate(cell,u)) for u in _box_vertices(a,b)):
                raise ValueError("nontrivial singular ear chart")
        elif all(0<=t<=1 for t in tags):
            corners = _evaluate(cell,tags)
            squared = _side_sq(corners)
            if squared:
                if min(tags)==0 or max(tags)==1:
                    raise ValueError("nontransverse ear contact")
                squares[tuple(sorted(corners))] = (squared,ids,tags)
    return sorted(squares.values())


def polygon(upper=Q(1),lower=Q(1)):
    def rim(scale):
        # Traverse the top side from right to left.
        return [pt(4,1)]+[add(pt(1,1),mul(p,scale)) for p in profile()[1:-1][::-1]]+[pt(-4,1)]
    return validate_jordan_polygon(tuple(rim(upper)+[mul(p,-1) for p in rim(lower)]))


def cut_crossings():
    scales = (Q(9,10),Q(1))
    target = Q(2,5)**2
    sizes = [s for s,_,_ in ear_census()]
    above = [sum(t*t*s>target for s in sizes) for t in scales]
    if len(sizes)!=2 or above!=[0,1] or any(t*t*s==target for t in scales for s in sizes):
        raise ArithmeticError("wrong cutoff census")
    return {"normalized_squared_sides":list(map(str,sizes)),
            "ear_counts_above_cutoff":above,"event_scale_squared":str(target/sizes[-1]),
            "symmetric_degrees":[1,1],"one_ear_moved_degrees":[1,0]}


def bounds(a,b):
    return (min(a.x,b.x),max(a.x,b.x),min(a.y,b.y),max(a.y,b.y))


def sub_range(a,b):
    return a[0]-b[1],a[1]-b[0]


def norm_bounds(x,y):
    return (sum(max(0,a,-b)**2 for a,b in (x,y)),
            sum(max(abs(a),abs(b))**2 for a,b in (x,y)))


def possible_contact(x,y,edges,pad):
    left,right = x[0]-pad,x[1]+pad
    low,high = y[0]-pad,y[1]+pad
    for (x0,x1,y0,y1),nx,ny,constant in edges:
        if x1<left or right<x0 or y1<low or high<y0:
            continue
        ax,bx = nx*left,nx*right
        ay,by = ny*low,ny*high
        if min(ax,bx)+min(ay,by)<=constant<=max(ax,bx)+max(ay,by):
            return True
    return False


def tube_certificate(vertices, side=Q(2,5), epsilon=Q(1,8192), max_depth=38, upper_side=None):
    """Exclude four vertices in an L-infinity tube at the specified side length."""
    upper_side = side if upper_side is None else upper_side
    if side<=0 or upper_side<side or epsilon<0 or max_depth<0:
        raise ValueError("invalid cutoff, tube width or depth")
    segments = list(zip(vertices,vertices[1:]+vertices[:1]))
    supports = []
    for a,b in segments:
        nx,ny = b.y-a.y,a.x-b.x
        supports.append((bounds(a,b),nx,ny,nx*a.x+ny*a.y))
    todo = [(p,q,0) for p,q in product(segments,repeat=2)]
    leaves = {"side":0,"first_completion":0,"second_completion":0}
    deepest = 0
    while todo:
        first,second,depth = todo.pop()
        deepest = max(depth,deepest)
        a,b = bounds(*first),bounds(*second)
        dx,dy = sub_range(b[:2],a[:2]),sub_range(b[2:],a[2:])
        e = 2*epsilon
        lo,hi = norm_bounds((dx[0]-e,dx[1]+e),(dy[0]-e,dy[1]+e))
        if hi<side*side or upper_side*upper_side<lo:
            leaves["side"]+=1
            continue
        for box,name in ((a,"first_completion"),(b,"second_completion")):
            x = sub_range(box[:2],dy)
            y = (box[2]+dx[0],box[3]+dx[1])
            if not possible_contact(x,y,supports,4*epsilon):
                leaves[name]+=1
                break
        else:
            if depth==max_depth:
                raise ValueError("unresolved cutoff box")
            split_first = max(a[1]-a[0],a[3]-a[2])>=max(b[1]-b[0],b[3]-b[2])
            start,end = first if split_first else second
            mid = mul(add(start,end),Q(1,2))
            for half in ((start,mid),(mid,end)):
                todo.append((half,second,depth+1) if split_first else (first,half,depth+1))
    return {"side_interval":[str(side),str(upper_side)],"tube_halfwidth":str(epsilon),"max_depth_used":deepest,
            "excluded_leaves":leaves}


def winding_witness():
    size,ids,tags = ear_census()[-1]
    points=profile()
    windows=[]
    for i,t in zip(ids[:2],tags[:2]):
        start,end = points[i:i+2]
        center = add(pt(1,1), add(start,mul(sub(end,start),t)))
        windows.append([add(center,mul(sub(end,start),h)) for h in (-Q(1,512),Q(1,512))])
    result = whole_curve_winding(polygon(),*windows)
    if abs(result["winding"])!=1:
        raise ArithmeticError("lost transverse positive-scale square")
    return dict(result,side_squared=str(size),
                windows=[[[str(p.x),str(p.y)] for p in arc] for arc in windows])


def grid_witnesses():
    p=polygon()
    j,k=(dyadic_arc(p,10,i) for i in (343,101))
    sep=separation_squared(j,k)
    local=window_alternative(p,j,k)
    coarse=scan_level(p,3,Q(7,5))
    if sep<=Q(4,25) or local.get("winding")!=1 or not coarse or coarse["outcome"]!="boundary square":
        raise ArithmeticError("missing dyadic witness")
    return [coarse,dict(local,level=10,cells=[343,101],separation_squared=str(sep))]


def report():
    return {"dyadic_witnesses":grid_witnesses(),"ear_motion":cut_crossings(),"cutoff_clearance":tube_certificate(polygon()),
            "large_gap":tube_certificate(polygon(),Q(7,5),Q(1,32),upper_side=Q(8,5)),
            "whole_curve_window":winding_witness()}


if __name__=="__main__":
    print(json.dumps(report(),indent=2))
