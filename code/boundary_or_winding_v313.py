"""Exact side-window perimeter test, including tied crossings and collinear overlaps."""
from fractions import Fraction as Q
from itertools import product

from polygon_orbit_census_v223 import (
    _intersect, _on_seg, add, cross, dot, mul, pt, sub, validate_jordan_polygon,
)
from rough_arc_winding_v312 import complete, edges, validate_arcs


def lerp(a,b,t):
    return add(a,mul(sub(b,a),t))


def contact_interval(a,b,c,d):
    """Closed time interval on a moving point [a,b] meeting the segment [c,d]."""
    v,w = sub(b,a),sub(d,c)
    den = cross(v,w)
    diff=sub(c,a)
    if den:
        s,t = cross(diff,w)/den,cross(diff,v)/den
        return (s,s) if 0<=s<=1 and 0<=t<=1 else None
    vv=dot(v,v)
    if not vv:
        return (Q(0),Q(1)) if _on_seg(a,c,d) else None
    if cross(diff,v):
        return None
    lo,hi=sorted((dot(diff,v)/vv,dot(sub(d,a),v)/vv))
    lo,hi=max(Q(0),lo),min(Q(1),hi)
    return (lo,hi) if lo<=hi else None


def point_sign(p,boundary):
    inside=False
    for a,b in boundary:
        if _on_seg(p,a,b):
            return 0
        if (a.y>p.y)!=(b.y>p.y) and p.x<a.x+(p.y-a.y)*(b.x-a.x)/(b.y-a.y):
            inside=not inside
    return 1 if inside else -1


def window_alternative(polygon,j,k):
    polygon=validate_jordan_polygon(tuple(polygon))
    boundary=list(zip(polygon,polygon[1:]+polygon[:1]))
    validate_arcs((j,k))
    if any(not any(_on_seg(a,c,d) and _on_seg(b,c,d) for c,d in boundary)
           for arc in (j,k) for a,b in edges(arc)):
        raise ValueError("split the marked subarcs at polygon vertices")
    faces=([(q,k[0]) for q in j],[(j[-1],b) for b in k],
           [(q,k[-1]) for q in reversed(j)],[(j[0],b) for b in reversed(k)])
    walk=[]
    for face in faces:
        for (q0,b0),(q1,b1) in edges(face):
            p0,d0=complete(q0,b0)
            p1,d1=complete(q1,b1)
            knots={Q(0),Q(1)}
            for a,b in ((p0,p1),(d0,d1)):
                for c,d in boundary:
                    interval=contact_interval(a,b,c,d)
                    if interval is not None:
                        knots.update(interval)
            ordered=sorted(knots)
            probes=sorted(knots|{(x+y)/2 for x,y in zip(ordered,ordered[1:])})
            for t in probes:
                pair=pt(point_sign(lerp(p0,p1,t),boundary),point_sign(lerp(d0,d1,t),boundary))
                if pair==pt(0,0):
                    q,b=lerp(q0,q1,t),lerp(b0,b1,t)
                    side=sub(b,q)
                    return {"outcome":"boundary square","side_squared":str(dot(side,side)),
                            "ordered_side":[[str(z.x),str(z.y)] for z in (q,b)]}
                walk.append(pair)
    winding=0
    for a,b in edges(walk+walk[:1]):
        if a.x*b.x<=0 and a.y*b.y<=0:
            raise ArithmeticError("sign cells failed to share an open half-plane")
        det=cross(a,b)
        if a.y<=0<b.y and det>0:
            winding+=1
        elif b.y<=0<a.y and det<0:
            winding-=1
    return {"outcome":"winding","winding":winding}


def dyadic_arc(p,level,index):
    count=2**level
    if level<1 or not 0<=index<count:
        raise ValueError("invalid dyadic cell")
    n=len(p)
    lo,hi=Q(n*index,count),Q(n*(index+1),count)
    def at(t):
        k=int(t)
        return lerp(p[k%n],p[(k+1)%n],t-k)
    tags=[lo]+[Q(i) for i in range(int(lo)+1,int(hi)+1) if i<hi]+[hi]
    return [at(t) for t in tags]


def point_segment_squared(a,c,d):
    v=sub(d,c)
    t=max(Q(0),min(Q(1),dot(sub(a,c),v)/dot(v,v)))
    gap=sub(a,lerp(c,d,t))
    return dot(gap,gap)


def separation_squared(j,k):
    def distance(first,second):
        a,b=first
        c,d=second
        if _intersect(a,b,c,d):
            return Q(0)
        return min(point_segment_squared(p,c,d) for p in (a,b))
    return min(min(distance(first,second),distance(second,first))
               for first,second in product(edges(j),edges(k)))


def scan_level(p,level,side):
    """Search one mesh. An empty result certifies no conclusion on this mesh."""
    if side<=0:
        raise ValueError("positive side threshold required")
    arcs=[dyadic_arc(p,level,i) for i in range(2**level)]
    for i,j in product(range(len(arcs)),repeat=2):
        sep=separation_squared(arcs[i],arcs[j])
        if sep<=side**2:
            continue
        result=window_alternative(p,arcs[i],arcs[j])
        if result["outcome"]=="boundary square" or result["winding"]:
            return dict(result,level=level,cells=[i,j],separation_squared=str(sep))
    return None
