"""Exact boundary crossing and independent marked-square census for four arcs."""
from fractions import Fraction as Q
from itertools import combinations, product
import json

from polygon_orbit_census_v223 import (
    _affine_system, _box_vertices, _evaluate, _intersect, _on_seg,
    _side_sq, _solve, add, cross, mul, orientation_area, pt, sub, validate_jordan_polygon,
)


def edges(arc):
    return list(zip(arc, arc[1:]))


def rot(v):
    return pt(-v.y, v.x)


def complete(q, b):
    offset = rot(sub(b, q))
    return add(q, offset), add(b, offset)


def cut_parameter(a, b, c, d):
    """Parameter of a transverse interior crossing. Reject ambiguous ray events."""
    if not _intersect(a, b, c, d):
        return None
    v, w = sub(b, a), sub(d, c)
    denominator = cross(v, w)
    if not denominator:
        raise ValueError("overlap in boundary sweep")
    t = cross(sub(c, a), w) / denominator
    u = cross(sub(c, a), v) / denominator
    if not (0 < t < 1 and 0 < u < 1):
        raise ValueError("boundary crossing at a segment endpoint")
    return t


def above(d, last):
    # Vertical ray from below, with the half-open vertex convention.
    result = False
    for a,b in edges(last):
        if _on_seg(d,a,b):
            raise ValueError("square on the parameter boundary")
        if (a.x > d.x) != (b.x > d.x):
            height = a.y + (d.x-a.x)*(b.y-a.y)/(b.x-a.x)
            if height < d.y:
                result = not result
    return result


def validate_arcs(arcs):
    for arc in arcs:
        if len(arc)<2 or len(set(arc))!=len(arc):
            raise ValueError("need simple arcs with distinct vertices")
        segments=edges(arc)
        for index,(a,b) in enumerate(segments):
            if index and _on_seg(arc[index-1],a,b):
                raise ValueError("adjacent backtracking")
            if index and _on_seg(b,arc[index-1],a):
                raise ValueError("adjacent backtracking")
            if any(_intersect(a,b,c,d) for c,d in segments[index+2:]):
                raise ValueError("self-intersecting arc")
    for first,second in combinations(arcs,2):
        if any(_intersect(a,b,c,d) for (a,b),(c,d) in product(edges(first),edges(second))):
            raise ValueError("marked arcs intersect")


def boundary_bit(arcs):
    first, bottom, top, last = arcs
    validate_arcs(arcs)
    # Infinite vertical / horizontal tails extend the separator arcs.
    first = sorted_ends(first, "y")
    last = sorted_ends(last, "x")
    low, high = first[0].y, first[-1].y
    left, right = last[0].x, last[-1].x
    if (any(not low<p.y<high for p in first[1:-1])
            or any(not left<p.x<right for p in last[1:-1])):
        raise ValueError("tails may meet their separator")
    # Coordinate extrema of the completions occur at input vertex pairs.
    # Use separate one-dimensional extrema to avoid an interior pair census.
    def extent(arc, fn):
        values = [fn(p) for p in arc]
        return min(values),max(values)
    by = extent(bottom,lambda p:p.y-p.x)
    tx = extent(top,lambda p:p.x)
    bx = extent(bottom,lambda p:p.y)
    ty = extent(top,lambda p:p.x-p.y)
    if not (low < by[0]+tx[0] <= by[1]+tx[1] < high
            and left < bx[0]+ty[0] <= bx[1]+ty[1] < right):
        raise ValueError("completion can reach an artificial tail")

    return ray_winding(first, bottom, top, last, positive_left=False)


def ray_winding(first, bottom, top, last, positive_left):
    faces = [[(q,top[0]) for q in bottom],[(bottom[-1],b) for b in top],
             [(q,top[-1]) for q in reversed(bottom)],[(bottom[0],b) for b in reversed(top)]]
    counts = []
    winding = 0
    for face in faces:
        count = 0
        path = [complete(q,b) for q,b in face]
        for (p,d),(p1,d1) in edges(path):
            for a,b in edges(first):
                t = cut_parameter(p,p1,a,b)
                if t is not None and above(add(d,mul(sub(d1,d),t)),last):
                    count += 1
                    # On the positive vertical ray, a right-to-left crossing is positive.
                    direction = 1 if cross(sub(p1,p),sub(b,a)) < 0 else -1
                    winding += -direction if positive_left else direction
        counts.append(count)
    return {"positive_ray_hits_by_face":counts,"winding":winding,"bit":sum(counts)%2}


def sorted_ends(arc, coord):
    return arc if getattr(arc[0],coord)<getattr(arc[-1],coord) else arc[::-1]


def whole_curve_winding(polygon, first_side_arc, second_side_arc):
    """Inputs are polygonal subarcs on distinct parameter windows of this polygon."""
    polygon = validate_jordan_polygon(tuple(polygon))
    if orientation_area(polygon)<0:
        polygon = polygon[::-1]
    closed = list(polygon)+[polygon[0]]
    validate_arcs((first_side_arc,second_side_arc))
    # Allow partial edges and subdivisions, with each supplied segment
    # lying on one input edge. Split an arc at any intervening polygon vertex.
    for arc in (first_side_arc,second_side_arc):
        for a,b in edges(arc):
            if not any(_on_seg(a,c,d) and _on_seg(b,c,d) for c,d in edges(closed)):
                raise ValueError("side arc is not on the polygon")
    return ray_winding(closed,first_side_arc,second_side_arc,closed,positive_left=True)


def triangle_example():
    triangle = [pt(0,0),pt(4,0),pt(1,3)]
    j,k = [[pt(x+t,0) for t in (-Q(1,16),Q(1,16))] for x in (Q(4,7),Q(16,7))]
    result = whole_curve_winding(triangle,j,k)
    clearance,eta = Q(1,40),Q(1,160)
    if Q(7,80)**2 <= 10*clearance**2 or (clearance-2*eta)**2<=2*eta**2:
        raise ArithmeticError("triangle clearance exhausted")
    return dict(result,wall_clearance_lower=str(clearance),registration=str(eta),
                side_floor=str(k[0].x-j[-1].x-2*eta))


def marked_census(arcs):
    """Solve every product of four segments. Singular feasible families fail closed."""
    found = {}
    for cell in product(*(edges(arc) for arc in arcs)):
        a,b = _affine_system(cell,1)
        tags = _solve(a,b)
        if tags is None:
            if tuple(_box_vertices(a,b)):
                raise ValueError("singular feasible marked cell")
            continue
        if all(0<=t<=1 for t in tags):
            if any(t in (0,1) for t in tags):
                raise ValueError("marked square on an edge endpoint")
            points = _evaluate(cell,tags)
            found[points] = _side_sq(points)
    return sorted(found.values())


def full_side_window_census(polygon, bottom, top):
    """Exact independent checks on both remaining vertices anywhere on the loop."""
    boundary = edges(list(polygon)+[polygon[0]])
    found = []
    def candidates(points):
        xmin,xmax = min(p.x for p in points),max(p.x for p in points)
        ymin,ymax = min(p.y for p in points),max(p.y for p in points)
        return [(a,b) for a,b in boundary if min(a.x,b.x)<=xmax and max(a.x,b.x)>=xmin
                and min(a.y,b.y)<=ymax and max(a.y,b.y)>=ymin]
    for j,k in product(edges(bottom),edges(top)):
        ps,ds = zip(*(complete(q,b) for q,b in product(j,k)))
        for i,l in product(candidates(ps),candidates(ds)):
            found.extend(marked_census((i,j,k,l)))
    return sorted(found)


def corridor(h=Q(1,4), e=Q(1,16), a=Q(5,8)):
    h,e,a = map(Q,(h,e,a))
    if not (0<=e<h and 2*h+e<a<1-h and a<1-e):
        raise ValueError("corridors lack end clearance or separation")
    return {"face_margin":str(h-3*e),"tail_margin":str(a-2*h-e),
            "side_floor":str(1-2*e),
            "status":"odd degree" if h>3*e else
                     "closed-arc existence" if h==3*e else "no uniform guarantee"}


def folded_arcs():
    # (longitudinal fraction, transverse displacement in units of 1/128).
    # Each input arc traverses a portion of its long-coordinate range backwards.
    patterns = [
        [(-1,0),(Q(-1,5),-3),(Q(1,7),-2),(Q(-1,8),4),(Q(2,5),5),(1,0)],
        [(-1,0),(Q(-1,4),4),(Q(1,6),5),(Q(-1,7),-3),(Q(2,5),-2),(1,0)],
        [(-1,0),(Q(-1,5),-5),(Q(1,9),-4),(Q(-1,6),3),(Q(1,3),4),(1,0)],
        [(-1,0),(Q(-1,4),3),(Q(1,8),5),(Q(-1,7),-4),(Q(2,5),-3),(1,0)]]
    arcs=[]
    for i,pattern in enumerate(patterns):
        width = Q(5,8) if i in (0,3) else Q(1,4)
        local = [(s*width,Q(t,128)) for s,t in pattern]
        arc = ([pt(t,s) for s,t in local] if i==0 else
               [pt(1+t,s) for s,t in local] if i==1 else
               [pt(1+s,1+t) for s,t in local] if i==2 else
               [pt(s,1+t) for s,t in local])
        arcs.append(arc)
    return arcs


def loop(arcs):
    i,j,k,l = arcs
    return validate_jordan_polygon(tuple(i[::-1] + [pt(0,Q(-3,4)),pt(1,Q(-3,4))]
        + j + [pt(Q(3,2),Q(1,3)),pt(Q(3,2),1)] + k[::-1] + l[::-1]
        + [pt(Q(-3,4),1),pt(Q(-3,4),Q(5,8))]))


def model_arcs(name):
    a,h = Q(5,8),Q(1,4)
    i = [pt(0,-a),pt(0,a)]
    l = [pt(-a,1),pt(a,1)]
    j = [pt(1,-h),pt(1,h)]
    k = [pt(1-h,1),pt(1+h,1)]
    if name=="folded":
        return folded_arcs()
    if name=="winding_only":
        j[0] = pt(1,Q(-1,32))
        k[1:1] = [pt(Q(13,16),Q(15,16)),pt(Q(15,16),1)]
    elif name=="two":
        j = [pt(Q(15,16),-h),pt(1,h),pt(Q(17,16),-h)]
    elif name=="empty":
        j = [pt(Q(15,16),-h),pt(Q(17,16),-Q(1,8))]
    else:
        raise ValueError("unknown model")
    return i,j,k,l


def notch_strip(e=Q(1,80)):
    e=Q(e)
    margins={"end of first arc":Q(1,32)-2*e, "bottom ends":Q(1,8)-2*e,
             "top ends":Q(1,8)-3*e, "last arc ends":Q(1,16)-2*e}
    if e<0 or min(margins.values())<=0:
        raise ValueError("notch strip face clearance exhausted")
    return {"halfwidth":str(e),"side_floor":str(1-2*e),
            "face_margins":{k:str(v) for k,v in margins.items()}}


def roughening_diamonds(arcs, alpha=Q(1,4096)):
    """Enclose all infinite series on marked edges and check the connecting polygon."""
    boundary = loop(arcs)
    marked = {frozenset(edge) for arc in arcs for edge in edges(arc)}
    cells = []
    for a,b in edges(list(boundary)+[boundary[0]]):
        if frozenset((a,b)) in marked:
            mid = mul(add(a,b), Q(1,2))
            offset = mul(rot(sub(b,a)),alpha/2)
            cell = [a,sub(mid,offset),b,add(mid,offset)]
        else:
            cell = [a,b]
        cells.append(cell)
    for n,cell in enumerate(cells):
        a,b = boundary[n],boundary[(n+1)%len(cells)]
        for m in range(n):
            other = cells[m]
            c,d = boundary[m],boundary[(m+1)%len(cells)]
            common = {a,b} & {c,d}
            for x,y in edges(cell+cell[:1]):
                for u,v in edges(other+other[:1]):
                    if _intersect(x,y,u,v):
                        shared = {x,y} & {u,v} & common
                        if len(shared)!=1:
                            raise ArithmeticError("overlapping roughening cells")
                        w, = shared
                        if any(_on_seg(z,u,v) for z in (x,y) if z!=w) or any(
                                _on_seg(z,x,y) for z in (u,v) if z!=w):
                            raise ArithmeticError("overlapping roughening edges")
            # Convex cells can also contain another without a proper crossing.
            for poly,vertices in ((cell,other),(other,cell)):
                if len(poly)==4 and any(all(cross(sub(y,x),sub(z,x))>=0
                        for x,y in edges(poly+poly[:1])) for z in vertices if z not in common):
                    raise ArithmeticError("nested roughening cells")
    for index,arc in enumerate(arcs):
        center = 0 if index==0 else 1
        long_center = 1 if index==2 else 0
        half_length = Q(5,8) if index in (0,3) else Q(1,4)
        for a,b in edges(arc):
            mid = mul(add(a,b),Q(1,2))
            offset = mul(rot(sub(b,a)),alpha/2)
            for point in (sub(mid,offset),add(mid,offset)):
                longitudinal,transverse = ((point.y,point.x) if index<2 else (point.x,point.y))
                if not (abs(longitudinal-long_center)<half_length
                        and abs(transverse-center)<Q(1,16)):
                    raise ArithmeticError("roughening exits corridor")
    return {"marked_edges":len(marked),"polygon_edges":len(cells),"coefficient":str(alpha)}


def report():
    data = {}
    for name in ("folded","winding_only","two","empty"):
        arcs = model_arcs(name)
        result = boundary_bit(arcs)
        sizes = marked_census(arcs)
        if len(sizes)%2!=result["bit"]:
            raise ArithmeticError("boundary disagrees with marked census")
        data[name] = dict(result,count=len(sizes),squared_sides=list(map(str,sizes)))
    arcs = folded_arcs()
    polygon = loop(arcs)
    whole = whole_curve_winding(polygon,arcs[1],arcs[2])
    side_sizes = full_side_window_census(polygon,arcs[1],arcs[2])
    whole["oriented_side_window_count"] = len(side_sizes)
    if len(side_sizes)%2 != whole["bit"]:
        raise ArithmeticError("whole-curve winding disagrees with side-window census")
    return {"whole_folded":whole,"whole_triangle":triangle_example(),"corridor":corridor(),"notch_strips":notch_strip(),
            "infinite_roughening":roughening_diamonds(folded_arcs()),"models":data}


if __name__=="__main__":
    print(json.dumps(report(),indent=2))
