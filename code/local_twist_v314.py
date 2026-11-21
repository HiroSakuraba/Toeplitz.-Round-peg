"""Half-turn fixed-circle degree and center-wall certificates in rational arithmetic."""
from fractions import Fraction as Q
import json

from polygon_orbit_census_v223 import (
    _on_seg, _solve, add, cross, dot, mul, pt, sub, validate_jordan_polygon,
)
from rough_arc_winding_v312 import edges, rot
from boundary_or_winding_v313 import point_sign
from half_turn_degree_v313 import polygon as ear_polygon, tube_certificate


def radial_polygon():
    rays = [pt(3,0),pt(3,1),pt(2,2),pt(1,3),pt(0,3)]
    first = [mul(p,r) for p,r in zip(rays,[3,1,1,3,1])]
    second = [rot(mul(p,r)) for p,r in zip(rays,[1,2,3,2,3])]
    half = first[:-1]+second[:-1]
    return validate_jordan_polygon(tuple(half+[mul(p,-1) for p in half]))


def direction(t):
    return pt((1-t*t)/(1+t*t),2*t/(1+t*t))


def box_clear(a,b,boundary):
    xmin,xmax=sorted((a.x,b.x))
    ymin,ymax=sorted((a.y,b.y))
    for c,d in boundary:
        if max(c.x,d.x)<xmin or min(c.x,d.x)>xmax or max(c.y,d.y)<ymin or min(c.y,d.y)>ymax:
            continue
        n=rot(sub(d,c))
        values=[n.x*x+n.y*y for x in (xmin,xmax) for y in (ymin,ymax)]
        if min(values)<=dot(n,c)<=max(values):
            return False
    return True


def twist(vertices,radius,max_depth=18):
    """Compute seam parity for a centered polygon about 0, when certified."""
    if radius<=0:
        raise ValueError("positive radius required")
    p=validate_jordan_polygon(tuple(vertices))
    boundary=list(zip(p,p[1:]+p[:1]))
    if {(mul(a,-1),mul(b,-1)) for a,b in boundary}!=set(boundary):
        raise ValueError("expected matching half-turn edges")
    def positions(t):
        u=mul(direction(t),radius)
        return u,rot(u)
    stack=[(Q(0),Q(1),0)]
    marks=[Q(0)]
    while stack:
        lo,hi,depth=stack.pop()
        a,b=positions(lo)
        c,d=positions(hi)
        if box_clear(a,c,boundary) or box_clear(b,d,boundary):
            marks.append(hi)
            continue
        if depth>=max_depth:
            raise ValueError("centered cutoff unresolved")
        mid=(lo+hi)/2
        stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
    # Work in evaluation coordinates (a,b), with the quarter seam swapping them.
    # point_sign is positive inside; use the exterior-positive convention here.
    path=[pt(*[-point_sign(z,boundary) for z in positions(t)]) for t in marks]
    if path[-1]!=pt(path[0].y,path[0].x):
        raise ArithmeticError("wrong seam")
    base=pt(1,1)
    connector=[base,path[0]]
    if _on_seg(pt(0,0),*connector):
        connector=[base,pt(-1,1),path[0]]
    walk=connector+path[1:]+[pt(v.y,v.x) for v in connector[-2::-1]]
    winding=0
    for a,b in edges(walk):
        if _on_seg(pt(0,0),a,b):
            raise ArithmeticError("unresolved angle")
        det=cross(a,b)
        if a.y<=0<b.y and det>0:
            winding+=1
        elif b.y<=0<a.y and det<0:
            winding-=1
    return {"radius":str(radius),"twist":winding%2,"subintervals":len(marks)-1}


def clip(poly,n,bound):
    result=[]
    for a,b in zip(poly,poly[1:]+poly[:1]):
        fa,fb=dot(n,a)-bound,dot(n,b)-bound
        if fa<=0:
            result.append(a)
        if (fa<0<fb) or (fb<0<fa):
            result.append(add(a,mul(sub(b,a),fa/(fa-fb))))
    return result


def pull(p,j):
    for _ in range((-j)%4):
        p=rot(p)
    return p


def center_charts(vertices,halfwidth,epsilon):
    if halfwidth<=0 or epsilon<0:
        raise ValueError("invalid center box")
    boundary=list(zip(vertices,vertices[1:]+vertices[:1]))
    margin=halfwidth+epsilon
    bound=max(max(abs(p.x),abs(p.y)) for p in vertices)+margin
    initial=[pt(-bound,-bound),pt(bound,-bound),pt(bound,bound),pt(-bound,bound)]
    stages=[((),initial)]
    for j in range(4):
        following=[]
        for i,(a,b) in enumerate(boundary):
            a,b=pull(a,j),pull(b,j)
            normal=rot(sub(b,a))
            normals=(pt(1,0),pt(-1,0),pt(0,1),pt(0,-1),normal,mul(normal,-1))
            tests=[(n,max(dot(n,a),dot(n,b))+margin*(abs(n.x)+abs(n.y))) for n in normals]
            for ids,poly in stages:
                for n,rhs in tests:
                    poly=clip(poly,n,rhs)
                    if not poly:
                        break
                if poly:
                    following.append((ids+(i,),poly))
        stages=following
    return [ids for ids,_ in stages]


def chart_system(boundary,ids):
    rows=[]; rhs=[]; weights=[]
    for j,i in enumerate(ids):
        a,b=boundary[i]
        n=rot(sub(b,a))
        u=pull(n,j)
        rows.append([n.x,n.y,u.x,u.y])
        rhs.append(dot(n,a))
        weights.append(abs(n.x)+abs(n.y))
    return rows,rhs,weights


def local_certificate(halfwidth=Q(1,8),epsilon=Q(1,128)):
    p=radial_polygon()
    boundary=list(zip(p,p[1:]+p[:1]))
    charts=center_charts(p,halfwidth,epsilon)
    sample=validate_jordan_polygon(tuple(pt(v.x+Q((3*i*i-i+1)%15-7,1024),
                                            v.y+Q((i*i+2*i+3)%15-7,1024))
                                         for i,v in enumerate(p)))
    sample_edges=list(zip(sample,sample[1:]+sample[:1]))
    shifted=set()
    radii={}
    max_center=Q(0)
    for ids in charts:
        rows,rhs,weights=chart_system(boundary,ids)
        root=_solve(rows,rhs)
        sample_rows,sample_rhs,_=chart_system(sample_edges,ids)
        new=_solve(sample_rows,sample_rhs)
        if new is None:
            raise ArithmeticError("singular sample chart")
        c,v=pt(*new[:2]),pt(*new[2:])
        if max(abs(c.x),abs(c.y))<=halfwidth:
            points=[]
            for i in ids:
                points.append(add(c,v))
                v=rot(v)
            if all(_on_seg(z,*sample_edges[i]) for z,i in zip(points,ids)):
                shifted.add((c.x,c.y,dot(v,v)))
        if root is None or tuple(root[:2])!=(0,0):
            raise ArithmeticError("noncentral or nonregular candidate",ids)
        v=pt(*root[2:])
        r2=dot(v,v)
        for j,i in enumerate(ids):
            a,b=boundary[i]
            if not _on_seg(v,a,b) or v in (a,b):
                raise ArithmeticError("invalid chart root",ids)
            v=rot(v)
        columns=[_solve(rows,[int(k==i) for k in range(4)]) for i in range(4)]
        errors=[epsilon*sum(abs(col[j])*w for col,w in zip(columns,weights)) for j in range(4)]
        max_center=max(max_center,*errors[:2])
        low=sum(max(Q(0),abs(x)-z)**2 for x,z in zip(root[2:],errors[2:]))
        high=sum((abs(x)+z)**2 for x,z in zip(root[2:],errors[2:]))
        radii.setdefault(r2,[]).append((low,high))
    if max_center>=halfwidth or not radii or any(len(group)!=4 for group in radii.values()):
        raise ArithmeticError("center wall or orbit multiplicity unresolved")
    bands=[(min(a for a,b in group),max(b for a,b in group)) for _,group in sorted(radii.items())]
    cuts=[Q(4),Q(5),Q(33,5),Q(7)]
    records=[twist(p,r) for r in cuts]
    for r,record in zip(cuts,records):
        if any(lo<=r*r<=hi for lo,hi in bands):
            raise ArithmeticError("tube reaches radius wall")
        if sum(k>r*r for k in radii)%2!=record["twist"]:
            raise ArithmeticError("circle and geometric census disagree")
    if len(shifted)!=3 or any(x==0 and y==0 for x,y,_ in shifted):
        raise ArithmeticError("unexpected displaced sample squares")
    return {"sample_centers_and_radii_squared":[list(map(str,z)) for z in sorted(shifted,key=lambda z:z[2])],
            "charts":len(charts),"radii_squared":list(map(str,sorted(radii))),
            "radius_squared_bands":[list(map(str,b)) for b in bands],
            "halfwidth":str(halfwidth),"epsilon":str(epsilon),
            "center_sup_bound":str(max_center),"circle_scan":records}


def cap_clearance(eta=Q(1,32)):
    H,X,W,s,h=Q(25,16),Q(3,2),Q(4),Q(2,5),Q(1,2)
    eta=Q(eta)
    if eta<0:
        raise ValueError("negative width")
    margins={"single_cap":1-eta-h,"opposing_caps":h-((H-1)/2+eta),
             "ear_to_wall":W-X-2*eta-2*(H+eta-h),
             "two_walls":2*W-2*eta-2*(H+eta),"two_rows":2-2*eta-2*(1+eta-h),
             "adjacent_strips":s-4*eta,"small_square":1-eta-h-s,
             "side_reach":W-X-2*eta-(2-2*eta)}
    if min(margins.values())<=0:
        raise ValueError("cap separation not certified")
    return {"eta":str(eta),"side_floor_squared":str((2-2*eta)**2/2),
            "slab_halfheight":str(h),"margins":{k:str(v) for k,v in margins.items()}}


def report():
    return {"star":local_certificate(),"caps":cap_clearance(),
            "split_cutoff_tube":tube_certificate(ear_polygon(Q(1),Q(9,10)))}


if __name__=="__main__":
    print(json.dumps(report(),indent=2))
