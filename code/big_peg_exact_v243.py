"""v243: exhaustive exact maximum inscribed-square computation on rational polygons.

All comparisons and the complete edge-assignment census use Fraction.
Degenerate positive-dimensional families are included: the squared side
length is convex on each compact feasible affine parameter polytope, so
its maximum occurs at a vertex. Unlike a generic-root finder, this method
cannot silently discard singular charts or vertex-on-edge configurations.

For a sharp Big Peg counterexample, supply an independently established
rational polygon AND a disk actually contained in its interior. This module
does not supply or pretend to reconstruct Pochuev's reported 10-gon.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import argparse
import json

from polygon_orbit_census_v223 import (
    Point, pt, add, sub, mul, dot, cross, validate_jordan_polygon,
)

def _turn(p, orientation):
    return pt(-orientation*p.y, orientation*p.x)

def _form_on(form, t, u):
    return form[0] + form[1]*t + form[2]*u

def _boundary_intersection(f, g):
    """Intersection of two nonparallel affine zero-lines in (t,u)."""
    det = f[1]*g[2] - g[1]*f[2]
    if not det:
        return None
    return ((f[2]*g[0] - f[0]*g[2])/det,
            (f[0]*g[1] - f[1]*g[0])/det)

def _candidate_vertices(equalities, nonnegative):
    """Vertices of a bounded feasible set in [0,1]^2, including rank 0/1/2."""
    active = []
    for f in equalities:
        if f[1] or f[2]:
            active.append(f)
        elif f[0]:
            return (), -1

    def permitted(t,u):
        return (all(_form_on(e,t,u)==0 for e in equalities)
                and all(_form_on(v,t,u)>=0 for v in nonnegative))

    if len(active)>1:
        hit = _boundary_intersection(active[0], active[1])
        if hit is not None:
            t,u=hit
            return ((hit,) if permitted(t,u) else ()), 2

    if active:
        f=active[0]
        if f[1]:
            start=(-f[0]/f[1],Q(0))
            direction=(-f[2]/f[1],Q(1))
        else:
            start=(Q(0),-f[0]/f[2])
            direction=(Q(1),Q(0))
        if any(_form_on(e,*start) for e in equalities):
            return (), -1
        lower=upper=None
        for v in nonnegative:
            value=_form_on(v,*start)
            slope=v[1]*direction[0]+v[2]*direction[1]
            if not slope:
                if value<0:
                    return (), 1
            else:
                endpoint=-value/slope
                if slope>0:
                    lower=endpoint if lower is None else max(lower,endpoint)
                else:
                    upper=endpoint if upper is None else min(upper,endpoint)
        if lower is None or upper is None:
            raise AssertionError("parameter box must bound every feasible line")
        if lower>upper:
            return (), 1
        hits=[]
        for value in (lower,upper):
            t=start[0]+value*direction[0]
            u=start[1]+value*direction[1]
            if permitted(t,u) and (t,u) not in hits:
                hits.append((t,u))
        return tuple(hits), 1

    hits=set()
    for f,g in combinations(nonnegative,2):
        hit=_boundary_intersection(f,g)
        if hit is not None and permitted(*hit):
            hits.add(hit)
    return tuple(sorted(hits)), 0

def _point_form(edge_origin, edge_direction, triple):
    a,b,c=triple
    return (
        (cross(edge_direction,sub(a,edge_origin)),
         cross(edge_direction,b), cross(edge_direction,c)),
        (dot(edge_direction,sub(a,edge_origin)),
         dot(edge_direction,b), dot(edge_direction,c)),
        dot(edge_direction,edge_direction),
    )

def _segment_constraints(edge_form):
    line, along, norm_sq=edge_form
    upper=(norm_sq-along[0],-along[1],-along[2])
    return line,(along,upper)

def _on_edge(x,a,b):
    v=sub(b,a)
    return (cross(sub(x,a),v)==0
            and min(a.x,b.x)<=x.x<=max(a.x,b.x)
            and min(a.y,b.y)<=x.y<=max(a.y,b.y))

def _real_square(p,q,r,s,orientation):
    d=sub(q,p)
    j=_turn(d,orientation)
    return (r==add(q,j) and s==add(p,j) and dot(d,d)>0)

@dataclass(frozen=True)
class SquareMaximum:
    side_squared: Q
    witness: tuple | None
    edge_assignment: tuple | None
    orientation: int | None
    parameters: tuple | None
    charts: int
    feasible_charts: int
    singular_feasible_charts: int

    @property
    def positive(self):
        return self.side_squared>0

def maximum_square(vertices):
    """Exhaustive rational upper AND attained lower bound on squared side.

    P=A+tE, Q=B+uH, R=Q+/-J(Q-P), S=P+/-J(Q-P).
    All segment conditions are affine in t,u, including singular families.
    Every maximum of the convex squared side length occurs at a vertex
    of the compact feasible parameter polytope.
    """
    P=validate_jordan_polygon(tuple(vertices));n=len(P)
    edges=tuple((P[i],P[(i+1)%n]) for i in range(n))
    best=Q(0); witness=None; assignment=None; direction=None; parameters=None
    charts=feasible=singular=0
    for i,(A,A2) in enumerate(edges):
        E=sub(A2,A)
        for j,(B,B2) in enumerate(edges):
            H=sub(B2,B)
            d0=sub(B,A)
            for sigma in (1,-1):
                j0=_turn(d0,sigma)
                jt=_turn(mul(E,-1),sigma)
                ju=_turn(H,sigma)
                R=(add(B,j0),jt,add(H,ju))
                S=(add(A,j0),add(E,jt),ju)
                r_forms=[_point_form(C,sub(D,C),R) for C,D in edges]
                s_forms=[_point_form(C,sub(D,C),S) for C,D in edges]
                for k in range(n):
                    line_r,bounds_r=_segment_constraints(r_forms[k])
                    for ell in range(n):
                        charts+=1
                        line_s,bounds_s=_segment_constraints(s_forms[ell])
                        constraints=(
                            (Q(0),Q(1),Q(0)),(Q(1),Q(-1),Q(0)),
                            (Q(0),Q(0),Q(1)),(Q(1),Q(0),Q(-1)),
                            *bounds_r,*bounds_s,
                        )
                        cand,rank=_candidate_vertices((line_r,line_s),constraints)
                        if not cand:
                            continue
                        feasible+=1
                        if rank<2:
                            singular+=1
                        for t,u in cand:
                            d=add(add(d0,mul(E,-t)),mul(H,u))
                            side_sq=dot(d,d)
                            if side_sq<=best:
                                continue
                            p=add(A,mul(E,t))
                            q=add(B,mul(H,u))
                            rot=_turn(d,sigma)
                            r=add(q,rot)
                            s=add(p,rot)
                            if not (_real_square(p,q,r,s,sigma)
                                    and _on_edge(p,*edges[i])
                                    and _on_edge(q,*edges[j])
                                    and _on_edge(r,*edges[k])
                                    and _on_edge(s,*edges[ell])):
                                raise AssertionError("invalid exact square witness")
                            best=side_sq
                            witness=(p,q,r,s)
                            assignment=(i,j,k,ell)
                            direction=sigma
                            parameters=(t,u)
    if witness is None:
        raise AssertionError("no positive square found; investigate census")
    return SquareMaximum(best,witness,assignment,direction,parameters,
                         charts,feasible,singular)

def _point_inside(P,x):
    hits=0
    for i,a in enumerate(P):
        b=P[(i+1)%len(P)]
        if _on_edge(x,a,b):
            return False
        if (a.y>x.y)!=(b.y>x.y):
            cut=a.x+(x.y-a.y)*(b.x-a.x)/(b.y-a.y)
            if cut>x.x:
                hits+=1
    return bool(hits%2)

def _distance_sq_to_segment(p,a,b):
    e=sub(b,a)
    w=sub(p,a)
    ee=dot(e,e)
    z=dot(w,e)
    if z<=0:
        return dot(w,w)
    if z>=ee:
        d=sub(p,b)
        return dot(d,d)
    return cross(w,e)**2/ee

def contained_disk(vertices,center,radius_squared):
    """Certify a given rational closed disk belongs to polygon's closure."""
    P=validate_jordan_polygon(tuple(vertices))
    if not isinstance(center,Point) or not isinstance(radius_squared,Q):
        raise TypeError("exact Point centre and Fraction radius_squared needed")
    if radius_squared<=0:
        raise ValueError("disk radius must be positive")
    if not _point_inside(P,center):
        return False
    return all(_distance_sq_to_segment(center,P[i],P[(i+1)%len(P)])
               >=radius_squared for i in range(len(P)))

def certify_big_peg(vertices,center,radius_squared):
    """If M^2 < 2*rho^2, certify violation of the sharp Big Peg constant."""
    P=validate_jordan_polygon(tuple(vertices))
    if not contained_disk(P,center,radius_squared):
        raise ValueError("the claimed positive disk is not certified contained")
    answer=maximum_square(P)
    return {
        "counterexample_certified":answer.side_squared<2*radius_squared,
        "maximum_side_squared":str(answer.side_squared),
        "certified_disk_radius_squared":str(radius_squared),
        "J_squared_upper_bound":str(answer.side_squared/radius_squared),
        "complete_edge_charts":answer.charts,
        "feasible_charts":answer.feasible_charts,
        "singular_feasible_charts":answer.singular_feasible_charts,
        "max_square_edge_assignment":answer.edge_assignment,
        "max_square_orientation":answer.orientation,
        "max_square_parameters":[str(v) for v in answer.parameters],
        "max_square_vertices":[[str(p.x),str(p.y)] for p in answer.witness],
    }

def read_polygon(path):
    """Read decimal/rational-coordinate JSON; reject rounded JSON floats."""
    data=json.loads(Path(path).read_text())
    def exact(value):
        if type(value) not in (int,str):
            raise ValueError("coordinates must be strings/integers, not floats")
        return Q(value)
    vertices=tuple(pt(exact(x),exact(y)) for x,y in data["vertices"])
    validate_jordan_polygon(vertices)
    disk=data.get("disk")
    if disk is None:
        return vertices,None,None
    centre=pt(*(exact(v) for v in disk["center"]))
    radius_squared=exact(disk["radius_squared"])
    return vertices,centre,radius_squared

def cli():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polygon-json",type=Path)
    args=parser.parse_args()
    if args.polygon_json is None:
        square=tuple(pt(*v) for v in ((-1,-1),(1,-1),(1,1),(-1,1)))
        answer=certify_big_peg(square,pt(0,0),Q(1))
        assert answer["maximum_side_squared"]=="4"
        assert not answer["counterexample_certified"]
        print(json.dumps({"self_test":"PASS","square":answer},indent=2))
    else:
        P,centre,r2=read_polygon(args.polygon_json)
        output=(certify_big_peg(P,centre,r2) if centre is not None else {
            "maximum_side_squared":str(maximum_square(P).side_squared),
            "disk":"not supplied; no Big Peg counterexample claim"})
        print(json.dumps(output,indent=2))

if __name__=="__main__":
    cli()
