"""v264: reusable exact moving-convex-polygon square-orbit event ledger.

Input: ANY two rational endpoint polygons with the same vertex count,
affinely interpolated over t in [0,1]. Certify strict convexity after
removing permanently collinear, forward-pointing subdivisions; preserve
the polygonal image exactly. Produce an exhaustive weakly-increasing
four-edge chart cover, with certified exclusions, stable interior/
boundary-locked branches, unresolved singular charts, and algebraic
event brackets. Rational chambers use v223's independent complete
geometric-orbit census; irrational event slices remain UNKNOWN unless
an independent event theorem is supplied.

An optional quantitative nonadjacent-edge separation proof certifies a
strict, uniform positive side-size lower bound, stable under exact
collinear refinement after compression. No claim for arbitrary C0
Jordan limits or nongraceful square orders on nonconvex curves.

Standard library Fraction and exact existing v223/v224/v249 machinery.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from math import comb
import json

from polygon_orbit_census_v223 import (
    pt, sub, dot, validate_jordan_polygon, census_polygon,
    _affine_system, PositiveDimensionalSquares,
)
from moving_polygon_events_v224 import (
    chart_polynomials, interpolated_polygon,
    moving_point, pairminus, pairdot, pairrot,
    minus, times, evaluate, iszero, isolate_roots,
)
from big_peg_multichart_v249 import (
    pol, psub, pscale, pneg, certify_sign,
)
from pentagon_complete_ledger_v263 import (
    singular_inconsistency_witness,
)
from near_square_transversality_v258 import det4


class InconclusiveAtlas(RuntimeError):
    """An unsupported case MUST NOT be converted to zero square orbits."""


def _strict_sign(p):
    """+/- iff strictly signed for EVERY t in the CLOSED [0,1].

    First use fast exact Bernstein certificates. If those fail, exact
    Sturm root isolation proves absence of zero and the sign at t=0.
    A root at an endpoint also prevents calling this strictly signed.
    """
    if iszero(p):
        return 0
    for sign,z in ((1,p),(-1,pneg(p))):
        try:
            certify_sign(z,Q(0),Q(1),strict=True,max_depth=12)
            return sign
        except ValueError:
            pass
    if not isolate_roots(p,Q(0),Q(1),width=Q(1,2**16)):
        val=evaluate(p,Q(0))
        if val==0:
            raise AssertionError("root-free polynomial vanished at boundary")
        return 1 if val>0 else -1
    return None


def _poly_cross(a,b):
    return psub(times(a[0],b[1]),times(a[1],b[0]))


def _valid_endpoints(P0,P1):
    A=tuple(P0);B=tuple(P1)
    if len(A)!=len(B) or len(A)<4:
        raise ValueError("matching polygons of at least four vertices required")
    validate_jordan_polygon(A)
    validate_jordan_polygon(B)
    if any(type(v.x) is not Q or type(v.y) is not Q for v in A+B):
        raise ValueError("exact Fraction vertex coordinates required")
    return A,B


def _straight_forward_at_all_times(P0,P1,i):
    n=len(P0)
    vertices=tuple(moving_point(a,b) for a,b in zip(P0,P1))
    prev,at,nxt=(vertices[(i-1)%n],vertices[i],vertices[(i+1)%n])
    e0=pairminus(at,prev)
    e1=pairminus(nxt,at)
    return iszero(_poly_cross(e0,e1)) and _strict_sign(pairdot(e0,e1))==1


def compress_collinear_refinement(P0,P1):
    """Collapse only ALL-TIME straight, forward edges, preserving geometry.

    The signed dot product is strictly positive for all times; the
    collinearity polynomial is identically zero. Thus each removed
    moving vertex stays strictly BETWEEN its neighbors at all t.
    The resulting polygon traces the same geometric boundary.
    """
    A,B=_valid_endpoints(P0,P1)
    a=list(A);b=list(B);labels=list(range(len(A)));removed=[]
    while len(a)>4:
        drop=next((i for i in range(len(a))
                   if _straight_forward_at_all_times(a,b,i)),None)
        if drop is None:
            break
        removed.append(labels.pop(drop))
        a.pop(drop);b.pop(drop)
    return tuple(a),tuple(b),tuple(removed)


def certify_strict_convex_motion(P0,P1):
    """Prove CCW strict convexity on [0,1] by ALL ordered triple signs."""
    A,B=_valid_endpoints(P0,P1)
    V=tuple(moving_point(a,b) for a,b in zip(A,B))
    signs=0
    for i,j,k in combinations(range(len(A)),3):
        e0=pairminus(V[j],V[i]);e1=pairminus(V[k],V[i])
        f=_poly_cross(e0,e1)
        if _strict_sign(f)!=1:
            raise InconclusiveAtlas(
                "CCW strict convexity not established for triple "+
                repr((i,j,k)))
        signs+=1
    return signs


def _vec(p,q):
    return pt(q.x-p.x,q.y-p.y)


def _projection_range(poly_at_a,poly_at_b,i,axis):
    n=len(poly_at_a)
    values=[dot(poly_at_a[j],axis) for j in (i,(i+1)%n)]
    values.extend(dot(poly_at_b[j],axis) for j in (i,(i+1)%n))
    return min(values),max(values)


def _separation_directions(A,B,i,j):
    n=len(A)
    a,b,c,d=A[i],A[(i+1)%n],A[j],A[(j+1)%n]
    e,f=_vec(a,b),_vec(c,d)
    ma=pt((a.x+b.x)/2,(a.y+b.y)/2)
    mb=pt((c.x+d.x)/2,(c.y+d.y)/2)
    dirs=[
        pt(1,0),pt(0,1),pt(1,1),pt(1,-1),
        e,f,pt(-e.y,e.x),pt(-f.y,f.x),
        _vec(ma,mb),
    ]
    for x in (a,b):
        for y in (c,d):
            dirs.append(_vec(x,y))
    return tuple(z for z in dict.fromkeys(dirs) if dot(z,z)>0)


def _one_pair_separation(P0,P1,i,j,max_depth):
    """Construct an exact rational squared-distance lower bound.

    In each closed dyadic time slab find a fixed rational projection
    separating the ENTIRE swept pair of line segments. Cauchy-Schwarz
    yields ||x-y||^2 >= gap^2 / ||axis||^2. The directions are heuristic
    witnesses, but each accepted inequality is exact and sufficient.
    """
    cells=0;depth_seen=0
    stack=[(Q(0),Q(1),0)]
    bounds=[]
    while stack:
        lo,hi,depth=stack.pop()
        depth_seen=max(depth_seen,depth)
        a=interpolated_polygon(P0,P1,lo)
        b=interpolated_polygon(P0,P1,hi)
        middle=interpolated_polygon(P0,P1,(lo+hi)/2)
        accepted=Q(0)
        for axis in _separation_directions(middle,middle,i,j):
            xlo,xhi=_projection_range(a,b,i,axis)
            ylo,yhi=_projection_range(a,b,j,axis)
            gap=max(xlo-yhi,ylo-xhi)
            if gap>0:
                accepted=max(accepted,gap*gap/dot(axis,axis))
        if accepted:
            cells+=1
            bounds.append(accepted)
        elif depth>=max_depth:
            raise InconclusiveAtlas(
                "no positive edge-pair separation witness "+
                repr((i,j,lo,hi)))
        else:
            mid=(lo+hi)/2
            stack.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
    return min(bounds),cells,depth_seen


def certify_uniform_side_floor(P0,P1,max_depth=8):
    """A positive square has corners on nonadjacent supporting edges.

    For a convex polygon with >=4 edges, two adjacent edges cannot
    alone carry 4 square corners, and any 3 different edges include a
    nonadjacent pair. Therefore some square corners have separation at
    least the uniform minimum distance between nonadjacent edges.
    Their separation <=sqrt(2)*square_side.
    """
    if type(max_depth) is not int or not 0<=max_depth<=14:
        raise ValueError("max_depth must be an integer in 0..14")
    A,B=_valid_endpoints(P0,P1)
    # A separate full-time convexity certificate is an essential premise.
    certify_strict_convex_motion(A,B)
    n=len(A)
    pairs=[(i,j) for i,j in combinations(range(n),2)
           if j!=i+1 and (i,j)!=(0,n-1)]
    if not pairs:
        raise InconclusiveAtlas("no nonadjacent edge pairs")
    bounds=[];cells=0;deepest=0
    for i,j in pairs:
        b,c,d=_one_pair_separation(A,B,i,j,max_depth)
        bounds.append(b);cells+=c;deepest=max(deepest,d)
    minimum=min(bounds)
    if minimum<=0:
        raise AssertionError("positive gap not established")
    return {
        "nonadjacent_edge_pairs":len(pairs),
        "dyadic_time_cells_certified":cells,
        "deepest_dyadic_refinement":deepest,
        "edge_pair_distance_squared_lower_bound":minimum,
        "every_positive_square_side_squared_lower_bound":minimum/2,
        "safe_side_squared_cutoff":minimum/4,
    }


def _signed_constraints(chart,sign):
    terms=[]
    for j in range(4):
        terms.append((f"corner-{j}-start",pscale(chart.numerators[j],sign)))
        terms.append((f"corner-{j}-end",pscale(chart.upper_differences[j],sign)))
    for j in range(3):
        if chart.edges[j]==chart.edges[j+1]:
            # Keep four corners in the chosen positive boundary order.
            terms.append((f"repeated-edge-order-{j}",
                          pscale(psub(chart.numerators[j+1],
                                      chart.numerators[j]),sign)))
    return tuple(terms)


def _chart_raw_factors(chart):
    factors=[("singular-det",chart.determinant)]
    for j in range(4):
        factors.append((f"corner-{j}-start",chart.numerators[j]))
        factors.append((f"corner-{j}-end",chart.upper_differences[j]))
    for j in range(3):
        if chart.edges[j]==chart.edges[j+1]:
            factors.append((f"repeated-edge-order-{j}",
                            psub(chart.numerators[j+1],
                                 chart.numerators[j])))
    factors.append(("side-cutoff",chart.cutoff))
    return factors


def classify_chart(P0,P1,ids,cutoff,width=Q(1,2**18)):
    """No genericity assumptions: strict exclusion or explicit event walls.

    The result is valid only once the caller has certified convexity and
    proved the geometric positive-square cutoff for this entire motion.
    """
    n=len(P0)
    if len(P1)!=n or len(ids)!=4 or any(not 0<=e<n for e in ids):
        raise ValueError("invalid moving polygon chart")
    if tuple(ids)!=tuple(sorted(ids)):
        raise ValueError("expected weakly increasing edge indices")
    if type(cutoff) is not Q or cutoff<=0:
        raise ValueError("strictly positive certified cutoff required")
    if max(ids.count(i) for i in set(ids))>=3:
        return {"class":"excluded_three_collinear"}
    if len(set(ids))==2:
        i,j=sorted(set(ids))
        if (j==i+1 or (i,j)==(0,n-1)):
            return {"class":"excluded_two_adjacent_edges"}
    chart=chart_polynomials(P0,P1,ids,orientation=1,
                            cutoff_squared=cutoff)
    D=chart.determinant
    if iszero(D):
        witness=singular_inconsistency_witness(P0,P1,ids)
        return ({"class":"excluded_singular_inconsistent","witness":witness}
                if witness is not None else
                {"class":"unresolved_singular","reason":"identically zero determinant, no universal augmented-rank contradiction"})
    sign=_strict_sign(D)
    if sign in (-1,1):
        factors=_signed_constraints(chart,sign)
        for key,p in factors:
            if _strict_sign(p)==-1:
                return {"class":"excluded_off_edge","witness":key,
                        "det_sign":sign}
        size_sign=_strict_sign(chart.cutoff)
        if size_sign in (-1,0):
            return {"class":"excluded_below_side_cutoff"}
        # A repeated corner cannot coincide with its neighbor in an
        # actual positive-size square.
        for key,p in factors:
            if key.startswith("repeated-edge-order-") and iszero(p):
                return {"class":"excluded_coincident_corners"}
        all_signs=tuple((key,_strict_sign(p)) for key,p in factors)
        if size_sign==1 and all(s in (0,1) for _,s in all_signs):
            locks=[key for key,s in all_signs if s==0]
            return {"class":("stable_boundary" if locks
                              else "stable_interior"),
                    "det_sign":sign,
                    "identically_locked_faces":locks}
    roots=[]
    for key,p in _chart_raw_factors(chart):
        if iszero(p):continue
        for bracket in isolate_roots(p,Q(0),Q(1),width=width):
            roots.append({
                "kind":key,"polynomial":tuple(p),
                "lower":bracket.lower,"upper":bracket.upper,
            })
    if not roots:
        # For any nonzero polynomial, signs cannot change in a closed
        # interval without a root. Reaching here means an unhandled
        # identically locked degeneracy; do NOT silently call it stable.
        return {"class":"unresolved_rootless",
                "reason":"symbolic rank/contact/side classification incomplete"}
    return {"class":"possible_event","roots":roots}


def _event_clusters(raw):
    roots=sorted(raw,key=lambda r:(r["lower"],r["upper"]))
    clusters=[]
    for root in roots:
        if clusters and root["lower"]<=clusters[-1]["upper"]:
            c=clusters[-1]
            c["upper"]=max(c["upper"],root["upper"])
            c["factors"].append((root["edges"],root["kind"]))
        else:
            clusters.append({"lower":root["lower"],"upper":root["upper"],
                             "factors":[(root["edges"],root["kind"])]})
    return clusters


def _sample_rational_chamber(P0,P1,low,high,cutoff,side_floor):
    t=(low+high)/2
    polygon=interpolated_polygon(P0,P1,t)
    try:
        result=census_polygon(polygon,cutoff_squared=cutoff)
    except PositiveDimensionalSquares as exc:
        raise InconclusiveAtlas(
            "positive-dimensional positive-square family in chamber "+
            repr((low,high))) from exc
    for sq in result.squares:
        if sq.side_squared<side_floor:
            raise AssertionError("square violates global side floor")
        if not sq.graceful:
            raise AssertionError("strict convex polygon has nongraceful square")
    signed=[]
    if all(not sq.has_polygon_vertex for sq in result.squares):
        edges=tuple((polygon[i],polygon[(i+1)%len(polygon)])
                    for i in range(len(polygon)))
        for sq in result.squares:
            positions=sq.curve_positions
            k=min(range(4),key=lambda j:positions[j])
            pos=positions[k:]+positions[:k]
            ids=tuple(int(v) for v in pos)
            if ids!=tuple(sorted(ids)):
                raise AssertionError("convex curve square corners out of order")
            M,b=_affine_system(tuple(edges[e] for e in ids),1)
            det=det4(M)
            if det==0:
                signed=[]
                break
            signed.append(1 if det>0 else -1)
    return {
        "left":str(low),"right":str(high),"sample":str(t),
        "geometric_square_orbits":result.count,"quotient_parity":result.parity,
        "sample_square_side_squared":[str(s.side_squared)
                                       for s in result.squares],
        "regular_oriented_index_sum":(
            sum(signed) if len(signed)==result.count else None),
        "all_sample_squares_regular_in_open_positive_sector":
            len(signed)==result.count,
    }


def audit_moving_convex_polygon(P0,P1,width=Q(1,2**18),max_depth=8,
                                allow_inconclusive=False):
    """Exhaust all positive-square charts; certify chambers, not guesses.

    If any singular or polynomial case cannot be excluded/certified,
    raise InconclusiveAtlas by default, rather than reporting a count.
    Irrational event brackets are explicitly NOT full-time root counts.
    """
    A,B=_valid_endpoints(P0,P1)
    C,D,removed=compress_collinear_refinement(A,B)
    triple_checks=certify_strict_convex_motion(C,D)
    gap=certify_uniform_side_floor(C,D,max_depth=max_depth)
    cutoff=gap["safe_side_squared_cutoff"]
    all_charts=tuple(combinations_with_replacement(range(len(C)),4))
    ledger={}
    classes=Counter()
    unresolved={}
    raw_roots=[]
    for ids in all_charts:
        key=",".join(map(str,ids))
        rec=classify_chart(C,D,ids,cutoff,width)
        classes[rec["class"]]+=1
        if rec["class"]=="possible_event":
            for r in rec["roots"]:
                raw_roots.append(dict(r,edges=key))
            ledger[key]={"class":rec["class"],
                         "candidate_factor_count":len(rec["roots"])}
        else:
            ledger[key]=rec
        if rec["class"].startswith("unresolved"):
            unresolved[key]=rec
    report={
        "initial_vertex_count":len(A),
        "compressed_vertex_count":len(C),
        "all_time_exact_collinear_refinement_vertices_removed":
            list(removed),
        "all_time_CCW_strict_convexity_triple_checks":triple_checks,
        "uniform_positive_square_size_certificate":{
            k:str(v) if isinstance(v,Q) else v for k,v in gap.items()},
        "all_weakly_increasing_support_charts":len(all_charts),
        "chart_classification":dict(sorted(classes.items())),
        "unresolved_charts":unresolved,
        "full_chart_ledger":ledger,
        "status":None,
    }
    if unresolved:
        report["status"]="INCONCLUSIVE: unresolved singular or algebraic charts"
        if not allow_inconclusive:
            raise InconclusiveAtlas("unresolved charts: "+repr(unresolved))
        return report
    clusters=_event_clusters(raw_roots)
    chambers=[]
    cursor=Q(0)
    for e in clusters:
        if cursor<e["lower"]:
            chambers.append(_sample_rational_chamber(
                C,D,cursor,e["lower"],cutoff,
                gap["every_positive_square_side_squared_lower_bound"]))
        cursor=max(cursor,e["upper"])
    if cursor<Q(1):
        chambers.append(_sample_rational_chamber(
            C,D,cursor,Q(1),cutoff,
            gap["every_positive_square_side_squared_lower_bound"]))
    if not clusters and len(chambers)!=1:
        raise AssertionError("root-free interval should be a single chamber")
    events=[]
    for e in clusters:
        if e["lower"]<Q(0) or e["upper"]>Q(1):
            raise AssertionError("event root outside requested interval")
        rational=e["lower"]==e["upper"]
        item={
            "lower":str(e["lower"]),"upper":str(e["upper"]),
            "exact_rational_time":str(e["lower"]) if rational else None,
            "irrational_or_merged_event_slice_unresolved":not rational,
            "candidate_factors":sorted(set(e["factors"])),
        }
        if rational:
            result=census_polygon(interpolated_polygon(C,D,e["lower"]),
                                  cutoff_squared=cutoff)
            item["exact_event_geometric_square_orbits"]=result.count
        else:
            item["exact_event_geometric_square_orbits"]=None
        events.append(item)
    report["possible_event_clusters"]=events
    report["rational_event_free_chambers"]=chambers
    # Endpoint counts independently check the polygon census. In a
    # fully event-free situation they also check the global constancy.
    endpoints=[]
    for t in (Q(0),Q(1)):
        res=census_polygon(interpolated_polygon(C,D,t),cutoff_squared=cutoff)
        endpoints.append(res.count)
    report["endpoint_geometric_square_orbits"]=endpoints
    if not events:
        if endpoints!=[chambers[0]["geometric_square_orbits"]]*2:
            raise AssertionError("endpoint count contradicts root-free event cover")
        report["status"]="PASS: complete all-real-time geometric square count"
        report["complete_all_real_time_geometric_square_orbits"]=endpoints[0]
        report["complete_all_real_time_quotient_parity"]=endpoints[0]%2
    else:
        report["status"]="PASS: complete event-free chamber cover; algebraic event slices separately flagged"
        report["complete_all_real_time_geometric_square_orbits"]=None
        report["complete_all_real_time_quotient_parity"]=None
    report["scope"]="Rational affine strictly convex polygon motions (after certified exact straight-edge compression), with uniform nonadjacent-edge gap. Complete only away from isolated candidate-event brackets unless every event slice is separately certified. No nonconvex/all-Jordan theorem."
    return report


def demo():
    from big_peg_pair_death_v254 import (
        family_vertices,ROOT_LEFT,ROOT_RIGHT,
    )
    from big_peg_multichart_v249 import ppoint_at
    P=family_vertices()
    a=tuple(ppoint_at(v,ROOT_LEFT) for v in P)
    b=tuple(ppoint_at(v,ROOT_RIGHT) for v in P)
    return audit_moving_convex_polygon(a,b)


if __name__=="__main__":
    print(json.dumps(demo(),sort_keys=True,indent=2))
