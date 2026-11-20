"""v263: exact all-time, unfiltered positive-square census on v254's pentagon.

Enumerate EVERY cyclically ordered weakly increasing edge quadruple, including
repeated edges; classify them by exact rational Bernstein inequalities for
ALL REAL times in [64722/100000,64723/100000]. Do not infer exhaustiveness
from sampled rational slices. If ANY chart cannot be disposed of, FAIL CLOSED.

Strict convexity of v254 implies that a nondegenerate square has its four
corners in cyclic polygon order, so these edge quadruples cover every
positive geometric square (up to rotation of labels).

A positive uniform side-length cutoff is derived from v254's exact
nonadjacent-edge rectangle separation: every square has corners on
nonadjacent edges and its diameter is sqrt(2) times its side. Every
square has s^2 >= gap^2/2, so eta=gap^2/4 excludes NO real square.

The code never treats a singular chart as empty just because det=0.
It uses an exact augmented-minor infeasibility certificate or explicitly
reports that it remains unresolved. Singular event analysis must also
certify isolated determinant roots; it never samples them away.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, product
import json

from polygon_orbit_census_v223 import census_polygon
from big_peg_pair_death_v254 import (
    family_vertices, certify_jordan_family, certify_branch,
    certify_collision, PERSISTENT, COLLIDING_LEFT, COLLIDING_RIGHT,
    ROOT_LEFT, ROOT_RIGHT, QPOLY,
)
from big_peg_multichart_v249 import (
    pol, pneg, psub, pmul, pscale, peval, certify_sign,
    ppoint_at,
)
from moving_polygon_events_v224 import (
    chart_polynomials as build_chart, determinant_poly,
    moving_point, pairminus, pairplus, residual_poly,
    neg, iszero,
)
from big_peg_cutoff_ledger_v252 import cutoff_numerator

START=ROOT_LEFT
END=ROOT_RIGHT
TRACKED=(PERSISTENT,COLLIDING_LEFT,COLLIDING_RIGHT)
ALL_EDGE_CHARTS=tuple(combinations_with_replacement(range(5),4))


def sign_constant(p, strict=True):
    """Return +-1 only when a Bernstein proof works for ALL real tau."""
    try:
        certify_sign(p,Q(0),Q(1),strict=strict,max_depth=12)
        return 1
    except ValueError:
        try:
            certify_sign(pneg(p),Q(0),Q(1),strict=strict,max_depth=12)
            return -1
        except ValueError:
            return None


def edge_parameter_sign_polynomials(chart,sign):
    D=chart.determinant
    return tuple(pscale(z,sign) for j in range(4)
                 for z in (chart.numerators[j],chart.upper_differences[j]))


def _affine_polynomial_system(start,end,ids):
    """Polynomial 4x4 residual matrix M and constant RHS b=-F(starts)."""
    n=len(start)
    verts=tuple(moving_point(a,b) for a,b in zip(start,end))
    origins=tuple(verts[i] for i in ids)
    directions=tuple(pairminus(verts[(i+1)%n],verts[i]) for i in ids)
    base=residual_poly(origins,orientation=1)
    cols=[]
    for j in range(4):
        r=list(origins)
        r[j]=pairplus(r[j],directions[j])
        f=residual_poly(r,orientation=1)
        cols.append(tuple(psub(f[i],base[i]) for i in range(4)))
    matrix=tuple(tuple(cols[j][i] for j in range(4))
                 for i in range(4))
    return matrix,tuple(neg(x) for x in base)


def minor(matrix,rows,cols):
    return determinant_poly(tuple(
        tuple(matrix[i][j] for j in cols)
        for i in rows))


def singular_inconsistency_witness(start,end,ids):
    """For det M identically zero, certify rank [M|b] > rank M.

    For any r=3,2,1,0, if ALL (r+1)-minors of M vanish
    identically yet one (r+1)-minor of [M|b] is uniformly
    nonzero, the system M u=b has NO solution at any real time.
    The witness is a rational polynomial strictly signed by Bernstein.
    """
    A,b=_affine_polynomial_system(start,end,ids)
    aug=tuple(row+(b[i],) for i,row in enumerate(A))
    for r in (3,2,1,0):
        size=r+1
        # Prove rank(M) <=r at every real time, by polynomial identities.
        if any(not iszero(minor(A,ri,ci))
               for ri in combinations(range(4),size)
               for ci in combinations(range(4),size)):
            continue
        for ri in combinations(range(4),size):
            for ci in combinations(range(5),size):
                if 4 not in ci:
                    continue
                witness=minor(aug,ri,ci)
                sign=sign_constant(witness)
                if sign is not None:
                    return {"rank_upper_bound":r,
                            "rows":list(ri),"augmented_columns":list(ci),
                            "sign":sign}
        return None
    return None


def adjacent_edges(i,j):
    return (i-j)%5 in (1,4)


def certified_side_cutoff(P):
    """Derive *geometric* strictly positive cutoff from v254 separation.

    Uses the proved global axis rectangle separation of nonadjacent
    edge trajectories. A positive square on a simple polygon has two
    corners supported on nonadjacent polygon edges (v250 lemma).
    Their separation is at most the square diagonal sqrt(2)*side.
    Therefore side^2 >= gap^2/2 > eta=gap^2/4.
    """
    jordan=certify_jordan_family(P)
    if not jordan["convex_Jordan_for_all_real_times"]:
        raise AssertionError("strictly convex Jordan family needed")
    gap=jordan["global_nonadjacent_edge_rect_gap"]
    if type(gap) is not Q or gap<=0:
        raise AssertionError("no uniform nonadjacent-edge separation")
    eta=gap*gap/Q(4)
    if not 0<eta<Q(1):
        raise AssertionError("unexpected cutoff scale")
    return eta,gap


def classify_chart(P,ids,eta,start,end):
    """Every chart is either certified excluded, tracked, or unresolved."""
    if tuple(P)!=family_vertices() or (start,end)!=(START,END):
        raise ValueError("only v254 pentagon in exact event bracket")
    if type(eta) is not Q or eta<=0:
        raise ValueError("positive exact size cutoff required")
    if tuple(ids) not in ALL_EDGE_CHARTS:
        raise ValueError("canonical increasing edge quadruple required")
    if max(ids.count(i) for i in set(ids))>=3:
        # Three distinct square vertices cannot lie on one line.
        return {"class":"three_collinear","witness":"geometric"}
    if len(set(ids))==2:
        # If two polygon edges are adjacent, a nondegenerate square
        # cannot have all its corners supported there: if pairs
        # form the two square diagonals the center would lie at
        # their intersection, which is on the boundary; two pairs
        # of opposite sides would require parallel supporting lines.
        i,j=sorted(set(ids))
        if adjacent_edges(i,j):
            return {"class":"adjacent_two_edges","witness":"convex_geometry"}

    P0=tuple(ppoint_at(v,start) for v in P)
    P1=tuple(ppoint_at(v,end) for v in P)
    chart=build_chart(P0,P1,ids,orientation=1,cutoff_squared=eta)
    D=chart.determinant
    if iszero(D):
        witness=singular_inconsistency_witness(P0,P1,ids)
        return ({"class":"singular_inconsistent","witness":witness}
                if witness is not None else
                {"class":"unresolved_singular","witness":None})

    detsign=sign_constant(D)
    if detsign is None:
        return {"class":"unresolved_determinant_event","witness":None}
    # A positive square on all four closed edges requires N_j/D
    # in [0,1], eight weak inequalities. A STRICT negative proof for
    # any one polynomial excludes that chart at EVERY real t.
    contacts=edge_parameter_sign_polynomials(chart,detsign)
    for k,q in enumerate(contacts):
        if sign_constant(q)==-1:
            return {"class":"off_edge","witness":k}
    # Any algebraic square with side²<=eta can be safely omitted:
    # the global geometric gap theorem proves no *positive* square
    # on the real polygon can have such a short side.
    size=cutoff_numerator(D,chart.numerators,eta)
    if sign_constant(size)==-1:
        return {"class":"below_certified_cutoff","witness":str(eta)}
    if ids in TRACKED:
        # Re-use the full-interval exact proof for each real time,
        # not a sampled or convenient interior root.
        branch=certify_branch(P,ids)
        if branch.minimum_side_squared_exceeds!=Q(1):
            raise AssertionError("tracked branch lost positivity")
        return {"class":"tracked","witness":str(branch.determinant_sign)}
    return {"class":"unresolved_other_candidate","witness":None}


def diagnostic():
    P=family_vertices()
    eta,gap=certified_side_cutoff(P)
    results={}
    classes=Counter()
    for ids in ALL_EDGE_CHARTS:
        result=classify_chart(P,ids,eta,START,END)
        key=",".join(str(i) for i in ids)
        results[key]=result
        classes[result["class"]]+=1
    return {
        "time":[str(START),str(END)],
        "all_weakly_increasing_edge_quadruples":len(ALL_EDGE_CHARTS),
        "uniform_nonadjacent_edge_separation":str(gap),
        "globally_safe_squared_side_cutoff":str(eta),
        "classes":dict(sorted(classes.items())),
        "unresolved_charts":{
            key:value for key,value in results.items()
            if value["class"].startswith("unresolved")},
        "full_classification":results,
    }


def audit():
    report=diagnostic()
    unresolved=report["unresolved_charts"]
    if unresolved:
        raise AssertionError(
            "v263 cannot claim all-time completeness: "+repr(unresolved))
    if report["classes"].get("tracked",0)!=3:
        raise AssertionError("not exactly three tracked algebraic root charts")
    collision=certify_collision(family_vertices())
    if not collision["root_identity_divisible_by_event_in_four_components"]:
        raise AssertionError("v254 collision not rigorously certified")
    # Complete independent rational slice census, *not* used to replace
    # the algebraic all-time classification.
    samples=[]
    for t,wanted in ((START,3),(END,1)):
        polygon=tuple(ppoint_at(v,t) for v in family_vertices())
        census=census_polygon(polygon)
        if census.count!=wanted:
            raise AssertionError("independent complete geometric census mismatch")
        samples.append({"time":str(t),"geometric_squares":wanted})
    report["independent_complete_rational_slices"]=samples
    report["all_real_time_positive_geometric_square_count"]={
        "before_algebraic_collision":3,
        "at_algebraic_collision":2,
        "after_algebraic_collision":1,
    }
    report["unique_square_event_polynomial"]=[str(q) for q in QPOLY]
    report["all_real_time_exhaustiveness_proved"]=True
    report["scope"]="Complete all-positive-square census ONLY for v254's moving convex pentagon in its narrow exact algebraic event window. Does NOT prove an arbitrary-Jordan square existence invariant or a universal Big Peg size floor."
    return report


if __name__=="__main__":
    import sys
    result=diagnostic() if "--diagnostic" in sys.argv else audit()
    print(json.dumps(result,sort_keys=True,indent=2))
