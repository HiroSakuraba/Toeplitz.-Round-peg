"""v254: rigorously certified PAIR DEATH of positive-size square orbits.

This is a genuinely changing, strictly convex rational pentagon homotopy,
not a square-size cutoff event, not a collinear relabeling, and NOT a
singular 4x4 Cramer determinant. Two distinct positive-size geometric
square branches collide at a genuine polygon VERTEX at the unique
irrational/algebraic root alpha of a rational quadratic, then both
become infeasible. A third positive-size square persists.

Proved for ALL REAL t in [3/5,13/20] by exact rational Bernstein signs:
- every polygon is strictly convex/Jordan;
- all 3 selected Cramer determinants are nonzero;
- each candidate square has side_squared >1 throughout the interval;
- the two departing branches have seven strictly interior contact
  inequalities and a final contact inequality with exact sign of Q(t);
- these two branches coincide at Q(alpha)=0 in ALL FOUR square parameters,
  but are distinct for all other t in the interval;
- the collision corner is exactly vertex P_2(alpha);
- the third square stays strictly on all four selected edges and remains
  distinct at the collision.

Independent pre-existing v223 FULL, unfiltered geometric-square censuses
confirm exactly three orbits at rational times to the left and exactly
one on the right (including near-isolator samples). This does NOT
classify every possible orbit at every irrational time: the entire
all-time theorem concerns the named THREE branches; exact full counts
are asserted only at the rational censused slices. This is a genuine
2->1->0 branch-pair collision, not generic determinant rank loss.

No numerical proof decisions, interval libraries or a Toeplitz claim.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
import json

from polygon_orbit_census_v223 import (
    pt, cross, dot, sub, validate_jordan_polygon, census_polygon,
)
from big_peg_multichart_v249 import (
    pol,padd,psub,pneg,pmul,pscale,peval,
    ppoint_at,pvsub,pvcross,pvdot,corner_numerators,
    chart_polynomials,root_at,certify_sign,
)
from moving_polygon_events_v224 import divide, iszero
from singular_continuation_v229 import (
    ExactTime,root_sign,side_signs,
)
from big_peg_cutoff_ledger_v252 import (
    cutoff_numerator,
    _edge_parameter_inequalities,
)

START=F(3,5)
END=F(13,20)
ROOT_LEFT=F(64722,100000)
ROOT_RIGHT=F(64723,100000)
CUT=F(1)

# Real t moves only the first point's x-coordinate: x0 = 0.5355-t/5.
P0=((F(5355,10000),F(6134,10000)),
    (F(-5770,10000),F(10693,10000)),
    (F(-12571,10000),F(2233,10000)),
    (F(-8305,10000),F(-3562,10000)),
    (F(285,10000),F(-7637,10000)))

# Exactly three relevant selected geometric square branches:
PERSISTENT=(0,1,2,4)
COLLIDING_LEFT=(0,1,3,4)
COLLIDING_RIGHT=(0,2,3,4)

# The two disappearing branches' contact inequalities are respectively
# "second corner at or before end of edge 1" and
# "second corner at or after start of edge 2".
CRITICAL={COLLIDING_LEFT:3,COLLIDING_RIGHT:2}

# Zero Q(t) is an algebraic event where the two squares meet P_2(t).
QPOLY=pol(291143052416293,-496738335230000,72468080000000)
EVENT=ExactTime(QPOLY,ROOT_LEFT,ROOT_RIGHT)


def family_vertices():
    result=[(pol(x),pol(y)) for x,y in P0]
    result[0]=(psub(pol(P0[0][0]),pol(0,F(1,5))),pol(P0[0][1]))
    return tuple(result)


def _strict_sign(p,a=START,b=END):
    """Return a rigorously constant nonzero sign, or fail closed."""
    try:
        certify_sign(p,a,b,strict=True)
        return 1
    except ValueError:
        certify_sign(pneg(p),a,b,strict=True)
        return -1


def certify_jordan_family(P):
    """Prove strict convexity, simplicity and edge length >0 for ALL t.

    Exact axis-aligned edge enclosures of this AFFINE family are strictly
    separated for all NONADJACENT edge pairs. All consecutive turn
    cross products are strictly positive over the full closed interval.
    This is global, not a finite sampling of vertex positions.
    """
    if tuple(P)!=family_vertices():
        raise ValueError("v254 Jordan certification is scoped to exact pentagon")
    if len(P)!=5 or any(len(p)!=2 or any(len(c)>2 for c in p) for p in P):
        raise ValueError("five affine vertex trajectories required")
    turns=0
    for k in range(5):
        u=P[k];v=P[(k+1)%5];w=P[(k+2)%5]
        turns+=certify_sign(pvcross(pvsub(v,u),pvsub(w,v)),
                            START,END,strict=True)
        turns+=certify_sign(pvdot(pvsub(v,u),pvsub(v,u)),
                            START,END,strict=True)
    rect=[]
    for i in range(5):
        box=[]
        for axis in (0,1):
            coords=[peval(P[j][axis],t) for j in (i,(i+1)%5)
                    for t in (START,END)]
            box.append((min(coords),max(coords)))
        rect.append(tuple(box))
    gap=None
    for i,j in combinations(range(5),2):
        if j==i+1 or (i==0 and j==4):
            continue
        gaps=[max(rect[i][ax][0]-rect[j][ax][1],
                  rect[j][ax][0]-rect[i][ax][1]) for ax in (0,1)]
        distance=max(gaps)
        if distance<=0:
            raise AssertionError("nonadjacent moving edges not separated")
        gap=distance if gap is None else min(gap,distance)
    # Endpoint snapshots independently prove correct orientation.
    for t in (START,END):
        snap=tuple(ppoint_at(p,t) for p in P)
        validate_jordan_polygon(snap)
    return {"strict_turn_sign_checks":turns,
            "global_nonadjacent_edge_rect_gap":gap,
            "convex_Jordan_for_all_real_times":True}


def _contact_polynomials(P,ids):
    D,N=chart_polynomials(P,ids)
    sign=_strict_sign(D)
    checks=tuple(_edge_parameter_inequalities(P,ids,D,N,sign))
    assert len(checks)==8
    return D,N,checks


@dataclass(frozen=True)
class BranchAudit:
    edges: tuple
    determinant_sign: int
    other_strict_contact_checks: int
    minimum_side_squared_exceeds: F
    critical_contact_index: int|None
    event_factor: tuple|None
    event_factor_sign: int|None


def certify_branch(P,ids):
    """Recompute every rational contact and positive square-size condition."""
    if tuple(P)!=family_vertices() or ids not in (
        PERSISTENT,COLLIDING_LEFT,COLLIDING_RIGHT):
        raise ValueError("uncertified source family or square edge assignment")
    D,N,ineq=_contact_polynomials(P,ids)
    det_sign=_strict_sign(D)
    # Side²>1 for the ENTIRE real interval, including event time.
    side2_minus1=cutoff_numerator(D,N,CUT)
    certify_sign(side2_minus1,START,END,strict=True)
    critical=CRITICAL.get(ids)
    strict=0;factor=None;factor_sign=None
    for k,q in enumerate(ineq):
        if k==critical:
            quotient,remainder=divide(q,QPOLY)
            if not iszero(remainder):
                raise AssertionError("critical contact lacks exact event factor")
            factor_sign=_strict_sign(quotient)
            if factor_sign<=0:
                raise AssertionError("critical contact factor orientation is wrong")
            factor=quotient
            continue
        strict+=certify_sign(q,START,END,strict=True)
    if ids==PERSISTENT and strict!=8:
        raise AssertionError("persistent square does not have four interior contacts")
    if critical is not None and strict!=7:
        raise AssertionError("colliding square has another active contact")
    return BranchAudit(ids,det_sign,strict,CUT,critical,factor,factor_sign)


def certify_collision(P):
    """Prove two DISTINCT squares join at P_2(alpha), then leave edges.

    Nonzero Cramer determinants and a positive common side are checked.
    Exact polynomial divisibility by Q shows the roots agree in ALL FOUR
    square coordinates at alpha, not merely sharing one corner or size.
    """
    if tuple(P)!=family_vertices():
        raise ValueError("uncertified moving polygon")
    if not (_strict_sign(pneg(QPOLY))==-1 and
            peval(QPOLY,START)>0 and peval(QPOLY,END)<0):
        raise AssertionError("wrong sign of square-collision event polynomial")
    certify_sign(pneg(pol(-496738335230000,2*72468080000000)),
                 START,END,strict=True)  # -Q'(t)>0, unique zero
    if not (START<ROOT_LEFT<ROOT_RIGHT<END
            and peval(QPOLY,ROOT_LEFT)>0
            and peval(QPOLY,ROOT_RIGHT)<0
            and root_sign(QPOLY,EVENT)==0):
        raise AssertionError("invalid unique algebraic collision root")
    if side_signs(QPOLY,EVENT)!=(1,0,-1):
        raise AssertionError("unexpected event one-sided sign structure")
    a=certify_branch(P,COLLIDING_LEFT)
    b=certify_branch(P,COLLIDING_RIGHT)
    persistent=certify_branch(P,PERSISTENT)

    Da,Na,_=_contact_polynomials(P,COLLIDING_LEFT)
    Db,Nb,_=_contact_polynomials(P,COLLIDING_RIGHT)
    Dc,Nc,_=_contact_polynomials(P,PERSISTENT)
    root_quotients=[]
    for j in range(4):
        difference=psub(pmul(Na[j],Db),pmul(Nb[j],Da))
        quotient,remainder=divide(difference,QPOLY)
        if not iszero(remainder):
            raise AssertionError("two square roots differ at collision")
        root_quotients.append(quotient)
    # The FIRST coordinate differs for every t!=alpha, so the two
    # square configurations do not merely give duplicate edge labels
    # on the pre-collision interval.
    distinct_first=_strict_sign(root_quotients[0])
    if distinct_first==0:
        raise AssertionError("first-coordinate quotient vanishes")

    # At alpha, SECOND square corner sits precisely at polygon P2.
    vertex=P[2]
    x,y=corner_numerators(Na,1)
    for coord_num,coord_vertex in ((x,vertex[0]),(y,vertex[1])):
        delta=psub(coord_num,pmul(Da,coord_vertex))
        quotient,remainder=divide(delta,QPOLY)
        if not iszero(remainder):
            raise AssertionError("merged square misses actual pentagon vertex")

    # Persistent branch is a third DIFFERENT square even at alpha,
    # as certified by exact algebraic sign of first-corner difference.
    third=psub(pmul(Na[0],Dc),pmul(Nc[0],Da))
    if root_sign(third,EVENT)==0:
        raise AssertionError("persistent and dying square collide")
    if root_sign(Da,EVENT)==0 or root_sign(Db,EVENT)==0:
        raise AssertionError("this is NOT a determinant singularity")
    if root_sign(cutoff_numerator(Da,Na,CUT),EVENT)<=0:
        raise AssertionError("shared polygon vertex-square is not positive-size")
    return {
        "event":EVENT,
        "persistent":persistent,
        "branch_left":a,
        "branch_right":b,
        "strict_distinctness_quotient_sign":distinct_first,
        "root_identity_divisible_by_event_in_four_components":True,
        "shared_second_corner_at_pentagon_vertex_2":True,
        "determinants_nonzero_at_event":True,
    }


def _sample_snapshots(P):
    """Independent FULL v223 orbit census at four RATIONAL times.

    All positive squares are counted, including repeated-edge charts,
    BOTH orientations and degeneracy handling; if a positive-dimensional
    family appears, v223 raises rather than fabricate a finite count.
    """
    sample_times=(START,ROOT_LEFT,ROOT_RIGHT,END)
    snapshots=[]
    for t in sample_times:
        pts=tuple(ppoint_at(v,t) for v in P)
        result=census_polygon(pts)
        if any(s.side_squared<=CUT for s in result.squares):
            raise AssertionError("sampled positive square unexpectedly below scale 1")
        expected=3 if t<=ROOT_LEFT else 1
        if result.count!=expected:
            raise AssertionError("independent full orbit count changed at "+str(t))
        snapshots.append({
            "time":str(t),
            "all_positive_geometric_square_orbits":result.count,
            "minimum_square_side_squared":str(min(
                s.side_squared for s in result.squares)),
            "singular_charts_examined":result.singular_charts,
            "redundant_edge_labelings":result.redundant_labelings,
        })
    return snapshots


def audit():
    P=family_vertices()
    jordan=certify_jordan_family(P)
    event=certify_collision(P)
    snapshots=_sample_snapshots(P)
    # Explicit algebraic interval event: t<alpha: A and B are
    # feasible distinct squares, t=alpha: same geometric vertex square,
    # t>alpha: each fails one necessary segment membership condition.
    # The persistent square C is feasible for every real t.
    return {
        "status":"PASS: exact positive-size GEOMETRIC square pair death at a true convex pentagon corner",
        "pentagon_start_coordinates":[[str(x),str(y)] for x,y in P0],
        "moving_vertex_0_x":"5355/10000 - t/5",
        "certified_real_time_interval":[str(START),str(END)],
        "jordan":{k:str(v) for k,v in jordan.items()},
        "event_defining_quadratic":[str(x) for x in QPOLY],
        "unique_algebraic_event_bracket":[str(ROOT_LEFT),str(ROOT_RIGHT)],
        "event_second_corner":"polygon vertex index 2 (noncollinear corner)",
        "source_square_branches":[list(COLLIDING_LEFT),
                                  list(COLLIDING_RIGHT)],
        "third_persistent_square_branch":list(PERSISTENT),
        "dying_branch_side_squared_strictly_greater_than":"1 at every real time as supporting-line squares",
        "geometric_branch_count_near_event":"2 distinct positive-size squares before; 1 shared vertex square at event; neither dying branch feasible after",
        "total_positive_geometric_square_counts_at_rational_times":snapshots,
        "singularity_type":"two nonsingular Cramer charts glued at a polygon vertex; NO rank drop of their incidence matrices",
        "local_mod2_balance":"two geometric orbits disappear as a pair (2 mod 2 = 0)",
        "scope":"Exact local pair-death/one-sided branches, plus full unfiltered rational-slice censuses. NOT an all-time complete orbit census at irrational times, not general singular event cobordism, positive universal Big Peg or Toeplitz.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
