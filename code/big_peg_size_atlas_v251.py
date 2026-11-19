"""v251: exact SIZE-AWARE, BRANCH-NEUTRAL square-existence atlas.

Extends v249's continuous-root atlas: a proof that at least one square
exists above a side cutoff does NOT require the chosen square at chart
handoffs to be the same geometric square. Individual Cramer square
branches may compete, so the verifier checks a complete closed cover
of the parameter interval, validating every chart by rational Bernstein
polynomial signs, and never infers existence from sampled times.

On v250's genuinely moving convex non-collinear nine-gon, two distinct
squares exchange WHICH IS LARGER over [1/2,501/1000]. An exact high-size
atlas proves that for EVERY real time in that window there is an
inscribed square with side^2 > 94257/10000 = 9.4257, even though NEITHER
of the two selected branches individually satisfies this cutoff on
the entire window. The equality of the two squares' sizes occurs at
one rigorously isolated crossing, with distinct square vertices.

The full [0,1] family is separately certified with a discontinuously
selected witness and side^2>9; a continuous witness also exists via
v250, and no new generic homotopy parity theorem is claimed.

ALL proof decisions use exact fractions and proven polynomial signs.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
import json

from big_peg_multichart_v249 import (
    Z,O,pol,padd,psub,pneg,pmul,pscale,peval,certify_sign,
    chart_polynomials,root_at,certify_chart,corner_numerators,
    pvsub,pvdot,HALF,
)
from big_peg_corner_v250 import (
    family_vertices, prove_jordan_family, LEFT_CHART, RIGHT_CHART,
)

# The first half continues the v250 square. The second half additionally
# contains a different large square near the opposite corner of the 9-gon.
BRANCH_A=LEFT_CHART
BRANCH_B=(1,3,5,8)
WINDOW_START=HALF
WINDOW_END=Q(501,1000)
ATLAS_SWITCH=Q(10011,20000)  # 0.50055
HIGH_S2=Q(94257,10000)      # 9.4257; both branches > at the switch
FULL_S2=Q(9)
CROSS_LEFT=Q(500549,1000000)
CROSS_RIGHT=ATLAS_SWITCH


@dataclass(frozen=True)
class ExistentialAtlas:
    start: Q
    end: Q
    side_squared_strictly_above: Q
    certificates: tuple
    distinct_seam_squares: int
    fully_covered: bool


def certify_existence_atlas(P,intervals,cutoff,*,start=Q(0),end=Q(1)):
    """Certify: at EACH time, there EXISTS a boundary square side^2>cutoff.

    Intervals are CLOSED rational time spans and can overlap. Each chart
    independently guarantees a square throughout its interval, including
    both endpoints. Unlike continuous-branch continuation in v249, the
    square at a handoff can jump to a different geometric square.

    Failure conditions: no intervals, invalid cutoff, un-covered endpoints,
    any open time gap, or a chart whose determinant/edge membership or
    cutoff cannot be certified by the inherited exact polynomial verifier.
    """
    if type(start) is not Q or type(end) is not Q or not Q(0)<=start<end<=Q(1):
        raise ValueError("exact rational nonempty target interval needed")
    if type(cutoff) is not Q or cutoff<=0:
        raise ValueError("exact positive squared-side cutoff required")
    if not intervals:
        raise ValueError("an empty square atlas cannot prove existence")
    normalized=[]
    for entry in intervals:
        if len(entry)!=3:
            raise ValueError("each interval must be (left,right,assignment)")
        left,right,edges=entry
        if type(left) is not Q or type(right) is not Q or not start<=left<right<=end:
            raise ValueError("chart times must be exact rationals inside target")
        normalized.append((left,right,edges))
    normalized.sort(key=lambda item:(item[0],item[1]))
    if normalized[0][0]!=start:
        raise ValueError("square atlas misses its first time")
    proof=[]
    covered=start
    switches=0
    previous=None
    for a,b,edges in normalized:
        if a>covered:
            raise ValueError("there is an uncovered open time interval")
        chart=certify_chart(P,edges,a,b,side2_limit=cutoff)
        if previous is not None and a==previous.right:
            # This is metadata about the chosen witness, NOT a new
            # continuity condition. Changing witnesses is legitimate.
            if chart.left_root!=previous.right_root:
                switches+=1
        previous=chart
        proof.append(chart)
        covered=max(covered,b)
    if covered!=end:
        raise ValueError("square atlas misses its final time")
    return ExistentialAtlas(start,end,cutoff,tuple(proof),switches,True)


def square_side2(root):
    return root[2]*root[2]+root[3]*root[3]


def side_comparison_polynomial(P,first,second):
    """Unnormalized exact polynomial with sign = side_A^2-side_B^2.

    For det(M_A)*det(M_B)!=0, the actual difference of squared sides is
    G/(D_A^2 D_B^2); this denominator is STRICTLY positive. Hence all
    comparisons reduce to the sign of a rational polynomial.
    """
    da,na=chart_polynomials(P,first)
    db,nb=chart_polynomials(P,second)
    sa=padd(pmul(na[2],na[2]),pmul(na[3],na[3]))
    sb=padd(pmul(nb[2],nb[2]),pmul(nb[3],nb[3]))
    return psub(pmul(sa,pmul(db,db)),pmul(sb,pmul(da,da)))


def derivative(p):
    return pol(*(Q(i)*p[i] for i in range(1,len(p)))) if len(p)>1 else Z


def exact_switch_proof(P):
    """Certify a unique equal-size CROSSING between two DISTINCT squares.

    The positive-square existence at every time uses a *separate*
    high-cutoff atlas. The size-comparison proof adds a rigorous
    monotonicity/IVT event statement for two chosen competing branches.
    It never claims they are globally maximizing over all squares.
    """
    # Both squares exist throughout the WHOLE observation window.
    a=certify_chart(P,BRANCH_A,WINDOW_START,WINDOW_END,side2_limit=FULL_S2)
    b=certify_chart(P,BRANCH_B,WINDOW_START,WINDOW_END,side2_limit=FULL_S2)
    # Even the higher cutoff is satisfied on a cover, though neither
    # individual branch meets it at both ends.
    atlas=certify_existence_atlas(
        P,
        ((WINDOW_START,ATLAS_SWITCH,BRANCH_A),
         (ATLAS_SWITCH,WINDOW_END,BRANCH_B)),
        HIGH_S2,start=WINDOW_START,end=WINDOW_END,
    )
    if atlas.distinct_seam_squares!=1:
        raise AssertionError("high-size atlas should switch between different squares")

    # These strict endpoint failures show the atlas' change of witness
    # is useful: neither selected branch alone proves this high cutoff.
    if not (square_side2(a.right_root)<HIGH_S2 and
            square_side2(b.left_root)<HIGH_S2):
        raise AssertionError("selected witness branches do not exchange high-cutoff roles")
    if not (square_side2(a.left_root)>square_side2(b.left_root) and
            square_side2(a.right_root)<square_side2(b.right_root)):
        raise AssertionError("no endpoint ordering reversal")

    G=side_comparison_polynomial(P,BRANCH_A,BRANCH_B)
    if not (peval(G,WINDOW_START)>0 and peval(G,WINDOW_END)<0):
        raise AssertionError("incorrect exact polynomial crossing signs")
    # Strict derivative sign makes the equal-size crossover UNIQUE
    # on the full observation interval, not just existent by IVT.
    derivative_checks=certify_sign(pneg(derivative(G)),
                                  WINDOW_START,WINDOW_END,strict=True)
    if not (peval(G,CROSS_LEFT)>0 and peval(G,CROSS_RIGHT)<0):
        raise AssertionError("exact crossover isolating rational bracket fails")

    # The two square configurations are geometrically distinct: at the
    # atlas switch they have *different* first corner coordinates.
    ad,an=chart_polynomials(P,BRANCH_A)
    bd,bn=chart_polynomials(P,BRANCH_B)
    root_a=root_at(ad,an,ATLAS_SWITCH)
    root_b=root_at(bd,bn,ATLAS_SWITCH)
    if root_a[:2]==root_b[:2]:
        raise AssertionError("two supposed square orbits have same first corner")
    if not (square_side2(root_a)>HIGH_S2 and
            square_side2(root_b)>HIGH_S2):
        raise AssertionError("both square branches must exceed cutoff at switch")
    return {
        "atlas":atlas,
        "crossing_polynomial":G,
        "strict_derivative_sign_checks":derivative_checks,
        "candidate_A_side2_at_start":square_side2(a.left_root),
        "candidate_B_side2_at_start":square_side2(b.left_root),
        "candidate_A_side2_at_end":square_side2(a.right_root),
        "candidate_B_side2_at_end":square_side2(b.right_root),
        "unique_crossover_bracket":(CROSS_LEFT,CROSS_RIGHT),
        "A_and_B_distinct_at_seam":True,
    }


def audit():
    P=family_vertices()
    jordan=prove_jordan_family(P)
    narrow=exact_switch_proof(P)

    # A deliberately DISCONTINUOUS selection of large-square witnesses.
    # The geometric polygon family already has a continuous square
    # branch by v250, so we claim only that one does NOT NEED a globally
    # continuous selected witness to establish a uniform size cutoff.
    whole=certify_existence_atlas(
        P,
        ((Q(0),HALF,RIGHT_CHART),
         (HALF,ATLAS_SWITCH,BRANCH_A),
         (ATLAS_SWITCH,WINDOW_END,BRANCH_B),
         (WINDOW_END,Q(1),BRANCH_A)),
        FULL_S2)
    if whole.distinct_seam_squares!=2:
        raise AssertionError("selected witness should jump at two handoffs")
    # The family lies in radius-3 disk, hence true inradius<=3. A side
    # strictly above 3 gives J>1 over the full family regardless of
    # which of its numerous square configurations is chosen.
    if not jordan["enclosed_in_radius_3_disk"]:
        raise AssertionError("no family-wide inradius upper bound")
    return {
        "status":"PASS: exact size-aware multi-branch square-existence atlas",
        "target_short_interval":[str(WINDOW_START),str(WINDOW_END)],
        "short_interval_side_squared_strictly_above":str(HIGH_S2),
        "short_interval_chart_count":len(narrow["atlas"].certificates),
        "short_interval_different_square_handoffs":narrow["atlas"].distinct_seam_squares,
        "candidates_reverse_relative_size":True,
        "unique_equal_size_crossing_bracket":[str(x) for x in narrow["unique_crossover_bracket"]],
        "polynomial_strict_derivative_checks":narrow["strict_derivative_sign_checks"],
        "candidate_sides_squared_at_window_endpoints":{
            key:str(narrow[key])
            for key in ("candidate_A_side2_at_start","candidate_B_side2_at_start",
                        "candidate_A_side2_at_end","candidate_B_side2_at_end")
        },
        "whole_family_interval":[str(whole.start),str(whole.end)],
        "whole_family_side_squared_strictly_above":str(FULL_S2),
        "whole_family_chart_count":len(whole.certificates),
        "whole_family_different_square_handoffs":whole.distinct_seam_squares,
        "restricted_family_true_big_peg_J_strictly_above":"1",
        "scope":"Finite exact existential size atlas for a specified Jordan polygon homotopy; no proof of global maximal-orbit switching, generic branch extinction, parity invariance, positive universal Big Peg c, or Toeplitz.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
