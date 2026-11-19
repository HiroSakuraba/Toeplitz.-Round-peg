"""v253: boundary-corrected square-count ledger + EXACT GLOBAL square-size bottleneck.

Only the certified v250 affine convex nine-gon, interval [1/2,501/1000],
and strict side² cutoff 94257/10000 are in scope.

v252 exhaustively classified ALL high-size square orbits at EVERY real
time: two and only two branches A and B can exceed the cutoff. v251
proves their maximum ALWAYS exceeds it. Thus the actual maximum squared
side over ALL inscribed squares, not merely selected witnesses, is
    M(t) = max(side_A²(t),side_B²(t)).
Exact rational Bernstein signs prove A strictly decreasing and B
strictly increasing over the whole interval. Their unique equal-size
crossing is therefore the UNIQUE GLOBAL MINIMIZER of M(t) on that
interval, with a rigorous algebraic root isolating bracket and narrow
rational squared-side enclosure.

Separately, v252 proves a strict size-filtered orbit-count sequence
1 -> 2 -> 1, with one enter (+1) and one exit (-1) event at distinct
algebraic cutoff roots. The INTEGER relative boundary-corrected count
    I(t) = N_>(t) - #entered_before(t) + #exited_by(t) = 1
is constant. For strict >, equality at a cutoff entry does not count
as entered, whereas equality at an exit is already excluded. Exact
cutoff-polynomial signs decide the event accounting for every RATIONAL
time without approximating the algebraic event times. The same
bookkeeping holds analytically at the event roots by monotonicity.

This is finite-branch event bookkeeping, NOT an unfiltered degree,
generic birth/death theorem, universal positive Big Peg constant, or
Toeplitz proof. All certificate comparisons use Fraction, not float.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
import json

from big_peg_multichart_v249 import (
    pol,padd,psub,pneg,pmul,pscale,peval,certify_sign,
    chart_polynomials,
)
from big_peg_corner_v250 import family_vertices,prove_jordan_family
from big_peg_size_atlas_v251 import (
    BRANCH_A,BRANCH_B,WINDOW_START,WINDOW_END,HIGH_S2,
    side_comparison_polynomial,certify_existence_atlas,
)
from big_peg_cutoff_ledger_v252 import (
    EXPECTED_CLASSIFICATION,B_ENTER_BRACKET,A_EXIT_BRACKET,
    cutoff_numerator,poly_derivative,unique_cutoff_event,
    exhaustive_cutoff_classification,exact_endpoint_census,
    audit as audit_v252,
)

# Tightened from v251's width 1e-6 bracket to width exactly 1e-11.
# Their rational signs are verified again, independently.
CROSS_LO=Q(50054929969,10**11)
CROSS_HI=Q(50054929970,10**11)
# Rigorous bound on MIN OVER TIME OF THE GLOBAL maximum square side²:
# 9.4257027585 < min_t M(t) < 9.4257027586
BOTTLE_S2_LO=Q(94257027585,10**10)
BOTTLE_S2_HI=Q(94257027586,10**10)


def _validate_scope(P,cutoff):
    if tuple(P)!=family_vertices() or type(cutoff) is not Q or cutoff!=HIGH_S2:
        raise ValueError("only exact certified v250 family and v252 cutoff allowed")


def _squared_side(P,assignment,t):
    """Rational square side² from a nonsingular exact Cramer chart."""
    D,N=chart_polynomials(P,assignment)
    d=peval(D,t)
    if d==0:
        raise ValueError("singular Cramer chart")
    return (peval(N[2],t)**2+peval(N[3],t)**2)/d**2


def _side_derivative_sign_polynomial(P,assignment):
    """Polynomial with the exact sign of d/dt(square side squared).

    Write s²=S/D², S=N_u²+N_v². Its derivative is
       (S' D - 2 S D')/D³.
    Multiplication by D gives (S'D-2SD')D; denominator D^4 >0.
    Hence the polynomial's sign is precisely the derivative's sign.
    """
    D,N=chart_polynomials(P,assignment)
    S=padd(pmul(N[2],N[2]),pmul(N[3],N[3]))
    return pmul(D,psub(
        pmul(poly_derivative(S),D),
        pscale(pmul(S,poly_derivative(D)),Q(2)),
    ))


@dataclass(frozen=True)
class GlobalEnvelopeCertificate:
    root_lower: Q
    root_upper: Q
    global_bottleneck_side2_lower: Q
    global_bottleneck_side2_upper: Q
    derivative_sign_checks_A: int
    derivative_sign_checks_B: int
    exact_root_bracket_signs: tuple
    witness_lower_at_crossing: tuple
    witness_upper_at_crossing: tuple


def certify_global_envelope(P):
    """Exactly locate the UNIQUE weakest time of the actual largest square.

    The global-maximum identity is conditional on the *separate* v252
    full 126-chart exclusion certificate and v251 positive-size atlas.
    The caller audit() first replays those independent prerequisites.
    """
    if tuple(P)!=family_vertices():
        raise ValueError("uncertified polygon homotopy")
    a,b=WINDOW_START,WINDOW_END
    # Strict A down / B up on the entire closed real interval.
    da=_side_derivative_sign_polynomial(P,BRANCH_A)
    db=_side_derivative_sign_polynomial(P,BRANCH_B)
    a_checks=certify_sign(pneg(da),a,b,strict=True)
    b_checks=certify_sign(db,a,b,strict=True)

    # Unique equal-size time: comparison difference has strict
    # negative derivative (v251, independently replayed here).
    G=side_comparison_polynomial(P,BRANCH_A,BRANCH_B)
    certify_sign(pneg(poly_derivative(G)),a,b,strict=True)
    if not (a<CROSS_LO<CROSS_HI<b
            and peval(G,CROSS_LO)>0
            and peval(G,CROSS_HI)<0):
        raise AssertionError("invalid exact root isolating bracket")
    if not (_squared_side(P,BRANCH_A,a)>_squared_side(P,BRANCH_B,a)
            and _squared_side(P,BRANCH_A,b)<_squared_side(P,BRANCH_B,b)):
        raise AssertionError("no monotone global maximum envelope switch")

    # At the root t*, A(t*)=B(t*), A decreases and B increases.
    # Since t* in (L,U):
    #  A(U)<A(t*)<A(L) and B(L)<B(t*)<B(U).
    # These are exact RATIONAL bounds on the irrational/algebraic value.
    lower_a=_squared_side(P,BRANCH_A,CROSS_HI)
    lower_b=_squared_side(P,BRANCH_B,CROSS_LO)
    upper_a=_squared_side(P,BRANCH_A,CROSS_LO)
    upper_b=_squared_side(P,BRANCH_B,CROSS_HI)
    if not (BOTTLE_S2_LO < max(lower_a,lower_b)
            <= min(upper_a,upper_b) < BOTTLE_S2_HI):
        raise AssertionError("sharp bottleneck enclosure not established")
    if not (BOTTLE_S2_LO>HIGH_S2>0):
        raise AssertionError("bottleneck must be strictly above ledger cutoff")
    return GlobalEnvelopeCertificate(
        CROSS_LO,CROSS_HI,BOTTLE_S2_LO,BOTTLE_S2_HI,
        a_checks,b_checks,(peval(G,CROSS_LO),peval(G,CROSS_HI)),
        (lower_a,lower_b),(upper_a,upper_b))


@dataclass(frozen=True)
class RelativeCountState:
    time: Q
    raw_large_square_count: int
    entered_cutoff_count: int
    exited_cutoff_count: int
    net_boundary_flux: int
    corrected_integer_index: int
    corrected_mod2_index: int
    largest_square_side_squared: Q
    maximizing_chart: tuple | None


def rational_time_state(P,time,cutoff=HIGH_S2):
    """Exact state for any rational time, including inside event brackets.

    The proof that the raw count is complete and that the two event
    polynomials have unique transverse crossings is supplied by audit().
    Because the roots are algebraic, this method uses exact polynomial
    signs rather than substituting rounded event times.

    For the strict > filter, +1 entry is counted ONLY when F_B>0;
    -1 exit is counted AS SOON AS F_A<=0 (including at equality).
    """
    _validate_scope(P,cutoff)
    if type(time) is not Q or not WINDOW_START<=time<=WINDOW_END:
        raise ValueError("time must be exact Fraction in certified window")
    Da,Na=chart_polynomials(P,BRANCH_A)
    Db,Nb=chart_polynomials(P,BRANCH_B)
    ga=peval(cutoff_numerator(Da,Na,cutoff),time)
    gb=peval(cutoff_numerator(Db,Nb,cutoff),time)
    if peval(Da,time)==0 or peval(Db,time)==0:
        raise ValueError("singular candidate at requested time")
    active_a=int(ga>0)
    active_b=int(gb>0)
    raw_count=active_a+active_b
    entered=active_b
    exited=int(ga<=0)
    flux=entered-exited
    corrected=raw_count-flux
    if corrected!=1:
        raise AssertionError("boundary-corrected index failed")
    sides=(_squared_side(P,BRANCH_A,time),
           _squared_side(P,BRANCH_B,time))
    best=max(sides)
    return RelativeCountState(
        time,raw_count,entered,exited,flux,corrected,corrected%2,best,
        BRANCH_A if sides[0]>sides[1] else
        BRANCH_B if sides[1]>sides[0] else None
    )


def _event_proof(P):
    """Two rigorously signed, unique algebraic boundary-crossing events."""
    a,b=WINDOW_START,WINDOW_END
    enter=unique_cutoff_event(P,BRANCH_B,HIGH_S2,a,b,
                               B_ENTER_BRACKET,direction=+1)
    exit=unique_cutoff_event(P,BRANCH_A,HIGH_S2,a,b,
                              A_EXIT_BRACKET,direction=-1)
    if not enter["bracket"][1]<exit["bracket"][0]:
        raise AssertionError("entry event does not precede exit event")
    return enter,exit


@lru_cache(maxsize=1)
def audit():
    """Recompute prerequisites BEFORE asserting a global or relative theorem."""
    P=family_vertices()
    geometry=prove_jordan_family(P)
    if not (geometry["strictly_convex_and_jordan_for_all_times"]
            and geometry["enclosed_in_radius_3_disk"]):
        raise AssertionError("changing polygons are not certified convex/Jordan")
    census=audit_v252()  # complete 126 charts and v223 endpoint counts
    if census["classifications"]!=EXPECTED_CLASSIFICATION:
        raise AssertionError("v252 exhaustive high-square ledger changed")
    if census["geometric_orbit_counts_strictly_above_cutoff_in_time_order"]!=[1,2,1]:
        raise AssertionError("v252 complete cutoff counts changed")
    # Independently use the v251 chart atlas to prove there is ALWAYS
    # at least one A/B square strictly above cutoff. Consequently the
    # v252 two-orbit exhaustion gives the TRUE GLOBAL s_max² as max(A,B).
    from big_peg_size_atlas_v251 import ATLAS_SWITCH
    existence=certify_existence_atlas(P,(
        (WINDOW_START,ATLAS_SWITCH,BRANCH_A),
        (ATLAS_SWITCH,WINDOW_END,BRANCH_B)),
        HIGH_S2,start=WINDOW_START,end=WINDOW_END)
    if not existence.fully_covered or len(existence.certificates)!=2:
        raise AssertionError("some real parameter time lacks high square")
    envelope=certify_global_envelope(P)
    enter,exit=_event_proof(P)

    # These exact rational times lie in the three REAL event chambers,
    # with additional two samples straddling the size crossover.
    samples=(
        WINDOW_START,
        Q(500545,10**6),
        Q(10011,20000),
        WINDOW_END,
    )
    states=tuple(rational_time_state(P,time) for time in samples)
    if tuple(s.raw_large_square_count for s in states)!=(1,2,2,1):
        raise AssertionError("incorrect chamber counts")
    if tuple(s.net_boundary_flux for s in states)!=(0,1,1,0):
        raise AssertionError("incorrect boundary-event accounting")
    if any(s.corrected_integer_index!=1 for s in states):
        raise AssertionError("relative integer count failed")

    return {
        "status":"PASS: exact boundary-corrected parity AND global maximum-square bottleneck",
        "interval":[str(WINDOW_START),str(WINDOW_END)],
        "cutoff_square_side_squared":str(HIGH_S2),
        "complete_high_square_charts_categorized":census["all_distinct_edge_charts_checked"],
        "number_of_possible_high_square_branches":2,
        "global_maximum_square_side_squared_identity":"smax²(t) = max(sA²(t),sB²(t)) for EVERY real t in interval",
        "global_maximum_unique_minimizer":"unique algebraic A=B side-length crossover",
        "global_minimizer_time_bracket":[str(CROSS_LO),str(CROSS_HI)],
        "minimum_of_global_max_square_side_squared_rigorous_interval":[
            str(envelope.global_bottleneck_side2_lower),
            str(envelope.global_bottleneck_side2_upper)],
        "A_strictly_decreasing_full_interval":envelope.derivative_sign_checks_A>0,
        "B_strictly_increasing_full_interval":envelope.derivative_sign_checks_B>0,
        "exact_bottleneck_lower_witnesses":[str(x) for x in envelope.witness_lower_at_crossing],
        "exact_bottleneck_upper_witnesses":[str(x) for x in envelope.witness_upper_at_crossing],
        "integer_signed_boundary_event_coefficients":{
            "B_enters_strict_cutoff":1,
            "A_exits_strict_cutoff":-1},
        "entry_root_bracket":[str(x) for x in enter["bracket"]],
        "exit_root_bracket":[str(x) for x in exit["bracket"]],
        "exact_rational_chamber_samples":[{
            "time":str(s.time),
            "strict_large_square_count":s.raw_large_square_count,
            "entered":s.entered_cutoff_count,
            "exited":s.exited_cutoff_count,
            "net_flux":s.net_boundary_flux,
            "corrected_integer_index":s.corrected_integer_index,
            "corrected_mod2_index":s.corrected_mod2_index,
            "largest_square_side_squared":str(s.largest_square_side_squared),
        } for s in states],
        "corrected_integer_count_for_ALL_real_times":"N_>(t) - entry_before(t) + exit_by(t) = 1 (strict-cutoff boundary convention)",
        "corrected_mod2_count_for_ALL_real_times":"(N_>(t)+total_crossings_before_or_at(t)) mod 2 = 1, with strict cutoff conventions",
        "restricted_Big_Peg_consequence":"J(t)>sqrt(94257027585/10^10)/3 > 1.02337 for every polygon in this interval",
        "scope":"Only a complete fixed-cutoff two-branch regular event ledger and algebraic global bottleneck in this affine convex nine-gon interval; no invariant across arbitrary orbit births/deaths, no universal c, no Toeplitz proof.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
