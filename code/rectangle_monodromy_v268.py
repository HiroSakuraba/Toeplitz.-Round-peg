"""v268: exact cyclic-monodromy audit on the actual rectangle locus.

The positive square-peg obstruction from diagonal "spin" is not the
individual chord-line w1. It is the sign representation of cyclic
relabeling on the diagonal dot-product along a CONTINUOUS family of
nondegenerate inscribed rectangles.

For any continuous lifted rectangle loop on a Jordan curve whose
cyclically ordered corners return with monodromy T^k, the diagonal
dot f=d_02.dot(d_13) obeys f(end)=(-1)^k f(start).
If k odd, a positive-size square is forced. If k even, this
sign argument is silent. Crucially, the odd-monodromy loop is NOT
known to exist for every Jordan curve.

We supply *exact rational-function* continuous rectangle loops on
a fixed circle radius 5:
  (1) NONSQUARE loop: fixed a=3,b=4, rotate 180 deg; k=2,
      individual diagonal lines both reverse, tensor monodromy 0,
      f=28 identically and NO square along this loop.
  (2) SQUARE-FORCING loop: vary a:3->4 and b:4->3 on
      a^2+b^2=25 while rotating 90 deg; k=1,
      f goes 28 to -28 and has one unique irrational zero
      u=9-6*sqrt(2). Both diagonal lengths are fixed at 10,
      so the square cannot collapse. The change in aspect
      is required to close the loop after an odd shift.

Identities are proved as exact rational polynomial congruences on
whole real intervals, not by finite numeric sampling. No claim
of universal Jordan square forcing or Floer theory.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from moving_polygon_events_v224 import (
    polynomial,plus,minus,neg,times,evaluate,divide,poly_gcd,
    iszero,isolate_roots,
)
from big_peg_multichart_v249 import certify_sign
import json


def P(*z):
    return polynomial(z)


@dataclass(frozen=True,eq=False)
class Rat:
    n:tuple
    d:tuple=(Q(1),)
    def __post_init__(self):
        p,q=polynomial(self.n),polynomial(self.d)
        if iszero(q):
            raise ZeroDivisionError("rational function zero denominator")
        if iszero(p):
            object.__setattr__(self,"n",P(0))
            object.__setattr__(self,"d",P(1))
            return
        g=poly_gcd(p,q)
        if len(g)>1:
            p,rem=divide(p,g)
            q,rem2=divide(q,g)
            if not iszero(rem) or not iszero(rem2):
                raise AssertionError("nonexact polynomial gcd division")
        lead=q[-1]
        object.__setattr__(self,"n",polynomial(x/lead for x in p))
        object.__setattr__(self,"d",polynomial(x/lead for x in q))
    @staticmethod
    def of(x):
        return x if isinstance(x,Rat) else Rat(P(x))
    def __add__(self,other):
        o=Rat.of(other)
        return Rat(plus(times(self.n,o.d),times(o.n,self.d)),
                   times(self.d,o.d))
    def __radd__(self,other):return self.__add__(other)
    def __neg__(self):return Rat(neg(self.n),self.d)
    def __sub__(self,other):return self+(-Rat.of(other))
    def __rsub__(self,other):return Rat.of(other)-self
    def __mul__(self,other):
        o=Rat.of(other)
        return Rat(times(self.n,o.n),times(self.d,o.d))
    def __rmul__(self,other):return self.__mul__(other)
    def __truediv__(self,other):
        o=Rat.of(other)
        return Rat(times(self.n,o.d),times(self.d,o.n))
    def __rtruediv__(self,other):return Rat.of(other)/self
    def __eq__(self,other):
        try:o=Rat.of(other)
        except (TypeError,ValueError):return NotImplemented
        return iszero(minus(times(self.n,o.d),times(o.n,self.d)))
    def at(self,u):
        den=evaluate(self.d,u)
        if not den:raise ZeroDivisionError("pole at evaluation point")
        return evaluate(self.n,u)/den


ZERO=Rat.of(0)
ONE=Rat.of(1)
U=Rat(P(0,1))


def _vadd(p,q):return tuple(a+b for a,b in zip(p,q))
def _vsub(p,q):return tuple(a-b for a,b in zip(p,q))
def _vneg(p):return tuple(-a for a in p)
def _vdot(p,q):return sum((a*b for a,b in zip(p,q)),ZERO)
def _eval_v(p,u):return tuple(x.at(u) for x in p)


def rectangle(a,b):
    """CCW four corners on centered axis-parallel rectangle."""
    a,b=Rat.of(a),Rat.of(b)
    return ((a,b),(-a,b),(-a,-b),(a,-b))


def rotate(corners,c,s):
    """Exact rational rotation in plane when c^2+s^2=1."""
    if c*c+s*s!=ONE:
        raise ValueError("rotation fractions do not lie on unit circle")
    return tuple((c*x-s*y,s*x+c*y) for x,y in corners)


def unit_quarter(u=U):
    """Rational stereographic quarter-arc (0<=u<=1), 0 to 90 degrees."""
    return ((ONE-u*u)/(ONE+u*u),2*u/(ONE+u*u))


def rectangle_defect(corners):
    """(midpoint-x,midpoint-y,diagonal lengths difference, dot)."""
    if len(corners)!=4:
        raise ValueError("four rectangle corners required")
    p0,p1,p2,p3=corners
    d0=_vsub(p2,p0);d1=_vsub(p3,p1)
    center=_vsub(_vadd(p0,p2),_vadd(p1,p3))
    return (center[0]/2,center[1]/2,
            _vdot(d0,d0)-_vdot(d1,d1),_vdot(d0,d1))


def diagonal_lengths_squared(corners):
    p0,p1,p2,p3=corners
    return (_vdot(_vsub(p2,p0),_vsub(p2,p0)),
            _vdot(_vsub(p3,p1),_vsub(p3,p1)))


def corner_radius_squared(corners):
    return tuple(_vdot(p,p) for p in corners)


def _strictly_positive(p):
    """Prove rational polynomial p>0 on entire CLOSED [0,1]."""
    certify_sign(polynomial(p),Q(0),Q(1),strict=True,max_depth=12)
    return True


def _positive_rational_function(r):
    """Certified >0 over ALL real u in [0,1]."""
    if not isinstance(r,Rat):
        raise TypeError("exact rational function needed")
    if evaluate(r.d,0)>0:
        return _strictly_positive(r.n) and _strictly_positive(r.d)
    return _strictly_positive(neg(r.n)) and _strictly_positive(neg(r.d))


def _assert_rectangle_identically(corners,expect_dot=None):
    defect=rectangle_defect(corners)
    if any(z!=ZERO for z in defect[:3]):
        raise AssertionError("rectangle equations not identically satisfied")
    if any(z!=Rat.of(100) for z in diagonal_lengths_squared(corners)):
        raise AssertionError("diagonal length changed from 10")
    if any(z!=Rat.of(25) for z in corner_radius_squared(corners)):
        raise AssertionError("corner left fixed radius-5 circle")
    if expect_dot is not None and defect[3]!=Rat.of(expect_dot):
        raise AssertionError("incorrect perpendicularity sign/constant")
    return defect[3]


def _monodromy(start,end):
    """Return T^k if an exact same-set loop closes with that cyclic lift."""
    if len(start)!=4 or len(end)!=4:
        raise ValueError("four corner endpoints required")
    if len(set(start))!=4 or len(set(end))!=4:
        raise ValueError("only nondegenerate corner loops admissible")
    possible=[k for k in range(4) if all(
        end[j]==start[(j+k)%4] for j in range(4))]
    if len(possible)!=1:
        raise ValueError("end not same geometric rectangle with C4 lift")
    return possible[0]


def diagonal_dot_character(cyclic_power):
    if type(cyclic_power) is not int:
        raise ValueError("C4 monodromy power must be integer")
    return -1 if cyclic_power%2 else 1


def _pointset(P,u):return tuple(_eval_v(p,u) for p in P)


def even_moebius_loop():
    """Non-square 180-degree loop, independent chord twists cancel."""
    a,b=Rat.of(3),Rat.of(4)
    base=rectangle(a,b)
    c0,s0=unit_quarter(U)
    path0=rotate(base,c0,s0)
    path1=rotate(base,-s0,c0) # second 90-degree arc
    g0=_assert_rectangle_identically(path0,28)
    g1=_assert_rectangle_identically(path1,28)
    if _pointset(path0,1)!=_pointset(path1,0):
        raise AssertionError("two rational arcs not glued continuously")
    start,end=_pointset(path0,0),_pointset(path1,1)
    k=_monodromy(start,end)
    if k!=2 or diagonal_dot_character(k)!=1:
        raise AssertionError("even loop monodromy mismatched")
    # Edges have lengths 6 and 8 throughout both arcs.
    min_side_sq=Q(36)
    return {
        "status":"CERTIFIED: non-square rectangle loop with two individually Möbius-twisted diagonals and TRIVIAL tensor twist",
        "base_rectangle_semiwidth":3,
        "base_rectangle_semiheight":4,
        "fixed_circumcircle_radius_squared":"25",
        "two_rational_circle_arcs":2,
        "cyclic_monodromy_T_power":k,
        "first_diagonal_reverses":True,
        "second_diagonal_reverses":True,
        "individual_chord_w1_values":[1,1],
        "tensor_w1":0,
        "diagonal_dot_constant":str(g0.at(Q(0))),
        "minimum_rectangle_side_squared":str(min_side_sq),
        "nondegenerate_square_count_on_loop":0,
        "all_real_parameter_identities_checked":True,
    }


def _nonnegative_rational_function(r):
    """Certify r>=0 throughout [0,1], including EXACT endpoint zeros."""
    if not isinstance(r,Rat):
        raise TypeError("exact rational function needed")
    den_positive=evaluate(r.d,Q(0))>0
    sign=1 if den_positive else -1
    certify_sign(r.d if sign==1 else neg(r.d),
                 Q(0),Q(1),strict=True,max_depth=12)
    certify_sign(r.n if sign==1 else neg(r.n),
                 Q(0),Q(1),strict=False,max_depth=12)
    return True


def odd_twisted_loop():
    """A genuine rational circle rectangle loop with ODD T monodromy.

    Exact algebraic square event: u=9-6sqrt(2), unique root of
    u^2-18u+9 in (0,1). Positively bounded sides prevent collapse.
    """
    # Stereographic p runs from 1/2 to 1/3. The semiwidth a grows
    # from 3 to 4 and semiheight b decreases from 4 to 3,
    # maintaining a^2+b^2=25 on the ENTIRE real interval.
    p=(3-U)/6
    a=5*(1-p*p)/(1+p*p)
    b=10*p/(1+p*p)
    if a*a+b*b!=Rat.of(25):
        raise AssertionError("radius-5 circle constraint changed")
    for r in (a,b):
        _positive_rational_function(r)
        _nonnegative_rational_function(r-3)
        _nonnegative_rational_function(4-r)
    c,s=unit_quarter(U)
    vertices=rotate(rectangle(a,b),c,s)
    f=_assert_rectangle_identically(vertices)
    expected=4*(b*b-a*a)
    if f!=expected:
        raise AssertionError("diagonal dot algebraic identity failed")
    start,end=_pointset(vertices,0),_pointset(vertices,1)
    k=_monodromy(start,end)
    if k!=1 or diagonal_dot_character(k)!=-1:
        raise AssertionError("odd rectangle loop did not return with T")
    if f.at(0)!=Q(28) or f.at(1)!=Q(-28):
        raise AssertionError("wrong certified end signs")
    # Denominator remains >0 on entire interval, so every zero of
    # f is captured by the numerator's complete Sturm root census.
    _positive_rational_function(Rat(f.d))
    event_poly=P(9,-18,1)
    quotient,remainder=divide(f.n,event_poly)
    if not iszero(remainder):
        raise AssertionError("expected irrational algebraic event factor missing")
    if len(quotient)<=0:
        raise AssertionError("missing event quotient")
    roots=isolate_roots(f.n,Q(0),Q(1),width=Q(1,2**22))
    exact_roots=isolate_roots(event_poly,Q(0),Q(1),width=Q(1,2**22))
    if (len(roots)!=1 or roots[0].rational or
            len(exact_roots)!=1 or exact_roots[0].rational):
        raise AssertionError("not exactly one irrational perpendicularity root")
    x,y=roots[0],exact_roots[0]
    if not (x.lower<y.upper and y.lower<x.upper):
        raise AssertionError("candidate root and derived minimal polynomial disagree")
    # At that instant, the 10-long diagonals are perpendicular,
    # giving a square with side^2=100/2=50.
    return {
        "status":"CERTIFIED: one odd-T geometric rectangle loop with exactly one positive algebraic square",
        "semiwidth_endpoints":[str(a.at(0)),str(a.at(1))],
        "semiheight_endpoints":[str(b.at(0)),str(b.at(1))],
        "fixed_circumcircle_radius_squared":"25",
        "cyclic_monodromy_T_power":k,
        "tensor_w1":1,
        "individual_diagonal_lines_globally_labeled":False,
        "diagonal_dot_endpoints":[str(f.at(0)),str(f.at(1))],
        "event_minimal_polynomial":[str(z) for z in event_poly],
        "event_algebraic_expression":"9 - 6*sqrt(2)",
        "event_isolation_bracket":[str(x.lower),str(x.upper)],
        "number_of_distinct_square_crossings":1,
        "square_side_squared_at_event":"50",
        "minimum_rectangle_side_squared":"36",
        "all_real_parameter_identities_checked":True,
        "unrestricted_toeplitz_proved":False,
    }


def monodromy_no_go():
    """Pure line-bundle implication, independent of numerical models.

    For a continuous loop of genuine nondegenerate rectangles,
    f(end)=(-1)^k f(start). If k is odd and f continuous, IVT
    gives an exact positive square because both diagonals stay
    nonzero. If all rectangles remain non-square, k must be even.

    IMPORTANT: The physical existence of an ODD loop on EVERY
    Jordan curve is NOT established by this abstract implication.
    """
    for k in range(4):
        if diagonal_dot_character(k)!=(-1 if k%2 else 1):
            raise AssertionError("C4 sign character failed")
    return {
        "sign_representation_C4":{"0":1,"1":-1,"2":1,"3":-1},
        "necessary_for_nontrivial_tensor_monodromy":"odd cyclic permutation of continuously lifted geometric rectangle",
        "odd_monodromy_forces_nonzero_size_square":True,
        "nonsquare_closed_loop_must_have_even_monodromy":True,
        "existence_of_odd_rectangle_loop_for_arbitrary_Jordan_curve_proved":False,
        "selected_halfturn_counterloop_has_both_individual_w1_equal_one":True,
    }


def audit():
    return {
        "even_no_square_rectangle_loop":even_moebius_loop(),
        "odd_algebraic_square_forcing_rectangle_loop":odd_twisted_loop(),
        "general_monodromy_character":monodromy_no_go(),
        "comparison_to_latest_Greene_Lobb_2026":"arXiv:2604.17116 proves positive-measure rectangle angles for arbitrary Jordan curves, but not a nondegenerate odd-monodromy rectangle loop, nor the required angle pi/2",
        "not_an_unrestricted_Toeplitz_proof":True,
    }


if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
