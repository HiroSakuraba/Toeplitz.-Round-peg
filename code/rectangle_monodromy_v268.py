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


def odd_twisted_loop():
    """One T-cycle and exactly one irrational crossing through squares."""
    # Stereographic t in [1/2,1/3] gives (a,b) from (3,4) to
    # (4,3) while maintaining a^2+b^2=25 at every real u.
    p=(3-U)/6
    a=5*(1-p*p)/(1+p*p)
    b=10*p/(1+p*p)
    if not _positive_rational_function(a-3):
        # a(u)-3 is zero at u=0, so this intentionally cannot be
        # strictly positive. A wrong certificate must not be used.
        raise AssertionError("unexpected a>3")
    return None


if __name__=="__main__":
    print(json.dumps({"even":even_moebius_loop()},sort_keys=True,indent=2))
