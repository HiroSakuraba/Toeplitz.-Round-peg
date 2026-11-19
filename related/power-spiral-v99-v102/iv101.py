"""
iv101: outward-rounded interval layer for the power-spiral face certificates.

Rewrite of the v100 transcendental layer.  Design rules:

1. Every arithmetic operation that touches an enclosure endpoint is performed
   under an explicit directed rounding context.  No step relies on
   "one ulp of half-even error is covered by next_minus".
2. Every series is a partial sum accumulated in interval arithmetic plus an
   explicitly bounded remainder interval.  No enclosure is obtained from an
   alternating-bracket argument, so no evaluator depends on an unstated
   term-monotonicity hypothesis.
3. Every evaluator asserts its validated argument domain.  Leaving the domain
   raises rather than returning silently.
4. pi is certified in-module by Machin's formula from the same arctan
   evaluator, not taken as a rational literal.
"""

from decimal import Decimal as D, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING

PREC = 40
getcontext().prec = PREC


class DomainError(Exception):
    pass


def _flo(f, *a):
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        return f(*a)


def _fhi(f, *a):
    with localcontext() as c:
        c.rounding = ROUND_CEILING
        return f(*a)


class IV:
    """Closed interval of Decimals, outward rounded."""

    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = D(lo)
        self.hi = D(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError((self.lo, self.hi))

    @staticmethod
    def _c(o):
        return o if isinstance(o, IV) else IV(o)

    def __add__(self, o):
        o = IV._c(o)
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            lo = self.lo + o.lo
            c.rounding = ROUND_CEILING
            hi = self.hi + o.hi
        return IV(lo, hi)

    __radd__ = __add__

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-IV._c(o))

    def __rsub__(self, o):
        return IV._c(o) - self

    def __mul__(self, o):
        o = IV._c(o)
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            lo = min(self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
            c.rounding = ROUND_CEILING
            hi = max(self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return IV(lo, hi)

    __rmul__ = __mul__

    def recip(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError((self.lo, self.hi))
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            a = D(1) / self.hi
            c.rounding = ROUND_CEILING
            b = D(1) / self.lo
        return IV(min(a, b), max(a, b))

    def __truediv__(self, o):
        return self * IV._c(o).recip()

    def __rtruediv__(self, o):
        return IV._c(o) / self

    def sqrt(self):
        if self.lo < 0:
            raise DomainError(('sqrt of negative', self.lo))
        with localcontext() as c:
            c.rounding = ROUND_FLOOR
            lo = self.lo.sqrt()
            c.rounding = ROUND_CEILING
            hi = self.hi.sqrt()
        return IV(lo, hi)

    def absmax(self):
        return max(abs(self.lo), abs(self.hi))

    def width(self):
        return _fhi(lambda: self.hi - self.lo)

    def __repr__(self):
        return f'[{self.lo},{self.hi}]'


def sym(bound):
    """Symmetric remainder interval [-bound, bound], bound assumed >= 0."""
    b = D(bound)
    if b < 0:
        raise ValueError(b)
    return IV(-b, b)


def _geom_tail(next_ub, ratio_ub):
    """Upper bound for next_ub * (1 + r + r^2 + ...) with r <= ratio_ub < 1."""
    if ratio_ub >= 1:
        raise DomainError(('tail ratio >= 1', ratio_ub))
    with localcontext() as c:
        c.rounding = ROUND_FLOOR
        den = D(1) - ratio_ub
        c.rounding = ROUND_CEILING
        return next_ub / den


# ----------------------------------------------------------------------
# Validated argument domains.  Measured ranges over the eight v100 faces
# are quoted; the asserted domain leaves margin on each side.
# ----------------------------------------------------------------------
DOM = {
    'exp':   (D(0),           D('1.6')),    # measured [0, 1.492538]
    'sin':   (D(0),           D('0.8')),    # measured [0, 0.716418]
    'cos':   (D(0),           D('0.8')),    # measured [0, 0.716418]
    'sinc':  (D(0),           D('0.8')),    # measured [0, 0.716418]
    'coshC': (D(0),           D('1.6')),    # measured [0, 1.492538]
    'L1':    (D('-0.8'),      D('2.5')),    # measured [-0.620527, 2.050882]
    'Efun':  (D('-1.0'),      D('3.0')),    # measured [-0.140899, 2.224946]
    'atanh': (D('-0.68'),     D('0.68')),   # induced by the L1 domain
}


def _chk(name, x):
    lo, hi = DOM[name]
    if not (lo <= x <= hi):
        raise DomainError((name, x, lo, hi))


# ----------------------------------------------------------------------
# exp(x) on [0, 1.6]
# ----------------------------------------------------------------------
N_EXP = 34


def exp_pt(x):
    _chk('exp', x)
    X = IV(x)
    term = IV(1)
    s = IV(1)
    for n in range(1, N_EXP + 1):
        term = term * X / n
        s = s + term
    # term now encloses x^N/N!.  Remainder is positive, ratio <= x/(N+2).
    nxt = _fhi(lambda: term.hi * x / D(N_EXP + 1))
    ratio = _fhi(lambda: x / D(N_EXP + 2))
    tail = _geom_tail(nxt, ratio)
    return s + IV(D(0), tail)


def exp_iv(v):
    return IV(exp_pt(v.lo).lo, exp_pt(v.hi).hi)


# ----------------------------------------------------------------------
# sin, cos, sinc on [0, 0.8]
# ----------------------------------------------------------------------
N_TRIG = 10


def sin_pt(x):
    _chk('sin', x)
    X = IV(x)
    X2 = X * X
    term = X
    s = X
    for k in range(1, N_TRIG + 1):
        term = term * X2 / D((2 * k) * (2 * k + 1))
        s = s + (term if k % 2 == 0 else -term)
    nxt = _fhi(lambda: term.absmax() * x * x / D((2 * N_TRIG + 2) * (2 * N_TRIG + 3)))
    ratio = _fhi(lambda: x * x / D((2 * N_TRIG + 4) * (2 * N_TRIG + 5)))
    return s + sym(_geom_tail(nxt, ratio))


def cos_pt(x):
    _chk('cos', x)
    X2 = IV(x) * IV(x)
    term = IV(1)
    s = IV(1)
    for k in range(1, N_TRIG + 1):
        term = term * X2 / D((2 * k - 1) * (2 * k))
        s = s + (term if k % 2 == 0 else -term)
    nxt = _fhi(lambda: term.absmax() * x * x / D((2 * N_TRIG + 1) * (2 * N_TRIG + 2)))
    ratio = _fhi(lambda: x * x / D((2 * N_TRIG + 3) * (2 * N_TRIG + 4)))
    return s + sym(_geom_tail(nxt, ratio))


def sinc_pt(x):
    """sin(x)/x by its own power series, so x = 0 needs no special case."""
    _chk('sinc', x)
    X2 = IV(x) * IV(x)
    term = IV(1)
    s = IV(1)
    for k in range(1, N_TRIG + 1):
        term = term * X2 / D((2 * k) * (2 * k + 1))
        s = s + (term if k % 2 == 0 else -term)
    nxt = _fhi(lambda: term.absmax() * x * x / D((2 * N_TRIG + 2) * (2 * N_TRIG + 3)))
    ratio = _fhi(lambda: x * x / D((2 * N_TRIG + 4) * (2 * N_TRIG + 5)))
    return s + sym(_geom_tail(nxt, ratio))


def sin_iv(v):
    return IV(sin_pt(v.lo).lo, sin_pt(v.hi).hi)          # increasing on [0,0.8]


def cos_iv(v):
    return IV(cos_pt(v.hi).lo, cos_pt(v.lo).hi)          # decreasing on [0,0.8]


def sinc_iv(v):
    return IV(sinc_pt(v.hi).lo, sinc_pt(v.lo).hi)        # decreasing on [0,0.8]


# ----------------------------------------------------------------------
# coshC(x) = (cosh x - 1)/x^2 on [0, 1.6]
# ----------------------------------------------------------------------
N_COSHC = 12


def coshC_pt(x):
    _chk('coshC', x)
    X2 = IV(x) * IV(x)
    term = IV(1) / 2
    s = term
    for k in range(N_COSHC):
        term = term * X2 / D((2 * k + 3) * (2 * k + 4))
        s = s + term
    m = 2 * N_COSHC + 3
    nxt = _fhi(lambda: term.hi * x * x / D(m * (m + 1)))
    ratio = _fhi(lambda: x * x / D((m + 2) * (m + 3)))
    tail = _geom_tail(nxt, ratio)
    return s + IV(D(0), tail)


def coshC_iv(v):
    return IV(coshC_pt(v.lo).lo, coshC_pt(v.hi).hi)      # increasing


# ----------------------------------------------------------------------
# atanh series and L1(z) = log(1+z)/z on [-0.8, 2.5]
# ----------------------------------------------------------------------
N_ATANH = 90


def atanh_pt(t):
    _chk('atanh', t)
    T = IV(t)
    T2 = T * T
    term = T
    s = T
    for n in range(1, N_ATANH + 1):
        term = term * T2
        s = s + term / D(2 * n + 1)
    nxt = _fhi(lambda: term.absmax() * t * t / D(2 * N_ATANH + 3))
    ratio = _fhi(lambda: t * t)
    tail = _geom_tail(nxt, ratio)
    return s + sym(tail)


def L1_pt(z):
    """log(1+z)/z via log(1+z) = 2 atanh(z/(2+z))."""
    _chk('L1', z)
    if z == 0:
        return IV(1)
    Z = IV(z)
    t = Z / (IV(2) + Z)
    if t.lo != t.hi:
        # t is an enclosure, not a point; bound atanh over it by monotonicity
        v = IV(atanh_pt(t.lo).lo, atanh_pt(t.hi).hi)
    else:
        v = atanh_pt(t.lo)
    return IV(2) * v / Z


def L1_iv(v):
    return IV(L1_pt(v.hi).lo, L1_pt(v.lo).hi)            # decreasing on z > -1


# ----------------------------------------------------------------------
# Efun(x) = (1 - e^{-x})/x on [-1, 3].  Valid for both signs: the
# remainder is bounded by magnitude, not by an alternation argument.
# ----------------------------------------------------------------------
N_EFUN = 42


def Efun_pt(x):
    _chk('Efun', x)
    X = IV(x)
    term = IV(1)          # (-x)^n/(n+1)! at n = 0
    s = IV(1)
    for n in range(1, N_EFUN + 1):
        term = term * (-X) / D(n + 1)
        s = s + term
    u = abs(x)
    nxt = _fhi(lambda: term.absmax() * u / D(N_EFUN + 2))
    ratio = _fhi(lambda: u / D(N_EFUN + 3))
    tail = _geom_tail(nxt, ratio)
    return s + sym(tail)


def Efun_iv(v):
    return IV(Efun_pt(v.hi).lo, Efun_pt(v.lo).hi)        # decreasing on all of R


# ----------------------------------------------------------------------
# arctan with four half-angle reductions, then a fast series.
#   atan(x) = 2 atan( x / (1 + sqrt(1 + x^2)) )
# Four reductions send any x >= 0 into [0, tan(pi/32)] < 0.0985.
# ----------------------------------------------------------------------
N_ATAN = 22
ATAN_REDUCTIONS = 4


def _atan_small(V):
    """arctan of an enclosure V with |V| <= 0.1."""
    if V.absmax() > D('0.1'):
        raise DomainError(('atan series argument too large', V))
    V2 = V * V
    term = V
    s = V
    for n in range(1, N_ATAN + 1):
        term = term * V2
        s = s + (term / D(2 * n + 1) if n % 2 == 0 else -(term / D(2 * n + 1)))
    b = V.absmax()
    nxt = _fhi(lambda: term.absmax() * b * b / D(2 * N_ATAN + 3))
    ratio = _fhi(lambda: b * b)
    return s + sym(_geom_tail(nxt, ratio))


def atan_pt(x):
    if x < 0:
        return -atan_pt(-x)
    V = IV(x)
    for _ in range(ATAN_REDUCTIONS):
        V = V / (IV(1) + (IV(1) + V * V).sqrt())
    r = _atan_small(V)
    for _ in range(ATAN_REDUCTIONS):
        r = IV(2) * r
    return r


def atan_iv(v):
    return IV(atan_pt(v.lo).lo, atan_pt(v.hi).hi)        # increasing


# ----------------------------------------------------------------------
# pi by Machin: pi = 16 atan(1/5) - 4 atan(1/239)
# ----------------------------------------------------------------------
def _machin_pi():
    a = atan_iv(IV(1) / IV(5))
    b = atan_iv(IV(1) / IV(239))
    return IV(16) * a - IV(4) * b


PI = _machin_pi()
