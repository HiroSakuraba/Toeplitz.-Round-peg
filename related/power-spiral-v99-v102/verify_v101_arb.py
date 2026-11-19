"""
Independent Arb reproduction of the eight v100/v101 face certificates.

This shares no code with the Decimal layer.  Differences on purpose:

  - cell endpoints are exact Fractions, subdivided exactly;
  - all enclosures are Arb balls at a fixed working precision;
  - every transcendental comes from Arb's own implementation
    (exp, sin, cos, sinh, log1p, expm1, atan, pi), not from a series
    written for this project;
  - removable singularities are handled by closed forms that do not
    cancel:  coshC(x) = (sinh(x/2)/(x/2))^2 / 2.

The purpose is to check the eight face signs, not to be fast.
"""

import sys
import time
from fractions import Fraction as F
from flint import arb, ctx

PREC = 160
ctx.prec = PREC


def fr(x):
    """Exact Fraction -> Arb ball containing it."""
    return arb(x.numerator) / arb(x.denominator)


def cell(a, b):
    """Arb ball containing the closed rational interval [a, b]."""
    return fr(a).union(fr(b))


def lo(X):
    return X.lower()


def hi(X):
    return X.upper()


def sq(X):
    """Square of a ball, via endpoints, avoiding ball-arithmetic blowup."""
    l, h = X.lower(), X.upper()
    if l >= 0:
        return (l * l).union(h * h)
    if h <= 0:
        return (h * h).union(l * l)
    return arb(0).union((l * l).union(h * h))


class Unevaluable(Exception):
    """Enclosure too coarse to evaluate on this cell; caller should subdivide."""


def clamp0(X):
    """Enclosure of X intersected with [0, inf) when X is known nonnegative."""
    l = X.lower()
    return l.union(X.upper()) if l >= 0 else arb(0).union(X.upper())


# ---- point evaluations, removable singularities handled explicitly ----

def sinc_pt(t):
    return arb(1) if t.is_zero() else t.sin() / t


def coshC_pt(x):
    if x.is_zero():
        return arb(1) / 2
    u = x / 2
    return (u.sinh() / u) ** 2 / 2


def L1_pt(z):
    return arb(1) if z.is_zero() else z.log1p() / z


def Efun_pt(x):
    return arb(1) if x.is_zero() else -((-x).expm1()) / x


# ---- monotone interval lifts ----

def mono(f, X, increasing):
    a = f(X.lower())
    b = f(X.upper())
    return a.union(b) if increasing else b.union(a)


def sinc_iv(X):
    return mono(sinc_pt, X, False)          # decreasing on [0, pi]


def coshC_iv(X):
    return mono(coshC_pt, X, True)


def L1_iv(X):
    return mono(L1_pt, X, False)            # decreasing on (-1, inf)


def Efun_iv(X):
    return mono(Efun_pt, X, False)          # decreasing on R


def exp_iv(X):
    return mono(lambda t: t.exp(), X, True)


def sin_iv(X):
    return mono(lambda t: t.sin(), X, True)      # valid on [0, pi/2]


def cos_iv(X):
    return mono(lambda t: t.cos(), X, False)     # valid on [0, pi]


def atan_iv(X):
    return mono(lambda t: t.atan(), X, True)


PI = arb.pi()


def evalE(cI, dI, qI):
    x = clamp0(cI * qI)
    y = clamp0(dI * qI)

    S = sinc_iv(y)
    C = coshC_iv(x)
    c2 = sq(cI)
    wp = c2 * qI * C + dI * S
    wm = c2 * qI * C - dI * S
    zp = qI * wp
    zm = qI * wm
    if zm.lower() <= -1 or zp.lower() <= -1:
        raise Unevaluable(('log1p domain', zp, zm))
    Lp = L1_iv(zp)
    Lm = L1_iv(zm)

    sp = (1 + wp * Lp / cI) / 2
    sm = (1 + wm * Lm / cI) / 2

    a = cI * (1 - qI)
    Ea = Efun_iv(a)
    gp = sp * Efun_iv(a * sp) / Ea
    gm = sm * Efun_iv(a * sm) / Ea

    rho = exp_iv(x)
    sy = sin_iv(y)
    cy = cos_iv(y)

    rap = (rho * (cy - sy) - 1) / (1 + rho * (cy + sy))
    ram = (rho * (cy + sy) - 1) / (1 + rho * (cy - sy))
    ap = atan_iv(rap)
    am = atan_iv(ram)

    Delta = 3 * PI - y
    Ep = arb(5) / 2 * PI + ap - Delta * gp
    Em = arb(3) / 2 * PI - am - Delta * gm
    return Ep, Em


def adaptive_face(box, kind, val, target, sign, alpha=F(0), maxdepth=22):
    cL, cU, dL, dU, qL, qU = box
    V = fr(val)
    stack = [((dL, dU) if kind == 'c' else (cL, cU)) + (qL, qU, 0)]
    leaves = 0
    bad = []
    margin = None
    nodes = 0
    maxdep = 0
    t0 = time.time()
    while stack:
        alo, ahi, qlo, qhi, dep = stack.pop()
        nodes += 1
        maxdep = max(maxdep, dep)
        A = cell(alo, ahi)
        cI, dI = (V, A) if kind == 'c' else (A, V)
        try:
            Ep, Em = evalE(cI, dI, cell(qlo, qhi))
            T = Ep if target == 'Ep' else Em if target == 'Em' else Ep + fr(alpha) * Em
            ok = (T > 0) if sign > 0 else (T < 0)
        except (Unevaluable, ZeroDivisionError):
            T = None
            ok = False
        if ok:
            leaves += 1
            mm = T.lower() if sign > 0 else -T.upper()
            if margin is None or mm < margin:
                margin = mm
            continue
        if dep >= maxdepth:
            bad.append((alo, ahi, qlo, qhi, dep, 'unevaluable' if T is None else 'sign'))
            continue
        aw = ahi - alo
        qw = qhi - qlo
        fulla = (dU - dL) if kind == 'c' else (cU - cL)
        fullq = qU - qL
        if aw / fulla >= qw / fullq:
            mid = (alo + ahi) / 2
            stack.append((alo, mid, qlo, qhi, dep + 1))
            stack.append((mid, ahi, qlo, qhi, dep + 1))
        else:
            mid = (qlo + qhi) / 2
            stack.append((alo, ahi, qlo, mid, dep + 1))
            stack.append((alo, ahi, mid, qhi, dep + 1))
    return {'nodes': nodes, 'leaves': leaves, 'bad': len(bad), 'margin': margin,
            'maxdepth': maxdep, 'secs': time.time() - t0, 'bad_examples': bad[:3]}


# Exact rational boxes.  Unlike v100 these are the true rational endpoints,
# with no decimal over-cover: q in [0, 1/3] and q in [1/3, 200/201].
BOX_A = (F(1), F(19, 10), F(9, 20), F(18, 25), F(0), F(1, 3))
BOX_B = (F(3, 4), F(3, 2), F(2, 5), F(18, 25), F(1, 3), F(200, 201))
ALPHA = F(57, 400)

TESTS = [
    (BOX_A, 'A_c_low',  'c', F(1),      'Em', +1, F(0)),
    (BOX_A, 'A_c_high', 'c', F(19, 10), 'Em', -1, F(0)),
    (BOX_A, 'A_d_low',  'd', F(9, 20),  'Ep', +1, F(0)),
    (BOX_A, 'A_d_high', 'd', F(18, 25), 'Ep', -1, F(0)),
    (BOX_B, 'B_c_low',  'c', F(3, 4),   'Em', +1, F(0)),
    (BOX_B, 'B_c_high', 'c', F(3, 2),   'Em', -1, F(0)),
    (BOX_B, 'B_d_low',  'd', F(2, 5),   'G',  +1, ALPHA),
    (BOX_B, 'B_d_high', 'd', F(18, 25), 'G',  -1, ALPHA),
]

if __name__ == '__main__':
    only = int(sys.argv[1]) if len(sys.argv) > 1 else None
    print(f'Arb reproduction at prec={PREC}, exact rational cells')
    tot_bad = 0
    for i, (box, name, kind, val, target, sign, alpha) in enumerate(TESTS):
        if only is not None and i != only:
            continue
        r = adaptive_face(box, kind, val, target, sign, alpha)
        tot_bad += r['bad']
        print(f'  {name:9s} nodes={r["nodes"]:6d} leaves={r["leaves"]:6d} '
              f'maxdepth={r["maxdepth"]:3d} bad={r["bad"]:3d} '
              f'margin={float(r["margin"]) if r["margin"] is not None else float("nan"):.6e} '
              f'secs={r["secs"]:.1f}')
    print()
    print('Arb result:', 'all requested faces certified' if tot_bad == 0
          else f'FAILED with {tot_bad} unresolved cells')
