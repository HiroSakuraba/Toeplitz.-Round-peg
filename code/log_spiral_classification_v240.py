"""v240: classify inscribed squares on open logarithmic spirals for a >= 3/50.

S_a = {exp((-a+i)t): t>=0} u {0}.  Proof outline (full text in
research/v240-log-spiral-classification.md):

  N  Normalise: rotate-dilate so the farthest vertex is v0 = 1 = z(0).  The
     other vertices then lie on S_a itself.  Label the square v0 v1 v2 v3 with
     v1 adjacent to v0 and orientation sg in {+1,-1}:
         v3 = 1 + sg*i*(v1 - 1),   v2 = v1 + sg*i*(v1 - 1).
  R  |v1|, |v3| >= sqrt(2)-1, so t1 <= ln(1+sqrt2)/a < 0.8814/a.
  L  Tiny squares (t1 <= eps) are impossible: all four vertices would lie on
     z([0, 6/5]), a strict 1-Lipschitz graph after rotation.
  B  Arb branch-and-bound over (t1, a) boxes, both orientations: a box is
     excluded when v2 or v3 is certainly outside the unit disk or certainly
     off the spiral.
  K  Around the one place B cannot exclude, a Krawczyk test proves a unique
     zero, i.e. exactly one square.

Result: for real a in [3/50, 9/10], S_a has an inscribed square if and only if
a = a* (about 0.1002746561), and then exactly one up to the spiral's
self-similarity.  With v239 (a >= 9/10) this classifies all a >= 3/50.

Requires python-flint (Arb ball arithmetic).  Every decision is a certain
Arb comparison; floats only choose split points and the Krawczyk matrix Y.
"""
import json
import sys
import time
from fractions import Fraction as F

import flint
from flint import arb

PREC = 128
A_MIN, A_MAX, BAND = F(3, 50), F(9, 10), F(1, 200)
TAU = F(6, 5)                 # tangent span for lemma L; 6/5 < pi/2
T1C = F(8814, 10000)          # certified upper bound for ln(1 + sqrt 2)
SQRT2_UP = F(1415, 1000)      # sqrt 2 < 1415/1000
DEPTH_MAX = 96

# The square at a*: high-precision seeds (t1, a) for each orientation,
# from Newton's method; the Krawczyk test below is what certifies them.
RESONANCE = {
    +1: ("1.060706413775613939922", "0.1002746561179786808441"),
    -1: ("5.078914619119260416439", "0.1002746561179786808441"),
}
K_RADIUS_T = F(1, 10 ** 9)
K_RADIUS_A = F(1, 10 ** 10)


def ctx():
    flint.ctx.prec = PREC


def ball(lo, hi=None):
    lo = F(lo)
    hi = lo if hi is None else F(hi)
    x = arb(lo.numerator) / arb(lo.denominator)
    if hi == lo:
        return x
    return x.union(arb(hi.numerator) / arb(hi.denominator))


# --------------------------------------------------------------- vertices --
def vertices(t, a, sg):
    """(x1,y1),(x2,y2),(x3,y3) for v1=z(t) and orientation sg; works on arb or Dual."""
    r1 = (-(a * t)).exp()
    x1, y1 = r1 * t.cos(), r1 * t.sin()
    dx, dy = x1 - 1, y1
    v3 = (1 - sg * dy, sg * dx)
    v2 = (x1 - sg * dy, y1 + sg * dx)
    return (x1, y1), v2, v3


def spiral_status(x, y, a):
    """'out' if |w|>1 for sure, 'off' if w is surely not on the spiral, else 'maybe'.

    w = |w| e^{i phi} lies on the spiral iff u = w*exp(i*ln|w|/a) is a positive
    real.  So Im u != 0, or Re u < 0, rules it out.
    """
    n2 = x * x + y * y
    if n2 > 1:
        return 'out'
    if not n2 > 0:
        return 'maybe'
    th = n2.log() / (2 * a)
    c, s = th.cos(), th.sin()
    im = x * s + y * c
    if im > 0 or im < 0:
        return 'off'
    if x * c - y * s < 0:
        return 'off'
    return 'maybe'


def box_excluded(tl, tu, al, au, sg):
    t, a = ball(tl, tu), ball(al, au)
    _, v2, v3 = vertices(t, a, sg)
    if spiral_status(*v3, a) != 'maybe':
        return True
    return spiral_status(*v2, a) != 'maybe'


# ---------------------------------------------------------------- lemma L --
def margin(al):
    """Rational M <= min(3/5, (4/5)(1 - exp(-2 pi al)))."""
    v = (1 - (-(2 * arb.pi() * ball(al))).exp()) * 4 / 5
    m = F(int(float(v.lower()) * 1000), 1000)
    assert ball(m) < v
    return min(F(3, 5), m)


def lemma_L(al, au):
    """Certify inf_{t >= TAU} |z(t) - 1| > M for all a in [al, au]; return eps."""
    M = margin(al)
    t_big = F(1)
    while not (-(ball(al) * ball(t_big))).exp() < 1 - ball(M):
        t_big *= F(11, 10)
    stack, cells = [(TAU, t_big, al, au, 0)], 0
    while stack:
        tl, tu, bl, bu, d = stack.pop()
        cells += 1
        t, a = ball(tl, tu), ball(bl, bu)
        r = (-(a * t)).exp()
        x, y = r * t.cos() - 1, r * t.sin()
        if x * x + y * y > ball(M) ** 2:
            continue
        if d > 50:
            raise AssertionError(('lemma L fails', float(tl), float(bl)))
        if tu - tl > bu - bl:
            m = (tl + tu) / 2
            stack += [(tl, m, bl, bu, d + 1), (m, tu, bl, bu, d + 1)]
        else:
            m = (bl + bu) / 2
            stack += [(tl, tu, bl, m, d + 1), (tl, tu, m, bu, d + 1)]
    # need sqrt2 * sqrt(1+au^2) * eps < M; use sqrt2 < SQRT2_UP, sqrt(1+x) <= 1+x/2
    eps = M / (SQRT2_UP * (1 + au * au / 2))
    eps = F(int(eps * 10000) - 1, 10000)
    assert SQRT2_UP ** 2 > 2 and 0 < eps < TAU
    assert SQRT2_UP * (1 + au * au / 2) * eps < M
    return eps, M, cells


# ------------------------------------------------------- Krawczyk (dual) --
class Dual:
    """Value plus gradient in (t1, a), all arb balls."""
    __slots__ = ('v', 'g')

    def __init__(self, v, g=(0, 0)):
        self.v, self.g = v, tuple(g)

    @staticmethod
    def _c(o):
        return o if isinstance(o, Dual) else Dual(arb(o) if not isinstance(o, arb) else o)

    def __add__(s, o):
        o = Dual._c(o); return Dual(s.v + o.v, (s.g[0] + o.g[0], s.g[1] + o.g[1]))
    __radd__ = __add__

    def __sub__(s, o):
        o = Dual._c(o); return Dual(s.v - o.v, (s.g[0] - o.g[0], s.g[1] - o.g[1]))

    def __rsub__(s, o):
        return Dual._c(o) - s

    def __neg__(s):
        return Dual(-s.v, (-s.g[0], -s.g[1]))

    def __mul__(s, o):
        o = Dual._c(o)
        return Dual(s.v * o.v, (s.g[0] * o.v + s.v * o.g[0], s.g[1] * o.v + s.v * o.g[1]))
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = Dual._c(o)
        q = s.v / o.v
        return Dual(q, ((s.g[0] - q * o.g[0]) / o.v, (s.g[1] - q * o.g[1]) / o.v))

    def _chain(s, f, df):
        d = df
        return Dual(f, (d * s.g[0], d * s.g[1]))

    def exp(s):
        e = s.v.exp(); return s._chain(e, e)

    def log(s):
        return s._chain(s.v.log(), 1 / s.v)

    def sin(s):
        return s._chain(s.v.sin(), s.v.cos())

    def cos(s):
        return s._chain(s.v.cos(), -s.v.sin())


def residual(t, a, sg):
    """(Im u3, Im u2): both vanish exactly when v3 and v2 lie on the spiral."""
    _, v2, v3 = vertices(t, a, sg)
    out = []
    for x, y in (v3, v2):
        th = (x * x + y * y).log() / (2 * a)
        out.append(x * th.sin() + y * th.cos())
    return out


def krawczyk(sg, tc, ac, rt=K_RADIUS_T, ra=K_RADIUS_A):
    """Krawczyk test on the box [tc+-rt] x [ac+-ra].  True means: exactly one
    zero of the residual there, and at it v2, v3 are genuinely on S_a."""
    tc, ac = F(tc), F(ac)
    T, A = ball(tc - rt, tc + rt), ball(ac - ra, ac + ra)
    # conditions that make a zero a genuine square, checked on the whole box
    _, v2, v3 = vertices(T, A, sg)
    for x, y in (v2, v3):
        n2 = x * x + y * y
        th = n2.log() / (2 * A)
        if not (n2 < 1 and n2 > 0 and x * th.cos() - y * th.sin() > 0):
            return False
    fm = residual(ball(tc), ball(ac), sg)
    J = residual(Dual(T, (1, 0)), Dual(A, (0, 1)), sg)
    Jm = residual(Dual(ball(tc), (1, 0)), Dual(ball(ac), (0, 1)), sg)
    m = [[float(Jm[i].g[j].mid()) for j in range(2)] for i in range(2)]
    det = m[0][0] * m[1][1] - m[0][1] * m[1][0]
    Y = [[F(m[1][1] / det), F(-m[0][1] / det)], [F(-m[1][0] / det), F(m[0][0] / det)]]
    Yb = [[ball(y) for y in row] for row in Y]
    dX = [ball(-rt, rt), ball(-ra, ra)]
    center = [ball(tc), ball(ac)]
    K = []
    for i in range(2):
        k = center[i] - Yb[i][0] * fm[0] - Yb[i][1] * fm[1]
        for j in range(2):
            c = (1 if i == j else 0) - (Yb[i][0] * J[0].g[j] + Yb[i][1] * J[1].g[j])
            k = k + c * dX[j]
        K.append(k)
    lo = [tc - rt, ac - ra]
    hi = [tc + rt, ac + ra]
    return all(K[i] > ball(lo[i]) and K[i] < ball(hi[i]) for i in range(2))


def resonance_boxes():
    boxes = {}
    for sg, (t, a) in RESONANCE.items():
        tc, ac = F(t), F(a)
        if not krawczyk(sg, tc, ac):
            raise AssertionError(('Krawczyk failed', sg))
        boxes[sg] = (tc - K_RADIUS_T, tc + K_RADIUS_T, ac - K_RADIUS_A, ac + K_RADIUS_A)
    return boxes


# -------------------------------------------------------- branch and bound --
def certify_band(al, au, kboxes=None, depth_max=DEPTH_MAX):
    """Exclude every (t1, a) box of the band except those inside a Krawczyk box.

    Returns dict(boxes=..., eps=..., stuck=[...]); stuck is empty on success.
    """
    eps, M, lcells = lemma_L(al, au)
    t1max = T1C / al
    kboxes = kboxes or {}
    boxes, stuck = 0, []
    for sg in (1, -1):
        kb = kboxes.get(sg)
        stack = [(eps, t1max, al, au, 0)]
        while stack:
            tl, tu, bl, bu, d = stack.pop()
            boxes += 1
            if kb and kb[0] <= tl and tu <= kb[1] and kb[2] <= bl and bu <= kb[3]:
                continue
            if box_excluded(tl, tu, bl, bu, sg):
                continue
            if d >= depth_max:
                stuck.append((sg, float(tl), float(bl)))
                continue
            amid = (bl + bu) / 2
            if (tu - tl) >= (bu - bl) * (1 + 1 / amid) * 6:
                m = (tl + tu) / 2
                stack += [(tl, m, bl, bu, d + 1), (m, tu, bl, bu, d + 1)]
            else:
                m = (bl + bu) / 2
                stack += [(tl, tu, bl, m, d + 1), (tl, tu, m, bu, d + 1)]
    return {'boxes': boxes, 'eps': eps, 'M': M, 'lemma_L_cells': lcells, 'stuck': stuck}


def classify(a_min=A_MIN, a_max=A_MAX, band=BAND):
    ctx()
    if not (1 + arb(2).sqrt()).log() < ball(T1C):
        raise AssertionError('ln(1+sqrt2) bound')
    if not ball(TAU) < arb.pi() / 2:
        raise AssertionError('tau bound')
    kboxes = resonance_boxes()
    bands, total, a = [], 0, F(a_min)
    while a < a_max:
        b = min(a + band, a_max)
        r = certify_band(a, b, kboxes)
        if r['stuck']:
            raise AssertionError(('band not certified', float(a), float(b), r['stuck'][:3]))
        total += r['boxes']
        bands.append((a, b, r['eps'], r['boxes']))
        a = b
    return {'bands': bands, 'total_boxes': total, 'kboxes': kboxes}


def report():
    t0 = time.time()
    res = classify()
    return {
        'theorem': ('For real a in [3/50, 9/10], S_a has an inscribed square iff a = a*, '
                    'and then exactly one up to self-similarity'),
        'a_star_enclosure': [str(res['kboxes'][1][2]), str(res['kboxes'][1][3])],
        'a_star_approx': RESONANCE[1][1],
        'square_parameters_at_a_star': {'v1 (sigma=+1)': RESONANCE[1][0],
                                         'v1 (sigma=-1)': RESONANCE[-1][0]},
        'bands': len(res['bands']),
        'total_boxes': res['total_boxes'],
        'min_eps': str(min(b[2] for b in res['bands'])),
        'precision_bits': PREC,
        'seconds': round(time.time() - t0, 1),
        'scope': 'open arc S_a, not a Jordan curve; combine with v239 for a >= 9/10',
    }


if __name__ == '__main__':
    print(json.dumps(report(), indent=2))
    sys.exit(0)
