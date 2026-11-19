"""v230 audit: machine checks behind Conditional Theorem 2 and the spiral identity.

Three independent checks of statements in research/v230-irregular-frontier-double-spiral.md.

1. certify_box_a()  [needs python-flint]
   Arb ball arithmetic over the WHOLE v100 Box A, i.e. c in [1,19/10],
   d in [9/20,18/25], q in [0,1/3] (q=1/(1+p), so every p>=2 at once), of the
   three inequalities that Theorem 2 uses:

     I1  (1-exp(-c(1-q))) / (3*pi - d*q)  >  1/25
           => r0^p > 1/25  =>  r0 > 25^(-1/p) >= 1/5
     I2  (exp(c(1-q)) - 1) / (3*pi - d*q)  <  1
           => outer vertex t3^p = rho^p r0^p < 1   (the v100 section 9 condition)
     I3  1 + rho^2 + 2*rho*cos(delta)  >  25/8,   rho = exp(c*q), delta = d*q
           => side = r0*|1+rho e^(i delta)|/sqrt(2) > (1/5)*sqrt(25/16) = 1/4

   These hold on the box, not just at the (unknown) Miranda zero, so they need
   no information about where the zero is.  They are still CONDITIONAL on the
   Box-A zero existing (the eight v101 face certificates) and on the v100
   analytic reduction.

2. y_dx_identity_residual(t, p)  [standard library]
   Complex-step derivatives check
       y * dx/dt  ==  d/dt[(t^2/4) sin(2 t^-p)] + (p/2) t^(1-p)
   for gamma_p(t) = t exp(i t^-p), without reusing the hand algebra.

3. track_box_a_square(p, ...)  [needs mpmath]  -- NUMERICAL, NOT A PROOF
   Newton continuation of the actual inscribed square with one diagonal on
   gamma_p(r0), gamma_p(r2).  Reports residuals, vertex radii, side, whether
   all four vertices are on the first arm, the recovered Box-A coordinates
   (c, d), and the angular coordinate theta = t^-p of the outermost vertex.
   That last number decides whether the square survives the theta >= e
   truncation used by Asano-Ike (arXiv:2609.11598v3, section 6.1).
"""
import cmath
import math
import os
import sys
from fractions import Fraction as F

BOX_A = (F(1), F(19, 10), F(9, 20), F(18, 25), F(0), F(1, 3))
I3_THRESHOLD = F(25, 8)


# ---------------------------------------------------------------------------
# 1. Arb certification over Box A
# ---------------------------------------------------------------------------
def _require_flint():
    try:
        import flint  # noqa: F401
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError('python-flint is required for certify_box_a') from exc
    return flint


def _box_a_margins(cell, flint):
    arb = flint.arb
    def ball(a, b):
        lo = arb(a.numerator) / arb(a.denominator)
        hi = arb(b.numerator) / arb(b.denominator)
        return lo.union(hi)
    (cl, cu), (dl, du), (ql, qu) = cell
    c, d, q = ball(cl, cu), ball(dl, du), ball(ql, qu)
    pi = arb.pi()
    den = 3 * pi - d * q
    a = c * (1 - q)
    i1 = (1 - (-a).exp()) / den - arb(1) / 25
    i2 = 1 - (a.exp() - 1) / den
    rho = (c * q).exp()
    i3 = 1 + rho * rho + 2 * rho * (d * q).cos() - arb(25) / 8
    return {'I1': i1, 'I2': i2, 'I3': i3}


def certify_box_a(prec=200, max_depth=30):
    """Certify I1, I2, I3 > 0 on every point of Box A by adaptive bisection.

    Returns {name: (cells_used, certified_lower_bound_of_margin_as_float)}.
    Raises AssertionError if any cell stays unresolved at max_depth.
    """
    flint = _require_flint()
    flint.ctx.prec = prec
    cl, cu, dl, du, ql, qu = BOX_A
    out = {}
    for name in ('I1', 'I2', 'I3'):
        stack = [(((cl, cu), (dl, du), (ql, qu)), 0)]
        cells = 0
        worst = None
        while stack:
            cell, depth = stack.pop()
            m = _box_a_margins(cell, flint)[name]
            if m > 0:
                cells += 1
                lb = m.lower()
                worst = lb if worst is None or lb < worst else worst
                continue
            if depth >= max_depth:
                raise AssertionError((name, 'unresolved cell', cell))
            widths = [(iv[1] - iv[0]) / (full[1] - full[0])
                      for iv, full in zip(cell, ((cl, cu), (dl, du), (ql, qu)))]
            k = max(range(3), key=lambda i: widths[i])
            lo, hi = cell[k]
            mid = (lo + hi) / 2
            for half in ((lo, mid), (mid, hi)):
                new = list(cell)
                new[k] = half
                stack.append((tuple(new), depth + 1))
        out[name] = (cells, float(worst.mid()))
        assert worst > 0
    return out


# ---------------------------------------------------------------------------
# 2. Complex-step check of the exact y dx identity
# ---------------------------------------------------------------------------
def _cstep(f, t, h=1e-30):
    return f(complex(t, h)).imag / h


def y_dx_identity_residual(t, p):
    """Relative residual of y x' - [F'(t) + (p/2) t^(1-p)], F=(t^2/4) sin(2 t^-p).

    x' and F' come from complex-step differentiation (no hand-derived
    derivative is reused), so this independently tests the identity.
    """
    x = lambda s: s * cmath.cos(s ** (-p))
    Fn = lambda s: s * s / 4 * cmath.sin(2 * s ** (-p))
    y = t * math.sin(t ** (-p))
    lhs = y * _cstep(x, t)
    rhs = _cstep(Fn, t) + p / 2 * t ** (1 - p)
    scale = max(1.0, abs(lhs), abs(rhs), p * t ** (1 - p))
    return abs(lhs - rhs) / scale


def liouville_primitive_tail(p, eps):
    """Exact value of the non-oscillating part  int_eps^1 (p/2) t^(1-p) dt."""
    if p == 2:
        return -math.log(eps)
    return p / (2 * (2 - p)) * (1 - eps ** (2 - p))


# ---------------------------------------------------------------------------
# 3. Numerical tracking of the Box-A square (NON-PROOF)
# ---------------------------------------------------------------------------
def _mp():
    try:
        import mpmath
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError('mpmath is required for square tracking') from exc
    return mpmath


def _square_from_diagonal(r0, r2, p, mp):
    g = lambda r: r * mp.exp(1j * r ** (-p))
    z0, z2 = g(r0), g(r2)
    m = (z0 + z2) / 2
    return (z0, m + 1j * (z0 - m), z2, m - 1j * (z0 - m))


def _arm_residual(w, p, mp):
    r = abs(w)
    return mp.im(mp.exp(-1j * r ** (-p)) * w) / r      # 0 on either arm


def _on_first_arm(w, p, mp):
    r = abs(w)
    return mp.re(mp.exp(-1j * r ** (-p)) * w) / r > 0


def track_box_a_square(p, c_guess=1.29, d_guess=0.587, dps=60):
    """Locate the Box-A square for exponent p by Newton's method (numerical)."""
    mp = _mp()
    mp.mp.dps = dps
    p = mp.mpf(p)
    q = 1 / (1 + p)
    c, d = mp.mpf(c_guess), mp.mpf(d_guess)
    rho = mp.e ** (c * q)
    delta_phase = 3 * mp.pi - d * q
    r0 = ((1 - rho ** (-p)) / delta_phase) ** (1 / p)
    f = lambda a, b: [_arm_residual(w, p, mp)
                      for w in _square_from_diagonal(a, b, p, mp)[1::2]]
    a, b = mp.findroot(f, (r0, rho * r0), solver='mdnewton', verify=False,
                       maxsteps=80)
    Z = _square_from_diagonal(a, b, p, mp)
    radii = [abs(z) for z in Z]
    sides = [abs(Z[(i + 1) % 4] - Z[i]) for i in range(4)]
    c_found = mp.log(b / a) / q
    d_found = (3 * mp.pi - (a ** (-p) - b ** (-p))) / q
    return {
        'p': p,
        'residual': max(abs(v) for v in f(a, b)),
        'side_spread': max(sides) - min(sides),
        'min_radius': min(radii),
        'max_radius': max(radii),
        'side': sides[0],
        'all_on_first_arm': all(_on_first_arm(z, p, mp) for z in Z),
        'c': c_found,
        'd': d_found,
        'in_box_a': (1 <= c_found <= mp.mpf(19) / 10
                     and mp.mpf(9) / 20 <= d_found <= mp.mpf(18) / 25),
        'outer_theta': max(radii) ** (-p),
    }


def track_family(ps, dps=60):
    """Continuation in p: each solve starts from the previous (c, d)."""
    c, d, out = 1.29, 0.587, []
    for p in ps:
        r = track_box_a_square(p, c, d, dps)
        out.append(r)
        c, d = r['c'], r['d']
    return out


def asano_ike_cut_crossover(lo=20, hi=30, steps=40, dps=60):
    """Bisect for the p where the outer vertex reaches theta = e (numerical)."""
    mp = _mp()
    start = track_family([2, 5, 10, lo], dps)[-1]
    c, d = start['c'], start['d']
    f = lambda p: track_box_a_square(p, c, d, dps)['outer_theta'] - mp.e
    a, b = mp.mpf(lo), mp.mpf(hi)
    fa = f(a)
    if fa * f(b) > 0:
        raise ValueError('no sign change in bracket')
    for _ in range(steps):
        m = (a + b) / 2
        fm = f(m)
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
    return (a + b) / 2


def main():
    print('I. Arb certification of the Theorem 2 inequalities on all of Box A')
    for name, (cells, margin) in certify_box_a().items():
        print(f'   {name}: certified on {cells} cells, margin >= {margin:.6g}')
    print('II. y dx identity, max complex-step residual over a (t,p) grid:',
          '%.2e' % max(y_dx_identity_residual(t, p)
                       for p in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0)
                       for t in (0.05, 0.1, 0.3, 0.5, 0.9)))
    try:
        rows = track_family([2, 2.5, 3, 5, 10, 30, 100])
    except RuntimeError as exc:
        print('III. skipped:', exc)
        return 0
    print('III. Box-A square, numerical tracking (not a proof):')
    for r in rows:
        print('   p=%-6g min r=%.5f max r=%.5f side=%.5f arm1=%s inBoxA=%s outer theta=%.4f'
              % (float(r['p']), float(r['min_radius']), float(r['max_radius']),
                 float(r['side']), r['all_on_first_arm'], r['in_box_a'],
                 float(r['outer_theta'])))
    print('   outer vertex crosses the Asano-Ike theta>=e cut near p = %.3f'
          % float(asano_ike_cut_crossover()))
    return 0


if __name__ == '__main__':
    sys.exit(main())
