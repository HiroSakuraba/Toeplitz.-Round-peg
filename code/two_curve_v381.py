"""v381: squares standing between two curves. Existence at stationary gaps, and localization.

Two graphs y = f1(x) (lower) and y = f2(x) (upper), gap g = f2 - f1 > 0. A square with bottom corners
P1 = (x1, f1(x1)), P2 = (x1 + u, f1(x1 + u)) and top corners P4 = P1 + J(P2 - P1), P3 = P2 + J(P2 - P1)
(J = quarter turn) is inscribed iff
    F1:  f1(x1) + u - f2(x1 - s) = 0,          s = f1(x1 + u) - f1(x1)          (P4 on the upper curve)
    H:   s - [f2(x1 + u - s) - f2(x1 - s)] = 0                                      (then P3 is too).
F1 has a unique root u(x1) > 0 when |f1'|, |f2'| <= sigma < 1, so the squares are the zeros of
    H(x1) = -[g(x1 + u) - g(x1)] + R,   |R| <= |s| * TV(f2') over the window [x1, x1 + u] widened by |s| on the side s points to  <=  sigma u V2.
"""
import numpy as np
from math import pi


class Graph:
    """A graph y = f(x) from sample points (piecewise linear: exact for polygons, dense for smooth curves)."""

    def __init__(self, X, Y):
        o = np.argsort(X)
        self.X, self.Y = np.asarray(X)[o], np.asarray(Y)[o]

    def __call__(self, x):
        return np.interp(x, self.X, self.Y)

    def slopes(self, lo, hi):
        """Slopes of every segment meeting [lo, hi] (partial end segments included, so TV is not underestimated)."""
        i0 = max(int(np.searchsorted(self.X, lo, side='right')) - 1, 0)
        i1 = min(int(np.searchsorted(self.X, hi, side='left')), len(self.X) - 1)
        X, Y = self.X[i0:i1 + 1], self.Y[i0:i1 + 1]
        return np.diff(Y) / np.diff(X)


def solve_u(f1, f2, x1, umax):
    lo, hi = 0.0, umax
    F = lambda u: f1(x1) + u - f2(x1 - (f1(x1 + u) - f1(x1)))
    if F(hi) < 0:
        return None
    for _ in range(80):
        mid = (lo + hi) / 2
        if F(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _inside(f, *xs):
    return all(f.X[0] <= x <= f.X[-1] for x in xs)


def H(f1, f2, x1, umax):
    """H at bottom-left position x1, or (None, None) if some corner would leave the sampled stretch of a curve."""
    if not (_inside(f1, x1, x1 + umax) and _inside(f2, x1 - umax, x1 + umax)):
        return None, None
    u = solve_u(f1, f2, x1, umax)
    if u is None:
        return None, None
    s = f1(x1 + u) - f1(x1)
    if not (_inside(f1, x1 + u) and _inside(f2, x1 - s, x1 + u - s)):
        return None, None
    return s - (f2(x1 + u - s) - f2(x1 - s)), u


def square(f1, f2, x1, u):
    P1 = np.array([x1, f1(x1)])
    P2 = np.array([x1 + u, f1(x1 + u)])
    d = P2 - P1
    Jd = np.array([-d[1], d[0]])
    return [P1, P2, P2 + Jd, P1 + Jd]


def zeros(f1, f2, xs, umax):
    """Zeros of H on a grid of bottom-left positions, refined by bisection. Returns (x1, u)."""
    vals = [H(f1, f2, x, umax) for x in xs]
    out = []
    for (x0, (h0, _)), (x1, (h1, _)) in zip(zip(xs, vals), zip(xs[1:], vals[1:])):
        if h0 is None or h1 is None or h0 * h1 > 0:
            continue
        lo, hi, hlo = x0, x1, h0
        for _ in range(80):
            mid = (lo + hi) / 2
            hm, _ = H(f1, f2, mid, umax)
            if hm * hlo > 0:
                lo, hlo = mid, hm
            else:
                hi = mid
        x = (lo + hi) / 2
        out.append((x, solve_u(f1, f2, x, umax)))
    return out


def hypothesis_margin(f1, f2, x, umax):
    """(gap change across the chord, error bound sigma u V2) at bottom-left position x."""
    u = solve_u(f1, f2, x, umax)
    s = f1(x + u) - f1(x)
    dg = (f2(x + u) - f1(x + u)) - (f2(x) - f1(x))
    lo, hi = x - max(s, 0.0) - 1e-12, x + u + max(-s, 0.0) + 1e-12   # every x - tau, tau between 0 and s
    sl2 = f2.slopes(lo - 1e-9, hi + 1e-9)
    V2 = float(np.sum(np.abs(np.diff(sl2)))) if len(sl2) > 1 else 0.0
    return dg, abs(s) * V2


# ---------------------------------------------------------------- spiral passes in a local frame
def pass_frame(rfun, psi_c, L_outer, half=0.6, n=20001, L_inner=None):
    """Lower (inner) and upper (outer) passes near angle psi_c as graphs in the frame with origin on the inner
    pass, x tangential and y radial. Pass radii: outer r(psi + L_outer), inner r(psi + L_inner), by default the
    neighbouring pass L_inner = L_outer + pi."""
    if L_inner is None:
        L_inner = L_outer + pi
    psi = np.linspace(psi_c - half, psi_c + half, n)
    ex, ey = 1j * np.exp(1j * psi_c), np.exp(1j * psi_c)
    O = rfun(psi_c + L_inner) * np.exp(1j * psi_c)
    graphs = []
    for L in (L_inner, L_outer):
        z = rfun(psi + L) * np.exp(1j * psi) - O
        X = (z * np.conj(ex)).real
        Y = (z * np.conj(ey)).real
        graphs.append(Graph(X, Y))
    to_plane = lambda P: O + P[0] * ex + P[1] * ey
    return graphs[0], graphs[1], to_plane


def polygon_pass_frame(inner_pts, outer_pts, psi_c):
    """Polygonal passes (vertex lists in the plane) near angle psi_c as exact piecewise-linear graphs in the frame
    x = <z, i e^{i psi_c}>, y = <z, e^{i psi_c}> (origin at the wild point)."""
    ex, ey = 1j * np.exp(1j * psi_c), np.exp(1j * psi_c)
    graphs = []
    for pts in (inner_pts, outer_pts):
        z = np.array(pts)
        graphs.append(Graph((z * np.conj(ex)).real, (z * np.conj(ey)).real))
    to_plane = lambda P: P[0] * ex + P[1] * ey
    return graphs[0], graphs[1], to_plane
