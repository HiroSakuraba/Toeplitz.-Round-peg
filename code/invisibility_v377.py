"""v377: near-squares with a free corner at the wild point.

Near a two-arm spiral point p, choose a side s and an angle theta. The free corner A is placed in the
comb cell nearest p so that corners B = A + s e^{i theta} and D = A + i s e^{i theta} lie exactly on the arc;
then theta is moved until C = B + D - A lies on the arc too. The result is a square of side s with three
corners on the arc and the fourth within about one turn spacing d(s) = a g s^(1+1/a) of p, which bounds the
invisibility radius of Lemma 374.2 from above.
"""
import numpy as np
from math import pi


def arc_function(z, lam, a=0.5):
    """Smooth function vanishing exactly on the two arms (for radii inside the arc)."""
    r = abs(z)
    psi = np.angle(z)
    Th = r ** (-1 / a)
    if lam == 0:
        return np.sin(psi - Th)
    ph = lam * 0.4 * np.sin(5 * np.log(1 / r))
    return np.sin((psi - Th) / 2) * np.sin((psi - Th - pi - ph) / 2)


def turn_spacing(s, lam, a=0.5):
    return a * (pi - 0.4 * lam) * s ** (1 + 1 / a)


def place_free_corner(A, th, s, lam, d0, a=0.5):
    e = np.exp(1j * th)
    for _ in range(50):
        B, D = A + s * e, A + 1j * s * e
        y = np.array([arc_function(B, lam, a), arc_function(D, lam, a)])
        if np.max(np.abs(y)) < 1e-14:
            return A
        h = d0 * 1e-6
        J = np.zeros((2, 2))
        for k, dz in enumerate((h, 1j * h)):
            J[:, k] = (np.array([arc_function(B + dz, lam, a), arc_function(D + dz, lam, a)]) - y) / h
        try:
            st = np.linalg.solve(J, -y)
        except np.linalg.LinAlgError:
            return None
        st = np.clip(st, -0.2 * d0, 0.2 * d0)
        A = A + st[0] + 1j * st[1]
    e = np.exp(1j * th)
    ok = max(abs(arc_function(A + s * e, lam, a)), abs(arc_function(A + 1j * s * e, lam, a))) < 1e-11
    return A if ok else None


def near_squares(s, lam, n_theta=20000, a=0.5, start=0j):
    """Squares of side s with B, C, D on the arc and free corner A in the comb cell nearest p + start.

    Returns a list of (A, theta, max residual of the arc function at B, C, D)."""
    d0 = turn_spacing(s, lam, a)
    gC = lambda A, th: arc_function(A + (1 + 1j) * s * np.exp(1j * th), lam, a)
    out, prev = [], None
    for th in np.linspace(0, 2 * pi, n_theta, endpoint=False):
        A = place_free_corner(start, th, s, lam, d0, a)
        if A is None:
            prev = None
            continue
        v = gC(A, th)
        if prev is not None and np.sign(v) != np.sign(prev[1]) and abs(A - prev[2]) < 0.5 * d0:
            lo, hi, Alo, Am, mid = prev[0], th, prev[2], None, None
            for _ in range(60):
                mid = (lo + hi) / 2
                Am = place_free_corner(Alo, mid, s, lam, d0, a)
                if Am is None:
                    break
                if np.sign(gC(Am, mid)) == np.sign(prev[1]):
                    lo, Alo = mid, Am
                else:
                    hi = mid
            if Am is not None:
                e = np.exp(1j * mid)
                res = max(abs(arc_function(Am + s * e, lam, a)), abs(gC(Am, mid)),
                          abs(arc_function(Am + 1j * s * e, lam, a)))
                if res < 1e-9:
                    out.append((Am, mid, res))
        prev = (th, v, A)
    return out
