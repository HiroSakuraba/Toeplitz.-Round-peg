"""v379: inscribed squares on a smooth half-turn-symmetric spiral, and the selection rule.

Arc: both arms of r = Theta^(-a) and its half-turn image, i.e. {z : Phi(z) = |z|^(-1/a) - arg z in pi Z}.
For an inscribed square A, B = A + p, C = A + p + q, D = A + q, write Phi(Z) = pi k_Z and
n1 = k_B - k_A, n3 = k_D - k_A, j = k_A - k_B + k_C - k_D. The exact identity
    pi j = Phi(A) - Phi(B) + Phi(C) - Phi(D) = int_0^1 int_0^1 p^T Hess(Phi)(A + u p + v q) q du dv
together with Hess(Phi) = (1/(K rho^2)) [-e_psi e_psi^T + (1 + 1/a) e_r e_r^T] + (1/rho^2)(e_r e_psi^T + e_psi e_r^T)
gives, to leading order, j = (2 + 1/a) pi K(r) n1 n3: tilted squares sit only at the radii this allows.
"""
import numpy as np
from math import atan2, hypot, pi

A_EXP = 0.5
RMAX = 3.0 ** -0.5


def arc_f(z, a=A_EXP):
    return np.sin(np.angle(z) - abs(z) ** (-1 / a))


def phase_rel(Z, a=A_EXP):
    """Phi at each corner, with arguments lifted continuously from the first corner."""
    base = np.angle(Z[0])
    out = []
    for z in Z:
        arg = base + np.angle(z / Z[0])
        out.append(abs(z) ** (-1 / a) - arg)
    return np.array(out)


def corners(v):
    A = v[0] + 1j * v[1]
    e = np.exp(1j * v[2])
    t = v[3]
    return [A, A + t * e, A + (1 + 1j) * t * e, A + 1j * t * e]   # A, B, C, D counterclockwise


def solve_square(v, it=80, a=A_EXP):
    F = lambda w: np.array([arc_f(z, a) for z in corners(w)])
    for _ in range(it):
        y = F(v)
        if np.max(np.abs(y)) < 1e-13:
            break
        J = np.zeros((4, 4))
        t = abs(v[3])
        h = 1e-7 * t
        for k in range(4):
            dv = np.zeros(4)
            dv[k] = h if k != 2 else h / t
            J[:, k] = (F(v + dv) - y) / dv[k]
        try:
            st = np.linalg.solve(J, -y)
        except np.linalg.LinAlgError:
            return None
        sc = max(1, np.max(np.abs(st[[0, 1, 3]])) / (0.1 * t), abs(st[2]) / 0.1)
        v = v + st / sc
        if v[3] <= 0:
            return None
    return v if np.max(np.abs(F(v))) < 1e-11 else None


def snap(z, a=A_EXP):
    r, u = abs(z), z / abs(z)
    for _ in range(40):
        f = arc_f(r * u, a)
        h = 1e-9 * r
        df = (arc_f((r + h) * u, a) - f) / h
        if df == 0 or abs(f) < 1e-14:
            break
        r += max(min(-f / df, 0.15 * a * r ** 3), -0.15 * a * r ** 3)
    return r * u


def classify(v, a=A_EXP):
    """(n1, n3, j, side, centre radius, min corner radius) of an inscribed square."""
    Z = corners(v)
    k = np.rint(phase_rel(Z, a) / pi).astype(int)
    n1, n3 = k[1] - k[0], k[3] - k[0]
    j = k[0] - k[1] + k[2] - k[3]
    return int(n1), int(n3), int(j), float(v[3]), float(abs(sum(Z) / 4)), float(min(abs(z) for z in Z))


def predicted_radius(n1, n3, j, a=A_EXP):
    """Radius where (2 + 1/a) pi K(r) n1 n3 = j, K = a r^(1/a)."""
    K = j / ((2 + 1 / a) * pi * n1 * n3)
    return (K / a) ** a if K > 0 else None


def min_side(r, a=A_EXP):
    """Smallest side a tilted square can have at radius r: r sqrt(2 pi K / (2 + 1/a))."""
    return r * (2 * pi * a * r ** (1 / a) / (2 + 1 / a)) ** 0.5


def outer_radius(t, a=A_EXP):
    """Inverse of min_side: radius beyond which no tilted square of side <= t sits."""
    return (t * t * (2 + 1 / a) / (2 * pi * a)) ** (a / (2 * a + 1))


def search(seeds_rng=0, n_random=20000, a=A_EXP):
    """Collect inscribed squares from random seeds and from seeds at the predicted radii."""
    rng = np.random.default_rng(seeds_rng)
    found = {}

    def add(v):
        if v is None:
            return
        Z = corners(v)
        if v[3] > 0.08 or max(abs(z) for z in Z) > RMAX:
            return
        key = tuple(np.round(sorted((z.real, z.imag) for z in Z), 7).ravel())
        found.setdefault(key, v)

    for _ in range(n_random):
        r = rng.uniform(0.09, 0.36)
        D = pi * a * r ** (1 + 1 / a)
        A = snap(r * np.exp(1j * rng.uniform(0, 2 * pi)), a)
        add(solve_square(np.array([A.real, A.imag, rng.uniform(0, 2 * pi), D * rng.uniform(0.8, 6)]), a=a))
    for n1 in range(1, 8):
        for n3 in range(1, 8):
            for j in (1, 2, 3):
                r = predicted_radius(n1, n3, j, a)
                if r is None or not 0.07 < r < 0.42:
                    continue
                D = pi * a * r ** (1 + 1 / a)
                t = D * hypot(n1, n3)
                for sgn in (1, -1):
                    beta = sgn * atan2(n1, n3)
                    for psi in np.linspace(0, 2 * pi, 48, endpoint=False):
                        A = snap(r * np.exp(1j * psi), a)
                        er = A / abs(A)
                        d = np.cos(beta) * 1j * er + np.sin(beta) * er
                        add(solve_square(np.array([A.real, A.imag, np.angle(d), t]), a=a))
    return list(found.values())


def aligned_near_squares(Theta_A, n3, sgn_b=1, sgn_d=1, a=A_EXP):
    """Aligned near-squares: A and B on the same pass (B = arm point at Theta_A + sgn_b dT, dT < pi/2),
    D = A + i sgn_d (B - A) exactly on the pass n3 levels from A, C free.

    Returns (side, m) for each solution, where m = (Phi(A) - Phi(B) + Phi(C) - Phi(D)) / pi. A genuine square
    needs m to be an integer; the selection rule predicts m close to 0 but not equal to it."""
    rA = Theta_A ** (-a)
    A = rA * np.exp(1j * Theta_A)

    def level(dT):
        TB = Theta_A + sgn_b * dT
        B = TB ** (-a) * np.exp(1j * TB)
        D = A + 1j * sgn_d * (B - A)
        ph = phase_rel([A, B, B + D - A, D], a)
        return (ph[3] - ph[0]) / pi, ph, abs(B - A)

    out = []
    grid = np.linspace(1e-4, pi / 2, 4000)
    vals = [level(x)[0] for x in grid]
    for x0, x1, v0, v1 in zip(grid, grid[1:], vals, vals[1:]):
        for target in (n3, -n3):
            if (v0 - target) * (v1 - target) < 0:
                lo, hi = x0, x1
                for _ in range(80):
                    mid = (lo + hi) / 2
                    if (level(mid)[0] - target) * (v0 - target) > 0:
                        lo = mid
                    else:
                        hi = mid
                _, ph, t = level((lo + hi) / 2)
                out.append((t, (ph[0] - ph[1] + ph[2] - ph[3]) / pi))
    return out
