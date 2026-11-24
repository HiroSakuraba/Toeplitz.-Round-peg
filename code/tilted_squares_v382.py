"""v382: tilted squares exist at every resonance.

Half-turn spiral of v379: the arc is {z : Phi(z) = |z|^(-1/a) - arg z in pi Z}. Put corner A on arm one at lift
Theta_A, and place B = A + t e^{i theta}, D = A + i t e^{i theta} exactly on the passes n1 and n3 levels from A:
    Phi(B) - Phi(A) = pi n1,   Phi(D) - Phi(A) = pi n3      (two equations for theta, t).
The opposite corner C = B + D - A is then on the arc iff
    h(Theta_A) := (Phi(C) - Phi(A)) / pi - n1 - n3
is an integer j, and that j is the alternating sum of v379. To leading order h = (2 + 1/a) pi K(r) n1 n3, which
falls steadily to 0 as A moves in, so h crosses every integer j between 0 and its starting value: one square each.
"""
import numpy as np
from math import pi

A_EXP = 0.5


def phi_rel(z, A, a=A_EXP):
    """Phi at z with the argument lifted continuously from A (Phi(A) uses the lift Theta_A itself)."""
    return abs(z) ** (-1 / a) - (np.angle(A) + np.angle(z / A))


def arm_point(Th, a=A_EXP):
    return Th ** (-a) * np.exp(1j * Th)


def three_corners(Th_A, n1, n3, guess=None, a=A_EXP, it=60):
    """Solve for (theta, t) with B, D on their passes; returns (theta, t) or None."""
    A = arm_point(Th_A, a)
    PA = phi_rel(A, A, a)
    r = abs(A)
    K = a * r ** (1 / a)
    if guess is None:
        # straight-comb guess: radial components of AB and AD are -pi K r n1 and -pi K r n3 (Phi grows inward)
        er, ep = A / r, 1j * A / r
        beta = np.arctan2(-n1, n3)
        t = pi * K * r * np.hypot(n1, n3)
        theta = np.angle(np.cos(beta) * ep + np.sin(beta) * er)
        x = np.array([theta, t])
    else:
        x = np.array(guess, dtype=float)

    def F(v):
        e = np.exp(1j * v[0])
        B, D = A + v[1] * e, A + 1j * v[1] * e
        return np.array([phi_rel(B, A, a) - PA - pi * n1, phi_rel(D, A, a) - PA - pi * n3])

    for _ in range(it):
        y = F(x)
        if np.max(np.abs(y)) < 1e-12:
            break
        J = np.zeros((2, 2))
        h = 1e-8
        for k in range(2):
            dx = np.zeros(2)
            dx[k] = h if k == 0 else h * x[1]
            J[:, k] = (F(x + dx) - y) / dx[k]
        try:
            st = np.linalg.solve(J, -y)
        except np.linalg.LinAlgError:
            return None
        sc = max(1.0, abs(st[0]) / 0.05, abs(st[1]) / (0.1 * x[1]))
        x = x + st / sc
        if x[1] <= 0:
            return None
    return (x[0], x[1]) if np.max(np.abs(F(x))) < 1e-10 else None


def h_value(Th_A, n1, n3, guess=None, a=A_EXP):
    sol = three_corners(Th_A, n1, n3, guess, a)
    if sol is None:
        return None, None
    theta, t = sol
    A = arm_point(Th_A, a)
    e = np.exp(1j * theta)
    C = A + (1 + 1j) * t * e
    return (phi_rel(C, A, a) - phi_rel(A, A, a)) / pi - n1 - n3, sol


def scan(n1, n3, Th_lo, Th_hi, steps=4000, a=A_EXP, rmax=None, max_ratio=0.3):
    """h along arm one, continued from the inner end outward, and the squares where h crosses a nonzero integer.

    Only the regime of the theorem is used: side at most max_ratio |A|, every corner within a quarter turn of A
    (beyond that the phase lift around A, and so the level count, is no longer meaningful)."""
    Ths = np.linspace(Th_hi, Th_lo, steps)       # start deep (small h), move outward
    vals, sols, guess = [], [], None
    for T in Ths:
        v, sol = h_value(T, n1, n3, guess, a)
        if v is not None:
            A = arm_point(T, a)
            e = np.exp(1j * sol[0])
            Zt = [A + sol[1] * e, A + (1 + 1j) * sol[1] * e, A + 1j * sol[1] * e]
            if sol[1] > max_ratio * abs(A) or max(abs(np.angle(z / A)) for z in Zt) > pi / 2:
                v = None
        if v is None:
            vals.append(None)
            sols.append(None)
            guess = None
            continue
        vals.append(v)
        sols.append(sol)
        guess = sol
    squares = []
    for i in range(len(Ths) - 1):
        v0, v1 = vals[i], vals[i + 1]
        if v0 is None or v1 is None:
            continue
        for j in range(int(np.floor(min(v0, v1))) + 1, int(np.floor(max(v0, v1))) + 1):
            if j == 0:
                continue
            lo, hi, glo = Ths[i], Ths[i + 1], sols[i]
            for _ in range(60):
                mid = (lo + hi) / 2
                vm, sm = h_value(mid, n1, n3, glo, a)
                if vm is None:
                    break
                if (vm - j) * (v0 - j) > 0:
                    lo, glo = mid, sm
                else:
                    hi = mid
            vm, sm = h_value((lo + hi) / 2, n1, n3, glo, a)
            if vm is None or abs(vm - j) > 1e-8:
                continue
            A = arm_point((lo + hi) / 2, a)
            e = np.exp(1j * sm[0])
            Z = [A, A + sm[1] * e, A + (1 + 1j) * sm[1] * e, A + 1j * sm[1] * e]
            if rmax is not None and max(abs(z) for z in Z) > rmax:
                continue
            squares.append(dict(j=j, Theta_A=(lo + hi) / 2, side=sm[1], corners=Z))
    return Ths, vals, squares


def arc_residual(Z, a=A_EXP):
    """Largest |sin(arg z - |z|^(-1/a))| over the corners: zero exactly on the arc."""
    return max(abs(np.sin(np.angle(z) - abs(z) ** (-1 / a))) for z in Z)


def predicted_radius(n1, n3, j, a=A_EXP):
    K = j / ((2 + 1 / a) * pi * n1 * n3)
    return (K / a) ** a if K > 0 else None
