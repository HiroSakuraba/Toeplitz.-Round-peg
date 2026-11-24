"""v380: comb laws for arbitrary arm shapes, tested on wobbling and backing-up spirals.

General comb law (Theorem 380.1): for any C^1 or polygonal arc and any inscribed square, at least two sides join
corners whose connecting subarc turns by at least pi/4. So the side is at least the distance from some corner to
the nearest point of the arc that is a quarter-turn's worth of turning away along the arc.

Test family (half-turn symmetric, two arms): arm one is r(Theta) e^{i Theta} with
    r(Theta) = Theta^(-a) (1 + (kappa / Theta) sin(omega Theta)),   Theta >= theta0,
arm two is its negative. Jordan when kappa < a pi / 2; the radius backs up (non-monotone) when kappa omega > a.
"""
import numpy as np
from math import pi

A_EXP = 0.5
THETA0 = 3.0


def radius(Th, kappa, omega, a=A_EXP):
    return Th ** (-a) * (1 + kappa / Th * np.sin(omega * Th))


def point(Th, arm, kappa, omega, a=A_EXP):
    z = radius(Th, kappa, omega, a) * np.exp(1j * Th)
    return z if arm == 1 else -z


def tangent_angle(Th, kappa, omega, a=A_EXP):
    """Direction angle of d/dTheta [r(Theta) e^{i Theta}] (the same for arm two, which is a half-turn copy)."""
    h = 1e-6
    r = radius(Th, kappa, omega, a)
    dr = (radius(Th + h, kappa, omega, a) - radius(Th - h, kappa, omega, a)) / (2 * h)
    return Th + np.arctan2(r, dr)


def turning(arm1, Th1, arm2, Th2, kappa, omega, a=A_EXP, n=400):
    """Range of the tangent angle along the subarc between two arc points (infinite through the wild point)."""
    if arm1 != arm2:
        return np.inf
    lo, hi = sorted((Th1, Th2))
    ang = tangent_angle(np.linspace(lo, hi, n), kappa, omega, a)
    return float(ang.max() - ang.min())


def pass_of(z, kappa, omega, a=A_EXP, kmax=400):
    """Nearest pass to z: (Theta on arm one's lift, arm). Passes at angle psi are Theta = psi + pi k."""
    psi = np.angle(z)
    best = None
    for k in range(kmax):
        Th = psi + pi * k
        if Th < THETA0:
            continue
        d = abs(abs(z) - radius(Th, kappa, omega, a))
        if best is None or d < best[0]:
            best = (d, Th, 1 if k % 2 == 0 else 2)
        if radius(Th, kappa, omega, a) < abs(z) * 0.5:
            break
    return best[1], best[2]


def corners(v):
    A = v[0] + 1j * v[1]
    e = np.exp(1j * v[2])
    t = v[3]
    return [A, A + t * e, A + (1 + 1j) * t * e, A + 1j * t * e]


def solve_square(v, kappa, omega, a=A_EXP, it=80):
    """Newton with each corner assigned to the pass nearest the seed (smooth equations |z| = r(psi + pi k))."""
    Z0 = corners(v)
    lifts = []
    for z in Z0:
        Th, _ = pass_of(z, kappa, omega, a)
        lifts.append(Th - np.angle(z))          # = pi k for the chosen pass

    def F(w):
        out = []
        for z, L in zip(corners(w), lifts):
            T = np.angle(z) + L
            out.append(abs(z) - radius(T, kappa, omega, a) if T > 1.0 else 1e3)
        return np.array(out)

    for _ in range(it):
        y = F(v)
        if np.max(np.abs(y)) < 1e-14:
            break
        J = np.zeros((4, 4))
        t = abs(v[3])
        h = 1e-8 * t
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
    if np.max(np.abs(F(v))) > 1e-12:
        return None
    # every corner must sit on the arc (Theta >= theta0) and the angle jump must not have changed the pass
    for z, L in zip(corners(v), lifts):
        if np.angle(z) + L < THETA0:
            return None
    return v, [np.angle(z) + L for z, L in zip(corners(v), lifts)]


def check_theorem(v, Ths, kappa, omega, a=A_EXP):
    """Number of square sides whose corners are joined by a subarc turning by at least pi/4 (theorem: >= 2)."""
    arms = [1 if round((T - np.angle(z)) / pi) % 2 == 0 else 2 for z, T in zip(corners(v), Ths)]
    sides = [(0, 1), (1, 2), (2, 3), (3, 0)]
    return sum(1 for i, k in sides if turning(arms[i], Ths[i], arms[k], Ths[k], kappa, omega, a) >= pi / 4)


def quarter_gap(z_Th, arm, kappa, omega, a=A_EXP, n=6000):
    """G(z): distance from the arc point to the nearest arc point a quarter-turn of turning away (or on the other arm)."""
    z = point(z_Th, arm, kappa, omega, a)
    Th = np.linspace(THETA0, z_Th + 40 * pi, n * 8)
    other = point(Th, 3 - arm, kappa, omega, a)
    g = np.min(np.abs(other - z))
    same = point(Th, arm, kappa, omega, a)
    ang = tangent_angle(Th, kappa, omega, a)
    a0 = tangent_angle(z_Th, kappa, omega, a)
    # running range of tangent angle from z_Th outward and inward
    idx = np.searchsorted(Th, z_Th)
    up = np.maximum.accumulate(ang[idx:]) - np.minimum.accumulate(ang[idx:])
    up = np.maximum(up, np.abs(ang[idx:] - a0))
    far = idx + np.nonzero(up >= pi / 4)[0]
    if len(far):
        g = min(g, np.min(np.abs(same[far] - z)))
    dn = (np.maximum.accumulate(ang[:idx][::-1]) - np.minimum.accumulate(ang[:idx][::-1]))
    far = idx - 1 - np.nonzero(dn >= pi / 4)[0]
    if len(far):
        g = min(g, np.min(np.abs(same[far] - z)))
    return float(g)


def search(kappa, omega, n_random=6000, seed=0, a=A_EXP):
    rng = np.random.default_rng(seed)
    found = {}
    for _ in range(n_random):
        r = rng.uniform(0.12, 0.38)
        D = pi * a * r ** (1 + 1 / a)
        psi = rng.uniform(0, 2 * pi)
        Th, arm = pass_of(r * np.exp(1j * psi), kappa, omega, a)
        A = point(Th, arm, kappa, omega, a)
        out = solve_square(np.array([A.real, A.imag, rng.uniform(0, 2 * pi), D * rng.uniform(0.6, 6)]),
                           kappa, omega, a)
        if out is None:
            continue
        v, Ths = out
        Z = corners(v)
        if v[3] > 0.08 or max(abs(z) for z in Z) > THETA0 ** -a * 0.98:
            continue
        key = tuple(np.round(sorted((z.real, z.imag) for z in Z), 7).ravel())
        found.setdefault(key, (v, Ths))
    return list(found.values())


def gap_criticality(v, Ths, kappa, omega, a=A_EXP):
    """For an aligned square (two short sides), the pair of neighbouring passes it straddles and where it sits
    relative to the gap function g(psi) = r(psi + L1) - r(psi + L2) between them.

    Returns (side / g(psi_c), |g'(psi_c)| / max |g'| over the surrounding turn)."""
    Z = corners(v)
    lifts = sorted(set(round(T - np.angle(z), 6) for z, T in zip(Z, Ths)))
    if len(lifts) != 2:
        return None
    L1, L2 = lifts
    psi = np.angle(sum(Z) / 4)
    g = lambda s: radius(s + L1, kappa, omega, a) - radius(s + L2, kappa, omega, a)
    h = 1e-5
    dg = lambda s: (g(s + h) - g(s - h)) / (2 * h)
    ss = np.linspace(psi - pi, psi + pi, 2001)
    amp = np.max(np.abs(dg(ss)))
    return float(v[3] / abs(g(psi))), float(abs(dg(psi)) / amp)


def polygon_far_sides(V, entry):
    """Theorem 380.1 on a polygon: cut the closed polygon at its last vertex (an outer closure vertex) to get an arc,
    and count square sides whose corners are joined along it by a subarc turning by at least pi/4."""
    n = len(V)
    P = np.array([[float(x), float(y)] for x, y in V])
    d = np.roll(P, -1, axis=0) - P
    ang = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    pos = [float(q) for q in entry[2]]
    def turn(p1, p2):
        lo, hi = sorted((p1, p2))
        i0, i1 = int(lo), int(hi)
        seg = ang[i0:i1 + 1]
        return float(seg.max() - seg.min())
    sides = [(0, 1), (1, 2), (2, 3), (3, 0)]
    return sum(1 for i, k in sides if turn(pos[i], pos[k]) >= pi / 4)
