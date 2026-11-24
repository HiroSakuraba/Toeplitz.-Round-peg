"""v376: the comb law for polygons inscribed in two spiral arms.

The v374 drawings: arm one has vertices at Theta_k = theta0 + k pi/12 (24 per turn, mesh mu = pi/12), rounded to
2^-24; arm two turns each vertex by pi + lam phi(r). For a square on the polygon, each corner x on an arm edge gets
a star x* on the smooth arc (the point with the same interpolation parameter in Theta), with |x - x*| <= h_e, where
h_e = dTheta^2/8 * sup|gamma''| + nu (linear interpolation error plus vertex rounding nu).
"""
from cmath import exp
from math import asin, atan, log, pi, sin, sqrt

from comb_law_v375 import spacing_bound


def phi(r, lam):
    return lam * 0.4 * sin(5 * log(1 / r))


def gamma(arm, Theta, lam, a=0.5):
    r = Theta ** (-a)
    if arm == 1:
        return r * exp(1j * Theta)
    return r * exp(1j * (Theta + pi + phi(r, lam)))


def curvature_factor(K, L=0.0, L2=0.0, a=0.5):
    """sup|gamma''(Theta)| / r on an arm, in terms of K = a r^(1/a) and |r phi'| <= L, |r^2 phi''| <= L2."""
    e = K * L
    return (1 + 1 / a) * K * K + (1 + e) ** 2 + 2 * K * (1 + e) + K * K * (L2 + (1 + 1 / a) * L)


def mesh_for_relative_error(r, theta, a=0.5, g=pi):
    """Largest Theta-mesh with 2h <= theta d_g(r) (to leading order): mu^2 r / 4 <= theta g K r."""
    return sqrt(4 * theta * g * a * r ** (1 / a))


def edge_data(V, depth, lam, theta0=3.0, a=0.5):
    """Map each arm edge index to (arm, Theta_hi_radius_end, Theta_other_end); also the vertex rounding nu."""
    steps = 12 * depth
    Th = [theta0 + depth * pi * k / steps for k in range(steps + 1)]
    n = steps + 1
    nu = 0.0
    for k in range(n):
        nu = max(nu, abs(complex(float(V[k][0]), float(V[k][1])) - gamma(1, Th[k], lam, a)))
        b = V[2 * n - 1 - k]
        nu = max(nu, abs(complex(float(b[0]), float(b[1])) - gamma(2, Th[k], lam, a)))
    edges = {}
    for i in range(n - 1):
        edges[i] = (1, Th[i], Th[i + 1])                       # from Theta_i to Theta_{i+1}
    for i in range(n, 2 * n - 1):
        edges[i] = (2, Th[2 * n - 1 - i], Th[2 * n - 2 - i])   # B runs back outward... start at larger Theta
    return edges, nu


def analyse_square(entry, edges, nu, V, lam, a=0.5, s_mesh=pi / 12):
    """Theorem 376.1 for one square: returns None if a corner is off the arm polygon, else a dict."""
    g, L, L2 = pi - 0.4 * lam, 2.0 * lam, lam * sqrt(104)
    side_sq, pts, pos, _ = entry
    t = sqrt(float(side_sq))
    stars, hs, Rp, om = [], [], 0.0, 0.0
    for q in pos:
        e = int(q)
        if e not in edges:
            return None
        f = float(q - e)
        arm, T0, T1 = edges[e]
        T = T0 + f * (T1 - T0)
        stars.append(gamma(arm, T, lam, a))
        Tl = min(T0, T1)
        r_hi = Tl ** (-a)
        K_hi = a * r_hi ** (1 / a)
        hs.append((T1 - T0) ** 2 / 8 * r_hi * curvature_factor(K_hi, L, L2, a) + nu)
        Rp = max(Rp, r_hi)
        P0, P1 = V[e], V[(e + 1) % len(V)]
        el = sqrt(float(P1[0] - P0[0]) ** 2 + float(P1[1] - P0[1]) ** 2)
        om = max(om, asin(min(1.0, 2 * nu / el)))
    h = max(hs)
    rho_star = min(abs(z) for z in stars)
    tp = t + sqrt(2) * h
    KR = a * Rp ** (1 / a)
    eps1 = KR * L
    tau = atan(KR / (1 - eps1))
    eta = (pi / sqrt(2) + sqrt(2) * L) * tp / rho_star
    ok = KR <= 0.25 and eps1 < 1 and eta <= pi / 4 and (eta + 2 * s_mesh) * (1 + eps1) + tau + om < pi / 2
    bound = spacing_bound(rho_star, a, g) - 2 * h
    return dict(t=t, h=h, rho_star=rho_star, bound=bound, eta=eta, ok=ok, holds=t >= bound,
                min_corner=min(sqrt(float(p[0]) ** 2 + float(p[1]) ** 2) for p in pts))


def refined_thetas(depth, theta0=3.0, mu0=pi / 12):
    """Theta-mesh mu(Theta) = mu0 sqrt(theta0 / Theta): proportional to r for a = 1/2, as Corollary 376.4 asks."""
    end = theta0 + depth * pi
    T = [theta0]
    while T[-1] < end:
        T.append(min(end, T[-1] + mu0 * sqrt(theta0 / T[-1])))
    return T


def drawing_from_thetas(T, lam=1, a=0.5, outer=None):
    """The v374 construction on an arbitrary Theta list (exact rational vertices)."""
    from fractions import Fraction as Q
    from math import hypot, tan
    from wild_point_parity_v370 import DEN, _unit, radius
    outer = Q(2) if outer is None else outer
    A = []
    for th in T:
        r = radius(th, a)
        ux, uy = _unit(th)
        A.append((r * ux, r * uy))
    B = []
    for x, y in reversed(A):
        rr = hypot(float(x), float(y))
        u = Q(round(tan(phi(rr, float(lam)) / 2) * DEN), DEN)
        c, sn = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
        B.append((-(c * x - sn * y), -(sn * x + c * y)))
    nrm = Q(round(sqrt(float(A[0][0] ** 2 + A[0][1] ** 2)) * DEN), DEN)
    far = lambda p: (p[0] * outer / nrm, p[1] * outer / nrm)
    ux, uy = _unit(T[0] + pi / 2)
    top = (ux * outer * 2, uy * outer * 2)
    return A + B + [far(B[-1]), top, far(A[0])]
