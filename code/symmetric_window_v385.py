"""v385: the path condition holds automatically when the asymmetry stays below the invisibility radius.

Setting of v374-v376: a two-arm point p = 0 with arms Theta^(-a) e^{i Theta} and Theta^(-a) e^{i(Theta + pi + phi(r))}.
A point z near radius r lies on arm j iff the comb phase
    h_j(z) = |z|^(-1/a) - arg z + [j = 2] (pi + phi(|z|))
is in 2 pi Z. For a polygon inscribed in the arms, the drawn points miss this by at most a phase defect nu.

Lemma 385.1. A near-square of side s (three corners on the arc) whose free corner A is within eta of p puts its two
neighbours B = A + s e^{i theta}, D = A + i s e^{i theta} on the arms. Their phase difference is pi/2 plus the arm
offsets plus errors controlled by eta. The offsets are 0, +-(pi + phi) or a difference of phi values, so the
difference can reach 2 pi Z only if
    2 arcsin(eta/s) + (eta sigma / a) (s - eta)^(-1/a - 1) + 2 nu >= pi/2 - Phi*,
    sigma = 2 sin(pi/4 + arcsin(eta/s)),  Phi* = max(sup |phi|, oscillation of phi) on [s - eta, s + eta].
The largest eta violating this, eta*(s), is a lower bound for the distance from p to every near-square free corner,
hence (Lemma 384.1) for the invisibility radius. To leading order eta* = (pi/2 - Phi* - 2 nu) a s^(1+1/a) / sqrt 2.
"""
import numpy as np
from math import asin, log, pi, sin, sqrt


def lemma_lhs(eta, s, a=0.5, nu=0.0):
    sig = 2 * sin(pi / 4 + asin(eta / s))
    return 2 * asin(eta / s) + eta * sig / a * (s - eta) ** (-1 / a - 1) + 2 * nu


def eta_star(s, a=0.5, Phi=0.0, nu=0.0):
    """Largest eta < s/2 with lemma_lhs(eta) < pi/2 - Phi (bisection; the left side increases with eta)."""
    m = pi / 2 - Phi
    if m <= 2 * nu:
        return 0.0
    lo, hi = 0.0, s / 2
    if lemma_lhs(hi * (1 - 1e-12), s, a, nu) < m:
        return hi
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if lemma_lhs(mid, s, a, nu) < m else (lo, mid)
    return lo


def eta_star_leading(s, a=0.5, Phi=0.0, nu=0.0):
    return (pi / 2 - Phi - 2 * nu) * a * s ** (1 + 1 / a) / sqrt(2)


def twist_profile(lam):
    return lambda r: lam * 0.4 * sin(5 * log(1 / r))


def Phi_star(s, eta, phi, n=200):
    """max(sup |phi|, oscillation of phi) on [s - eta, s + eta]."""
    rr = np.linspace(s - eta, s + eta, n)
    v = np.array([phi(r) for r in rr])
    return float(max(np.max(np.abs(v)), np.max(v) - np.min(v)))


def eta_star_for(s, phi, a=0.5, nu=0.0):
    """eta* with Phi* evaluated self-consistently on the window it uses."""
    eta = eta_star(s, a, 0.0, nu)
    for _ in range(20):
        P = Phi_star(s, eta, phi)
        new = eta_star(s, a, P, nu)
        if abs(new - eta) <= 1e-15 * s:
            break
        eta = min(eta, new) if new < eta else new
    return eta, P


def phase_defect(h, r, a=0.5, L=0.0):
    """Bound for the comb-phase error of a point within h of the arm point at radius r (Lemma 385.1, polygons)."""
    return (1 / a) * (r - h) ** (-1 / a - 1) * h + asin(min(1.0, h / (r - h))) + L * log(r / (r - h))


# ---------------------------------------------------------------- the window theorem

def d_g(r, a=0.5, g=pi):
    """v375 turn spacing at radius r."""
    return g * a * r ** (1 + 1 / a) / (1 + pi / 2 * a * r ** (1 / a))


def best_cutoff(eps, a=0.5, g=pi):
    """Largest cutoff for which condition (a) holds on the symmetrized smooth arc at radius eps: s = d_g(eps/8)."""
    return d_g(eps / 8, a, g)


def window_inner_radius(eps, a=0.5, g=pi, nu=0.0):
    """Largest rho such that symmetry on [rho, eps] gives a square (Theorem 385.3): min(s/2, eta*(s)), s = d_g(eps/8)."""
    s = best_cutoff(eps, a, g)
    return min(s / 2, eta_star(s, a, 0.0, nu)), s


# ---------------------------------------------------------------- checks

def lattice_check(s_values, lam, M=4, a=0.5):
    """Nearest near-square free corner (all 8 arm triples, v378 lattice) against eta*(s)."""
    from surgery_v378 import lattice_count
    phi = twist_profile(lam)
    out = []
    for s in s_values:
        found = lattice_count(s, M, lam, a)
        Ks = a * s ** (1 + 1 / a)
        dmin = min(np.hypot(*z) for _, z in found) * Ks
        es, P = eta_star_for(s, phi, a)
        out.append(dict(s=float(s), nearest=float(dmin), eta_star=float(es), ratio=float(dmin / es), Phi=P,
                        count=len(found)))
    return out


def solver_check(s_values, lam, a=0.5, n_theta=6000):
    """Same comparison with v377's independent near-square solver (free corner in the comb cell nearest p)."""
    from invisibility_v377 import near_squares
    phi = twist_profile(lam)
    out = []
    for s in s_values:
        sq = near_squares(s, lam, n_theta=n_theta, a=a)
        if not sq:
            continue
        dmin = min(abs(A) for A, _, _ in sq)
        es, P = eta_star_for(s, phi, a)
        out.append(dict(s=float(s), nearest=float(dmin), eta_star=float(es), ratio=float(dmin / es)))
    return out


def run_checks(path=None):
    import json
    res = {}
    S = np.geomspace(0.003, 0.03, 80)
    for lam in (0, 1):
        rows = lattice_check(S, lam)
        res[f"lattice_lam{lam}"] = dict(rows=rows, min_ratio=min(r["ratio"] for r in rows),
                                         median_ratio=float(np.median([r["ratio"] for r in rows])))
    for lam in (0, 1):
        rows = solver_check([0.02, 0.01, 0.005], lam)
        res[f"solver_lam{lam}"] = dict(rows=rows, min_ratio=min(r["ratio"] for r in rows))
    res["leading_order"] = [dict(s=s, exact=eta_star(s), leading=eta_star_leading(s)) for s in (0.05, 0.02, 0.01, 0.005)]
    res["window"] = []
    for eps in (0.3, 0.1, 0.03, 0.01):
        rho, s = window_inner_radius(eps)
        res["window"].append(dict(eps=eps, s=s, rho_max=rho, log_ratio=log(rho) / log(eps)))
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    from pathlib import Path
    import json
    r = run_checks(Path(__file__).resolve().parents[1] / "results" / "v385" / "window_check.json")
    print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "rows"})
                      for k, v in r.items()}, indent=1))
