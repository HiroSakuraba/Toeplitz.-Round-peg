"""v388: protection for symmetric windows on smooth arms, with a slowly turned-on twist.

Combines v385 (asymmetry below eta*(s) is invisible), v386 (smooth parity theorems; graceful squares on twisted smooth
arms are confined by the twist's size through c2 sqrt(1 - Phi/pi) and by its gradient L through c_L / L) and v387 (the
classification those bounds use holds for moderate L).

Graceful confinement radius (leading order, Proposition 386.5):
    r_g(t) = max( r_smooth(t) (1 - Phi/pi)^(-a/(2a+1)),  (L t / c_L)^(a/(1+a)) ),   c_L = (2a+1) pi^2 / (2 sqrt 2).
Theorem 388.1: if phi = 0 on [rho, eps], the twist outside has |phi| <= Phi and |r phi'| <= L, and for some cutoff s
    8 r_g(s) <= eps   and   rho < min(s/2, eta*(s)),
then the point inscribes a square of side greater than s.
"""
import numpy as np
from math import pi, sqrt

from smooth_cores_v386 import c2
from symmetric_window_v385 import eta_star, window_inner_radius


def c_L(a=0.5):
    return (2 * a + 1) * pi ** 2 / (2 * sqrt(2))


def r_graceful(t, a=0.5, Phi=0.0, L=0.0):
    rs = (t / (c2(a) * sqrt(1 - Phi / pi))) ** (2 * a / (2 * a + 1))
    rl = (L * t / c_L(a)) ** (a / (1 + a)) if L > 0 else 0.0
    return max(rs, rl)


def best_cutoff(eps, a=0.5, Phi=0.0, L=0.0):
    """Largest s with 8 r_g(s) <= eps (r_g increases with s)."""
    lo, hi = 0.0, eps
    for _ in range(200):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if 8 * r_graceful(m, a, Phi, L) <= eps else (lo, m)
    return lo


def rho_max(eps, a=0.5, Phi=0.0, L=0.0):
    s = best_cutoff(eps, a, Phi, L)
    return min(s / 2, eta_star(s, a)), s


def exponent_table(eps_list=(0.1, 0.01, 1e-3, 1e-5), L_list=(0.0, 0.01, 0.1, 1.0, 2.0), a=0.5, Phi=0.4):
    rows = []
    for eps in eps_list:
        poly, s_poly = window_inner_radius(eps, a)
        for L in L_list:
            r, s = rho_max(eps, a, Phi, L)
            rows.append(dict(eps=eps, L=L, s=s, rho_max=r, polygon_rho_max=poly, gain=r / poly,
                             log_ratio=float(np.log(r) / np.log(eps))))
    return rows


def graceful_check(lams=(1.0, 2.0, 3.0), n_seeds=3000, a=0.5):
    """Graceful squares on twisted smooth arms (twist lam * 0.4 sin(5 log 1/r)): is r_min <= r_g(t)?"""
    from smooth_cores_v386 import sector_two_arm
    from twisted_classification_v387 import default_phi_derivs, twisted_random_search
    phi, _, _ = default_phi_derivs()
    out = {}
    for lam in lams:
        f = lambda r, lam=lam: lam * phi(r)
        found = twisted_random_search(1.0, f, n_seeds=n_seeds, seed=3, a=a)
        Phi, L = 0.4 * lam, 2.0 * lam
        rows = []
        for Z, ks in found:
            sec = sector_two_arm(ks, Z, a)
            t = abs(Z[1] - Z[0])
            rmin = min(abs(z) for z in Z)
            rows.append(dict(sector=sec, t=t, rmin=rmin, ratio=rmin / r_graceful(t, a, Phi, L),
                             ratio_plain=rmin / r_graceful(t, a, 0.0, 0.0)))
        g = [r for r in rows if r["sector"] == "graceful"]
        out[str(lam)] = dict(squares=len(rows), graceful=len(g),
                             max_ratio_graceful=max((r["ratio"] for r in g), default=None),
                             max_plain_ratio_graceful=max((r["ratio_plain"] for r in g), default=None),
                             max_plain_ratio_crossed=max((r["ratio_plain"] for r in rows if r["sector"] == "crossed"),
                                                         default=None))
    return out


def run_checks(path=None):
    import json
    res = dict(exponents=exponent_table(), graceful=graceful_check())
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import json
    from pathlib import Path
    print(json.dumps(run_checks(Path(__file__).resolve().parents[1] / "results" / "v388" / "slow_twist_check.json"),
                     indent=1))
