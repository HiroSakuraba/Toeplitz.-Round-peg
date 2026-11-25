"""v386: smooth cores. An exact selection identity, and smooth confinement made rigorous.

Exact flag selection identity (Theorem 386.3). Corners A, B, C, D (cyclic) of a square on passes of the half-turn
spiral rho = (theta + pi k)^(-a), lifts Theta_X = theta_X + pi k_X, squared radii F(Theta) = Theta^(-2a). With
delta_X = Theta_X - Theta_A, j = k_A + k_C - k_B - k_D and Dtheta = theta_A + theta_C - theta_B - theta_D,
    pi j + Dtheta = (2a + 1) zeta^(2a+1) xi^(-2a-2) delta_B delta_D        (exact)
for some xi, zeta in the hull [2 Theta_lo - Theta_hi, 2 Theta_hi - Theta_lo]. This is the British flag identity
F(Theta_A) + F(Theta_C) = F(Theta_B) + F(Theta_D) written as a mixed second difference plus one first difference.
It is v379's selection rule with no error terms.

With a twist, a corner on arm two has Theta_X = theta_X + pi k_X - phi(rho_X) (pi k_X - pi plus pi), so the left side
gains the alternating sum of phi over the arm-two corners: the effective level j + Sigma/pi.
"""
import numpy as np
from math import asin, pi, sqrt

from two_turn_squares_v383 import window, split_type


def c2(a=0.5):
    """Smooth comb constant: t >= c2 r^(1 + 1/(2a)) (v379 Corollary 379.4)."""
    return a * sqrt(2 * pi / (2 * a + 1))


def r_smooth(t, a=0.5):
    return (t / c2(a)) ** (2 * a / (2 * a + 1))


def choose_A(k):
    """Index of a corner whose two neighbours are on passes other than its own (None if there is none)."""
    for i in range(4):
        if k[(i + 1) % 4] != k[i] and k[(i + 3) % 4] != k[i]:
            return i
    return None


def selection_data(Z, lifts, a=0.5):
    """Exact selection identity data for a square on the plain half-turn spiral."""
    th, k, W = window(Z, lifts)
    i = choose_A(k)
    if i is None:
        return None
    o = [(i + m) % 4 for m in range(4)]                    # A, B, C, D
    T = [th[x] + pi * k[x] for x in o]
    rho = [abs(Z[x]) for x in o]
    j = k[o[0]] + k[o[2]] - k[o[1]] - k[o[3]]
    dth = th[o[0]] + th[o[2]] - th[o[1]] - th[o[3]]
    dB, dD = T[1] - T[0], T[3] - T[0]
    lo, hi = min(rho), max(rho)
    Tlo, Thi = hi ** (-1 / a), lo ** (-1 / a)
    L1, L2 = 2 * Tlo - Thi, 2 * Thi - Tlo
    # implied value of zeta^(2a+1) / xi^(2a+2), and its admissible range
    implied = (pi * j + dth) / ((2 * a + 1) * dB * dD)
    rng = (L1 ** (2 * a + 1) / L2 ** (2 * a + 2), L2 ** (2 * a + 1) / L1 ** (2 * a + 2))
    t = abs(Z[1] - Z[0])
    Zc = sum(Z) / 4
    flag = abs(Z[o[0]]) ** 2 + abs(Z[o[2]]) ** 2 - abs(Z[o[1]]) ** 2 - abs(Z[o[3]]) ** 2
    return dict(j=int(j), dtheta=float(dth), dB=float(dB), dD=float(dD), implied=float(implied),
                lo=float(rng[0]), hi=float(rng[1]), inside=bool(rng[0] * (1 - 1e-9) <= implied <= rng[1] * (1 + 1e-9)),
                t=float(t), rmin=float(lo), rmax=float(hi), absZ=float(abs(Zc)), W=float(W),
                pattern=split_type(k), flag=float(flag), dtheta_bound=float(2 * asin(min(1, t * t / (2 * abs(Zc) ** 2)))))


def confinement_terms(d, a=0.5):
    """The two alternatives of Theorem 386.4 evaluated on one square: returns (case, lhs, rhs, holds)."""
    t, r, R = d["t"], d["rmin"], d["rmax"]
    Tlo, Thi = R ** (-1 / a), r ** (-1 / a)
    L1, L2 = 2 * Tlo - Thi, 2 * Thi - Tlo
    W = d["W"]
    db = d["dtheta_bound"]
    if d["j"] == 0:
        lhs, rhs = db, (2 * a + 1) * (pi - W) ** 2 * L1 ** (2 * a + 1) / L2 ** (2 * a + 2)
        return "j=0", lhs, rhs, lhs >= rhs
    lhs = pi - db
    rhs = (2 * a + 1) * L2 ** (2 * a + 1) / L1 ** (2 * a + 2) * t * t / (2 * a * a) * r ** (-2 / a - 2) \
        * (1 + t / (sqrt(2) * r)) ** 2
    return "j!=0", lhs, rhs, lhs <= rhs


# ---------------------------------------------------------------- twisted arms

def twisted_phase(z, theta_unwrapped, k, lam, a=0.5, phi=None):
    """Comb phase residual of a corner meant to lie on pass k (arm one for even k, arm two for odd k):
    |z|^(-1/a) - (theta + pi k) + [k odd] (pi - pi + phi) ... written as Theta - theta - pi k + [odd] phi(|z|)."""
    r = abs(z)
    val = r ** (-1 / a) - theta_unwrapped - pi * k
    if k % 2:
        val = val + lam * phi(r)
    return val


def solve_twisted(v, ks, lam, phi, a=0.5, it=80):
    """Newton for a square (A, angle, side) with corner X on pass ks[X] of the twisted spiral."""
    def corners(w):
        A = w[0] + 1j * w[1]
        e = np.exp(1j * w[2])
        return [A, A + w[3] * e, A + (1 + 1j) * w[3] * e, A + 1j * w[3] * e]

    def F(w):
        Z = corners(w)
        th, _, _ = window(Z, [0, 0, 0, 0])
        return np.array([twisted_phase(z, t, kk, lam, a, phi) for z, t, kk in zip(Z, th, ks)])

    for _ in range(it):
        y = F(v)
        if np.max(np.abs(y)) < 1e-12:
            break
        J = np.zeros((4, 4))
        t = abs(v[3])
        h = 1e-8 * t
        for m in range(4):
            dv = np.zeros(4)
            dv[m] = h if m != 2 else h / t
            J[:, m] = (F(v + dv) - y) / dv[m]
        try:
            st = np.linalg.solve(J, -y)
        except np.linalg.LinAlgError:
            return None
        sc = max(1, np.max(np.abs(st[[0, 1, 3]])) / (0.1 * t), abs(st[2]) / 0.1)
        v = v + st / sc
        if v[3] <= 0:
            return None
    if np.max(np.abs(F(v))) > 1e-10:
        return None
    return v, corners(v)


def default_twist(r):
    from math import log, sin
    return 0.4 * sin(5 * log(1 / r))


def twisted_selection(Z, ks, phi, a=0.5):
    """Identity data on the twisted spiral: lifts from radii (Theta = rho^(-1/a) on both arms), effective level
    j_eff = (pi j - Sigma*) / pi with Sigma* the alternating sum of phi over arm-two (odd) corners."""
    th, _, W = window(Z, [0, 0, 0, 0])
    k = list(ks)
    i = choose_A(k)
    o = [(i + m) % 4 for m in range(4)]
    rho = [abs(Z[x]) for x in o]
    T = [r ** (-1 / a) for r in rho]
    kk = [k[x] for x in o]
    j = kk[0] + kk[2] - kk[1] - kk[3]
    ph = [phi(r) if kk[m] % 2 else 0.0 for m, r in enumerate(rho)]
    sig = ph[0] + ph[2] - ph[1] - ph[3]
    dth = th[o[0]] + th[o[2]] - th[o[1]] - th[o[3]]
    dB, dD = T[1] - T[0], T[3] - T[0]
    lo, hi = min(rho), max(rho)
    Tlo, Thi = hi ** (-1 / a), lo ** (-1 / a)
    L1, L2 = 2 * Tlo - Thi, 2 * Thi - Tlo
    implied = (pi * j - sig + dth) / ((2 * a + 1) * dB * dD)
    rng = (L1 ** (2 * a + 1) / L2 ** (2 * a + 2), L2 ** (2 * a + 1) / L1 ** (2 * a + 2))
    t = abs(Z[1] - Z[0])
    odd = [m for m in range(4) if kk[m] % 2]
    kind = {0: "none", 4: "all"}.get(len(odd), "diagonal" if odd in ([0, 2], [1, 3]) else "adjacent" if len(odd) == 2 else "odd")
    return dict(j=int(j), j_eff=float((pi * j - sig) / pi), sigma=float(sig), implied=float(implied), lo=float(rng[0]),
                hi=float(rng[1]), inside=bool(rng[0] * (1 - 1e-9) <= implied <= rng[1] * (1 + 1e-9)), t=float(t),
                rmin=float(lo), ratio=float(t / (c2(a) * lo ** (1 + 1 / (2 * a)))), arm_two=kind, W=float(W),
                t_over_Z=float(t / abs(sum(Z) / 4)))


def twisted_search(radii, labels, lam=1.0, phi=default_twist, a=0.5, max_ratio=0.2):
    found = {}
    for r in radii:
        Th = r ** (-1 / a)
        A = r * np.exp(1j * Th)
        K = a * r ** (1 / a)
        k0 = int(round((Th - np.angle(A)) / pi))
        for n1, n3 in labels:
            for shift in (0, 1):                       # A on arm one (even pass) or arm two (odd pass)
                kA = k0 + shift if (k0 + shift) % 2 == shift else k0 + shift
                if kA % 2 != shift:
                    continue
                Ar = r * np.exp(1j * (Th + (pi + lam * phi(r) if shift else 0)))
                er, ep = Ar / r, 1j * Ar / r
                beta = np.arctan2(-n1, n3)
                t = pi * K * r * np.hypot(n1, n3)
                th = np.angle(np.cos(beta) * ep + np.sin(beta) * er)
                kA = int(round((Th + (pi + lam * phi(r) if shift else 0) - np.angle(Ar)
                                - (lam * phi(r) if shift else 0)) / pi))
                ks = [kA, kA + n1, kA + n1 + n3, kA + n3]
                res = solve_twisted(np.array([Ar.real, Ar.imag, th, t]), ks, lam, phi, a)
                if res is None:
                    continue
                v, Z = res
                if v[3] > max_ratio * abs(sum(Z) / 4):
                    continue
                key = tuple(np.round(sorted((z.real, z.imag) for z in Z), 8).ravel())
                found.setdefault(key, (Z, ks))
    return list(found.values())


def effective_bound_holds(d_tw, Z, a=0.5):
    """Twisted version of alternative (ii): pi |j_eff| - |Dtheta| <= (2a+1) L2^(2a+1)/L1^(2a+2) t^2/(2a^2) r^(-2/a-2)
    (1 + t/(sqrt2 r))^2, with Dtheta bounded by 2 arcsin(t^2 / (2|Z|^2))."""
    t, r = d_tw["t"], d_tw["rmin"]
    R = max(abs(z) for z in Z)
    Tlo, Thi = R ** (-1 / a), r ** (-1 / a)
    L1, L2 = 2 * Tlo - Thi, 2 * Thi - Tlo
    Zc = abs(sum(Z) / 4)
    lhs = pi * abs(d_tw["j_eff"]) - 2 * asin(min(1, t * t / (2 * Zc * Zc)))
    rhs = (2 * a + 1) * L2 ** (2 * a + 1) / L1 ** (2 * a + 2) * t * t / (2 * a * a) * r ** (-2 / a - 2) \
        * (1 + t / (sqrt(2) * r)) ** 2
    return lhs <= rhs


def sector_two_arm(ks, Z, a=0.5):
    """'graceful' or 'crossed' on the two-arm arc: arm one (even passes) runs inward with Theta increasing, then the
    wild point, then arm two (odd passes) outward. Theta = rho^(-1/a) on both arms."""
    keys = []
    for k, z in zip(ks, Z):
        Th = abs(z) ** (-1 / a)
        # (arm, position): arm one first with Theta increasing, then arm two with Theta decreasing.
        # (A float key big - Theta would round to big and tie for Theta << big.)
        keys.append((0, Th) if k % 2 == 0 else (1, -Th))
    order = sorted(range(4), key=lambda i: keys[i])
    m = order.index(0)
    cyc = tuple(order[m:] + order[:m])
    return "graceful" if cyc in ((0, 1, 2, 3), (0, 3, 2, 1)) else "crossed"


def run_checks(path=None):
    import json
    from collections import Counter
    from general_comb_v380 import corners as c380, search
    from tilted_squares_v382 import arm_point, three_corners
    res = {}
    # plain spiral: v382's labelled squares
    d382 = json.load(open(Path(__file__).resolve().parents[1] / "results" / "v382" / "tilted_check.json"))
    rows = []
    for s in d382["squares"]:
        sol = three_corners(s["Theta_A"], s["n1"], s["n3"])
        if sol is None:
            continue
        th, t = sol
        A = arm_point(s["Theta_A"])
        e = np.exp(1j * th)
        Z = [A, A + t * e, A + (1 + 1j) * t * e, A + 1j * t * e]
        d = selection_data(Z, [abs(z) ** -2 for z in Z])
        if d is None:
            continue
        case, l, r, ok = confinement_terms(d)
        rows.append(dict(d, case=case, ok=bool(ok)))
    small = [r for r in rows if r["t"] / r["absZ"] <= 0.2]
    res["v382"] = dict(squares=len(rows), small_ratio=len(small),
                       identity_inside=sum(r["inside"] for r in small), alternative_holds=sum(r["ok"] for r in small),
                       min_t_over_c2r2=min(r["t"] / (c2() * r["rmin"] ** 2) for r in rows),
                       flag_max=max(abs(r["flag"]) for r in rows), cases=dict(Counter(r["case"] for r in small)))
    # plain spiral: v380 random search
    sq = search(0.0, 1.0, n_random=3000, seed=2)
    rows = []
    for v, Ths in sq:
        Z = c380(v)
        d = selection_data(Z, Ths)
        case, l, r, ok = confinement_terms(d)
        rows.append(dict(d, ok=bool(ok)))
    small = [r for r in rows if r["t"] / r["absZ"] <= 0.2]
    res["v380_plain"] = dict(squares=len(rows), small_ratio=len(small),
                             identity_inside=sum(r["inside"] for r in small),
                             alternative_holds=sum(r["ok"] for r in small),
                             min_t_over_c2r2=min(r["t"] / (c2() * r["rmin"] ** 2) for r in rows),
                             max_t_over_Z_of_rest=max((r["t"] / r["absZ"] for r in rows if r["t"] / r["absZ"] > 0.2),
                                                      default=None))
    # twisted spiral: crossed j = 0 squares, broad and near the twist's zeros
    labels = [(n1, n3) for n1 in range(-5, 6) for n3 in range(-5, 6) if n1 and n3]
    found = twisted_search(np.geomspace(0.06, 0.3, 120), labels)
    near = [(1, 1), (1, -1), (-1, 1), (-1, -1), (1, 2), (2, 1), (-1, 2), (2, -1), (1, -2), (-2, 1), (-1, -2), (-2, -1)]
    for r0 in (np.exp(-2 * pi / 5), np.exp(-3 * pi / 5), np.exp(-4 * pi / 5)):
        found += twisted_search(r0 * np.exp(np.linspace(-0.12, 0.12, 161)), near)
    tw = {}
    for Z, ks in found:
        key = tuple(np.round(sorted((z.real, z.imag) for z in Z), 8).ravel())
        if key in tw:
            continue
        d = twisted_selection(Z, ks, default_twist)
        d["sector"] = sector_two_arm(ks, Z)
        d["bound"] = bool(effective_bound_holds(d, Z))
        d["t_over_gap"] = float(d["t"] / (pi * 0.5 * d["rmin"] ** 3))
        tw[key] = d
    tw = list(tw.values())
    res["twisted"] = dict(squares=len(tw), identity_inside=sum(d["inside"] for d in tw),
                          bound_holds=sum(d["bound"] for d in tw),
                          sectors=dict(Counter(d["sector"] for d in tw)),
                          kinds=dict(Counter((d["j"], d["arm_two"]) for d in tw).most_common()) if False else
                          dict(Counter(f"j={d['j']},{d['arm_two']}" for d in tw)),
                          rows=sorted([dict(ratio=d["ratio"], j_eff=d["j_eff"], rmin=d["rmin"], t=d["t"],
                                            t_over_gap=d["t_over_gap"], sector=d["sector"]) for d in tw],
                                      key=lambda x: x["ratio"]))
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


from pathlib import Path  # noqa: E402

if __name__ == "__main__":
    import json
    out = run_checks(Path(__file__).resolve().parents[1] / "results" / "v386" / "smooth_check.json")
    print(json.dumps({k: ({kk: vv for kk, vv in v.items() if kk != "rows"} if isinstance(v, dict) else v)
                      for k, v in out.items()}, indent=1))
    print(out["twisted"]["rows"][:5])
