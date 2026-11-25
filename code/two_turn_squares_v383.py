"""v383: no small square on the plain spiral touches only one or two passes.

Every rectangle P1P2P3P4 and every point O satisfy the British flag identity
    |OP1|^2 + |OP3|^2 = |OP2|^2 + |OP4|^2.
With O at the wild point, the corners' squared radii are F(Theta) = Theta^(-2a) at their lifts. Near the point the
two-arm spiral is the family of polar graphs rho = R_k(theta) = (theta + pi k)^(-a), one per pass k. Splitting the
four corners among at most two passes:

* one pass: impossible on a convex arc that turns by less than pi (three sides of the square would need pi);
* 3-1, or 2-2 with the pairs diagonal: the flag sums cannot balance, because F jumps between passes;
* 2-2 with the pairs adjacent (aligned): the flag identity becomes a slope identity,
      R_in'(eta1) (theta2 - theta1)(rho1 + rho2) = R_out'(eta2) (theta3 - theta4)(rho3 + rho4),
  and the equal chords P1P2 = P4P3 force the ratio Q = R_in'(eta1) / R_out'(eta2) to lie within second-order terms
  of 1. On the plain spiral Q = (Theta_out / Theta_in)^(a+1) is first-order below 1. Contradiction.

This module computes the explicit thresholds, the window data of a square, and the checks.
"""
import numpy as np
from math import pi


# ---------------------------------------------------------------- window data of a square

def window(Z, lifts):
    """Unwrapped polar angles theta_X (all within one short arc around the centre's direction), pass indices k_X with
    lift = theta_X + pi k_X, and the angular width W of the window."""
    Zc = sum(Z) / 4
    th = np.array([np.angle(Zc) + np.angle(z / Zc) for z in Z])
    k = np.array([int(round((L - t) / pi)) for L, t in zip(lifts, th)])
    return th, k, float(th.max() - th.min())


def split_type(k):
    """Pattern of the pass indices of corners A, B, C, D (cyclic order)."""
    k = list(k)
    vals = sorted(set(k))
    if len(vals) == 1:
        return "4"
    if len(vals) == 2:
        c = [k.count(v) for v in vals]
        if sorted(c) == [1, 3]:
            return "3-1"
        return "2-2 diagonal" if k[0] == k[2] else "2-2 adjacent"
    if len(vals) == 3:
        return "2-1-1 diagonal" if (k[0] == k[2] or k[1] == k[3]) else "2-1-1 adjacent"
    return "1-1-1-1"


def flag_residual(Z):
    A, B, C, D = Z
    return abs(A) ** 2 + abs(C) ** 2 - abs(B) ** 2 - abs(D) ** 2


def slope_ratio(Z, th, k):
    """For a 2-2 adjacent square: the inner and outer secant slopes of R and the chord data.

    Labels are rotated so that P1, P2 are on the inner pass (larger k) and P4P3 is the opposite side.
    Returns dict with Q = (inner secant slope) / (outer secant slope) and the chord terms."""
    idx = [0, 1, 2, 3]
    for r in range(4):
        j = idx[r:] + idx[:r]
        if k[j[0]] == k[j[1]] and k[j[2]] == k[j[3]] and k[j[0]] > k[j[2]]:
            break
    else:
        raise ValueError("not a 2-2 adjacent square")
    P = [Z[i] for i in j]
    T = [th[i] for i in j]
    rho = [abs(z) for z in P]
    phi_in, phi_out = T[1] - T[0], T[2] - T[3]
    m_in = (rho[1] - rho[0]) / phi_in
    m_out = (rho[2] - rho[3]) / phi_out
    t = abs(P[1] - P[0])
    return dict(Q=m_in / m_out, m_in=m_in, m_out=m_out, phi_in=phi_in, phi_out=phi_out, rho=rho, t=t)


def q_bounds_measured(d, W):
    """Two-sided bounds for Q from the chord identity alone, using the measured window (valid for any polar graphs)."""
    g, T = consts(W)
    r1, r2, r3, r4 = d["rho"]
    t = d["t"]
    p_out = abs(r3 - r4) / t
    p_in = abs(r2 - r1) / t
    s = np.sqrt(1 - p_out ** 2) * np.sqrt(r1 * r2 / (r3 * r4))      # lower bound, sin ratio out/in
    s_rev = np.sqrt(1 - p_in ** 2) * np.sqrt(r3 * r4 / (r1 * r2))   # lower bound, sin ratio in/out
    fac = (r3 + r4) / (r1 + r2)
    return x_lower(s, T, g) * fac, fac / x_lower(s_rev, T, g)


def x_lower(s, T, g):
    """Lower bound for x = |phi'| / |phi| from sin(|phi'|/2) >= s sin(|phi|/2), both angles at most W <= pi.

    s >= 1: x >= 1 + (s - 1) tan(u)/u >= s.  s < 1: either x >= 1, or the mean value theorem gives
    1 - x <= (1 - s) tan(u)/u <= T (1 - s).  Always x >= s/g, since xu >= sin(xu) >= s sin u >= s u / g."""
    if s >= 1:
        return s
    return max(1 - T * (1 - s), s / g)


# ---------------------------------------------------------------- explicit thresholds (Theorem 383.3)

def consts(W):
    """g(W) = (W/2)/sin(W/2) and T(W) = tan(W/2)/(W/2)."""
    return (W / 2) / np.sin(W / 2), np.tan(W / 2) / (W / 2)


def p_bound(L, a, W):
    """Upper bound for |rho3 - rho4| / t on the outer pass with window floor L."""
    g, _ = consts(W)
    return a * g * (1 + W / L) ** a / L


def sigma(L, a, W):
    p = p_bound(L, a, W)
    return np.sqrt(max(0.0, 1 - p * p)) * (1 + W / L) ** (-2 * a)


def phi(kappa, a, W):
    """phi(kappa) = T - (T - 1) kappa^(-a) - kappa^(a+1): concave, phi(1) = 0."""
    _, T = consts(W)
    return T - (T - 1) * kappa ** (-a) - kappa ** (a + 1)


def cond_small_kappa(L, a, W, k1):
    """(A): for kappa <= k1 the bound Q >= sigma/g beats Q <= kappa^(a+1)."""
    g, _ = consts(W)
    return p_bound(L, a, W) < 1 and sigma(L, a, W) / g > k1 ** (a + 1)


def cond_large_kappa(L, a, W, k1):
    """(B): for k1 <= kappa <= kappa_max, concavity of phi gives phi(kappa) >= phi(k1)(1-kappa)/(1-k1), and
    1 - kappa >= (pi - W)/(L + pi) beats T(1 - sigma) <= T(p^2 + 2aW/L)."""
    _, T = consts(W)
    p = p_bound(L, a, W)
    return p < 1 and (pi - W) / (L + pi) * phi(k1, a, W) / (1 - k1) > T * (p * p + 2 * a * W / L)


def _first_true(f, lo=0.6, hi=1e9):
    """Smallest L in [lo, hi] (to bisection accuracy) beyond which the monotone condition f holds."""
    if not f(hi):
        return np.inf
    if f(lo):
        return lo
    for _ in range(200):
        m = np.sqrt(lo * hi)
        lo, hi = (lo, m) if f(m) else (m, hi)
    return hi


def threshold_adjacent(a, W=1.0, grid=500):
    best = (np.inf, None)
    for k1 in np.linspace(0.5, 0.999, grid):
        if phi(k1, a, W) <= 0:
            continue
        L = max(_first_true(lambda L: cond_small_kappa(L, a, W, k1)),
                _first_true(lambda L: cond_large_kappa(L, a, W, k1)))
        if L < best[0]:
            best = (float(L), float(k1))
    return best


def threshold_three_one(a, W=1.0):
    """Lone corner inside the triple's pass: (pi - W)(L/(L+pi))^(2a+1) > W suffices (lone corner outside: W < pi/2)."""
    return float(_first_true(lambda L: (pi - W) * (L / (L + pi)) ** (2 * a + 1) > W))


def threshold(a, W=1.0):
    La, k1 = threshold_adjacent(a, W)
    L31 = threshold_three_one(a, W)
    return dict(adjacent=La, kappa1=k1, three_one=L31, overall=max(La, L31, 0.5))


def exact_q_gap(L, n, a, W=1.0):
    """The two sides of the adjacent-case contradiction at outer floor L and pass gap n:
    best geometric lower bound for Q, and the spiral's upper bound kappa^(a+1)."""
    g, T = consts(W)
    Lin = L + pi * n
    kappa = (L + W) / Lin
    s = sigma(L, a, W) * kappa ** a
    q_lo = max((1 - T * (1 - s)) * kappa ** (-a), sigma(L, a, W) / g)
    return float(q_lo), float(kappa ** (a + 1))


# ---------------------------------------------------------------- Proposition 383.5 (2-1-1 adjacent)

def ratio_limit_211(a, W=1.0):
    """Largest corner-radius ratio lambda = max|P|/min|P| for which (pi - W) > W lambda^((2a+1)/a)."""
    return float(((pi - W) / W) ** (a / (2 * a + 1)))


def side_limit_211(a, W=1.0):
    """Side bound t <= tau |Z| that forces lambda below ratio_limit_211: lambda <= (1 + t/(sqrt2|Z|))/(1 - t/(sqrt2|Z|)),
    capped by the window condition t <= sqrt2 sin(W/2) |Z|."""
    lam = ratio_limit_211(a, W)
    tau = np.sqrt(2) * (lam - 1) / (lam + 1)
    return float(min(tau, np.sqrt(2) * np.sin(W / 2)))


# ---------------------------------------------------------------- Newton with prescribed passes

def corners(v):
    A = v[0] + 1j * v[1]
    e = np.exp(1j * v[2])
    t = v[3]
    return [A, A + t * e, A + (1 + 1j) * t * e, A + 1j * t * e]


def solve_labelled(v, ks, rfun, it=80, theta_min=1.0):
    """Newton for a square whose corner X lies on pass ks[X]: |z| = rfun(theta_X + pi ks[X]), with theta_X the polar
    angle unwrapped around the centre's direction. Returns (v, lifts) or None."""
    def lifts_of(w):
        Z = corners(w)
        th, _, _ = window(Z, [0, 0, 0, 0])
        return [t + pi * k for t, k in zip(th, ks)]

    def F(w):
        out = []
        for z, L in zip(corners(w), lifts_of(w)):
            out.append(abs(z) - rfun(L) if L > theta_min else 1e3)
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
    return v, lifts_of(v)


PATTERNS = {
    "4": (0, 0, 0, 0),
    "3-1": (0, 0, 0, 1),
    "2-2 diagonal": (0, 1, 0, 1),
    "2-2 adjacent": (0, 0, 1, 1),
    "2-1-1 adjacent": (0, 0, 1, 2),
}


def labelled_search(rfun, a, patterns=PATTERNS, n_seeds=300, seed=0, Th_range=(8.0, 60.0), gaps=(1, 2, 3),
                    max_ratio=0.3):
    """Seed squares with prescribed pass labels near arm one and count Newton solutions per pattern."""
    rng = np.random.default_rng(seed)
    out = {}
    for name, base in patterns.items():
        found = {}
        tries = 0
        for _ in range(n_seeds):
            n = int(rng.choice(gaps)) * int(rng.choice([-1, 1]))
            Th = rng.uniform(*Th_range)
            A = rfun(Th) * np.exp(1j * Th)
            r = abs(A)
            gap = pi * a * r / Th
            tang = np.angle(1j * A) + rng.normal(0, 0.15) + (pi if rng.random() < 0.5 else 0)
            t = gap * abs(n) * rng.uniform(0.5, 3.0) if name != "4" else r * rng.uniform(0.05, 0.3)
            k0 = int(round((Th - np.angle(A)) / pi))
            ks = [k0 + n * b for b in base]
            tries += 1
            res = solve_labelled(np.array([A.real, A.imag, tang, t]), ks, rfun)
            if res is None:
                continue
            v, L = res
            Z = corners(v)
            if v[3] > max_ratio * abs(sum(Z) / 4):
                continue
            key = tuple(np.round(sorted((z.real, z.imag) for z in Z), 8).ravel())
            found.setdefault(key, (v, L, ks))
        out[name] = dict(seeds=tries, solutions=list(found.values()))
    return out


# ---------------------------------------------------------------- checks

def classify_found(sq, a, kappa, omega):
    """Pattern counts and hypothesis data for squares from the v380 random search."""
    from collections import Counter
    from general_comb_v380 import corners as c380
    counts = Counter()
    worst_flag = 0.0
    maxW = 0.0
    floor_min = np.inf
    aligned = []
    for v, Ths in sq:
        Z = c380(v)
        th, k, W = window(Z, Ths)
        counts[split_type(k)] += 1
        worst_flag = max(worst_flag, abs(flag_residual(Z)) / max(abs(z) for z in Z) ** 2)
        maxW = max(maxW, W)
        floor_min = min(floor_min, float(th.min() + pi * k.min()))
        if split_type(k) == "2-2 adjacent":
            d = slope_ratio(Z, th, k)
            lo, hi = q_bounds_measured(d, W)
            aligned.append(dict(Q=float(d["Q"]), lo=float(lo), hi=float(hi),
                                inside=bool(lo - 1e-12 <= d["Q"] <= hi + 1e-12)))
    return dict(counts=dict(counts), flag_residual_max=float(worst_flag), W_max=maxW, floor_min=floor_min,
                aligned=len(aligned), aligned_inside_bounds=sum(x["inside"] for x in aligned),
                aligned_one_minus_Q_max=max((abs(1 - x["Q"]) for x in aligned), default=None),
                aligned_band_max=max((x["hi"] - x["lo"] for x in aligned), default=None))


def run_checks(path=None, n_random=6000, n_labelled=600):
    import json
    from general_comb_v380 import radius, search
    a = 0.5
    res = {"note": "half-turn spiral r = Theta^(-a)(1 + (kappa/Theta) sin(omega Theta)), a = 1/2; plain spiral is "
                   "kappa = 0. Window W <= 1 throughout."}
    res["thresholds"] = {str(x): threshold(x) for x in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0)}
    res["side_limit_211"] = {str(x): side_limit_211(x) for x in (0.25, 0.5, 1.0, 2.0)}
    res["adjacent_margin"] = [dict(L=L, n=n, q_lower=ql, q_upper=qu)
                              for L in (6.78, 10.0, 30.0, 100.0) for n in (1, 2, 3)
                              for ql, qu in [exact_q_gap(L, n, a)]]
    fam = {}
    for kappa, omega in ((0.0, 1.0), (0.2, 2.0), (0.3, 1.0), (0.5, 3.0), (0.7, 2.0)):
        sq = search(kappa, omega, n_random=n_random, seed=0)
        fam[f"{kappa},{omega}"] = classify_found(sq, a, kappa, omega)
    res["random_search"] = fam
    lab = {}
    for kappa, omega in ((0.0, 1.0), (0.5, 3.0)):
        out = labelled_search(lambda T: radius(T, kappa, omega), a, n_seeds=n_labelled, seed=1)
        lab[f"{kappa},{omega}"] = {k: dict(seeds=v["seeds"], solutions=len(v["solutions"])) for k, v in out.items()}
    res["labelled_search"] = lab
    v382 = json.load(open(path.parent.parent / "v382" / "tilted_check.json")) if path else None
    if v382:
        from collections import Counter
        res["v382_patterns"] = dict(Counter(split_type((0, s["n1"], s["n1"] + s["n3"] + s["j"], s["n3"]))
                                            for s in v382["squares"]))
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    from pathlib import Path
    out = run_checks(Path(__file__).resolve().parents[1] / "results" / "v383" / "two_turn_check.json")
    import json
    print(json.dumps({k: v for k, v in out.items() if k != "adjacent_margin"}, indent=1))
