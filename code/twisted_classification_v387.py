"""v387: closing the gaps in the square classification.

Part 1. v383 on twisted arms. Arm two is arm one turned by pi + phi(r), with |phi| <= Phi and |r phi'| <= L. A corner on
pass k has lift Theta = theta + pi k - phi*(rho), phi* = phi on odd passes (arm two), 0 on even ones. So:
  * corners on one pass share a lift window of width W1 = W + (oscillation of phi) <= W + 2 L atanh(u),
    u = sin(W/2) bounding t / (sqrt2 |Z|);
  * neighbouring passes are at least G = pi - Phi apart in lift at equal angle.
v383's flag arguments then go through with (pi, W) replaced by (G, W1). The aligned case (2-2 adjacent) also picks up
the twist's slope: on arm two dR/dtheta = -a Theta^(-a-1) / (1 + mu), |mu| <= a L / Theta.

Part 2. Completeness of v382's families on the plain spiral: every small square found by an independent search has
a corner A with both neighbours off its pass, and the three-corner solve of v382 at that corner's lift returns exactly
its neighbours (so the square is the one on the branch with its label).
"""
import numpy as np
from math import atanh, log, pi, sin, sqrt

from two_turn_squares_v383 import consts, split_type, window, x_lower


def default_profile():
    """phi = 0.4 sin(5 log 1/r): Phi = 0.4, |r phi'| <= 2, |r^2 phi''| <= 12."""
    return dict(Phi=0.4, L=2.0, L2=12.0)


def W1_of(W, L):
    u = sin(W / 2)
    return W + 2 * L * atanh(u)


# ---------------------------------------------------------------- unbalanced splits (flag)

def diagonal_ok(W, Phi, L):
    """2-2 diagonal: every inner lift exceeds every outer lift."""
    return W1_of(W, L) < pi - Phi


def three_one_outside_ok(W, Phi, L):
    return pi - Phi > 2 * W1_of(W, L)


def three_one_inside_threshold(W, Phi, L, a=0.5):
    """Smallest triple floor Theta_0 with (G - W1)(Theta_0/(Theta_0 + G))^(2a+1) > W1 (monotone in Theta_0)."""
    G, W1 = pi - Phi, W1_of(W, L)
    if G <= W1:
        return np.inf
    lo, hi = 1e-6, 1e9
    f = lambda T: (G - W1) * (T / (T + G)) ** (2 * a + 1) > W1
    if not f(hi):
        return np.inf
    for _ in range(200):
        m = sqrt(lo * hi)
        lo, hi = (lo, m) if f(m) else (m, hi)
    return hi


def ratio_limit_211(W, Phi, L, a=0.5):
    """2-1-1 adjacent is impossible when the corner radii differ by less than this factor."""
    G, W1 = pi - Phi, W1_of(W, L)
    return ((G - W1) / W1) ** (a / (2 * a + 1)) if G > W1 else 1.0


# ---------------------------------------------------------------- the aligned case

def aligned_bounds(Lf, n, W, Phi, L, a=0.5, odd="outer"):
    """Chord lower bound and twisted-spiral upper bound for Q = R_in'(eta1)/R_out'(eta2) at window floor Lf.

    odd: which pass is on arm two ('outer', 'inner', 'both', 'none'); for n odd exactly one of the first two."""
    g, T = consts(W)
    osc = 2 * L * atanh(sin(W / 2))
    po, pi_ = (Phi, 0.0) if odd == "outer" else (0.0, Phi) if odd == "inner" else (osc, osc) if odd == "both" else (0, 0)
    out_lo, out_hi = Lf - po, Lf + W + po
    in_lo, in_hi = Lf + pi * n - pi_, Lf + pi * n + W + pi_
    if odd == "both":            # the common phi shift cancels up to the oscillation
        out_lo, out_hi, in_lo, in_hi = Lf, Lf + W + osc, Lf + pi * n - osc, Lf + pi * n + W
    if out_lo <= a * L or in_lo <= out_hi:
        return None
    mu_out = a * L / out_lo if odd in ("outer", "both") else 0.0
    mu_in = a * L / in_lo if odd in ("inner", "both") else 0.0
    q_up = (out_hi / in_lo) ** (a + 1) * (1 + mu_out) / (1 - mu_in)
    p = a * g * out_hi ** a * out_lo ** (-a - 1) / (1 - mu_out)
    if p >= 1:
        return None
    s = sqrt(1 - p * p) * (out_lo / in_hi) ** a
    fac = (in_lo / out_hi) ** a
    q_lo = x_lower(min(s, 1.0), T, g) * fac
    return q_lo, q_up


def aligned_margin_leading(n, W, Phi, L, a=0.5, odd="outer"):
    """Leading-order margin (times the floor): (2a+1)(gap - W) - T a (gap + W) - (slope terms)."""
    _, T = consts(W)
    if odd in ("outer", "inner"):
        gap_lo, gap_hi, slope = pi * n - Phi, pi * n + Phi, a * L
    elif odd == "both":
        osc = 2 * L * atanh(sin(W / 2))
        gap_lo, gap_hi, slope = pi * n - osc, pi * n + osc, 2 * a * L
    else:
        gap_lo, gap_hi, slope = pi * n, pi * n, 0.0
    return (2 * a + 1) * (gap_lo - W) - T * a * (gap_hi + W) - slope


def aligned_threshold(W, Phi, L, a=0.5, n_max=400, grid=2500):
    """Largest floor at which the aligned bound fails, over n <= n_max and all parity cases (grid check)."""
    worst = 0.0
    floors = np.geomspace(0.5, 1e9, grid)
    for n in range(1, n_max + 1):
        cases = ("outer", "inner") if n % 2 else ("both", "none")
        for c in cases:
            if aligned_margin_leading(n, W, Phi, L, a, c) <= 0:
                return np.inf
            for Lf in floors[::-1]:
                b = aligned_bounds(Lf, n, W, Phi, L, a, c)
                if b is None or b[0] <= b[1]:
                    worst = max(worst, Lf)
                    break
    return float(worst)


# ---------------------------------------------------------------- one pass: convexity of the twisted arm

def arm_two_curvature_ok(Theta, phi, dphi, d2phi, a=0.5):
    """Sign of the curvature of arm two r = Theta^-a at polar angle theta = Theta + pi + phi(r) (positive = convex)."""
    r = Theta ** (-a)
    dr = -a * Theta ** (-a - 1)
    d2r = a * (a + 1) * Theta ** (-a - 2)
    dth = 1 + dphi(r) * dr
    d2th = d2phi(r) * dr * dr + dphi(r) * d2r
    # r as a function of theta: r_t = dr/dth, r_tt = (d2r dth - dr d2th)/dth^3
    rt = dr / dth
    rtt = (d2r * dth - dr * d2th) / dth ** 3
    return r * r + 2 * rt * rt - r * rtt > 0 and dth > 0


def default_phi_derivs():
    from math import cos
    phi = lambda r: 0.4 * sin(5 * log(1 / r))
    dphi = lambda r: -2.0 * cos(5 * log(1 / r)) / r
    d2phi = lambda r: (2.0 * cos(5 * log(1 / r)) - 10.0 * sin(5 * log(1 / r))) / (r * r)
    return phi, dphi, d2phi


# ---------------------------------------------------------------- searches on twisted arms

def twisted_corners(v):
    A = v[0] + 1j * v[1]
    e = np.exp(1j * v[2])
    return [A, A + v[3] * e, A + (1 + 1j) * v[3] * e, A + 1j * v[3] * e]


def nearest_pass(z, theta, lam, phi, a=0.5):
    """Pass index whose comb phase is closest to zero at the corner z (theta = unwrapped polar angle)."""
    Th = abs(z) ** (-1 / a)
    best = None
    k0 = int(np.floor((Th - theta) / pi)) - 1
    for k in range(k0, k0 + 4):
        val = Th - theta - pi * k + (lam * phi(abs(z)) if k % 2 else 0.0)
        if best is None or abs(val) < abs(best[0]):
            best = (val, k)
    return best[1]


def twisted_random_search(lam, phi, n_seeds=2000, seed=0, r_range=(0.08, 0.3), a=0.5, max_ratio=0.3):
    """Random seeds near the twisted arms, passes assigned from the seed, Newton (v386 solve_twisted)."""
    from smooth_cores_v386 import solve_twisted
    rng = np.random.default_rng(seed)
    found = {}
    for _ in range(n_seeds):
        r = rng.uniform(*r_range)
        Th = r ** (-1 / a)
        arm2 = rng.random() < 0.5
        ang = Th + (pi + lam * phi(r) if arm2 else 0.0)
        A = r * np.exp(1j * ang)
        gap = pi * a * r ** (1 + 1 / a)
        v = np.array([A.real, A.imag, rng.uniform(0, 2 * pi), gap * rng.uniform(0.7, 8)])
        Z = twisted_corners(v)
        th, _, _ = window(Z, [0, 0, 0, 0])
        ks = [nearest_pass(z, t, lam, phi, a) for z, t in zip(Z, th)]
        res = solve_twisted(v, ks, lam, phi, a)
        if res is None:
            continue
        w, Zs = res
        if w[3] > max_ratio * abs(sum(Zs) / 4):
            continue
        key = tuple(np.round(sorted((z.real, z.imag) for z in Zs), 8).ravel())
        found.setdefault(key, (Zs, ks))
    return list(found.values())


def twisted_labelled_search(lam, phi, patterns, n_seeds=200, seed=0, Th_range=(12.0, 60.0), a=0.5, gaps=(1, 2, 3)):
    """Prescribed two-pass labels on the twisted arms (v383's experiment with a twist)."""
    from smooth_cores_v386 import solve_twisted
    rng = np.random.default_rng(seed)
    out = {}
    for name, base in patterns.items():
        found = {}
        for _ in range(n_seeds):
            n = int(rng.choice(gaps)) * int(rng.choice([-1, 1]))
            Th = rng.uniform(*Th_range)
            r = Th ** (-a)
            arm2 = rng.random() < 0.5
            ang = Th + (pi + lam * phi(r) if arm2 else 0.0)
            A = r * np.exp(1j * ang)
            gap = pi * a * r / Th
            tang = ang + pi / 2 + rng.normal(0, 0.15) + (pi if rng.random() < 0.5 else 0)
            t = gap * abs(n) * rng.uniform(0.5, 3.0) if name != "4" else r * rng.uniform(0.05, 0.3)
            Z0 = twisted_corners(np.array([A.real, A.imag, tang, t]))
            th, _, _ = window(Z0, [0, 0, 0, 0])
            k0 = nearest_pass(A, th[0], lam, phi, a)
            ks = [k0 + n * b for b in base]
            res = solve_twisted(np.array([A.real, A.imag, tang, t]), ks, lam, phi, a)
            if res is None:
                continue
            w, Zs = res
            if w[3] > 0.3 * abs(sum(Zs) / 4):
                continue
            key = tuple(np.round(sorted((z.real, z.imag) for z in Zs), 8).ravel())
            found.setdefault(key, (Zs, ks))
        out[name] = dict(seeds=n_seeds, solutions=list(found.values()))
    return out


# ---------------------------------------------------------------- completeness of v382's families

def on_branch(Z, a=0.5, tol=1e-7):
    """Is this plain-spiral square the v382 branch square of its label? Returns (label, ok)."""
    from smooth_cores_v386 import choose_A
    from tilted_squares_v382 import three_corners
    lifts = [abs(z) ** (-1 / a) for z in Z]
    th, k, W = window(Z, lifts)
    i = choose_A(list(k))
    if i is None:
        return None, False
    o = [(i + m) % 4 for m in range(4)]
    # v382 places A on arm one; by the half-turn symmetry an arm-two corner is the image of an arm-one corner
    sgn = 1 if k[o[0]] % 2 == 0 else -1
    P = [sgn * Z[x] for x in o]
    kk = [k[x] for x in o]
    n1, n3 = kk[1] - kk[0], kk[3] - kk[0]
    j = kk[0] + kk[2] - kk[1] - kk[3]
    Th_A = abs(P[0]) ** (-1 / a)
    best = False
    for orient in (1, -1):                    # square order A,B,C,D or A,D,C,B
        Q = P if orient == 1 else [P[0], P[3], P[2], P[1]]
        m1, m3 = (n1, n3) if orient == 1 else (n3, n1)
        t0 = abs(Q[1] - Q[0])
        th0 = np.angle(Q[1] - Q[0])
        if np.angle((Q[3] - Q[0]) / (Q[1] - Q[0])) < 0:        # need D = A + i (B - A)
            continue
        sol = three_corners(Th_A, m1, m3, guess=(th0, t0), a=a)
        if sol is None:
            continue
        # uniqueness: also solve from the generic straight-comb guess and require the same answer
        sol2 = three_corners(Th_A, m1, m3, a=a)
        if sol2 is not None and abs(sol[1] - t0) < tol * t0 and abs(sol2[1] - sol[1]) < tol * t0 \
                and abs(np.angle(np.exp(1j * (sol2[0] - sol[0])))) < 1e-6:
            best = True
    return (n1, n3, j), best


PATTERNS = {
    "4": (0, 0, 0, 0),
    "3-1": (0, 0, 0, 1),
    "2-2 diagonal": (0, 1, 0, 1),
    "2-2 adjacent": (0, 0, 1, 1),
    "2-1-1 adjacent": (0, 0, 1, 2),
}


def classify_twisted(found, a=0.5):
    from collections import Counter
    c = Counter()
    W_max, tZ_max = 0.0, 0.0
    floor_min = np.inf
    for Z, ks in found:
        th, _, W = window(Z, [0, 0, 0, 0])
        k = list(ks)
        c[split_type(k)] += 1
        W_max = max(W_max, W)
        tZ_max = max(tZ_max, abs(Z[1] - Z[0]) / abs(sum(Z) / 4))
        floor_min = min(floor_min, float(min(th) + pi * min(k)))
    return dict(counts=dict(c), W_max=float(W_max), t_over_Z_max=float(tZ_max), floor_min=float(floor_min))


def run_checks(path=None):
    import json
    from general_comb_v380 import corners as c380, search
    P = default_profile()
    phi, dphi, d2phi = default_phi_derivs()
    res = {"profile": P}
    # criteria
    crit = {}
    for W in (0.2, 0.3, 0.4):
        crit[str(W)] = dict(W1=W1_of(W, P["L"]), diagonal=diagonal_ok(W, P["Phi"], P["L"]),
                            three_one_outside=three_one_outside_ok(W, P["Phi"], P["L"]),
                            three_one_inside_from=three_one_inside_threshold(W, P["Phi"], P["L"]),
                            ratio_211=ratio_limit_211(W, P["Phi"], P["L"]),
                            aligned_margin_n1=aligned_margin_leading(1, W, P["Phi"], P["L"]),
                            aligned_from=aligned_threshold(W, P["Phi"], P["L"], n_max=200, grid=1500))
    res["criteria"] = crit
    # convexity of arm two
    Ths = np.geomspace(3.0, 1e8, 20000)
    bad = [float(T) for T in Ths if not arm_two_curvature_ok(T, phi, dphi, d2phi)]
    res["arm_two_convex_from"] = max(bad) if bad else 3.0
    # where the aligned criterion fails: steep twists
    res["aligned_margin_vs_L"] = {str(L): aligned_margin_leading(1, 0.3, P["Phi"], L) for L in (0, 1, 2, 3, 4, 5, 6)}
    # searches on the twisted arms
    lab = twisted_labelled_search(1.0, phi, PATTERNS, n_seeds=300, seed=1)
    res["labelled_lam1"] = {k: len(v["solutions"]) for k, v in lab.items()}
    steep = lambda r: 3 * phi(r)
    lab3 = twisted_labelled_search(1.0, steep, {"2-2 adjacent": PATTERNS["2-2 adjacent"]}, n_seeds=300, seed=1)
    sols = lab3["2-2 adjacent"]["solutions"]
    res["labelled_steep_aligned"] = dict(solutions=len(sols), **(classify_twisted(sols) if sols else {}))
    rand = twisted_random_search(1.0, phi, n_seeds=3000, seed=2)
    res["random_lam1"] = classify_twisted(rand)
    # completeness of v382's families on the plain spiral
    sq = search(0.0, 1.0, n_random=6000, seed=7)
    ok = [on_branch(c380(v))[1] for v, _ in sq]
    res["v382_completeness"] = dict(squares=len(ok), on_branch=int(sum(ok)))
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import json
    from pathlib import Path
    print(json.dumps(run_checks(Path(__file__).resolve().parents[1] / "results" / "v387" / "classification_check.json"),
                     indent=1))
