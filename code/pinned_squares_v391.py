"""v391: squares with an odd number of corners on arm two cannot collapse, and big ones are pinned.

Exact twisted identity (v386): Theta_C - C' = pi j - Sigma* + Dtheta = (2a+1) zeta^(2a+1) xi^(-2a-2) delta_B delta_D,
with C' = Theta_B + Theta_D - Theta_A, xi in the hull of Theta_A, Theta_B, Theta_D, C' and zeta between Theta_C and C'.
For a label with an odd number of corners on arm two, j is odd and |Sigma*| <= Phi + osc, so the left side has size at
least pi - Phi - osc - |Dtheta| > 0 for every twist of size Phi < pi. With a self-consistent bound on the side, this
gives, for each label and twist size:
  * an upper wall: every lift Theta_A above a finite hi is excluded, so the label's small squares cannot shrink into
    the wild point (no collapse), for every label checked;
  * a lower wall for big enough labels (both pass offsets >= 10-14 for Phi = 0.4-1.2): the band is then a region no
    square of the label can leave, so by v390 its degree is that of the symmetric spiral, +-1 under uniqueness (U0).
Lifts too small for the side bound are reported as 'unknown', never as excluded.
"""
import numpy as np
from math import asin, atanh, log, pi, sin, sqrt

from smooth_cores_v386 import choose_A, twisted_phase
from two_turn_squares_v383 import window

A_EXP = 0.5


def odd_arm_count(ks):
    return sum(1 for k in ks if k % 2) % 2 == 1


Q_MAX = 0.3   # label regions only contain squares with side <= Q_MAX |A|


def side_map(q, Theta, n1, n3, Phi, L, a=A_EXP):
    """G(q): the bound on side / |A| implied by assuming side / |A| <= q (increasing in q); None if out of range.

    AB and AD are perpendicular with length t, so t^2 = <AB, u>^2 + <AD, u>^2 exactly, u the radial unit vector at A.
    <AB, u> = (|B| - |A|) - |B| (1 - cos angle(A, B)), with ||B| - |A|| <= a dTheta Theta_min^(-a-1), the angle at
    most asin q and |B| <= (1 + q) |A|. B's lift step is pi n1 plus its polar-angle change (<= asin q) minus the twist
    difference (<= Phi if one of A, B is on arm one, <= L |log(|B| / |A|)| <= -L log(1 - q) if both are on arm two)."""
    if q >= 0.5:
        return None
    alpha = asin(q)                       # angle AB (or AD) subtends at the origin, |AB| = q |A|
    e = alpha + max(Phi, -L * log(1 - q))  # B's lift step deviates from pi n1 by its angle and the twist difference
    comp = []
    for n in (n1, n3):
        dT = pi * abs(n) + e
        T_min = Theta - dT
        if T_min <= 0:
            return None
        comp.append(a * dT * Theta ** a / T_min ** (a + 1) + (1 + q) * (1 - sqrt(1 - q * q)))
    return sqrt(comp[0] ** 2 + comp[1] ** 2)


def side_bound(Theta, n1, n3, Phi, L, a=A_EXP, q_max=Q_MAX):
    """Upper bound q on side / |A| for every square of the label with corner A at lift Theta and side <= q_max |A|.

    The true ratio q_t satisfies q_t <= G(q_t). Iterate G from 0 to its least fixed point, step slightly above it to
    q_hat with G(q_hat) <= q_hat, and check G(q) < q on [q_hat, q_max] (adaptive cells, using that G increases). Then
    q_t <= q_hat, and no square of the label has side exactly q_max |A|. Returns (q_hat, u) with
    u = (t / sqrt 2) / |centre| bounded through |centre| >= |A| (1 - q / sqrt 2); None if any step fails."""
    q = 0.0
    for _ in range(5000):
        g = side_map(q, Theta, n1, n3, Phi, L, a)
        if g is None or g >= q_max:
            return None
        if g - q < 1e-13:
            break
        q = g
    else:
        return None
    q_star = q
    q_hat = q_star * (1 + 1e-4) + 1e-12
    g = side_map(q_hat, Theta, n1, n3, Phi, L, a)
    if g is None or g > q_hat:
        return None
    # march up to q_max: on [lo, hi], G(q) <= G(hi) < lo <= q since G increases
    lo, h = q_hat, q_hat - q_star
    for _ in range(100000):
        if lo >= q_max:
            break
        hi = min(q_max, lo + h)
        g = side_map(hi, Theta, n1, n3, Phi, L, a)
        if g is not None and g < lo:
            lo, h = hi, 2 * h
        else:
            h /= 2
            if h < 1e-9 * q_hat:
                return None
    else:
        return None
    return q_hat, (q_hat / sqrt(2)) / (1 - q_hat / sqrt(2))


def band_status(Theta, n1, n3, j, Phi, L, a=A_EXP):
    """What the identity says about a square with this label (B, D are n1, n3 passes from A, level j) whose corner A
    has lift Theta: 'excl' (ruled out), 'feas' (allowed) or 'unknown' (outside the small-square estimates).

    Identity (v386), with C' = Theta_B + Theta_D - Theta_A:
        Theta_C - C' = pi j - Sigma* + Dtheta = (2a+1) zeta^(2a+1) xi^(-2a-2) delta_B delta_D,
    where xi lies in the hull of Theta_A, Theta_B, Theta_D, C' (second difference of Theta^(-2a)) and zeta between
    Theta_C and C' (first difference). Each corner's lift is Theta_A + pi o_X + eps_X with |eps_X| <= e.
    """
    sb = side_bound(Theta, n1, n3, Phi, L, a)
    if sb is None:
        return "unknown"
    q, u = sb
    W = 2 * asin(u)
    osc = 2 * L * atanh(u)
    ud = u / (1 - u)
    if ud >= 1:
        return "unknown"
    d = 2 * asin(ud * ud)
    e = W + max(Phi, osc)
    xs = [0.0, pi * n1 - e, pi * n1 + e, pi * n3 - e, pi * n3 + e, pi * (n1 + n3) - 2 * e, pi * (n1 + n3) + 2 * e]
    zs = [pi * (n1 + n3) - 2 * e, pi * (n1 + n3) + 2 * e, pi * (n1 + n3 + j) - e, pi * (n1 + n3 + j) + e]
    x_lo, x_hi = Theta + min(xs), Theta + max(xs)
    z_lo, z_hi = Theta + min(zs), Theta + max(zs)
    if x_lo <= 0 or z_lo <= 0:
        return "unknown"
    N_lo, N_hi = pi * abs(j) - Phi - osc - d, pi * abs(j) + Phi + osc + d
    if N_lo <= 0:
        return "unknown"
    D_lo = max(0.0, pi * abs(n1) - e) * max(0.0, pi * abs(n3) - e)
    D_hi = (pi * abs(n1) + e) * (pi * abs(n3) + e)
    X_lo = (2 * a + 1) * z_lo ** (2 * a + 1) / x_hi ** (2 * a + 2)
    X_hi = (2 * a + 1) * z_hi ** (2 * a + 1) / x_lo ** (2 * a + 2)
    ratio_hi = N_hi / D_lo if D_lo > 0 else float("inf")
    return "feas" if (N_lo / D_hi <= X_hi and ratio_hi >= X_lo) else "excl"


def band_feasible(Theta, n1, n3, j, Phi, L, a=A_EXP):
    return band_status(Theta, n1, n3, j, Phi, L, a) != "excl"


def pinned_band(n1, n3, j, Phi, L, a=A_EXP, grid=6000):
    """The band of lifts around the resonance that the identity allows, and whether each end is a wall.

    hi_wall: every lift above hi (out to 1000 x centre) is excluded, so squares cannot collapse into the point.
    lo_wall: the lift just below lo is excluded, so squares cannot drift outward past lo.
    closed: both walls, so the band is a region whose boundary no square of the label can touch.
    """
    centre = (2 * a + 1) * pi * abs(n1 * n3) / max(abs(j), 1)
    T = np.geomspace(centre / 50, centre * 1000, grid)
    st = [band_status(x, n1, n3, j, Phi, L, a) for x in T]
    ok = np.array([q != "excl" for q in st])
    if not ok.any():
        return None
    # the run of non-excluded lifts containing the one closest to the centre
    c = int(np.argmin(np.abs(np.log(T / centre)) + 1e9 * (~ok)))
    i0 = c
    while i0 > 0 and ok[i0 - 1]:
        i0 -= 1
    i1 = c
    while i1 < len(T) - 1 and ok[i1 + 1]:
        i1 += 1
    hi_wall = i1 < len(T) - 1 and all(q == "excl" for q in st[i1 + 1:])
    lo_wall = i0 > 0 and st[i0 - 1] == "excl"
    return dict(lo=float(T[i0]), hi=float(T[i1]), centre=float(centre), hi_wall=bool(hi_wall),
                lo_wall=bool(lo_wall), closed=bool(hi_wall and lo_wall),
                unknown_in_band=int(sum(1 for q in st[i0:i1 + 1] if q == "unknown")))


def corners_of(v):
    A = v[0] + 1j * v[1]
    e = np.exp(1j * v[2])
    return [A, A + v[3] * e, A + (1 + 1j) * v[3] * e, A + 1j * v[3] * e]


def label_residual(v, ks, lam, phi, ref, a=A_EXP):
    """Comb-phase residuals with the polar angle unwrapped continuously about the reference direction ref.

    (v383's window takes the angle of the centre from np.angle, which jumps by 2 pi when the centre crosses the
    negative real axis; inside a Newton step that shifts every residual by 2 pi.)"""
    Z = corners_of(v)
    th = [ref + np.angle(z * np.exp(-1j * ref)) for z in Z]
    return np.array([twisted_phase(z, t, k, lam, a, phi) for z, t, k in zip(Z, th, ks)])


def solve_label(v, ks, lam, phi, ref, a=A_EXP, it=60):
    for _ in range(it):
        y = label_residual(v, ks, lam, phi, ref, a)
        if np.max(np.abs(y)) < 1e-12:
            break
        t = abs(v[3])
        h = 1e-8 * t
        J = np.zeros((4, 4))
        for m in range(4):
            dv = np.zeros(4)
            dv[m] = h if m != 2 else h / t
            J[:, m] = (label_residual(v + dv, ks, lam, phi, ref, a) - y) / dv[m]
        try:
            stp = np.linalg.solve(J, -y)
        except np.linalg.LinAlgError:
            return None
        sc = max(1, np.max(np.abs(stp[[0, 1, 3]])) / (0.1 * t), abs(stp[2]) / 0.1)
        v = v + stp / sc
        if v[3] <= 0:
            return None
    if np.max(np.abs(label_residual(v, ks, lam, phi, ref, a))) > 1e-10:
        return None
    return v


def centre_angle(v, ref):
    Zc = sum(corners_of(v)) / 4
    return ref + np.angle(Zc * np.exp(-1j * ref))


def continue_in_twist(v0, ks, phi, m, steps=40, a=A_EXP):
    """Follow a square from twist 0 to twist m * phi with its label fixed (angles unwrapped along the way)."""
    lam, v = 0.0, v0.copy()
    ref = float(np.angle(sum(corners_of(v)) / 4))
    path = [(0.0, v.copy())]
    h = 1.0 / steps
    f = lambda r: m * phi(r)
    while lam < 1.0 - 1e-12:
        nl = 1.0 if lam + h > 1.0 - 1e-9 else lam + h
        res = solve_label(v.copy(), ks, nl, f, ref, a)
        if res is None or abs(res[3] - v[3]) > 0.2 * v[3]:
            if h < 1.0 / (steps * 128):
                break
            h /= 2
            continue
        lam, v = nl, res
        ref = centre_angle(v, ref)
        path.append((lam, v.copy()))
        h = min(1.0 / steps, 2 * h)
    return path


def side_estimate(Theta, n1, n3, Phi, L, a=A_EXP):
    sb = side_bound(Theta, n1, n3, Phi, L, a)
    return float("nan") if sb is None else sb[0]


def local_degree(v, ks, lam, phi, a=A_EXP, h=1e-7):
    """Sign of det of the label map's Jacobian at a zero (angles unwrapped about the square's own centre)."""
    ref = float(np.angle(sum(corners_of(v)) / 4))
    y = label_residual(v, ks, lam, phi, ref, a)
    J = np.zeros((4, 4))
    t = abs(v[3])
    for m in range(4):
        dv = np.zeros(4)
        dv[m] = h * t if m != 2 else h
        J[:, m] = (label_residual(v + dv, ks, lam, phi, ref, a) - y) / dv[m]
    d = np.linalg.det(J)
    return int(np.sign(d)) if abs(d) > 1e-9 * np.prod(np.linalg.norm(J, axis=0)) else 0


def lift_of_A(v, ks, a=A_EXP):
    Z = corners_of(v)
    i = choose_A(list(ks))
    return abs(Z[i]) ** (-1 / a), i


def survival_check(ms=(1, 2, 3), max_ratio=0.2, a=A_EXP):
    """v382's odd-arm squares (small relative to their radius), continued to twist m * 0.4 sin(5 log 1/r)."""
    import json
    from pathlib import Path
    from tilted_squares_v382 import arm_point, three_corners
    phi0 = lambda r: 0.4 * sin(5 * log(1 / r))
    d382 = json.loads((Path(__file__).resolve().parents[1] / "results" / "v382" / "tilted_check.json").read_text())
    starts, seen = [], set()
    for s in d382["squares"]:
        sol = three_corners(s["Theta_A"], s["n1"], s["n3"])
        if sol is None:
            continue
        th, t = sol
        A = arm_point(s["Theta_A"])
        Z = [A, A + t * np.exp(1j * th), A + (1 + 1j) * t * np.exp(1j * th), A + 1j * t * np.exp(1j * th)]
        if t > max_ratio * abs(sum(Z) / 4):
            continue
        key = tuple(np.round(sorted((z.real, z.imag) for z in Z), 8).ravel())
        if key in seen:
            continue
        seen.add(key)
        thw, _, _ = window(Z, [0, 0, 0, 0])
        ks = [int(round((abs(z) ** (-1 / a) - tt) / pi)) for z, tt in zip(Z, thw)]
        if not odd_arm_count(ks):
            continue
        starts.append((np.array([A.real, A.imag, th, t]), ks, (s["n1"], s["n3"], s["j"])))
    out = {}
    for m in ms:
        rows = []
        for v0, ks, lab in starts:
            i = choose_A(ks)
            o = [(i + q) % 4 for q in range(4)]
            n1, n3 = ks[o[1]] - ks[o[0]], ks[o[3]] - ks[o[0]]
            j = ks[o[0]] + ks[o[2]] - ks[o[1]] - ks[o[3]]
            band = pinned_band(n1, n3, j, 0.4 * m, 2.0 * m, a)
            path = continue_in_twist(v0, ks, phi0, m, a=a)
            lifts = [lift_of_A(v, ks, a)[0] for _, v in path]
            sides = [v[3] for _, v in path]
            ratios = [v[3] / abs(corners_of(v)[choose_A(ks)]) for _, v in path]
            est = [side_estimate(L_, n1, n3, 0.4 * m, 2.0 * m, a) for L_ in lifts]
            reached = path[-1][0] == 1.0
            sgn = local_degree(path[-1][1], ks, 1.0, lambda r: m * phi0(r), a) if reached else None
            below_hi = band is not None and max(lifts) <= band["hi"]
            above_lo = band is not None and min(lifts) >= band["lo"]
            rows.append(dict(label=list(lab), n1=n1, n3=n3, j=j, reached=bool(reached), lam_end=float(path[-1][0]),
                             lift_min=float(min(lifts)), lift_max=float(max(lifts)), band=band,
                             below_hi=bool(below_hi), above_lo=bool(above_lo),
                             in_band=bool(below_hi and above_lo),
                             side_ratio_over_estimate=float(max(r_ / e_ for r_, e_ in zip(ratios, est))),
                             side_change=float(max(sides) / min(sides) - 1), degree_start=local_degree(v0, ks, 0.0, phi0, a),
                             degree_end=sgn))
        out[str(m)] = dict(squares=len(rows), reached=sum(r["reached"] for r in rows),
                           below_hi=sum(r["below_hi"] for r in rows),
                           in_band=sum(r["in_band"] for r in rows),
                           closed_bands=sum(1 for r in rows if r["band"] and r["band"]["closed"]),
                           hi_walls=sum(1 for r in rows if r["band"] and r["band"]["hi_wall"]),
                           in_closed_band=sum(1 for r in rows if r["in_band"] and r["band"]["closed"]),
                           degree_kept=sum(1 for r in rows if r["reached"] and r["degree_end"] == r["degree_start"]),
                           max_side_ratio_over_estimate=max(r["side_ratio_over_estimate"] for r in rows),
                           max_side_change=max(r["side_change"] for r in rows),
                           rows=rows)
    return out


def closure_table(Phis=(0.4, 0.8, 1.2), nmax=20, a=A_EXP, grid=800):
    """For each twist size Phi (with L = 5 Phi, as for Phi sin(5 log 1/r)), which odd-arm labels (|j| = 1) have a
    closed band, and the smallest N such that every label with min(|n1|, |n3|) >= N (and both up to nmax) closes."""
    out = {}
    for P in Phis:
        closed, hi_walls, total = {}, 0, 0
        for n1 in range(1, nmax + 1):
            for n3 in range(1, nmax + 1):
                ok = True
                for s1, s3 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    for j in (1, -1):
                        b = pinned_band(s1 * n1, s3 * n3, j, P, 5 * P, a, grid)
                        total += 1
                        hi_walls += bool(b and b["hi_wall"])
                        ok = ok and bool(b and b["closed"])
                closed[(n1, n3)] = ok
        N = next((m for m in range(1, nmax + 1)
                  if all(closed[(x, y)] for x in range(m, nmax + 1) for y in range(m, nmax + 1))), None)
        out[str(P)] = dict(labels=total, hi_walls=hi_walls, closed_pairs=sum(closed.values()),
                           pairs=len(closed), N_all_close=N,
                           closed=[list(k) for k, v in closed.items() if v])
    return out


def large_label_checks(Phis=(0.4, 0.8, 1.2), a=A_EXP, grid=800):
    """Spot checks well beyond the closure table: do bigger labels (and |j| = 3) still close?"""
    labs = [(25, 25, 1), (30, 40, 1), (50, 50, 1), (40, 120, 1), (100, 100, 1), (25, 200, 1), (60, 60, 3), (90, 150, 3)]
    out = {}
    for P in Phis:
        rows = {}
        for n1, n3, j in labs:
            ok = True
            for s1, s3 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                for sj in (1, -1):
                    b = pinned_band(s1 * n1, s3 * n3, sj * j, P, 5 * P, a, grid)
                    ok = ok and bool(b and b["closed"])
            rows[f"{n1},{n3},{j}"] = ok
        out[str(P)] = rows
    return out


def run_checks(path=None):
    import json
    res = dict(survival=survival_check(), closure=closure_table(), large_labels=large_label_checks())
    res["example_bands"] = {f"{n1},{n3},{j}": {str(P): pinned_band(n1, n3, j, P, 5 * P) for P in (0.4, 0.8, 1.2, 2.0)}
                            for n1, n3, j in ((5, 6, 1), (-6, 7, -1), (8, 8, 1), (12, -13, -1))}
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import json
    from pathlib import Path
    r = run_checks(Path(__file__).resolve().parents[1] / "results" / "v391" / "pinned_check.json")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "rows"} for k, v in r["survival"].items()},
                     indent=1))
    print(json.dumps(r["example_bands"], indent=1))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "closed"} for k, v in r["closure"].items()}))
