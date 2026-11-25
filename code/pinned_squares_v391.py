"""v391: squares with an odd number of corners on arm two cannot collapse, and big ones are pinned.

Labels and coordinates are v390's anchored ones: A on its arm at lift Theta_A, the side vector w in A's frame,
offsets o_B, o_C, o_D. Write n1 = o_B, n3 = o_D, j = o_C - o_B - o_D (A's neighbours and the level).

Exact identity (v386), C' = Theta_B + Theta_D - Theta_A:
    Theta_C - C' = pi j - Sigma* + Dtheta = (2a+1) zeta^(2a+1) xi^(-2a-2) delta_B delta_D,
with xi in the hull of Theta_A, Theta_B, Theta_D, C' and zeta between Theta_C and C'. For an odd label j is odd and
|Sigma*| <= Phi + osc, so the left side has size at least pi - Phi - osc - |Dtheta| > 0 for a twist of size Phi < pi.

Two layers, kept apart:
  * sampled diagnostics (floats on a grid): band_status, pinned_band. They locate a band and candidate walls but
    prove nothing between or beyond the samples;
  * certificates (mpmath interval arithmetic, outward rounding; any failed step means 'not certified'):
      - certified_excluded(T): no square of the label has A's lift in the interval T;
      - certified_tail(T): no square of the label has A's lift >= T (monotonicity in the lift, proved in the note);
      - certify_region: both lift faces, the outer side face |w| = Q_MAX over the whole lift interval, and the inner
        side face |w| = q_in from a proven positive lower bound on the side. A certified region is one whose whole
        boundary no square of the label can touch, for every twist with |phi| <= Phi and |r phi'| <= L.
"""
import numpy as np
from math import asin, log, pi, sin, sqrt

from mpmath import iv

from label_degree_v390 import anchored_from_square, anchored_jacobian_sign, anchored_map
from smooth_cores_v386 import choose_A
from two_turn_squares_v383 import window

A_EXP = 0.5
Q_MAX = 0.3   # label regions only contain squares with side <= Q_MAX |A|


def arms_of(label):
    arm = label[0]
    return [arm] + [(arm + o) % 2 for o in label[1:]]


def odd_arm_count(label):
    """Odd number of corners on arm two (the twisted arm)."""
    return sum(arms_of(label)) % 2 == 1


def band_params(label):
    _, oB, oC, oD = label
    return oB, oD, oC - oB - oD


def compatible(n1, n3, j):
    """Both sides of the identity have definite signs (sign j and sign n1 n3) for small squares, so squares of the
    label need sign j = sign(n1 n3)."""
    return j != 0 and n1 * n3 != 0 and (j > 0) == (n1 * n3 > 0)


# ----------------------------------------------------------------------------------------------------------------
# Sampled diagnostics (floating point)
# ----------------------------------------------------------------------------------------------------------------

def side_map(q, Theta, n1, n3, Phi, L, a=A_EXP):
    """G(q): bound on side / |A| implied by side / |A| <= q (increasing in q, decreasing in Theta); None if out of range.

    AB and AD are perpendicular with length t, so t^2 = <AB, u>^2 + <AD, u>^2 exactly (u the radial unit vector at A).
    <AB, u> = (|B| - |A|) - |B| (1 - cos angle(A, B)), with ||B| - |A|| <= a |delta_B| min(Theta_A, Theta_B)^(-a-1),
    angle(A, B) <= asin q and |B| <= (1 + q) |A|. delta_B = pi n1 + arg(B / A) - (twist difference), the twist
    difference being <= Phi if one of A, B is on arm one and <= L |log(|B| / |A|)| <= -L log(1 - q) otherwise."""
    if q >= 0.5:
        return None
    e = asin(q) + max(Phi, -L * log(1 - q))
    comp = []
    for n in (n1, n3):
        dT = pi * abs(n) + e
        T_min = Theta - dT
        if T_min <= 0:
            return None
        comp.append(a * dT * Theta ** a / T_min ** (a + 1) + (1 + q) * (1 - sqrt(1 - q * q)))
    return sqrt(comp[0] ** 2 + comp[1] ** 2)


def side_bound(Theta, n1, n3, Phi, L, a=A_EXP, q_max=Q_MAX):
    """Float version of the side bound (diagnostic): least fixed point of G, slightly enlarged, with G(q) < q checked
    on adaptive cells up to q_max. See certified_side_bound for the certificate."""
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
    q_star, q_hat = q, q * (1 + 1e-4) + 1e-12
    g = side_map(q_hat, Theta, n1, n3, Phi, L, a)
    if g is None or g > q_hat:
        return None
    lo, h = q_hat, q_hat - q_star
    for _ in range(100000):
        if lo >= q_max:
            return q_hat
        hi = min(q_max, lo + h)
        g = side_map(hi, Theta, n1, n3, Phi, L, a)
        if g is not None and g < lo:
            lo, h = hi, 2 * h
        else:
            h /= 2
            if h < 1e-9 * q_hat:
                return None
    return None


def error_terms(q, Phi, L):
    """Upper bounds (floats) on: e, a corner's lift deviation from Theta_A + pi o_X; osc, the twist's variation over
    the square; d, |Dtheta|. q bounds side / |A|."""
    s2q = sqrt(2) * q
    e = asin(s2q) + max(Phi, -L * log(1 - s2q))          # worst corner is C, |C - A| = sqrt 2 t
    u = (q / sqrt(2)) / (1 - q / sqrt(2))                 # (t / sqrt 2) / |centre|
    osc = L * log((1 + u) / (1 - u))
    z = q * q / (1 - q) ** 2                              # Dtheta = arg(1 - bd / ((1 + b)(1 + d))), |b| = |d| = q
    d = asin(min(1.0, z))
    return e, osc, d


def band_status(Theta, n1, n3, j, Phi, L, a=A_EXP):
    """Sampled: 'excl' (identity rules out A at this lift), 'feas', or 'unknown' (outside the small-square bounds)."""
    q = side_bound(Theta, n1, n3, Phi, L, a)
    if q is None:
        return "unknown"
    e, osc, d = error_terms(q, Phi, L)
    xs = [0.0, pi * n1 - e, pi * n1 + e, pi * n3 - e, pi * n3 + e, pi * (n1 + n3) - 2 * e, pi * (n1 + n3) + 2 * e]
    zs = [pi * (n1 + n3) - 2 * e, pi * (n1 + n3) + 2 * e, pi * (n1 + n3 + j) - e, pi * (n1 + n3 + j) + e]
    x_lo, x_hi = Theta + min(xs), Theta + max(xs)
    z_lo, z_hi = Theta + min(zs), Theta + max(zs)
    if x_lo <= 0 or z_lo <= 0:
        return "unknown"
    slack = Phi + osc + d
    N_lo, N_hi = pi * abs(j) - slack, pi * abs(j) + slack
    if N_lo <= 0:
        return "unknown"
    D_lo = max(0.0, pi * abs(n1) - e) * max(0.0, pi * abs(n3) - e)
    D_hi = (pi * abs(n1) + e) * (pi * abs(n3) + e)
    RHS_lo = (2 * a + 1) * z_lo ** (2 * a + 1) / x_hi ** (2 * a + 2) * D_lo
    RHS_hi = (2 * a + 1) * z_hi ** (2 * a + 1) / x_lo ** (2 * a + 2) * D_hi
    return "excl" if (RHS_hi < N_lo or RHS_lo > N_hi) else "feas"


def pinned_band(n1, n3, j, Phi, L, a=A_EXP, grid=6000):
    """Sampled band: the run of non-excluded grid lifts containing the one nearest the resonance, and the candidate
    walls, the excluded samples just outside it. Nothing here is a certificate."""
    centre = (2 * a + 1) * pi * abs(n1 * n3) / max(abs(j), 1)
    T = np.geomspace(centre / 50, centre * 1000, grid)
    st = [band_status(x, n1, n3, j, Phi, L, a) for x in T]
    ok = np.array([s != "excl" for s in st])
    if not ok.any():
        return None
    c = int(np.argmin(np.abs(np.log(T / centre)) + 1e9 * (~ok)))
    i0 = c
    while i0 > 0 and ok[i0 - 1]:
        i0 -= 1
    i1 = c
    while i1 < len(T) - 1 and ok[i1 + 1]:
        i1 += 1
    lo_wall = i0 > 0 and st[i0 - 1] == "excl"
    hi_wall = i1 < len(T) - 1 and st[i1 + 1] == "excl"
    return dict(allowed_lo=float(T[i0]), allowed_hi=float(T[i1]), centre=float(centre),
                lo_wall_lift=float(T[i0 - 1]) if lo_wall else None,
                hi_wall_lift=float(T[i1 + 1]) if hi_wall else None,
                unknown_in_band=int(sum(1 for s in st[i0:i1 + 1] if s == "unknown")))


# ----------------------------------------------------------------------------------------------------------------
# Certificates (interval arithmetic)
# ----------------------------------------------------------------------------------------------------------------

def _I(x, y=None):
    return iv.mpf(x) if y is None else iv.mpf([x, y])


def _hull(cands):
    return iv.mpf([min(c.a for c in cands), max(c.b for c in cands)])


def _up(x):
    return float(x.b) if hasattr(x, "b") else float(x)


def _asin_up(x):
    """Interval whose upper end bounds asin over x (asin y <= y / sqrt(1 - y^2) on [0, 1))."""
    return x / iv.sqrt(1 - x * x)


def iv_side_map(q, T, n1, n3, Phi, L, a=A_EXP):
    """Enclosure of G(q, Theta) over q in the interval q and Theta in the interval T; None if out of range."""
    if q.b >= 0.5:
        return None
    PI = iv.pi
    lg = -L * iv.log(1 - q)
    e = _asin_up(q) + _I(max(Phi, float(lg.a)), max(Phi, float(lg.b)))
    curv = (1 + q) * (1 - iv.sqrt(1 - q * q))
    comp = []
    for n in (n1, n3):
        dT = PI * abs(n) + e
        T_min = T - dT
        if T_min.a <= 0:
            return None
        comp.append(a * dT * T ** a / T_min ** (a + 1) + curv)
    return iv.sqrt(comp[0] ** 2 + comp[1] ** 2)


def certified_side_bound(T, n1, n3, Phi, L, a=A_EXP, q_max=Q_MAX):
    """q_hat such that every square of the label with A's lift in the interval T and side <= q_max |A| has
    side / |A| <= q_hat, and none has side exactly q_max |A|; None if not certified.

    The true ratio q_t satisfies q_t <= G(q_t, Theta). Check G(q_hat) <= q_hat and G(cell) < lo(cell) on cells
    covering [q_hat, q_max]; then q_t in any cell would give q_t <= G(q_t) < lo <= q_t.

    G decreases in the lift (the radial term a dT Theta^a / (Theta - dT)^(a+1) has logarithmic derivative
    a / Theta - (a+1) / (Theta - dT) < 0), so the checks are made at the interval's lower end, which covers it all."""
    T = _I(T.a)
    qf = side_bound(float(T.a), n1, n3, Phi, L, a, q_max)
    if qf is None:
        return None
    q_hat = qf * (1 + 1e-3) + 1e-12
    g = iv_side_map(_I(q_hat), T, n1, n3, Phi, L, a)
    if g is None or g.b > q_hat:
        return None
    lo, h = q_hat, max(q_hat * 1e-3, 1e-9)
    for _ in range(20000):
        if lo >= q_max:
            return q_hat
        hi = min(q_max, lo + h)
        g = iv_side_map(_I(lo, hi), T, n1, n3, Phi, L, a)
        if g is not None and g.b < lo:
            lo, h = hi, 2 * h
        else:
            h /= 2
            if h < 1e-12:
                return None
    return None


def _iv_terms(q_hat, Phi, L):
    q = _I(q_hat)
    s2q = iv.sqrt(2) * q
    lg = -L * iv.log(1 - s2q)
    e = _asin_up(s2q) + _I(max(Phi, float(lg.a)), max(Phi, float(lg.b)))
    u = (q / iv.sqrt(2)) / (1 - q / iv.sqrt(2))
    osc = L * iv.log((1 + u) / (1 - u))
    z = q * q / (1 - q) ** 2
    if z.b >= 1:
        return None
    d = _asin_up(z)
    return e.b, osc.b, d.b


def _identity_sides(T, q_hat, n1, n3, j, Phi, L, a=A_EXP):
    """Enclosures of |left side| and of the right side of the identity over lifts in T (None if out of range)."""
    terms = _iv_terms(q_hat, Phi, L)
    if terms is None:
        return None
    e, osc, d = terms
    PI = iv.pi
    E = _I(-e, e)
    xi = _hull([T, T + PI * n1 + E, T + PI * n3 + E, T + PI * (n1 + n3) + 2 * E])
    zeta = _hull([T + PI * (n1 + n3) + 2 * E, T + PI * (n1 + n3 + j) + E])
    if xi.a <= 0 or zeta.a <= 0:
        return None
    slack = _I(Phi) + _I(osc) + _I(d)
    N = PI * abs(j) + _I(-1, 1) * slack
    if N.a <= 0:
        return None
    dB = PI * abs(n1) + E
    dD = PI * abs(n3) + E
    D = _I(max(0.0, float(dB.a)), float(dB.b)) * _I(max(0.0, float(dD.a)), float(dD.b))
    X = (2 * a + 1) * zeta ** (2 * a + 1) / xi ** (2 * a + 2)
    return N, X * D, (e, osc, d)


def certified_excluded(T, n1, n3, j, Phi, L, a=A_EXP):
    """True if no square of the label (side <= Q_MAX |A|) has A's lift in the interval T, for any admissible twist."""
    if not hasattr(T, "a"):
        T = _I(T)
    q_hat = certified_side_bound(T, n1, n3, Phi, L, a)
    if q_hat is None:
        return False
    sides = _identity_sides(T, q_hat, n1, n3, j, Phi, L, a)
    if sides is None:
        return False
    N, R, _ = sides
    return bool(R.b < N.a or R.a > N.b)


def certified_tail(T0, n1, n3, j, Phi, L, a=A_EXP):
    """True if no square of the label has A's lift >= T0 (any admissible twist).

    For Theta >= T0: G(q, Theta) <= G(q, T0) (the radial term a dT Theta^a / (Theta - dT)^(a+1) decreases in Theta),
    so q_hat(T0) bounds the side for all Theta >= T0 and the error terms are uniform; |left side| >= N_lo; and the
    right side is at most (2a+1) D_hi (Theta + P)^(2a+1) / (Theta + Q)^(2a+2), with P >= Q the largest zeta and
    smallest xi offsets, which decreases in Theta for Theta + Q > 0. So one check at T0 covers [T0, infinity)."""
    T = _I(T0)
    q_hat = certified_side_bound(T, n1, n3, Phi, L, a)
    if q_hat is None:
        return False
    terms = _iv_terms(q_hat, Phi, L)
    if terms is None:
        return False
    e, osc, d = terms
    PI = iv.pi
    E = _I(-e, e)
    xs = [_I(0), PI * n1 + E, PI * n3 + E, PI * (n1 + n3) + 2 * E]
    zs = [PI * (n1 + n3) + 2 * E, PI * (n1 + n3 + j) + E]
    Q = min(x.a for x in xs)
    P = max(z.b for z in zs)
    if T.a + Q <= 0 or P < Q:
        return False
    N_lo = PI * abs(j) - (_I(Phi) + _I(osc) + _I(d))
    D_hi = (PI * abs(n1) + _I(e)) * (PI * abs(n3) + _I(e))
    R_hi = (2 * a + 1) * D_hi * (T + _I(P)) ** (2 * a + 1) / (T + _I(Q)) ** (2 * a + 2)
    return bool(N_lo.a > 0 and R_hi.b < N_lo.a)


def side_floor(T, q_hat, n1, n3, Phi, L, a=A_EXP):
    """Certified lower bound on side / |A| over lifts in the interval T: t >= |<AB, u>| >= ||B| - |A|| - |B| (1 - cos),
    with ||B| - |A|| >= a |delta_B| max(Theta_A, Theta_B)^(-a-1) and |delta_B| >= pi |n1| - e_B, the true ratio being
    at most q_hat."""
    q = _I(q_hat)
    lg = -L * iv.log(1 - q)
    eB = _asin_up(q) + _I(max(Phi, float(lg.a)), max(Phi, float(lg.b)))
    curv = (1 + q) * (1 - iv.sqrt(1 - q * q))
    best = None
    for n in (n1, n3):
        dlo = iv.pi * abs(n) - eB
        c = iv.pi * abs(n) + eB
        if dlo.a <= 0 or T.a <= (a * c).b:
            continue
        # Theta^a / (Theta + c)^(a+1) decreases for Theta > a c, so its minimum over T is at the upper end
        Tb = _I(T.b)
        rel = a * dlo * Tb ** a / (Tb + c) ** (a + 1) - curv
        best = rel.a if best is None else max(best, rel.a)
    return float(best) if best is not None and best > 0 else None


def certify_region(n1, n3, j, Phi, L, T_lo, T_hi, a=A_EXP):
    """Certify every boundary component of Omega = {T_lo < Theta_A < T_hi, q_in < |w| < Q_MAX} for the label."""
    out = dict(lo_face=bool(certified_excluded(T_lo, n1, n3, j, Phi, L, a)),
               hi_face=bool(certified_excluded(T_hi, n1, n3, j, Phi, L, a)))
    TI = _I(T_lo, T_hi)
    q_hat = certified_side_bound(TI, n1, n3, Phi, L, a)
    out["outer_side_face"] = q_hat is not None
    out["q_hat"] = q_hat
    fl = side_floor(TI, q_hat, n1, n3, Phi, L, a) if q_hat is not None else None
    out["q_in"] = fl / 2 if fl else None
    out["inner_side_face"] = fl is not None
    out["certified"] = all(out[k] for k in ("lo_face", "hi_face", "outer_side_face", "inner_side_face"))
    return out


def certified_band(n1, n3, j, Phi, L, a=A_EXP, grid=6000):
    """Sampled band plus certificates: a certified region at the candidate walls, and the smallest certified tail
    start found (the collapse bound)."""
    b = pinned_band(n1, n3, j, Phi, L, a, grid)
    if b is None:
        return None
    b["region"] = None
    if b["lo_wall_lift"] and b["hi_wall_lift"]:
        # a sampled wall can sit right at the band's edge, with no margin for the interval check: slide each face
        # outward (while the sampled status stays 'excl') until it certifies
        faces = []
        for T0, f in ((b["lo_wall_lift"], 0.995), (b["hi_wall_lift"], 1.005)):
            T, chosen = T0, T0
            for _ in range(60):
                if band_status(T, n1, n3, j, Phi, L, a) != "excl":
                    break
                if certified_excluded(T, n1, n3, j, Phi, L, a):
                    chosen = T
                    break
                T *= f
            faces.append(chosen)
        b["region"] = certify_region(n1, n3, j, Phi, L, faces[0], faces[1], a)
        b["region"]["faces"] = [float(x) for x in faces]
    b["certified_closed"] = bool(b["region"] and b["region"]["certified"])
    b["tail_from"] = None
    T0 = b["hi_wall_lift"] or b["allowed_hi"]
    while T0 < 1000 * b["centre"]:
        if certified_tail(T0, n1, n3, j, Phi, L, a):
            b["tail_from"] = float(T0)
            break
        T0 *= 1.05
    return b


# ----------------------------------------------------------------------------------------------------------------
# Continuation in the twist (anchored map)
# ----------------------------------------------------------------------------------------------------------------

def solve_anchored(x, label, lam, phi, a=A_EXP, it=60):
    """Newton for a zero of v390's label map; x = (Theta_A, Re w, Im w)."""
    for _ in range(it):
        y = anchored_map(x, label, lam, phi, a)
        if np.max(np.abs(y)) < 1e-12:
            break
        q = abs(complex(x[1], x[2]))
        steps = (1e-8 * x[0], 1e-8 * q, 1e-8 * q)
        J = np.zeros((3, 3))
        for m in range(3):
            dx = np.zeros(3)
            dx[m] = steps[m]
            J[:, m] = (anchored_map(x + dx, label, lam, phi, a) - y) / steps[m]
        try:
            st = np.linalg.solve(J, -y)
        except np.linalg.LinAlgError:
            return None
        sc = max(1.0, abs(st[0]) / (0.02 * x[0]), np.hypot(st[1], st[2]) / (0.1 * q))
        x = x + st / sc
        if x[0] <= 0 or abs(complex(x[1], x[2])) >= 1 / sqrt(2):
            return None
    if np.max(np.abs(anchored_map(x, label, lam, phi, a))) > 1e-10:
        return None
    return x


def continue_in_twist(x0, label, phi, m, steps=40, a=A_EXP):
    """Follow a zero of the label map from twist 0 to twist m * phi."""
    lam, x = 0.0, x0.copy()
    path = [(0.0, x.copy())]
    h = 1.0 / steps
    f = lambda r: m * phi(r)
    while lam < 1.0:
        nl = 1.0 if lam + h > 1.0 - 1e-9 else lam + h
        res = solve_anchored(x.copy(), label, nl, f, a)
        if res is None or abs(abs(complex(res[1], res[2])) / abs(complex(x[1], x[2])) - 1) > 0.2:
            if h < 1.0 / (steps * 128):
                break
            h /= 2
            continue
        lam, x = nl, res
        path.append((lam, x.copy()))
        h = min(1.0 / steps, 2 * h)
    return path


def phi0(r):
    return 0.4 * sin(5 * log(1 / r))


def survival_starts(max_ratio=0.2, a=A_EXP):
    """v382's odd squares with side <= max_ratio |centre|, in anchored coordinates (A chosen as in v386)."""
    import json
    from pathlib import Path
    from tilted_squares_v382 import arm_point, three_corners
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
        x, label = anchored_from_square(Z, ks, choose_A(ks), a)
        if not odd_arm_count(label):
            continue
        starts.append((x, label, (s["n1"], s["n3"], s["j"])))
    return starts


def survival_row(x0, label, m, a=A_EXP, band=None):
    n1, n3, j = band_params(label)
    zero = lambda r: 0.0
    f = lambda r: m * phi0(r)
    path = continue_in_twist(x0, label, phi0, m, a=a)
    lifts = [x[0] for _, x in path]
    qs = [abs(complex(x[1], x[2])) for _, x in path]
    reached = path[-1][0] == 1.0
    xe = path[-1][1]
    q_diag = [side_bound(T, n1, n3, 0.4 * m, 2.0 * m, a) for T in lifts]
    row = dict(label=list(label), n1=n1, n3=n3, j=j, reached=bool(reached), lam_end=float(path[-1][0]),
               lift_min=float(min(lifts)), lift_max=float(max(lifts)), q_min=float(min(qs)), q_max=float(max(qs)),
               side_change=float(max(qs) / min(qs) - 1),
               side_over_bound=max((float(qv / qb) for qv, qb in zip(qs, q_diag) if qb), default=None),
               degree_start=anchored_jacobian_sign(x0, label, 0.0, zero, a)[0],
               degree_end=anchored_jacobian_sign(xe, label, 1.0, f, a)[0] if reached else None,
               residual_end=float(np.max(np.abs(anchored_map(xe, label, path[-1][0], f, a)))))
    if band is not None:
        row["tail_from"] = band["tail_from"]
        row["below_tail"] = bool(band["tail_from"] is not None and row["lift_max"] < band["tail_from"])
        reg = band["region"]
        row["certified_closed"] = band["certified_closed"]
        row["inside_certified_region"] = bool(
            band["certified_closed"] and reg["faces"][0] < row["lift_min"] and row["lift_max"] < reg["faces"][1]
            and reg["q_in"] < row["q_min"] and row["q_max"] < Q_MAX)
    return row


def survival_check(ms=(1, 2, 3), a=A_EXP):
    """v382's odd squares continued to twist m * 0.4 sin(5 log 1/r) (|phi| <= 0.4 m, |r phi'| <= 2 m)."""
    starts = survival_starts(a=a)
    out = {}
    for m in ms:
        rows = []
        for x0, label, _ in starts:
            band = certified_band(*band_params(label), 0.4 * m, 2.0 * m, a, grid=1500)
            rows.append(survival_row(x0, label, m, a, band))
        out[str(m)] = dict(squares=len(rows), reached=sum(r["reached"] for r in rows),
                           degree_kept=sum(1 for r in rows if r["reached"] and r["degree_end"] == r["degree_start"]
                                           and r["degree_start"] != 0),
                           certified_tails=sum(1 for r in rows if r["tail_from"] is not None),
                           below_tail=sum(r["below_tail"] for r in rows),
                           certified_closed=sum(r["certified_closed"] for r in rows),
                           inside_certified_region=sum(r["inside_certified_region"] for r in rows),
                           with_side_bound=sum(1 for r in rows if r["side_over_bound"] is not None),
                           max_side_over_bound=max(r["side_over_bound"] for r in rows
                                                   if r["side_over_bound"] is not None),
                           max_side_change=max(r["side_change"] for r in rows),
                           max_residual_end=max(r["residual_end"] for r in rows),
                           rows=rows)
    return out


def closure_table(Phis=(0.4, 0.8, 1.2), nmax=20, a=A_EXP, grid=800):
    """Sign-compatible odd labels with |j| = 1 and offsets 1..nmax, L = 5 Phi (as for Phi sin(5 log 1/r)): which have
    a certified region and a certified tail. Covers only the listed labels."""
    out = {}
    for P in Phis:
        closed, tails, total, sampled = {}, 0, 0, 0
        for n1 in range(1, nmax + 1):
            for n3 in range(1, nmax + 1):
                ok = True
                for s1, s3 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    j = 1 if s1 * s3 > 0 else -1
                    b = certified_band(s1 * n1, s3 * n3, j, P, 5 * P, a, grid)
                    total += 1
                    tails += bool(b and b["tail_from"] is not None)
                    sampled += bool(b and b["lo_wall_lift"] and b["hi_wall_lift"])
                    ok = ok and bool(b and b["certified_closed"])
                closed[(n1, n3)] = ok
        N = next((m for m in range(1, nmax + 1)
                  if all(closed[(x, y)] for x in range(m, nmax + 1) for y in range(m, nmax + 1))), None)
        out[str(P)] = dict(labels=total, certified_tails=tails, sampled_two_walls=sampled,
                           certified_closed_pairs=sum(closed.values()), pairs=len(closed), N_all_certified=N,
                           closed=[list(k) for k, v in closed.items() if v])
    return out


def large_label_checks(Phis=(0.4, 0.8, 1.2), a=A_EXP, grid=800):
    """Certified regions well beyond the table, including |j| = 3."""
    labs = [(25, 25, 1), (30, 40, 1), (50, 50, 1), (40, 120, 1), (100, 100, 1), (25, 200, 1), (60, 60, 3), (90, 150, 3)]
    out = {}
    for P in Phis:
        rows = {}
        for n1, n3, j in labs:
            ok = True
            for s1, s3 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                jj = j if s1 * s3 > 0 else -j
                b = certified_band(s1 * n1, s3 * n3, jj, P, 5 * P, a, grid)
                ok = ok and bool(b and b["certified_closed"] and b["tail_from"] is not None)
            rows[f"{n1},{n3},{j}"] = ok
        out[str(P)] = rows
    return out


def run_checks(path=None):
    import json
    res = dict(survival=survival_check(), closure=closure_table(), large_labels=large_label_checks())
    res["example_bands"] = {f"{n1},{n3},{j}": {str(P): certified_band(n1, n3, j, P, 5 * P, grid=1500)
                                               for P in (0.4, 0.8, 1.2)}
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
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "closed"} for k, v in r["closure"].items()}))
    print(json.dumps(r["large_labels"]))
