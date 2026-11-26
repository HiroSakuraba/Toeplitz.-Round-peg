"""v392: big odd squares are pinned uniformly, not just label by label.

v391 certified a closed label region one label at a time, for offsets up to 20. Here the same estimates are written
in scaled variables, so that one interval-arithmetic check covers infinitely many labels at once.

For an odd label (j odd) with offsets n1, n3 and level j, let m = min(|n1|, |n3|), M = max(|n1|, |n3|),
    eps = 1 / m,   mu = 1 / M <= eps,   Theta* = (2a+1) pi m M / |j|   (the resonance),   s = Theta_A / Theta*.
Then
    pi M / Theta_A = eps beta0,   pi m / Theta_A = mu beta0,   pi |j| / Theta_A = eps mu j^2 / ((2a+1) s),
    e / Theta_A <= eps mu e gamma0,   with beta0 = |j| / ((2a+1) s), gamma0 = beta0 / pi,
and the identity of v386 reads
    pi |j| +- (Phi + osc + d)  =  (pi |j| / s) F,
    F in [(1-K)^(2a+1) (1+K)^(-2a-2) (1 - e eps/pi)(1 - e mu/pi),  (1+K)^(2a+1) (1-K)^(-2a-2) (1 + e eps/pi)(1 + e mu/pi)],
    K = (eps + mu) beta0 + eps mu (j^2 / ((2a+1) s) + 2 e gamma0)            (hull of the intermediate lifts).
Every error term increases with eps, with mu and with the side bound q, and decreases with s. So:

  * family check (eps, mu): one check certifies every odd label of that level with 1/min offset <= eps and
    1/max offset <= mu, i.e. min offset >= 1/eps and max offset >= 1/mu, whatever the signs;
  * uniform check N: eps = mu = 1/N covers every odd label of that level with both offsets >= N.

The region is Omega = {s_lo Theta* < Theta_A < s_hi Theta*, q_in < |w| < Q_MAX} with
s_lo = pi|j| / (pi|j| + Phi) (1 - delta) and s_hi = pi|j| / (pi|j| - Phi) (1 + delta), and q_in = eps f / 2 for a certified
scaled side floor f. All checks are in mpmath interval arithmetic with outward rounding; a failed step is "not
certified".
"""
from mpmath import iv

A_EXP = 0.5
Q_MAX = 0.3


def _I(x, y=None):
    return iv.mpf(x) if y is None else iv.mpf([x, y])


def _asin_up(x):
    return x / iv.sqrt(1 - x * x)          # >= asin on [0, 1)


def _curv(q):
    return (1 + q) * (1 - iv.sqrt(1 - q * q))


def _maxI(c, x):
    return _I(max(float(c), float(x.a)), max(float(c), float(x.b)))


def walls(Phi, j, delta):
    pj = iv.pi * abs(j)
    return pj / (pj + Phi) * (1 - _I(delta)), pj / (pj - Phi) * (1 + _I(delta))


def _slopes(L):
    """L is the twist's slope bound |r phi'|, or (v395) a pair (L, L1) adding a common angular deviation e1 of both
    arms with |r e1'| <= L1. A common deviation enters only through differences between corners."""
    return (float(L[0]), float(L[1])) if isinstance(L, (tuple, list)) else (float(L), 0.0)


def e_side(q, Phi, L):
    """Lift deviation of B or D (|B - A| = q |A|): its angle from A, the twist difference (<= Phi if one of them is on
    arm one, <= L |log(|B|/|A|)| otherwise) and the common deviation's difference (<= L1 |log(|B|/|A|)|)."""
    Lt, L1 = _slopes(L)
    lg = -iv.log(1 - q)
    return _asin_up(q) + _maxI(Phi, Lt * lg) + L1 * lg


def e_all(q, Phi, L):
    """Lift deviation of any corner (worst is C, |C - A| = sqrt 2 q |A|)."""
    Lt, L1 = _slopes(L)
    s2q = iv.sqrt(2) * q
    lg = -iv.log(1 - s2q)
    return _asin_up(s2q) + _maxI(Phi, Lt * lg) + L1 * lg


def side_terms(q, Phi, L):
    """Upper enclosures of e (all corners), osc and d at side bound q. osc bounds the twist sum beyond Phi:
    L log((1+u)/(1-u)) for the twist, plus 2 L1 log((1+u)/(1-u)) for the common deviation, whose signed sum over the
    corners e1(A) + e1(C) - e1(B) - e1(D) is at most two differences."""
    Lt, L1 = _slopes(L)
    u = (q / iv.sqrt(2)) / (1 - q / iv.sqrt(2))
    osc = (Lt + 2 * L1) * iv.log((1 + u) / (1 - u))
    z = q * q / (1 - q) ** 2
    return e_all(q, Phi, L), osc, _asin_up(z)


def side_G(q, eps, mu, s, Phi, L, j, a=A_EXP):
    """Enclosure of an upper bound on G(q) for every label with 1/min offset <= eps, 1/max offset <= mu at scaled
    lift >= s: t^2 = <AB,u>^2 + <AD,u>^2, each <= a x / (1 - x)^(a+1) + curv(q) with x = pi |n| / Theta + e / Theta."""
    b0 = _I(abs(j)) / ((2 * a + 1) * s)
    g0 = b0 / iv.pi
    emax = e_side(_I(Q_MAX), Phi, L)                        # uniform over q <= Q_MAX
    xM = eps * b0 + eps * mu * emax * g0
    xm = mu * b0 + eps * mu * emax * g0
    if xM.b >= 1:
        return None
    k = _curv(q)
    cM = a * xM / (1 - xM) ** (a + 1) + k
    cm = a * xm / (1 - xm) ** (a + 1) + k
    return iv.sqrt(cM ** 2 + cm ** 2)


def scaled_side_bound(eps, mu, s, Phi, L, j, a=A_EXP, uniform=False):
    """Q_hat with side / |A| <= eps' Q_hat for every covered label (eps' its own 1/min offset), or None.

    Fixed eps: march G(q) < lo on cells covering [eps Q_hat, Q_MAX] (so q_t <= G(q_t) forces q_t < eps Q_hat).
    uniform=True (all eps' <= eps, mu' <= eps'): the march at eps covers q in [eps' Q, Q_MAX] for Q <= 1/eps... via
    H(Q; eps') = G(eps' Q; eps') / eps' increasing in eps'; for q in [Q_MAX eps'/eps, Q_MAX] it adds the scalar
    condition sqrt2 C1 / (Q_MAX / eps) + sqrt2 curv(Q_MAX) / Q_MAX < 1, C1 = (a x / (1 - x)^(a+1)) / eps at eps."""
    e_, m_, s_ = _I(eps), _I(mu), _I(s)
    # least fixed point of G (floats), then enlarge
    q = 0.0
    for _ in range(4000):
        g = side_G(_I(q), e_, m_, s_, Phi, L, j, a)
        if g is None or g.b >= Q_MAX:
            return None
        if g.b - q < 1e-15:
            break
        q = float(g.b)
    q_hat = q * (1 + 1e-3) + 1e-14
    g = side_G(_I(q_hat), e_, m_, s_, Phi, L, j, a)
    if g is None or g.b > q_hat:
        return None
    lo, h = q_hat, q_hat * 1e-3
    for _ in range(50000):
        if lo >= Q_MAX:
            break
        hi = min(Q_MAX, lo + h)
        g = side_G(_I(lo, hi), e_, m_, s_, Phi, L, j, a)
        if g is not None and g.b < lo:
            lo, h = hi, 2 * h
        else:
            h /= 2
            if h < 1e-15:
                return None
    else:
        return None
    if uniform:
        b0 = _I(abs(j)) / ((2 * a + 1) * s_)
        emax = e_side(_I(Q_MAX), Phi, L)
        xM = e_ * b0 + e_ * m_ * emax * b0 / iv.pi
        C1 = a * xM / (1 - xM) ** (a + 1) / e_
        cond = iv.sqrt(2) * C1 * e_ / Q_MAX + iv.sqrt(2) * _curv(_I(Q_MAX)) / Q_MAX
        if cond.b >= 1:
            return None
    return q_hat / eps


def _hullI(cands):
    return iv.mpf([min(c.a for c in cands), max(c.b for c in cands)])


SIGNS = ((1, 1), (1, -1), (-1, 1), (-1, -1))   # signs of the larger and the smaller offset


def identity_check(eps_I, mu_I, s, q_bar, Phi, L, j, signs, a=A_EXP):
    """Enclosures of |left side| and of the right side (pi|j|/s) zeta^(2a+1) xi^(-2a-2) D~ at scaled lift s, for one
    sign pattern, every 1/min offset in eps_I and 1/max offset in mu_I, and side bound q <= q_bar.

    Relative to Theta_A the larger offset sits at sM eps beta0, the smaller at sm mu beta0, the level at
    j |j| eps mu / ((2a+1) s), each corner's deviation at most e eps mu beta0 / pi (sign(j) = sM sm for squares)."""
    sM, sm = signs
    s_ = _I(s)
    q = _I(q_bar)
    if q.b >= Q_MAX:
        return None
    e, osc, d = side_terms(q, Phi, L)
    E = _I(0, e.b)
    b0 = _I(abs(j)) / ((2 * a + 1) * s_)
    uM = sM * eps_I * b0
    um = sm * mu_I * b0
    uj = (sM * sm) * eps_I * mu_I * _I(j * j) / ((2 * a + 1) * s_)
    eT = E * eps_I * mu_I * b0 / iv.pi
    pm = _I(-1, 1)
    xi = _hullI([_I(1), 1 + uM + pm * eT, 1 + um + pm * eT, 1 + uM + um + 2 * pm * eT])
    zeta = _hullI([1 + uM + um + 2 * pm * eT, 1 + uM + um + uj + pm * eT])
    if xi.a <= 0 or zeta.a <= 0:
        return None
    Dt = (1 + pm * E * eps_I / iv.pi) * (1 + pm * E * mu_I / iv.pi)
    if Dt.a <= 0:
        return None
    pj = iv.pi * abs(j)
    R = pj / s_ * zeta ** (2 * a + 1) / xi ** (2 * a + 2) * Dt
    slack = _I(Phi) + osc + d
    return dict(N_lo=(pj - slack).a, N_hi=(pj + slack).b, R_lo=R.a, R_hi=R.b, e=e.b, osc=osc.b, d=d.b)


def side_floor_scaled(eps, mu, s, Q_hat, Phi, L, j, a=A_EXP):
    """Certified f with side / |A| >= eps f at every scaled lift <= s (for 1/min offset <= eps, 1/max <= mu):
    t >= ||D| - |A|| - |D| (1 - cos), using the larger offset M (pi M / Theta = eps beta0)."""
    e_, m_, s_ = _I(eps), _I(mu), _I(s)
    q = e_ * Q_hat
    eB = e_side(q, Phi, L)
    b0 = _I(abs(j)) / ((2 * a + 1) * s_)
    g0 = b0 / iv.pi
    num = a * (b0 - m_ * eB * g0)
    den = (1 + e_ * b0 + e_ * m_ * eB * g0) ** (a + 1)
    f = num / den - _curv(q) / e_
    return float(f.a) if f.a > 0 else None


DELTAS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


def _mode(eps, mu, uniform, single):
    if uniform:
        return _I(0, eps), _I(0, eps)
    return _I(eps), (_I(mu) if single else _I(0, mu))


def lower_face(Phi, L, j, eps, mu, dlo, a=A_EXP, uniform=False, single=False):
    """Side bound valid for s >= s_lo, and the lower face at s_lo, for all sign patterns. Returns (Q(s_lo), ok)."""
    pj = iv.pi * abs(j)
    s_lo = float((pj / (pj + Phi) * (1 - _I(dlo))).a)
    Q = scaled_side_bound(eps, mu, s_lo, Phi, L, j, a, uniform)
    if Q is None:
        return s_lo, None, False
    eI, mI = _mode(eps, mu, uniform, single)
    for sg in SIGNS:
        r = identity_check(eI, mI, s_lo, Q * eps, Phi, L, j, sg, a)
        if not (r and r["N_lo"] > 0 and r["R_lo"] > r["N_hi"]):
            return s_lo, Q, False
    return s_lo, Q, True


def upper_face(Phi, L, j, eps, mu, dhi, a=A_EXP, uniform=False, single=False):
    """The upper face at s_hi and the whole tail beyond it (all sign patterns), with the side bound at s_hi."""
    pj = iv.pi * abs(j)
    s_hi = float((pj / (pj - Phi) * (1 + _I(dhi))).b)
    Q = scaled_side_bound(eps, mu, s_hi, Phi, L, j, a, uniform)
    if Q is None:
        return s_hi, False
    eI, mI = _mode(eps, mu, uniform, single)
    for sg in SIGNS:
        r = identity_check(eI, mI, s_hi, Q * eps, Phi, L, j, sg, a)
        if not (r and r["N_lo"] > 0 and r["R_hi"] < r["N_lo"]):
            return s_hi, False
    return s_hi, True


def certify_family(Phi, L, j, eps, mu, deltas=DELTAS, a=A_EXP, uniform=False, single=False):
    """Certify the regions of every odd label of level j covered by (eps, mu): min offset exactly 1/eps and max offset
    >= 1/mu (family); both offsets >= 1/eps (uniform); or the single offset pair (1/eps, 1/mu) (single).

    The two lift faces are independent checks, so each picks its own wall width from deltas; the inner side face
    couples them (side bound at s_lo, floor evaluated at s_hi)."""
    if j % 2 == 0 or iv.pi * abs(j) <= Phi:
        return dict(certified=False, reason="level")
    los = [lower_face(Phi, L, j, eps, mu, d, a, uniform, single) + (d,) for d in deltas]
    his = [upper_face(Phi, L, j, eps, mu, d, a, uniform, single) + (d,) for d in deltas]
    best = None
    for s_lo, Q, ok_lo, dlo in los:
        if not ok_lo or Q * eps >= Q_MAX:
            continue
        for s_hi, ok_hi, dhi in his:
            if not ok_hi:
                continue
            f = side_floor_scaled(eps, mu if not uniform else eps, s_hi, Q, Phi, L, j, a)
            if f is not None and (best is None or s_hi - s_lo < best["s_hi"] - best["s_lo"]):
                best = dict(s_lo=s_lo, s_hi=s_hi, delta_lo=dlo, delta_hi=dhi, Q_hat=Q, side_floor=f)
    out = dict(Phi=Phi, L=L, j=j, eps=eps, mu=mu, uniform=uniform, single=single, certified=best is not None)
    if best:
        out.update(best)
    return out


def best_uniform(Phi, L, j, N_max=400, a=A_EXP):
    """Smallest N with a certified uniform check at eps = mu = 1/N (certification is monotone in N: every error term
    decreases as eps does, and each face and wall width is tried)."""
    lo, hi = 2, N_max
    if not certify_family(Phi, L, j, 1 / hi, 1 / hi, a=a, uniform=True)["certified"]:
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if certify_family(Phi, L, j, 1 / mid, 1 / mid, a=a, uniform=True)["certified"]:
            hi = mid
        else:
            lo = mid
    return hi


def best_tail_in_M(Phi, L, j, m, M_max=5000, a=A_EXP):
    """For a fixed min offset m: the smallest M0 with a certified family check (every max offset >= M0)."""
    if not certify_family(Phi, L, j, 1 / m, 1 / M_max, a=a)["certified"]:
        return None
    lo, hi = m - 1, M_max
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if certify_family(Phi, L, j, 1 / m, 1 / mid, a=a)["certified"]:
            hi = mid
        else:
            lo = mid
    return hi


def label_ok(Phi, L, j, m, M, a=A_EXP):
    """Certified region for the offset pair (m, M), m <= M, in all sign patterns: v392's single-label check, or else
    v391's per-label certificate (certified_band, every sign pattern with sign j = sign(n1 n3))."""
    if certify_family(Phi, L, j, 1 / m, 1 / M, a=a, single=True)["certified"]:
        return "v392"
    from pinned_squares_v391 import certified_band
    for sM, sm in SIGNS:
        jj = abs(j) * sM * sm
        b = certified_band(sM * M, sm * m, jj, Phi, L, a, grid=1500)
        if not (b and b["certified_closed"]):
            return None
    return "v391"


def coverage(Phi, L, j, m_from=3, M_max=5000, a=A_EXP):
    """Smallest m* such that every odd label of level |j| with min offset >= m* has a certified region: m >= N (the
    uniform check) or, for m* <= m < N, a family check for max offsets >= M0(m) plus single-label checks for
    m <= M < M0(m)."""
    N = best_uniform(Phi, L, j, a=a)
    if N is None:
        return None
    per_m = {}
    for m in range(N - 1, m_from - 1, -1):
        M0 = best_tail_in_M(Phi, L, j, m, M_max=M_max, a=a)
        if M0 is None:
            per_m[m] = dict(family_M0=None)
            break
        how = {M: label_ok(Phi, L, j, m, M, a) for M in range(m, M0)}
        fails = [M for M, h in how.items() if h is None]
        per_m[m] = dict(family_M0=M0, singles=M0 - m, by_v392=sum(1 for h in how.values() if h == "v392"),
                        by_v391=sum(1 for h in how.values() if h == "v391"), single_failures=fails)
        if fails:
            break
    thr = N
    for m in range(N - 1, m_from - 1, -1):
        r = per_m.get(m)
        if r is None or r.get("family_M0") is None or r["single_failures"]:
            break
        thr = m
    return dict(Phi=Phi, L=L, j=j, uniform_N=N, threshold=thr, per_m={str(k): v for k, v in sorted(per_m.items())})


# ----------------------------------------------------------------------------------------------------------------
# All levels at once: et = |j| / m, mt = |j| / M <= et, kap = 1 / |j|
# ----------------------------------------------------------------------------------------------------------------

def side_G_levels(q, et, mt, kap, s, Phi, L, a=A_EXP):
    """side_G in level-free variables: pi M / Theta = et / ((2a+1) s), e / Theta <= et mt kap e / ((2a+1) pi s)."""
    c0 = 1 / ((2 * a + 1) * s)
    emax = e_side(_I(Q_MAX), Phi, L)
    xM = et * c0 + et * mt * kap * emax * c0 / iv.pi
    xm = mt * c0 + et * mt * kap * emax * c0 / iv.pi
    if xM.b >= 1:
        return None
    k = _curv(q)
    return iv.sqrt((a * xM / (1 - xM) ** (a + 1) + k) ** 2 + (a * xm / (1 - xm) ** (a + 1) + k) ** 2)


def side_bound_levels(et, kap, s, Phi, L, a=A_EXP):
    """Q with side / |A| <= et' Q for every label with |j| / min offset = et' <= et, |j| / max offset <= et',
    1 / |j| <= kap, at scaled lift >= s (same uniform argument as scaled_side_bound, in et)."""
    E, K, S = _I(et), _I(kap), _I(s)
    q = 0.0
    for _ in range(4000):
        g = side_G_levels(_I(q), E, E, K, S, Phi, L, a)
        if g is None or g.b >= Q_MAX:
            return None
        if g.b - q < 1e-15:
            break
        q = float(g.b)
    q_hat = q * (1 + 1e-3) + 1e-14
    g = side_G_levels(_I(q_hat), E, E, K, S, Phi, L, a)
    if g is None or g.b > q_hat:
        return None
    lo, h = q_hat, q_hat * 1e-3
    for _ in range(50000):
        if lo >= Q_MAX:
            break
        hi = min(Q_MAX, lo + h)
        g = side_G_levels(_I(lo, hi), E, E, K, S, Phi, L, a)
        if g is not None and g.b < lo:
            lo, h = hi, 2 * h
        else:
            h /= 2
            if h < 1e-15:
                return None
    else:
        return None
    c0 = 1 / ((2 * a + 1) * S)
    emax = e_side(_I(Q_MAX), Phi, L)
    xM = E * c0 + E * E * K * emax * c0 / iv.pi
    C1 = a * xM / (1 - xM) ** (a + 1) / E
    if (iv.sqrt(2) * C1 * E / Q_MAX + iv.sqrt(2) * _curv(_I(Q_MAX)) / Q_MAX).b >= 1:
        return None
    return q_hat / et


def identity_levels(et_I, mt_I, kap_I, s, q_bar, Phi, L, signs, a=A_EXP):
    """The identity divided by pi |j|: 1 +- slack kap / pi  vs  zeta^(2a+1) xi^(-2a-2) D~ / s."""
    sM, sm = signs
    S = _I(s)
    e, osc, d = side_terms(_I(q_bar), Phi, L)
    Ee = _I(0, e.b)
    c0 = 1 / ((2 * a + 1) * S)
    uM = sM * et_I * c0
    um = sm * mt_I * c0
    uj = (sM * sm) * et_I * mt_I * c0
    eT = Ee * et_I * mt_I * kap_I * c0 / iv.pi
    pm = _I(-1, 1)
    xi = _hullI([_I(1), 1 + uM + pm * eT, 1 + um + pm * eT, 1 + uM + um + 2 * pm * eT])
    zeta = _hullI([1 + uM + um + 2 * pm * eT, 1 + uM + um + uj + pm * eT])
    Dt = (1 + pm * Ee * et_I * kap_I / iv.pi) * (1 + pm * Ee * mt_I * kap_I / iv.pi)
    if xi.a <= 0 or zeta.a <= 0 or Dt.a <= 0:
        return None
    R = zeta ** (2 * a + 1) / xi ** (2 * a + 2) * Dt / S
    slack = (_I(Phi) + osc + d) * kap_I / iv.pi
    return dict(N_lo=(1 - slack).a, N_hi=(1 + slack).b, R_lo=R.a, R_hi=R.b)


def certify_levels(Phi, L, c, j0, deltas=DELTAS, a=A_EXP, deltas_lo=None):
    """Certify the region of every odd label with |j| >= j0 and min offset >= c |j| (any max offset, any signs).
    deltas_lo (optional) restricts the lower wall's widths (used by v394 to keep s_lo >= 1/4)."""
    et, kap = 1 / c, 1 / j0
    ET, KAP = _I(0, et), _I(0, kap)
    pj0 = iv.pi * j0

    def s_lo_of(d):   # the walls use the worst level j0 (pi|j|/(pi|j| +- Phi) moves toward 1 as |j| grows)
        return float((pj0 / (pj0 + Phi) * (1 - _I(d))).a)

    def s_hi_of(d):
        return float((pj0 / (pj0 - Phi) * (1 + _I(d))).b)

    los = []
    for d in (deltas_lo if deltas_lo is not None else deltas):
        s_lo = s_lo_of(d)
        Q = side_bound_levels(et, kap, s_lo, Phi, L, a)
        ok = Q is not None and all(
            (lambda r: bool(r and r["N_lo"] > 0 and r["R_lo"] > r["N_hi"]))(
                identity_levels(ET, ET, KAP, s_lo, Q * et, Phi, L, sg, a)) for sg in SIGNS)
        los.append((d, s_lo, Q, ok))
    his = []
    for d in deltas:
        s_hi = s_hi_of(d)
        Q = side_bound_levels(et, kap, s_hi, Phi, L, a)
        ok = Q is not None and all(
            (lambda r: bool(r and r["N_lo"] > 0 and r["R_hi"] < r["N_lo"]))(
                identity_levels(ET, ET, KAP, s_hi, Q * et, Phi, L, sg, a)) for sg in SIGNS)
        his.append((d, s_hi, ok))
    best = None
    for dlo, s_lo, Q, ok_lo in los:
        if not ok_lo:
            continue
        for dhi, s_hi, ok_hi in his:
            if not ok_hi:
                continue
            # scaled floor: side / |A| >= et f, larger offset at et / ((2a+1) s)
            q = _I(Q * et)
            eB = e_side(q, Phi, L)
            c0 = 1 / ((2 * a + 1) * _I(s_hi))
            f = a * (c0 - _I(et) * _I(kap) * eB * c0 / iv.pi) / (1 + _I(et) * c0 + _I(et) * _I(et) * _I(kap) * eB * c0 / iv.pi) ** (a + 1) - _curv(q) / et
            if f.a > 0 and (best is None or s_hi - s_lo < best["s_hi"] - best["s_lo"]):
                best = dict(s_lo=s_lo, s_hi=s_hi, delta_lo=dlo, delta_hi=dhi, Q_hat=Q, side_floor=float(f.a))
    out = dict(Phi=Phi, L=L, c=c, j0=j0, certified=best is not None)
    if best:
        out.update(best)
    return out


def best_c(Phi, L, j0, c_max=400.0, a=A_EXP):
    """Smallest integer c with certify_levels(Phi, L, c, j0) certified (searched by bisection; each reported c is
    certified by its own check)."""
    if not certify_levels(Phi, L, c_max, j0, a=a)["certified"]:
        return None
    lo, hi = 1, int(c_max)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if certify_levels(Phi, L, mid, j0, a=a)["certified"]:
            hi = mid
        else:
            lo = mid
    return hi


# ----------------------------------------------------------------------------------------------------------------
# Using the certificates on concrete labels
# ----------------------------------------------------------------------------------------------------------------

def resonance(n1, n3, j, a=A_EXP):
    return (2 * a + 1) * 3.141592653589793 * abs(n1 * n3) / abs(j)


def region_for(cert, n1, n3, j, a=A_EXP):
    """The concrete region (lift interval, side interval) that an all-levels certificate gives for one label."""
    m = min(abs(n1), abs(n3))
    T = resonance(n1, n3, j, a)
    return dict(T_lo=cert["s_lo"] * T, T_hi=cert["s_hi"] * T, q_in=abs(j) / m * cert["side_floor"] / 2, q_max=Q_MAX,
                covered=m >= cert["c"] * abs(j) and abs(j) >= cert["j0"] and j % 2 == 1)


def plain_square(n1, n3, j, a=A_EXP):
    """A zero of v390's label map on the plain spiral near the resonance (A on arm one, offsets n1, n3, n1 + n3 + j),
    by Newton from the leading-order guess Theta_A = Theta*, w = -a pi (n1 - i n3) / Theta*."""
    import numpy as np
    from label_degree_v390 import anchored_map
    from pinned_squares_v391 import solve_anchored
    T = resonance(n1, n3, j, a)
    label = (0, n1, n1 + n3 + j, n3)
    zero = lambda r: 0.0
    for scale in (1.0, 0.9, 1.1, 0.8, 1.25):
        w = -a * np.pi * complex(n1, -n3) / (T * scale)
        x = solve_anchored(np.array([T * scale, w.real, w.imag]), label, 0.0, zero, a, it=200)
        if x is not None and np.max(np.abs(anchored_map(x, label, 0.0, zero, a))) < 1e-10:
            return x, label
    return None, label


def cross_check(cert, labels, a=A_EXP):
    """v391's independent per-label certifier at this certificate's face lifts, every sign pattern."""
    from pinned_squares_v391 import certified_excluded
    rows = []
    for n1, n3, j in labels:
        ok = True
        for sM, sm in SIGNS:
            jj = abs(j) * sM * sm
            r = region_for(cert, n1, n3, j, a)
            ok = ok and certified_excluded(r["T_lo"], sM * n1, sm * n3, jj, cert["Phi"], cert["L"], a) \
                and certified_excluded(r["T_hi"], sM * n1, sm * n3, jj, cert["Phi"], cert["L"], a)
        rows.append(dict(label=[n1, n3, j], v391_agrees=bool(ok)))
    return rows


def big_label_survival(certs, labels, ms=(1, 2, 3), a=A_EXP):
    """Plain squares at big labels, continued to twist m 0.4 sin(5 log 1/r): do they stay strictly inside the region
    certified for Phi = 0.4 m, L = 2 m, with their local degree?"""
    import numpy as np
    from label_degree_v390 import anchored_jacobian_sign
    from pinned_squares_v391 import continue_in_twist, phi0
    rows = []
    for n1, n3, j in labels:
        x0, label = plain_square(n1, n3, j, a)
        if x0 is None:
            rows.append(dict(label=[n1, n3, j], found=False))
            continue
        for m in ms:
            cert = certs[str(round(0.4 * m, 1))]
            reg = region_for(cert, n1, n3, j, a)
            path = continue_in_twist(x0, label, phi0, m, a=a)
            lifts = [x[0] for _, x in path]
            qs = [abs(complex(x[1], x[2])) for _, x in path]
            f = lambda r: m * phi0(r)
            reached = path[-1][0] == 1.0
            rows.append(dict(label=[n1, n3, j], m=m, found=True, covered=reg["covered"], reached=bool(reached),
                             s_min=min(lifts) / resonance(n1, n3, j, a), s_max=max(lifts) / resonance(n1, n3, j, a),
                             s_lo=cert["s_lo"], s_hi=cert["s_hi"], q_min=min(qs), q_max=max(qs), q_in=reg["q_in"],
                             inside=bool(reg["T_lo"] < min(lifts) and max(lifts) < reg["T_hi"]
                                         and reg["q_in"] < min(qs) and max(qs) < Q_MAX),
                             degree_start=anchored_jacobian_sign(x0, label, 0.0, lambda r: 0.0, a)[0],
                             degree_end=anchored_jacobian_sign(path[-1][1], label, 1.0, f, a)[0] if reached else None))
    return rows


BIG_LABELS = [(11, 11, 1), (12, 13, 1), (16, 16, 1), (16, 23, 1), (20, 35, 1), (30, 31, 1), (40, 90, 1),
              (60, 60, 1), (12, -13, -1), (-20, -35, 1), (-16, 23, -1), (-40, -90, 1), (35, 40, 3), (50, 80, 3),
              (90, 100, 5)]


def run_checks(path=None):
    import json
    import random
    res = dict(levels={}, per_level={}, cross_check={}, survival=None)
    certs = {}
    for P in (0.4, 0.8, 1.2):
        L = 5 * P
        res["levels"][str(P)] = {}
        for j0 in (1, 3, 9):
            c = best_c(P, L, j0)
            cert = certify_levels(P, L, c, j0)
            res["levels"][str(P)][str(j0)] = cert
            if j0 == 1:
                certs[str(P)] = cert
        res["per_level"][str(P)] = {str(j): coverage(P, L, j) for j in (1, 3, 5, 7)}
        rng = random.Random(392)
        c1 = certs[str(P)]["c"]
        labels = []
        while len(labels) < 30:
            j = rng.choice([1, 1, 1, 3, 5])
            m = rng.randint(c1 * j, c1 * j + 40)
            labels.append((m, m + rng.randint(0, 200), j))
        res["cross_check"][str(P)] = cross_check(certs[str(P)], labels)
    res["survival"] = big_label_survival(certs, BIG_LABELS)
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import json
    from pathlib import Path
    r = run_checks(Path(__file__).resolve().parents[1] / "results" / "v392" / "uniform_check.json")
    for P, v in r["levels"].items():
        print(P, {k: (c["c"], round(c["s_lo"], 3), round(c["s_hi"], 3), round(c["Q_hat"], 3), round(c["side_floor"], 4))
                  for k, c in v.items()})
        print("  per level", {j: (c["uniform_N"], c["threshold"]) for j, c in r["per_level"][P].items()})
        print("  cross-check", sum(x["v391_agrees"] for x in r["cross_check"][P]), "/", len(r["cross_check"][P]))
    for row in r["survival"]:
        print(row)
