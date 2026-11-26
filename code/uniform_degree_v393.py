"""v393: the starting degree is +-1 for every big label with comparable offsets, so those squares exist outright.

v392 certified that a bounded twist can't push a big odd label's square out of its region. Existence then needed the
symmetric (plain) spiral's degree on that region to be nonzero, the uniqueness assumption U0. Here that degree is
computed, uniformly, for every label with min offset m >= N and max offset M <= rho m.

Scaled coordinates (a = 1/2). With eps = 1/m, r = m / M in [1/rho, 1], Theta* = (2a+1) pi m M / |j| and n1 = sig1 k1,
n3 = sig3 k3 (k1, k3 the two offsets, t1 = k1 / M, t3 = k3 / M):
    Theta_A = s Theta*,   Re w = -eps sig1 a |j| t1 U / (2a+1),   Im w = eps sig3 a |j| t3 V / (2a+1),
so w = eps omega with omega = O(1). On the plain spiral the label map's residuals are
R_X = Theta_A (|1 + w m_X|^(-1/a) - 1) - arg(1 + w m_X) - pi o_X (v390, twist 0, independent of A's arm), and
    Psi = (R_B / (sig1 pi k1),  R_D / (sig3 pi k3),  (R_C - R_B - R_D) / (pi |j|))
       = (sU - 1 + O(eps / t1),  sV - 1 + O(eps / t3),  sig1 sig3 s U V - sign(j) + O(eps / r)).
Psi is the label map composed with an invertible diagonal change of variables and an invertible linear recombination,
so its degree on a box is +-1 times the label map's degree on the corresponding region. Every remainder is written
in closed form with eps factored out (h(w) = |1 + w|^(-2) - 1 + 2 Re w = eps^2 hhat(omega, eps), and so on), so it
can be enclosed for eps in [0, 1/N] and r in an interval, by interval arithmetic with outward rounding.

Two certificates per sign pattern and parameter box:
  * Krawczyk: K(X) = y - C Psi(y) + (I - C J(X))(X - y) inside the interior of the box X around (1, 1, 1) proves,
    for every parameter value, a unique zero in X with every Jacobian in J(X) nonsingular and det J of the sign of
    det C; so deg(Psi, X) = sign det C;
  * localization: on v392's region at twist 0 (s in [s_lo, s_hi], side <= (|j| / m) Q_hat |A|), every zero of Psi
    lies in X (branch and bound: each sub-box outside X has a component of Psi bounded away from 0).
By excision deg(label map, v392 region) = +-1, so by v390 and v392 the label's square exists for every twist with the
certified bounds.
"""
from mpmath import iv

A_EXP = 0.5


def I(x, y=None):
    return iv.mpf(x) if y is None else iv.mpf([x, y])


# ----------------------------------------------------------------------------------------------------------------
# Forward-mode automatic differentiation over intervals (value and gradient in (s, U, V))
# ----------------------------------------------------------------------------------------------------------------

class D:
    """Interval with an optional gradient in (s, U, V). g is None for constants and in value-only evaluations
    (psi_box), which skips all derivative work; arithmetic treats None as a zero gradient."""
    __slots__ = ("v", "g")

    def __init__(self, v, g=None):
        self.v = v if isinstance(v, type(I(0))) else I(v)
        self.g = g

    @staticmethod
    def var(x, k):
        g = [I(0), I(0), I(0)]
        g[k] = I(1)
        return D(x, tuple(g))

    def _c(self, o):
        return o if isinstance(o, D) else D(o)

    def __add__(self, o):
        o = self._c(o)
        if self.g is None and o.g is None:
            return D(self.v + o.v)
        if o.g is None:
            return D(self.v + o.v, self.g)
        if self.g is None:
            return D(self.v + o.v, o.g)
        return D(self.v + o.v, tuple(a + b for a, b in zip(self.g, o.g)))

    __radd__ = __add__

    def __neg__(self):
        return D(-self.v, None if self.g is None else tuple(-a for a in self.g))

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        v = self.v * o.v
        if self.g is None and o.g is None:
            return D(v)
        if o.g is None:
            return D(v, tuple(a * o.v for a in self.g))
        if self.g is None:
            return D(v, tuple(self.v * b for b in o.g))
        return D(v, tuple(a * o.v + self.v * b for a, b in zip(self.g, o.g)))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        q = self.v / o.v
        if self.g is None and o.g is None:
            return D(q)
        ga = self.g if self.g is not None else (I(0), I(0), I(0))
        gb = o.g if o.g is not None else (I(0), I(0), I(0))
        return D(q, tuple((a - q * b) / o.v for a, b in zip(ga, gb)))

    def __rtruediv__(self, o):
        return self._c(o) / self


def d_sq(x):
    """x^2 with the tight interval square (never negative), unlike x * x."""
    v = x.v ** 2
    return D(v) if x.g is None else D(v, tuple(2 * x.v * a for a in x.g))


def d_atan(x):
    v = iv.atan2(x.v, I(1))
    return D(v) if x.g is None else D(v, tuple(a / (1 + x.v * x.v) for a in x.g))


# ----------------------------------------------------------------------------------------------------------------
# The scaled map
# ----------------------------------------------------------------------------------------------------------------

def hhat(x, y, eps):
    """h(w) / eps^2 for w = eps (x + i y), a = 1/2: h = |1 + w|^(-2) - 1 + 2 Re w
    = (-|w|^2 + 4 (Re w)^2 + 2 Re w |w|^2) / |1 + w|^2."""
    n2 = x * x + y * y
    return (-n2 + 4 * x * x + 2 * eps * x * n2) / (1 + 2 * eps * x + eps * eps * n2)


def hhat1(x, y, eps):
    """(hhat - (3x^2 - y^2)) / eps, exactly: hhat = (3x^2 - y^2) + eps hhat1, with
    hhat1 = (2x |w|^2 - (4x^2 - |w|^2)(2x + eps |w|^2)) / (1 + 2 eps x + eps^2 |w|^2)."""
    n2 = x * x + y * y
    return (2 * x * n2 - (4 * x * x - n2) * (2 * x + eps * n2)) / (1 + 2 * eps * x + eps * eps * n2)


_HULL_CACHE = {}


def _hull_pow(rho, k, p):
    """Enclosure of (1 + t rho)^(-(p + k)) over t in [0, 1] (rho an interval with 1 + rho > 0). Memoized on the
    endpoints (K, G and H ask for the same hull_4)."""
    lo = min(1.0, float((1 + rho).a))
    hi = max(1.0, float((1 + rho).b))
    if lo <= 0:
        raise ValueError("rho too negative")
    key = (lo, hi, k, float(p.a) if hasattr(p, "a") else float(p), float(p.b) if hasattr(p, "b") else float(p))
    out = _HULL_CACHE.get(key)
    if out is None:
        if len(_HULL_CACHE) > 200000:
            _HULL_CACHE.clear()
        out = _HULL_CACHE[key] = I(lo, hi) ** (-(p + k))
    return out


def _G_value(rho, p):
    """Enclosure of G(rho) = (K(rho) - 1) / rho to second order: G = int_0^1 K'(u rho) du with
    K'(sigma) = -(p+2)/3 + sigma int_0^1 K''(v sigma) dv and K'' in (p+2)(p+3)/6 hull_4, so
    G in -(p+2)/3 + rho (p+2)(p+3)/12 hull_4 (width O(rho^2) instead of O(rho))."""
    return -(p + 2) / 3 + rho * (p + 2) * (p + 3) / 12 * _hull_pow(rho, 4, p)


def d_K(rho, p):
    """K(rho) = 2 int_0^1 (1 - t)(1 + t rho)^(-p-2) dt, so that (1 + rho)^(-p) = 1 - p rho + p(p+1)/2 rho^2 K(rho)
    exactly (Taylor with integral remainder). Value 1 + rho G(rho); K'(rho) = -2(p+2) int (1-t) t (1 + t rho)^(-p-3) dt
    lies in -(p+2)/3 hull_3 (positive weights, mean value)."""
    v = 1 + rho.v * _G_value(rho.v, p)
    if rho.g is None:
        return D(v)
    dk = -(p + 2) / 3 * _hull_pow(rho.v, 3, p)
    return D(v, tuple(dk * g for g in rho.g))


def d_G(rho, p):
    """G(rho) = (K(rho) - 1) / rho = int_0^1 K'(u rho) du (value: _G_value), with
    G'(rho) = int_0^1 u K''(u rho) du in (p+2)(p+3)/12 hull_4 (K'' = 2(p+2)(p+3) int (1-t) t^2 (1+t rho)^(-p-4) dt)."""
    v = _G_value(rho.v, p)
    if rho.g is None:
        return D(v)
    dg = (p + 2) * (p + 3) / 12 * _hull_pow(rho.v, 4, p)
    return D(v, tuple(dg * g for g in rho.g))


def d_H(rho, p):
    """H(rho) = (G(rho) + (p+2)/3) / rho = int int u K2(u v rho) dv du (K2, K3: second and third derivatives of K):
    value in (p+2)(p+3)/12 hull_4, and H'(rho) = int int u^2 v K3(u v rho) in -(p+2)(p+3)(p+4)/60 hull_5
    (K3 = -2(p+2)(p+3)(p+4) int (1-t) t^3 (1 + t rho)^(-p-5) dt)."""
    v = (p + 2) * (p + 3) / 12 * _hull_pow(rho.v, 4, p)
    if rho.g is None:
        return D(v)
    dh = -(p + 2) * (p + 3) * (p + 4) / 60 * _hull_pow(rho.v, 5, p)
    return D(v, tuple(dh * g for g in rho.g))


def t_terms(u, v, ec, p):
    """hhat1(u, v, ec) = p(p+1) (T0 + ec T1) exactly, with T0 = 2 u v^2 - (4p+2)/3 u^3 (the ec = 0 part) and
    T1 = 2u^2 (2u lin H + G |w|^2) + 2u |w|^2 lin G + |w|^4 K / 2, lin = 2u + ec |w|^2, rho = ec lin
    (from G = -(p+2)/3 + rho H and K = 1 + rho G)."""
    n2 = u * u + v * v
    lin = 2 * u + ec * n2
    rho = ec * lin
    if float((1 + rho.v).a) <= 0:
        # ec appears twice in ec * lin, and with ec an interval the enclosure of 1 + rho can reach 0 while the true
        # value doesn't; the same quantity written as |1 + w|^2 = (1 + ec u)^2 + (ec v)^2 is never negative and much
        # tighter. (Used only as a fallback, so evaluations where the first form works are unchanged.)
        rho = d_sq(1 + ec * u) + d_sq(ec * v) - 1
    K, G, H = d_K(rho, p), d_G(rho, p), d_H(rho, p)
    T0 = 2 * u * v * v - D((4 * p + 2) / 3) * u * u * u
    T1 = 2 * u * u * (2 * u * lin * H + G * n2) + 2 * u * n2 * lin * G + n2 * n2 * K / 2
    return n2, T0, T1


def hhat_general(x, y, eps, p):
    """For w = eps (x + i y) and any exponent (p = 1/(2a)): h = |1 + w|^(-1/a) - 1 + (1/a) Re w = eps^2 hhat with
    hhat = Q + eps hhat1, Q = -p |omega|^2 + 2p(p+1) x^2, rho = eps (2x + eps |omega|^2),
    hhat1 = 2p(p+1) x^2 G(rho)(2x + eps |omega|^2) + p(p+1)/2 (4x |omega|^2 + eps |omega|^4) K(rho)."""
    n2 = x * x + y * y
    lin = 2 * x + eps * n2
    rho = eps * lin
    K = d_K(rho, p)
    G = d_G(rho, p)
    P = D(p)                       # mpmath intervals refuse dual operands on the left
    P1 = D(p * (p + 1))
    Q = -1 * P * n2 + 2 * P1 * x * x
    h1 = 2 * P1 * x * x * G * lin + P1 * (4 * x * n2 + eps * n2 * n2) * K / 2
    return Q, h1


def psi(s, U, V, eps, r, sig1, sig3, b_is_min, j, a=A_EXP):
    """Psi at (s, U, V) (D numbers or intervals) for parameters eps, r (intervals), the signs and which corner (B or D)
    carries the smaller offset. a = 1/2 (a float) uses the closed forms below; any other a (a float or an interval,
    v396) uses the integral-remainder forms of psi_general, whose limit system is the same."""
    if not (isinstance(a, float) and a == 0.5):
        return psi_general(s, U, V, eps, r, sig1, sig3, b_is_min, j, a)
    wrap = lambda z: z if isinstance(z, D) else D(z)   # mpmath intervals refuse foreign operands
    s, U, V, eps, r = wrap(s), wrap(U), wrap(V), wrap(eps), wrap(r)
    t1, t3 = (r, D(1)) if b_is_min else (D(1), r)
    aj = abs(j)
    c = a * aj / (2 * a + 1)
    x = -sig1 * c * t1 * U            # omega for corner B: w = eps (x + i y)
    y = sig3 * c * t3 * V
    PI = D(iv.pi)
    # corner B: m_B = 1
    eqB = s * U - 1 + s * (2 * a + 1) * eps * hhat(x, y, eps) / (sig1 * aj * t1) \
        - r * eps * d_atan(eps * y / (1 + eps * x)) / (PI * sig1 * t1)
    # corner D: m_D = i, i omega = (-y, x)
    eqD = s * V - 1 + s * (2 * a + 1) * eps * hhat(-y, x, eps) / (sig3 * aj * t3) \
        - r * eps * d_atan(eps * x / (1 - eps * y)) / (PI * sig3 * t3)
    # corner C: m_C = 1 + i, (1 + i) omega = (x - y, x + y). Linear parts cancel in R_C - R_B - R_D, and the quadratic
    # parts cancel to (3(x-y)^2 - (x+y)^2) - (3x^2 - y^2) - (3y^2 - x^2) = -8xy; -8xy (2a+1) s / (j^2 r) = sig1 sig3 s U V
    # exactly (t1 t3 = r), so only the eps hhat1 terms are left to enclose
    hd1 = hhat1(x - y, x + y, eps) - hhat1(x, y, eps) - hhat1(-y, x, eps)
    # Delta arg = arg(1 + z), z = -i w^2 / ((1 + w)(1 + i w)) = eps^2 zhat
    w2r, w2i = x * x - y * y, 2 * x * y                      # omega^2
    nr, ni = w2i, -w2r                                        # -i omega^2
    pr = (1 + eps * x) * (1 - eps * y) - (eps * y) * (eps * x)   # (1 + eps omega)(1 + i eps omega), real part
    pim = (1 + eps * x) * (eps * x) + (eps * y) * (1 - eps * y)  # imaginary part
    den = pr * pr + pim * pim
    zr = (nr * pr + ni * pim) / den
    zi = (ni * pr - nr * pim) / den
    darg = d_atan(eps * eps * zi / (1 + eps * eps * zr))
    sj = 1 if j > 0 else -1
    eq3 = (sig1 * sig3) * s * U * V + s * (2 * a + 1) * eps * hd1 / (j * j * r) - darg / (PI * aj) - sj
    return eqB, eqD, eq3


def psi_general(s, U, V, eps, r, sig1, sig3, b_is_min, j, a):
    """Psi for any exponent a (float or interval), arranged so the exponent appears without interval dependency:
      * the scale c = a|j|/(2a+1) is factored out: w = ec (xh + i yh), ec = eps c, xh = -sig1 t1 U, yh = sig3 t3 V,
        and eps^2 hhat(x, y, eps) = ec^2 hhat(xh, yh, ec);
      * hhat = Q + ec hhat1 with a Q = -|w|^2/2 + (1+p) xh^2 and a hhat1 = (1+p)/2 (T0 + ec T1) (t_terms), using
        a p = 1/2 and a p(p+1) = (2a+1)/(4a);
      * in R_C - R_B - R_D the quadratic parts cancel to sig1 sig3 s U V, and the cubic parts T0 cancel to the
        polynomial S = 2(xh^3 - yh^3) + (4p+6) xh yh (xh - yh), with prefactor (a^2/(2a+1)) p(p+1) = 1/4 exactly:
        eq3 = sig1 sig3 sUV + s ec [S + ec (T1(C) - T1(B) - T1(D))] / (4r) - Delta arg / (pi |j|) - sign j."""
    wrap = lambda z: z if isinstance(z, D) else D(z)
    s, U, V, eps, r = wrap(s), wrap(U), wrap(V), wrap(eps), wrap(r)
    aI = a if isinstance(a, type(I(0))) else I(a)
    p = 1 / (2 * aI)
    t1, t3 = (r, D(1)) if b_is_min else (D(1), r)
    aj = abs(j)
    ec = eps * D(aI * aj / (2 * aI + 1))
    xh = -sig1 * t1 * U
    yh = sig3 * t3 * V
    PI = D(iv.pi)
    P1h = D((1 + p) / 2)
    n2B, T0B, T1B = t_terms(xh, yh, ec, p)
    n2D, T0D, T1D = t_terms(-1 * yh, xh, ec, p)
    n2C, T0C, T1C = t_terms(xh - yh, xh + yh, ec, p)
    aQB = -1 * n2B / 2 + D(1 + p) * xh * xh
    aQD = -1 * n2D / 2 + D(1 + p) * yh * yh
    eqB = s * U - 1 + s * ec * (aQB + ec * P1h * (T0B + ec * T1B)) / (sig1 * t1) \
        - r * eps * d_atan(ec * yh / (1 + ec * xh)) / (PI * sig1 * t1)
    eqD = s * V - 1 + s * ec * (aQD + ec * P1h * (T0D + ec * T1D)) / (sig3 * t3) \
        - r * eps * d_atan(ec * xh / (1 - ec * yh)) / (PI * sig3 * t3)
    S = 2 * (xh * xh * xh - yh * yh * yh) + D(4 * p + 6) * xh * yh * (xh - yh)
    w2r, w2i = xh * xh - yh * yh, 2 * xh * yh
    nr, ni = w2i, -1 * w2r
    pr = (1 + ec * xh) * (1 - ec * yh) - (ec * yh) * (ec * xh)
    pim = (1 + ec * xh) * (ec * xh) + (ec * yh) * (1 - ec * yh)
    den = pr * pr + pim * pim
    zr = (nr * pr + ni * pim) / den
    zi = (ni * pr - nr * pim) / den
    darg = d_atan(ec * ec * zi / (1 + ec * ec * zr))
    sj = 1 if j > 0 else -1
    eq3 = (sig1 * sig3) * s * U * V + s * ec * (S + ec * (T1C - T1B - T1D)) / (4 * r) - darg / (PI * aj) - sj
    return eqB, eqD, eq3


def _val(z):
    return z.v if isinstance(z, D) else z


def psi_box(box, par):
    """Interval enclosure of Psi over a box (three intervals)."""
    s, U, V = (D(b) for b in box)
    return [_val(z) for z in psi(s, U, V, *par)]


def jac_box(box, par):
    s, U, V = D.var(box[0], 0), D.var(box[1], 1), D.var(box[2], 2)
    out = psi(s, U, V, *par)
    return [list(z.g) for z in out]


# ----------------------------------------------------------------------------------------------------------------
# Krawczyk
# ----------------------------------------------------------------------------------------------------------------

def _inv3(Mx):
    import numpy as np
    return np.linalg.inv(np.array(Mx, dtype=float))


def krawczyk(box, par):
    """True (and sign of det C) if K(box) lies in the interior of box for every parameter value in par."""
    import numpy as np
    y = [I(float(b.mid)) for b in box]
    J = jac_box(box, par)
    Jm = [[float(e.mid) for e in row] for row in J]
    C = _inv3(Jm)
    Py = psi_box(y, par)
    K = []
    for i in range(3):
        acc = y[i]
        for k in range(3):
            acc = acc - I(C[i, k]) * Py[k]
        for k in range(3):
            m_ik = (I(1) if i == k else I(0))
            for l in range(3):
                m_ik = m_ik - I(C[i, l]) * J[l][k]
            acc = acc + m_ik * (box[k] - y[k])
        K.append(acc)
    inside = all(K[i].a > box[i].a and K[i].b < box[i].b for i in range(3))
    return inside, int(np.sign(np.linalg.det(C))), K


# ----------------------------------------------------------------------------------------------------------------
# Localization by branch and bound
# ----------------------------------------------------------------------------------------------------------------

def _inside(sub, box):
    return all(sub[i].a >= box[i].a and sub[i].b <= box[i].b for i in range(3))


def _disjoint(sub, box):
    return any(sub[i].b < box[i].a or sub[i].a > box[i].b for i in range(3))


def localize(region, box, par, max_boxes=200000):
    """True if every zero of Psi in region lies in box: each piece of region outside box has a component of Psi
    excluding 0 (bisection on the widest coordinate)."""
    stack = [list(region)]
    n = 0
    while stack:
        sub = stack.pop()
        n += 1
        if n > max_boxes:
            return False, n
        if _inside(sub, box):
            continue
        try:
            P = psi_box(sub, par)
        except ValueError:
            # the general-exponent remainder needs 1 + rho > 0 on the whole box; a box too large for that
            # enclosure cannot be excluded, so it is split like any other undecided box
            P = None
        if P is not None and any(p.a > 0 or p.b < 0 for p in P):
            continue
        ends = [(float(sub[i].a), float(sub[i].b)) for i in range(3)]
        w = [b - a for a, b in ends]
        k = max(range(3), key=lambda i: w[i])
        if w[k] < 1e-9:
            return False, n
        a_, b_ = ends[k]
        # cut at a face of the target box when it lies inside this piece (pieces then fall inside or outside it),
        # otherwise at the midpoint
        cut = (a_ + b_) / 2
        for f in (float(box[k].a), float(box[k].b)):
            if a_ < f < b_ and min(f - a_, b_ - f) > 0.05 * w[k]:
                cut = f
                break
        lo, hi = list(sub), list(sub)
        lo[k] = I(a_, cut)
        hi[k] = I(cut, b_)
        stack.extend([lo, hi])
    return True, n


# ----------------------------------------------------------------------------------------------------------------
# Certifying families
# ----------------------------------------------------------------------------------------------------------------

def patterns(j):
    """Sign patterns (sig1, sig3) compatible with level j (sign j = sig1 sig3), and which corner has the min offset."""
    out = []
    for s1 in (1, -1):
        for s3 in (1, -1):
            if s1 * s3 == (1 if j > 0 else -1):
                for b_is_min in (True, False):
                    out.append((s1, s3, b_is_min))
    return out


def float_zero(par_mid):
    """Approximate zero of Psi for point parameters (floats), by Newton from (1, 1, 1)."""
    import numpy as np
    z = np.array([1.0, 1.0, 1.0])
    for _ in range(30):
        box = [I(float(t)) for t in z]
        F = np.array([float(v.mid) for v in psi_box(box, par_mid)])
        J = np.array([[float(e.mid) for e in row] for row in jac_box(box, par_mid)])
        z = z - np.linalg.solve(J, F)
        if np.max(np.abs(F)) < 1e-14:
            break
    return z


B_STAR = 0.15   # coarse box (1, 1, 1) +- B_STAR that all zeros are first localized into


def certify_piece(eps, r, sig1, sig3, b_is_min, j, outer, depth=0, max_depth=10,
                  halves=(0.06, 0.04, 0.025, 0.015), log=None):
    """Krawczyk on a small box around this piece's approximate zero, and localization of every zero in the box
    outer into it; split the piece (eps or r) when needed. Returns the certified leaves, or None."""
    par = (eps, r, sig1, sig3, b_is_min, j)
    mid = (I(float(eps.mid)), I(float(r.mid)), sig1, sig3, b_is_min, j)
    z = float_zero(mid)
    reason = "krawczyk"
    for h in halves:
        box = [I(z[k] - h, z[k] + h) for k in range(3)]
        if not all(outer[k].a < box[k].a and box[k].b < outer[k].b for k in range(3)):
            continue
        kin, sgn, _ = krawczyk(box, par)
        if not kin:
            continue
        loc, n = localize(outer, box, par, max_boxes=20000)
        if loc:
            return [dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], half=h,
                         centre=[float(c) for c in z], det_sign=sgn, boxes=n)]
        reason = "localize"
    if depth >= max_depth:
        if log is not None:
            log.append(dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], reason=reason))
        return None
    de, dr = float(eps.b) - float(eps.a), float(r.b) - float(r.a)
    if de * 8 >= dr:
        m = (float(eps.a) + float(eps.b)) / 2
        parts = [(I(float(eps.a), m), r), (I(m, float(eps.b)), r)]
    else:
        m = (float(r.a) + float(r.b)) / 2
        parts = [(eps, I(float(r.a), m)), (eps, I(m, float(r.b)))]
    out = []
    for e2, r2 in parts:
        sub = certify_piece(e2, r2, sig1, sig3, b_is_min, j, outer, depth + 1, max_depth, halves, log)
        if sub is None:
            return None
        out.extend(sub)
    return out


def coarse_localize(eps, r, sig1, sig3, b_is_min, j, cert, depth=0, max_depth=6, log=None):
    """Every zero of Psi in v392's twist-0 region lies in the box (1 +- B_STAR)^3 (splitting parameters as needed).
    The region: s in [s_lo, s_hi]; side <= (|j|/m) Q_hat |A| gives |omega| <= |j| Q_hat, so |U| <= (2a+1) Q_hat / (a t1)
    and |V| <= (2a+1) Q_hat / (a t3)."""
    par = (eps, r, sig1, sig3, b_is_min, j)
    Q = cert["Q_hat"]
    t1 = float(r.a) if b_is_min else 1.0
    t3 = 1.0 if b_is_min else float(r.a)
    region = [I(cert["s_lo"], cert["s_hi"]), I(-4 * Q / t1, 4 * Q / t1), I(-4 * Q / t3, 4 * Q / t3)]
    star = [I(1 - B_STAR, 1 + B_STAR)] * 3
    ok, n = localize(region, star, par, max_boxes=60000)
    if ok:
        return [dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], boxes=n)]
    if depth >= max_depth:
        if log is not None:
            log.append(dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], reason="coarse"))
        return None
    de, dr = float(eps.b) - float(eps.a), float(r.b) - float(r.a)
    if de * 8 >= dr:
        m = (float(eps.a) + float(eps.b)) / 2
        parts = [(I(float(eps.a), m), r), (I(m, float(eps.b)), r)]
    else:
        m = (float(r.a) + float(r.b)) / 2
        parts = [(eps, I(float(r.a), m)), (eps, I(m, float(r.b)))]
    out = []
    for e2, r2 in parts:
        sub = coarse_localize(e2, r2, sig1, sig3, b_is_min, j, cert, depth + 1, max_depth, log)
        if sub is None:
            return None
        out.extend(sub)
    return out


def certify_degree(N, rho, j, cert):
    """For every label of level j with min offset m >= N and max offset M <= rho m, in every compatible sign pattern
    and with either corner carrying the min offset: every zero of Psi in v392's twist-0 region lies in (1 +- B_STAR)^3,
    and there it is unique and nondegenerate, with a fixed sign of det."""
    pats = []
    ok = True
    star = [I(1 - B_STAR, 1 + B_STAR)] * 3
    for sig1, sig3, b_is_min in patterns(j):
        log = []
        coarse = coarse_localize(I(0, 1 / N), I(1 / rho, 1), sig1, sig3, b_is_min, j, cert, log=log)
        leaves = certify_piece(I(0, 1 / N), I(1 / rho, 1), sig1, sig3, b_is_min, j, star, log=log) \
            if coarse is not None else None
        signs = sorted({lf["det_sign"] for lf in leaves}) if leaves else []
        pats.append(dict(sig1=sig1, sig3=sig3, b_is_min=b_is_min,
                         certified=coarse is not None and leaves is not None and len(signs) == 1,
                         det_sign=signs[0] if len(signs) == 1 else None,
                         coarse_pieces=len(coarse) if coarse else 0, fine_pieces=len(leaves) if leaves else 0,
                         coarse_boxes=sum(c["boxes"] for c in coarse) if coarse else 0,
                         fine_boxes=sum(lf["boxes"] for lf in leaves) if leaves else 0,
                         failures=log[:5]))
        ok = ok and pats[-1]["certified"]
    return dict(N=N, rho=rho, j=j, certified=bool(ok), patterns=pats)


# ----------------------------------------------------------------------------------------------------------------
# Running everything
# ----------------------------------------------------------------------------------------------------------------

def scaled_coords(n1, n3, j, x, a=A_EXP):
    """(s, U, V) and the parameters of an actual square x = (Theta_A, Re w, Im w) with offsets n1, n3."""
    from uniform_closure_v392 import resonance
    m, M = min(abs(n1), abs(n3)), max(abs(n1), abs(n3))
    eps, r = 1 / m, m / M
    b_is_min = abs(n1) == m
    s1, s3 = (1 if n1 > 0 else -1), (1 if n3 > 0 else -1)
    t1, t3 = (r, 1.0) if b_is_min else (1.0, r)
    c = a * abs(j) / (2 * a + 1)
    s = x[0] / resonance(n1, n3, j, a)
    U = x[1] / (-eps * s1 * c * t1)
    V = x[2] / (eps * s3 * c * t3)
    return (s, U, V), dict(eps=eps, r=r, sig1=s1, sig3=s3, b_is_min=b_is_min)


def validate(labels, fine, a=A_EXP):
    """Plain squares at concrete labels: their scaled coordinates lie in a certified leaf box for their parameters,
    and the label map's local degree has one sign per pattern."""
    from label_degree_v390 import anchored_jacobian_sign
    from uniform_closure_v392 import plain_square
    rows = []
    for n1, n3, j in labels:
        x, label = plain_square(n1, n3, j, a)
        (s, U, V), p = scaled_coords(n1, n3, j, x, a)
        key = f"{p['sig1']},{p['sig3']},{p['b_is_min']}"
        leaf = next((lf for lf in fine[key]["leaves"] if lf["eps"][0] <= p["eps"] <= lf["eps"][1]
                     and lf["r"][0] <= p["r"] <= lf["r"][1]), None)
        in_box = leaf is not None and all(abs(v - c) <= leaf["half"] for v, c in zip((s, U, V), leaf["centre"]))
        rows.append(dict(label=[n1, n3, j], pattern=key, scaled=[s, U, V], in_certified_box=bool(in_box),
                         local_degree=anchored_jacobian_sign(x, label, 0.0, lambda r: 0.0, a)[0]))
    return rows


VALIDATION_LABELS = [(11, 11, 1), (11, 21, 1), (12, 13, 1), (13, 25, 1), (16, 23, 1), (20, 35, 1), (30, 31, 1),
                     (40, 79, 1), (60, 60, 1), (-11, -15, 1), (-20, -35, 1), (-50, -51, 1), (12, -13, -1),
                     (-16, 23, -1), (25, -40, -1), (-33, 60, -1), (21, -11, -1), (-17, 12, -1), (14, 27, 1)]


def star_inside_region(cert, rho, j, N, a=A_EXP):
    """Excision needs the coarse box, mapped back to (Theta_A, w), inside v392's region: s in (s_lo, s_hi), and
    q_in < |w| < Q_MAX, where |w| = eps c sqrt(t1^2 U^2 + t3^2 V^2), c = a|j| / (2a+1), q_in = (|j| / m) f / 2."""
    lo, hi = I(1 - B_STAR), I(1 + B_STAR)
    c = I(a * abs(j)) / (2 * a + 1)
    rmin = I(1 / rho)
    w_lo = c * lo * iv.sqrt(rmin * rmin + 1)                 # |w| / eps, smallest
    w_hi = c * hi * iv.sqrt(I(2))                            # |w| / eps, largest
    ok_s = cert["s_lo"] < float(lo.a) and float(hi.b) < cert["s_hi"]
    ok_in = (I(abs(j)) * cert["side_floor"] / 2).b < w_lo.a  # q_in / eps < |w| / eps
    ok_out = (w_hi / N).b < 0.3
    return dict(s=bool(ok_s), inner=bool(ok_in), outer=bool(ok_out), ok=bool(ok_s and ok_in and ok_out))


def _certs():
    from pathlib import Path
    import json
    return json.loads((Path(__file__).resolve().parents[1] / "results" / "v392" / "uniform_check.json").read_text())["levels"]


def run_level(j, N_fine=11, rho=2.0):
    """Fine (plain spiral, all twist sizes at once) and coarse (per twist size) certificates for one level."""
    import time
    certs = _certs()
    star = [I(1 - B_STAR, 1 + B_STAR)] * 3
    res = dict(N_fine=N_fine, rho=rho, j=j, B_star=B_STAR, fine={}, coarse={})
    for sig1, sig3, b_is_min in patterns(j):
        key = f"{sig1},{sig3},{b_is_min}"
        t0 = time.time()
        log = []
        leaves = certify_piece(I(0, 1 / N_fine), I(1 / rho, 1), sig1, sig3, b_is_min, j, star, log=log)
        signs = sorted({lf["det_sign"] for lf in leaves}) if leaves else []
        res["fine"][key] = dict(certified=leaves is not None and len(signs) == 1,
                                det_sign=signs[0] if len(signs) == 1 else None, pieces=len(leaves or []),
                                boxes=sum(lf["boxes"] for lf in leaves or []), seconds=round(time.time() - t0, 1),
                                failures=log[:5], leaves=leaves or [])
        for P in ("0.4", "0.8", "1.2"):
            cert = certs[P]["1"]
            log = []
            t0 = time.time()
            co = coarse_localize(I(0, 1 / cert["c"]), I(1 / rho, 1), sig1, sig3, b_is_min, j, cert, log=log)
            res["coarse"].setdefault(P, {})[key] = dict(N=cert["c"], certified=co is not None,
                                                        pieces=len(co or []), boxes=sum(c["boxes"] for c in co or []),
                                                        seconds=round(time.time() - t0, 1), failures=log[:5])
    return res


def merge(levels, rho=2.0):
    """Combine the per-level runs, check the excision containment and validate against plain squares."""
    certs = _certs()
    res = dict(rho=rho, B_star=B_STAR, levels={str(lv["j"]): lv for lv in levels})
    fine = {}
    for lv in levels:
        fine.update(lv["fine"])
    res["star_inside"] = {P: {str(lv["j"]): star_inside_region(certs[P]["1"], rho, lv["j"], certs[P]["1"]["c"])
                              for lv in levels} for P in ("0.4", "0.8", "1.2")}
    res["validation"] = validate(VALIDATION_LABELS, fine)
    res["certified"] = {P: all(res["star_inside"][P][str(lv["j"])]["ok"] and
                               all(lv["fine"][k]["certified"] and lv["coarse"][P][k]["certified"] for k in lv["fine"])
                               for lv in levels) for P in ("0.4", "0.8", "1.2")}
    return res


if __name__ == "__main__":
    import json
    import sys
    from pathlib import Path
    out = Path(__file__).resolve().parents[1] / "results" / "v393"
    out.mkdir(parents=True, exist_ok=True)
    if sys.argv[1:2] == ["level"]:
        j = int(sys.argv[2])
        (out / f"level_{j}.json").write_text(json.dumps(run_level(j), indent=1))
    else:
        levels = [json.loads((out / f"level_{j}.json").read_text()) for j in (1, -1)]
        r = merge(levels)
        (out / "degree_check.json").write_text(json.dumps(r, indent=1))
        for lv in levels:
            for k, v in lv["fine"].items():
                print("fine", lv["j"], k, {kk: vv for kk, vv in v.items() if kk != "leaves"})
            for P, d in lv["coarse"].items():
                print("coarse", lv["j"], P, d)
        print("star", r["star_inside"])
        print("certified", r["certified"])
        for row in r["validation"]:
            print(row)
