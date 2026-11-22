"""Exact angular-window enclosures, ellipse calibration and spiral carrier exclusion.

v328  certified enclosures of the normalized window average D_m(L) from interval
      data for L on a mesh, using the half-angle variable t = tan(theta/2);
v329  exact ellipse calibration: every subsequential spectral limit of C^1-convergent
      approximants of an ellipse has D_m between two explicit integrals, which are
      computed here as rigorous rational enclosures;
v330  spiral fixtures with turn count independent of scale, the oscillation of the
      boundary area primitive across each gadget, and the resulting exclusion of all
      small rectangles as spectral carriers near the square angle.

Everything is exact rational arithmetic. pi and arctan enter only through rigorous
alternating-series enclosures. Nothing here computes Floer homology.
"""
from fractions import Fraction as Q
from math import comb
import json

# ---------------------------------------------------------------- enclosures

def arctan_enclosure(x, terms):
    """Rigorous [lo, hi] for arctan(x), 0 <= x <= 1, from alternating partial sums."""
    x = Q(x)
    if not 0 <= x <= 1 or terms < 1:
        raise ValueError("need 0 <= x <= 1 and at least one term")
    total, power, x2 = Q(0), x, x*x
    sums = []
    for k in range(terms+1):
        total += (-1)**k * power / (2*k+1)
        sums.append(total)
        power *= x2
    a, b = sums[-2], sums[-1]
    return (min(a, b), max(a, b))


def pi_enclosure(terms):
    """Machin: pi = 16 arctan(1/5) - 4 arctan(1/239)."""
    a = arctan_enclosure(Q(1, 5), terms)
    b = arctan_enclosure(Q(1, 239), terms)
    return (16*a[0]-4*b[1], 16*a[1]-4*b[0])


def linear_enclosure(coeffs, intervals):
    """Enclose sum c_i * x_i for rational c_i and intervals x_i."""
    lo = hi = Q(0)
    for c, (a, b) in zip(coeffs, intervals):
        lo += c*a if c >= 0 else c*b
        hi += c*b if c >= 0 else c*a
    return lo, hi


def divide(num, den):
    """Interval quotient for intervals of one strict sign."""
    if den[0] <= 0 <= den[1]:
        raise ZeroDivisionError("denominator interval contains zero")
    values = [a/b for a in num for b in den]
    return min(values), max(values)


def wallis_full(m):
    """Z_m / pi, where Z_m = integral_0^pi sin^(2m) = pi C(2m,m) / 4^m."""
    return Q(comb(2*m, m), 4**m)


# ------------------------------------------------------------------ v328

INF = None  # mesh endpoint t = infinity, i.e. theta = pi


def kernel(t, m):
    """sin^(2m) theta at t = tan(theta/2); sin theta = 2t/(1+t^2)."""
    if t is INF:
        return Q(0)
    return (2*t/(1+t*t))**(2*m)


def theta_enclosure(t, terms=40):
    """theta = 2 arctan t as an interval, for rational t >= 0 or infinity."""
    p = pi_enclosure(terms)
    if t is INF:
        return p
    if t <= 1:
        a = arctan_enclosure(t, terms)
        return (2*a[0], 2*a[1])
    a = arctan_enclosure(1/Q(t), terms)          # arctan t = pi/2 - arctan(1/t)
    return (p[0]-2*a[1], p[1]-2*a[0])


def window_enclosure(mesh, values, m, lipschitz, terms=40):
    """Enclose D_m(L) = Z_m^{-1} integral sin^(2m) dL.

    mesh: increasing t-values from 0 to INF and containing 1 (theta = pi/2);
    values: certified intervals for L at the mesh points;
    lipschitz: a valid Lipschitz constant for L in theta.
    Only monotonicity, the Lipschitz bound and the interval data are used.
    """
    if mesh[0] != 0 or mesh[-1] is not INF or 1 not in mesh:
        raise ValueError("mesh must run from 0 to infinity and contain t = 1")
    if any(b is not INF and a >= b for a, b in zip(mesh, mesh[1:])):
        raise ValueError("mesh must increase")
    if len(values) != len(mesh) or any(lo > hi for lo, hi in values):
        raise ValueError("one nonempty interval per mesh point")
    thetas = [theta_enclosure(t, terms) for t in mesh]
    low = high = Q(0)
    for i in range(len(mesh)-1):
        w = (kernel(mesh[i], m), kernel(mesh[i+1], m))   # kernel is monotone on each side of t = 1
        wmin, wmax = min(w), max(w)
        inc_lo = max(Q(0), values[i+1][0]-values[i][1])
        inc_hi = min(values[i+1][1]-values[i][0], lipschitz*(thetas[i+1][1]-thetas[i][0]))
        if inc_hi < inc_lo:
            raise ValueError("interval data contradict monotonicity or the Lipschitz bound")
        low += wmin*inc_lo
        high += wmax*inc_hi
    p = pi_enclosure(terms)
    z = wallis_full(m)
    return divide((low, high), (z*p[0], z*p[1]))


def circle_values(mesh, a2, terms=40):
    """Certified values of the circle model L(theta) = a^2 theta."""
    return [(a2*lo, a2*hi) for lo, hi in (theta_enclosure(t, terms) for t in mesh)]


def plateau_values(mesh, K, b, terms=40):
    """Certified values of the v319 plateau L(theta) = pi K f_b(theta/pi), slope K."""
    a = Q(b)/2
    out = []
    p = pi_enclosure(terms)
    for t in mesh:
        th = theta_enclosure(t, terms)
        # f_b(x) = min(x,a) + max(0, x-(1-a)) is nondecreasing; evaluate on the
        # outer bounds x in [th_lo/p_hi, th_hi/p_lo] and scale by pi K in [p_lo K, p_hi K].
        xs = (th[0]/p[1], min(Q(1), th[1]/p[0]))
        f = [min(x, a)+max(Q(0), x-(1-a)) for x in xs]
        out.append((p[0]*K*f[0], p[1]*K*f[1]))
    return out


def standard_mesh(cells):
    """Uniform in t on [0,1] and in 1/t on [1,inf), always containing t = 1."""
    left = [Q(i, cells) for i in range(cells+1)]
    right = [Q(cells, i) for i in range(cells-1, 0, -1)]
    return left + right + [INF]


# ------------------------------------------------------------------ v329

def _poly_divide(m, alpha, beta):
    """(1-c^2)^m = Qpoly(c) (alpha - beta c) + R; returns (Qpoly coefficients low->high, R)."""
    num = [Q(0)]*(2*m+1)
    for j in range(m+1):
        num[2*j] = Q((-1)**j*comb(m, j))
    # synthetic division by (alpha - beta c) = -beta (c - alpha/beta)
    root = Q(alpha)/beta
    deg = 2*m
    if deg == 0:
        return [], num[0]
    quotient = [Q(0)]*deg
    carry = Q(0)
    for k in range(deg, 0, -1):
        carry = num[k] + carry*root if k < deg else num[k]
        quotient[k-1] = carry
    remainder = num[0] + quotient[0]*root
    quotient = [c/(-beta) for c in quotient]
    return quotient, remainder


def _cos_power_quarter(k):
    """integral_0^{pi/2} cos^k as (rational part, pi-coefficient): (k-1)!!/k!!, times pi/2 if k is even."""
    num = den = 1
    for j in range(k-1, 0, -2):
        num *= j
    for j in range(k, 0, -2):
        den *= j
    return (Q(0), Q(num, den)/2) if k % 2 == 0 else (Q(num, den), Q(0))


def ellipse_quarter_integral(m, a, b, sign):
    """integral_0^{pi/2} sin^(2m) / (alpha - sign*beta cos) as rational coefficients of
    (1, pi, arctan(b/a)), with alpha = a^2+b^2, beta = a^2-b^2."""
    a, b = Q(a), Q(b)
    alpha, beta = a*a+b*b, sign*(a*a-b*b)
    quotient, remainder = _poly_divide(m, alpha, beta)
    r0 = rpi = Q(0)
    for k, c in enumerate(quotient):
        x, y = _cos_power_quarter(k)
        r0 += c*x
        rpi += c*y
    # integral_0^{pi/2} d theta / (alpha - beta cos) = arctan(sqrt((alpha+beta)/(alpha-beta)))/(ab)
    # = arctan(a/b)/(ab) for sign=+1 and arctan(b/a)/(ab) for sign=-1.
    if sign > 0:      # arctan(a/b) = pi/2 - arctan(b/a)
        rpi += remainder/(2*a*b)
        rT = -remainder/(a*b)
    else:
        rT = remainder/(a*b)
    return (r0, rpi, rT)


def ellipse_full_integral(m, a, b):
    """integral_0^pi sin^(2m) / (alpha - beta cos) = pi * rational."""
    a, b = Q(a), Q(b)
    alpha, beta = a*a+b*b, a*a-b*b
    quotient, remainder = _poly_divide(m, alpha, beta)
    total = Q(0)
    for k, c in enumerate(quotient):
        if k % 2 == 0:
            total += c*Q(comb(k, k//2), 2**k)
    return total + remainder/(2*a*b)     # integral_0^pi d theta/(alpha - beta cos) = pi/sqrt(alpha^2-beta^2) = pi/(2ab)


def ellipse_window(m, a, b, tol=Q(1, 10**12)):
    """Return exact D_m for each single branch and enclosures for the branch-free bounds.

    D_single: Z_m^{-1} integral sin^(2m) q_i^2 (identical for both branches by symmetry);
    D_low/D_high: Z_m^{-1} integral sin^(2m) min/max(q_1^2, q_2^2).
    """
    a, b = Q(a), Q(b)
    if not a >= b > 0:
        raise ValueError("need semiaxes a >= b > 0")
    q2sq = 2*a*a*b*b/(a*a+b*b)
    if a == b:
        one = (a*a, a*a)
        return {"D_single": a*a, "D_low": one, "D_high": one, "square_radius_sq": q2sq}
    z = wallis_full(m)
    single = 2*a*a*b*b*ellipse_full_integral(m, a, b)/z
    terms = 40
    while True:
        p = pi_enclosure(terms)
        T = arctan_enclosure(b/a, terms)
        out = {}
        for key, sign in (("D_low", -1), ("D_high", 1)):
            r0, rpi, rT = ellipse_quarter_integral(m, a, b, sign)
            num = linear_enclosure((4*a*a*b*b*r0, 4*a*a*b*b*rpi, 4*a*a*b*b*rT),
                                   ((Q(1), Q(1)), p, T))
            out[key] = divide(num, (z*p[0], z*p[1]))
        if all(hi-lo < tol for lo, hi in out.values()):
            break
        terms *= 2
    out.update(D_single=single, square_radius_sq=q2sq, series_terms=terms)
    return out


# ------------------------------------------------------------------ v330

E = ((1, 0), (0, 1), (-1, 0), (0, -1))


def gadget(N):
    """The v324 interleaved double spiral with N turns (unscaled)."""
    if type(N) is not int or N < 1:
        raise ValueError("positive integer number of turns required")
    m, d = 4*N, Q(1, 16*N)
    r = lambda j: Q(3, 4)-j*d
    u = [(r(j)*E[j % 4][0], r(j)*E[j % 4][1]) for j in range(m+1)]
    v = [((r(j)-2*d)*E[j % 4][0], (r(j)-2*d)*E[j % 4][1]) for j in range(m-1, -2, -1)]
    return ([(Q(1), Q(0))]+u+[(Q(0), Q(0))]+v+
            [(Q(0), Q(-1)), (Q(-1), Q(-1)), (Q(-1), Q(1)), (Q(0), Q(1)), (Q(1), Q(1))])


def placed_gadget(N, eps, height):
    """Gadget scaled by eps, attached to x = 0 between heights height and height+eps."""
    return [((x-1)*eps, height+y*eps) for x, y in gadget(N)]


def square_with_gadgets(gadgets):
    """Clockwise boundary of [0,4]^2 with gadgets (N, eps, height) on the left side,
    listed by increasing height. Starts at (0,0)."""
    pts = [(Q(0), Q(0))]
    for N, eps, h in gadgets:
        g = placed_gadget(N, Q(eps), Q(h))
        if g[0] == pts[-1]:
            g = g[1:]
        pts.extend(g)
    pts.extend([(Q(0), Q(4)), (Q(4), Q(4)), (Q(4), Q(0))])
    return pts


def spiral_fixture(k, N):
    """One gadget of N turns at scale 1/k (the v324 fixture when N = k^2)."""
    return square_with_gadgets([(N, Q(1, k), 0)])


def retuned_turns(k):
    """Nearest integer to (192/19) k^2, so that the tagged return integral tends to 8."""
    x = Q(192, 19)*k*k
    return int(x + Q(1, 2))


def tagged_return_integral(k, N):
    """The v324 joint return integral eps^2 (19N/24 - 1/(384N) - 1/2) at scale eps = 1/k."""
    return (Q(19*N, 24)-Q(1, 384*N)-Q(1, 2))/(k*k)


def primitive(points, origin=(0, 0)):
    """Cumulative integral of ((x-ox)dy-(y-oy)dx)/2 along the polyline."""
    ox, oy = map(Q, origin)
    F = [Q(0)]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        F.append(F[-1]+((x0-ox)*(y1-oy)-(y0-oy)*(x1-ox))/2)
    return F


def corner_oscillation(points):
    """Oscillation of the primitive (origin (0,0)) on the arc of the boundary that runs from
    the bottom edge through every gadget to the top of the left side.

    The primitive is constant on the lines x = 0 and y = 0, so this equals the oscillation over
    the attached gadgets themselves."""
    top = points.index((Q(0), Q(4)))
    arc = [(Q(4), Q(0))] + points[:top+1]
    F = primitive(arc)
    return max(F)-min(F)


def area(points):
    return abs(sum(x0*y1-x1*y0 for (x0, y0), (x1, y1) in zip(points, points[1:]+points[:1])))/2


def _hull(points):
    pts = sorted(set(points))
    cross = lambda o, a, b: (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1]+upper[:-1]


def radius_squared(points):
    """R^2 = diam^2/4, exact (the diameter is attained at hull vertices)."""
    h = _hull(points)
    return max((a[0]-b[0])**2+(a[1]-b[1])**2 for i, a in enumerate(h) for b in h[i+1:])/4


def carrier_exclusion(points, terms=40):
    """Compare the small-rectangle action bound with the spectral gap at pi/2.

    Every rectangle of half-diagonal < 1 has its contacts in one corner neighbourhood, where the
    action residue is at most 2*Osc (Osc = 0 at the three plain corners). At the square angle,
    A/2 <= ell(pi/2) <= pi R^2/2, so dist(ell(pi/2), A Z) >= G = min(A/2, A - pi R^2/2).
    On |theta - pi/2| <= h the spectral value moves by at most R^2 h. The returned window_half
    is a rational h with R^2 h < G - 2 Osc, or None when the comparison fails."""
    A, R2, osc = area(points), radius_squared(points), corner_oscillation(points)
    p = pi_enclosure(terms)
    gap = min(A/2, A-p[1]*R2/2)          # a rigorous lower bound for G
    slack = gap-2*osc
    window = slack/R2 if slack > 0 else None
    return {"vertices": len(points), "area": A, "radius_sq": R2, "oscillation": osc,
            "small_rectangle_residue_bound": 2*osc, "spectral_gap_lower": gap,
            "window_half_width": window}


def original_family_bounds(k, terms=40):
    """Analytic bounds for the v324 family (N = k^2, eps = 1/k), valid for every k >= 1.

    * Osc <= 19/24 + (19/4) eps^2: with origin at the gadget centre the primitive rises
      monotonically along the inbound spiral to eps^2 (19N/24 - 1/(384N)) <= 19/24, never drops
      below -5 eps^2/2 afterwards, and moving the origin to (0,0) costs at most sqrt(5) eps^2/2
      at each end (sqrt 5 < 9/4).
    * A >= 16 and R^2 <= 8 + 6 eps + 5 eps^2/4 (bounding box [-2eps,4] x [-eps,4]).
    Returns rational bounds; the margin is decreasing in eps."""
    e = Q(1, k)
    osc_upper = Q(19, 24)+Q(19, 4)*e*e
    r2_upper = 8+6*e+Q(5, 4)*e*e
    p = pi_enclosure(terms)
    gap_lower = min(Q(8), 16-p[1]*r2_upper/2)
    return {"osc_upper": osc_upper, "radius_sq_upper": r2_upper, "gap_lower": gap_lower,
            "margin": gap_lower-2*osc_upper}


def fixed_target_bounds(terms=40):
    """Bound for every truncation of the fixed rough target (eps_j = 2^(-j-3), N_j = eps_j^-2,
    heights h_j = 2^-j): Osc <= max_j (19/24 + (19/4) eps_j^2 + (9/4) h_j eps_j) + sum_j 4 eps_j^2.

    Each gadget obeys the single-gadget bound about its own attachment point; moving the origin to
    (0,0) changes its oscillation by at most sqrt(5) h_j eps_j < (9/4) h_j eps_j; the primitive at the
    start of each gadget drifts by the enclosed areas of the earlier gadgets, each at most 4 eps^2.
    Together with A >= 16 and R^2 = 8 this gives the margin for every truncation."""
    worst = max(Q(19, 24)+Q(19, 4)*Q(1, 4**(j+3))+Q(9, 4)*Q(1, 2**j)*Q(1, 2**(j+3)) for j in range(1, 40))
    drift = 4*Q(1, 4**4)*Q(4, 3)      # sum_{j>=1} 4 eps_j^2
    osc_upper = worst+drift
    p = pi_enclosure(terms)
    gap_lower = 16-4*p[1]
    return {"osc_upper": osc_upper, "gap_lower": gap_lower, "margin": gap_lower-2*osc_upper,
            "window_half_width_lower": (gap_lower-2*osc_upper)/8}


# ------------------------------------------------------------------ report

def fmt(x):
    if isinstance(x, tuple):
        return [fmt(y) for y in x]
    if isinstance(x, Q):
        return {"exact": str(x), "approx": float(x)} if x.denominator < 10**30 else {"approx": float(x)}
    return x


def report():
    out = {}
    # v328 models
    mesh = standard_mesh(64)
    rows = []
    for m in (1, 4, 16):
        circ = window_enclosure(mesh, circle_values(mesh, Q(9)), m, Q(9))
        plat = window_enclosure(mesh, plateau_values(mesh, Q(1), Q(1, 2)), m, Q(1))
        coarse_mesh = [Q(0), Q(1), INF]
        coarse = window_enclosure(coarse_mesh, circle_values(coarse_mesh, Q(9)), m, Q(9))
        rows.append({"m": m, "circle_a2_9": [float(x) for x in circ],
                     "plateau_K1_b_half": [float(x) for x in plat],
                     "circle_three_point_mesh": [float(x) for x in coarse]})
    out["v328_models"] = {"mesh_points": len(mesh), "rows": rows}
    # v329 ellipses
    ellipses = []
    for scale in (1, 3):
        for ratio in (1, 2, 8, 32):
            a, b = Q(ratio*scale), Q(scale)
            row = {"a": str(a), "b": str(b), "inradius_sq": str(b*b),
                   "square_radius_sq": str(2*a*a*b*b/(a*a+b*b)),
                   "mean_slope_A_over_pi": str(a*b), "windows": []}
            for m in (0, 1, 4, 16, 64, 256):
                w = ellipse_window(m, a, b)
                row["windows"].append({"m": m, "D_single": float(w["D_single"]),
                                       "D_low": [float(x) for x in w["D_low"]],
                                       "D_high": [float(x) for x in w["D_high"]]})
            ellipses.append(row)
    out["v329_ellipses"] = ellipses
    # v330 spirals
    spirals = []
    for k in (2, 3, 4, 6, 7, 8, 10, 12):
        families = [("original", k*k)] + ([("retuned", retuned_turns(k))] if k in (2, 4, 8) else [])
        for label, N in families:
            P = spiral_fixture(k, N)
            c = carrier_exclusion(P)
            row = {"family": label, "k": k, "turns": N,
                   "tagged_return_integral": float(tagged_return_integral(k, N)),
                   **{key: (float(val) if isinstance(val, Q) else val) for key, val in c.items()}}
            if label == "original":
                bound = original_family_bounds(k)
                if not (c["oscillation"] <= bound["osc_upper"] and c["area"] >= 16
                        and c["radius_sq"] <= bound["radius_sq_upper"]):
                    raise ArithmeticError("analytic bound violated", k)
                row["analytic_margin"] = float(bound["margin"])
            spirals.append(row)
    out["v330_spirals"] = spirals
    if original_family_bounds(7)["margin"] <= 0:
        raise ArithmeticError("analytic margin fails at k = 7")
    out["v330_all_k_at_least_7"] = {key: float(val) for key, val in original_family_bounds(7).items()}
    fixed, exact = [], []
    for J in (1, 2, 3):
        gadgets = [(4**(j+3), Q(1, 2**(j+3)), Q(1, 2**j)) for j in range(J, 0, -1)]
        c = carrier_exclusion(square_with_gadgets(gadgets))
        exact.append(c)
        fixed.append({"gadgets": J, **{key: (float(val) if isinstance(val, Q) else val)
                                       for key, val in c.items()}})
    out["v330_fixed_rough_target_truncations"] = fixed
    bound = fixed_target_bounds()
    if bound["margin"] <= 0 or any(c["oscillation"] > bound["osc_upper"] or c["area"] < 16
                                   or c["radius_sq"] != 8 for c in exact):
        raise ArithmeticError("fixed-target bound fails")
    out["v330_fixed_target_all_truncations"] = {key: float(val) for key, val in bound.items()}
    return out


if __name__ == "__main__":
    print(json.dumps(report(), indent=2))
