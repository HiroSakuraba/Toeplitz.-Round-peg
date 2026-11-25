"""v395: twin arms, and stability under every slowly varying angular deviation.

1. Twin-arm lemma. Let p lie on a Jordan curve, with arcs alpha, beta ending at p, each radially monotone on
   0 < |z - p| < r0 (|z - p| strictly monotone along it, so it meets each small circle once). Take continuous angle
   functions theta_alpha(r), theta_beta(r). The arcs are disjoint, so theta_beta(r) - theta_alpha(r) is never a
   multiple of 2 pi, and by continuity on the connected interval (0, r0) it stays in one interval (2 pi k, 2 pi (k+1)).
   Hence: if one arm winds infinitely often around p so does the other, in the same direction, and with
   theta_beta = theta_alpha + pi + phi(r), |phi| < pi after choosing the lift. A radially monotone wild point is
   therefore a two-arm spiral with twist below contact; only the common winding profile theta_alpha(r) is free.

2. Stability. Write the arms as theta_1(r) = r^(-1/a) + e1(r) and theta_2(r) = r^(-1/a) - pi + e1(r) + phi(r): any
   radially monotone two-arm curve has this form. v390's anchored label map sees angles only through differences
   between corners, so e1 enters only through its variation: with |r e1'| <= L1 it adds L1 |log(|X|/|A|)| to each
   corner's lift deviation and 2 L1 log((1+u)/(1-u)) to the identity's twist sum. The size of e1 is irrelevant, so
   the arms may rotate together without bound, as long as they do it slowly. v392's estimates take the pair (L, L1)
   (uniform_closure_v392._slopes); v393/v394's degree step involves only the plain spiral and is unchanged.
"""
from math import log, pi, sin, sqrt

import numpy as np

from contact_twist_v394 import DELTAS_LO, annulus, entry_ok
from uniform_closure_v392 import Q_MAX, certify_levels, plain_square, resonance

A_EXP = 0.5
M_CORNER = (1.0, 1 + 1j, 1j)


# ----------------------------------------------------------------------------------------------------------------
# Twin-arm lemma, numerically
# ----------------------------------------------------------------------------------------------------------------

def first_meeting(theta1, theta2, r_hi, r_lo, n=200000):
    """Largest r in [r_lo, r_hi] where two radially monotone arms with angle functions theta1, theta2 meet
    (theta2 - theta1 a multiple of 2 pi), or None. Sign changes of sin((theta2 - theta1) / 2) on a log grid."""
    rs = np.geomspace(r_hi, r_lo, n)
    g = np.sin((np.array([theta2(r) for r in rs]) - np.array([theta1(r) for r in rs])) / 2)
    idx = np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]
    return float(rs[idx[0]]) if len(idx) else None


def relative_twist(theta1, theta2, rs):
    """phi(r) = theta2 - theta1 - pi, with the lift chosen so that phi(rs[0]) is in (-pi, pi]."""
    d = np.array([theta2(r) - theta1(r) - pi for r in rs])
    k = np.round(d[0] / (2 * pi))
    return d - 2 * pi * k


# ----------------------------------------------------------------------------------------------------------------
# Label map with both arms deformed
# ----------------------------------------------------------------------------------------------------------------

def arm_angle(r, arm, e1, phi, lam=1.0, a=A_EXP):
    """Polar angle (unreduced) of arm `arm` at radius r on the curve deformed by lam (e1, phi)."""
    return r ** (-1 / a) - pi * arm + lam * (e1(r) + arm * phi(r))


def corners2(x, arm_A, lam, e1, phi, a=A_EXP):
    T = x[0]
    r = T ** (-a)
    A = r * np.exp(1j * arm_angle(r, arm_A, e1, phi, lam, a))
    w = complex(x[1], x[2])
    return [A] + [A * (1 + w * m) for m in M_CORNER]


def label_map2(x, label, lam, e1, phi, a=A_EXP):
    """v390's anchored residuals on the curve deformed by lam (e1, phi):
    |X|^(-1/a) - Theta_A - arg(X/A) - pi o_X + lam (E_X - E_A), E = e1 + arm phi."""
    arm_A, offs = label[0], label[1:]
    Z = corners2(x, arm_A, lam, e1, phi, a)
    A = Z[0]
    EA = e1(abs(A)) + arm_A * phi(abs(A))
    out = []
    for X, o in zip(Z[1:], offs):
        arm_X = (arm_A + o) % 2
        EX = e1(abs(X)) + arm_X * phi(abs(X))
        out.append(abs(X) ** (-1 / a) - x[0] - np.angle(X / A) - pi * o + lam * (EX - EA))
    return np.array(out)


def on_curve_error(x, label, lam, e1, phi, a=A_EXP):
    """Largest angular distance (mod 2 pi) between a corner and the deformed curve's arm at the corner's radius:
    checks that zeros of label_map2 really are squares inscribed in the deformed curve."""
    Z = corners2(x, label[0], lam, e1, phi, a)
    arms = [label[0]] + [(label[0] + o) % 2 for o in label[1:]]
    err = 0.0
    for X, arm in zip(Z, arms):
        dth = np.angle(X) - arm_angle(abs(X), arm, e1, phi, lam, a)
        err = max(err, abs((dth + pi) % (2 * pi) - pi))
    return err


def solve2(x, label, lam, e1, phi, a=A_EXP, it=80):
    for _ in range(it):
        y = label_map2(x, label, lam, e1, phi, a)
        if np.max(np.abs(y)) < 1e-12:
            break
        q = abs(complex(x[1], x[2]))
        steps = (1e-8 * x[0], 1e-8 * q, 1e-8 * q)
        J = np.zeros((3, 3))
        for m in range(3):
            dx = np.zeros(3)
            dx[m] = steps[m]
            J[:, m] = (label_map2(x + dx, label, lam, e1, phi, a) - y) / steps[m]
        try:
            st = np.linalg.solve(J, -y)
        except np.linalg.LinAlgError:
            return None
        sc = max(1.0, abs(st[0]) / (0.02 * x[0]), np.hypot(st[1], st[2]) / (0.1 * q))
        x = x + st / sc
        if x[0] <= 0 or abs(complex(x[1], x[2])) >= 1 / sqrt(2):
            return None
    return x if np.max(np.abs(label_map2(x, label, lam, e1, phi, a))) <= 1e-9 else None


def continue2(x0, label, e1, phi, steps=400, a=A_EXP):
    lam, x = 0.0, x0.copy()
    path = [(0.0, x.copy())]
    h = 1.0 / steps
    while lam < 1.0:
        nl = 1.0 if lam + h > 1.0 - 1e-9 else lam + h
        res = solve2(x.copy(), label, nl, e1, phi, a)
        if res is None or abs(abs(complex(res[1], res[2])) / abs(complex(x[1], x[2])) - 1) > 0.2:
            if h < 1.0 / (steps * 256):
                break
            h /= 2
            continue
        lam, x = nl, res
        path.append((lam, x.copy()))
        h = min(1.0 / steps, 2 * h)
    return path


def degree_sign(x, label, lam, e1, phi, a=A_EXP, h=1e-7):
    y = label_map2(x, label, lam, e1, phi, a)
    q = abs(complex(x[1], x[2]))
    steps = (h * x[0], h * q, h * q)
    J = np.zeros((3, 3))
    for m in range(3):
        dx = np.zeros(3)
        dx[m] = steps[m]
        J[:, m] = (label_map2(x + dx, label, lam, e1, phi, a) - y) / steps[m]
    return int(np.sign(np.linalg.det(J)))


# ----------------------------------------------------------------------------------------------------------------
# Certificates with a common deviation
# ----------------------------------------------------------------------------------------------------------------

def certify_entry(Phi, L, L1, c):
    return certify_levels(Phi, (L, L1), c, 1, deltas_lo=DELTAS_LO)


def best_c(Phi, L, L1, c0=11, c_cap=200000):
    c = c0
    while not certify_entry(Phi, L, L1, c)["certified"]:
        c *= 2
        if c > c_cap:
            return None
    lo, hi = max(c0 - 1, c // 2), c
    if hi == c0:
        return c0
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if certify_entry(Phi, L, L1, mid)["certified"]:
            hi = mid
        else:
            lo = mid
    return hi


DELTA_GRID = (1.0, 0.3, 0.1, 0.03)
L_GRID = (0, 5, 20)
L1_GRID = (0, 1, 2, 5, 10)


def threshold_table():
    rows = []
    for d in DELTA_GRID:
        for L in L_GRID:
            for L1 in L1_GRID:
                c = best_c(pi - d, L, L1)
                cert = certify_entry(pi - d, L, L1, c) if c else dict(certified=False)
                if c:
                    cert["c"] = c
                rows.append(dict(delta=d, Phi=pi - d, L=L, L1=L1, c=c, cert=cert,
                                 degree_applies=bool(c and entry_ok(cert))))
    return rows


# ----------------------------------------------------------------------------------------------------------------
# An example whose arms rotate together without bound
# ----------------------------------------------------------------------------------------------------------------

def example_e1(r, k=5.0, b=1.0, f=3.0):
    """Common rotation k log log(1/r) + b sin(f log(1/r)): unbounded, |r e1'| <= k / log(1/r) + b f."""
    return k * log(log(1 / r)) + b * sin(f * log(1 / r))


def example_phi(r, amp=2.0, f=4.0):
    return amp * sin(f * log(1 / r))


def example_bounds(r_out, k=5.0, b=1.0, f=3.0, amp=2.0, fphi=4.0):
    """Exact bounds on (0, r_out]: |phi| <= amp, |r phi'| <= amp fphi, |r e1'| <= k / log(1/r_out) + b f."""
    return amp, amp * fphi, k / log(1 / r_out) + b * f


def run_example(label, table):
    """Follow the plain square of label through the deformation lam (e1, phi), lam: 0 -> 1, and check it against the
    region of the smallest table entry whose bounds hold on the label's annulus."""
    n1, n3, j = label
    m = min(abs(n1), abs(n3))
    hit = None
    for row in sorted((r for r in table if r["degree_applies"]), key=lambda r: r["c"]):
        if m < row["c"]:
            continue
        r_in, r_out = annulus(row["cert"], n1, n3, j)
        if r_out >= 1 / np.e:
            continue
        Phi, L, L1 = example_bounds(r_out)
        if Phi <= row["Phi"] and L <= row["L"] and L1 <= row["L1"]:
            hit = row
            break
    if hit is None:
        return dict(label=list(label), covered=False)
    cert = hit["cert"]
    x0, lab = plain_square(*label)
    path = continue2(x0, lab, example_e1, example_phi)
    T = resonance(*label)
    lifts = [x[0] / T for _, x in path]
    qs = [abs(complex(x[1], x[2])) for _, x in path]
    q_in = abs(j) / m * cert["side_floor"] / 2
    xe = path[-1][1]
    reached = path[-1][0] == 1.0
    r_in, r_out = annulus(cert, n1, n3, j)
    return dict(label=list(label), covered=True, entry=dict(delta=hit["delta"], L=hit["L"], L1=hit["L1"], c=hit["c"]),
                annulus=[r_in, r_out], rotation_at_annulus=[example_e1(r_out), example_e1(r_in)],
                reached=bool(reached), s_range=[min(lifts), max(lifts)], walls=[cert["s_lo"], cert["s_hi"]],
                q_range=[min(qs), max(qs)], q_in=q_in,
                inside=bool(cert["s_lo"] < min(lifts) and max(lifts) < cert["s_hi"] and q_in < min(qs)
                            and max(qs) < Q_MAX),
                on_curve_error=float(on_curve_error(xe, lab, 1.0, example_e1, example_phi)) if reached else None,
                degree_start=degree_sign(x0, lab, 0.0, example_e1, example_phi),
                degree_end=degree_sign(xe, lab, 1.0, example_e1, example_phi) if reached else None)


def twin_checks():
    """Numerical illustrations of the twin-arm lemma."""
    a = A_EXP
    base = lambda r: r ** (-1 / a)
    out = {}
    # opposite winding, and one arm winding while the other doesn't: the arms meet
    out["opposite_winding_meets_at"] = first_meeting(base, lambda r: -base(r) + 1.0, 0.5, 1e-3)
    out["winding_vs_straight_meets_at"] = first_meeting(base, lambda r: 1.0, 0.5, 1e-3)
    # unequal rates (exponents 1/2 and 0.45): the relative angle grows, so they meet
    out["unequal_rates_meet_at"] = first_meeting(base, lambda r: r ** (-1 / 0.45) + pi, 0.5, 1e-3)
    # twins: common profile plus a twist below contact never meet, and the relative angle stays in (-pi, pi)
    tw = lambda r: base(r) + example_e1(r) + pi + example_phi(r)
    th = lambda r: base(r) + example_e1(r)
    out["twins_meet_at"] = first_meeting(th, tw, 0.3, 1e-4)
    rs = np.geomspace(0.3, 1e-4, 20000)
    ph = relative_twist(th, tw, rs)
    out["twins_relative_twist_range"] = [float(ph.min()), float(ph.max())]
    return out


EXAMPLE_LABELS = [(150, 160, 1), (-200, -220, 1), (250, -260, -1), (300, 310, 1)]


if __name__ == "__main__":
    import json
    from pathlib import Path
    out = Path(__file__).resolve().parents[1] / "results" / "v395"
    table = json.loads((out / "table.json").read_text())
    res = dict(twin=twin_checks(),
               table=[{k: v for k, v in r.items() if k != "cert"} | dict(s_lo=r["cert"].get("s_lo"),
                                                                          s_hi=r["cert"].get("s_hi")) for r in table],
               examples=[run_example(lab, table) for lab in EXAMPLE_LABELS])
    (out / "twin_check.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res["twin"], indent=1))
    for r in res["table"]:
        print(r["delta"], r["L"], r["L1"], r["c"], r["degree_applies"])
    for e in res["examples"]:
        print(e)
