"""v390: a square's label is a conserved quantity, and its degree survives any deformation that cannot push the
square out of its label region.

Coordinates (no polar-angle chart anywhere). Corner A sits on arm arm_A at its own lift Theta_A:
    A = Theta_A^(-a) exp(i (Theta_A - pi arm_A + arm_A lam phi(Theta_A^(-a)))),
and the other corners are X = A (1 + w m_X), m_B = 1, m_C = 1 + i, m_D = i, with w = (side / |A|) e^(i orientation)
measured in A's own frame. The label is (arm_A, o_B, o_C, o_D): A's arm and the other corners' pass offsets from A.
Corner X's residual is
    |X|^(-1/a) - Theta_A - arg(X / A) - pi o_X + lam (arm_X phi(|X|) - arm_A phi(|A|)),   arm_X = arm_A + o_X mod 2,
with the principal argument, which is continuous because |w m_X| < 1 keeps X / A off the negative real axis. So the
label map F_lam : (Theta_A, Re w, Im w) -> R^3 is continuous on the closure of
    Omega = {Theta_lo < Theta_A < Theta_hi, q_in < |w| < q_max},   q_max < 1 / sqrt 2,
a bounded open subset of R^3, and its Brouwer degree there is defined whenever F_lam has no zero on the boundary.

(An absolute label, pass indices k_X = (lift - theta_X) / pi with theta_X a polar angle, is chart-dependent: when
the angle crosses its branch cut every k_X changes by 2 although nothing geometric happens. phase_map below keeps
that older form only to show the jump.)

Consequences:
  * a square moving continuously keeps its label: A cannot change arm without passing through the wild point, and
    the offsets are integers that move continuously while |w| < 1 / sqrt 2;
  * the degree on a label region is constant in lam while no square of the label reaches the region's boundary;
  * parity bookkeeping alone is tautological (events at s = L(gamma_0, s) + L(gamma_1, s) mod 2);
  * a cascade through half-turn partners is impossible: the half-turn swaps arm_A and keeps the offsets, so
    partners carry different labels;
  * so a square of gamma_0 can disappear only by leaving its label region: collapsing into the wild point, or
    drifting out of its band of lifts or sides.
"""
import numpy as np
from math import pi

from smooth_cores_v386 import solve_twisted, sector_two_arm
from twisted_classification_v387 import nearest_pass, twisted_corners
from two_turn_squares_v383 import window

A_EXP = 0.5


M_CORNER = (1.0, 1 + 1j, 1j)   # B, C, D relative to A, in units of w


def anchor_point(Theta_A, arm_A, lam, phi, a=A_EXP):
    r = Theta_A ** (-a)
    return r * np.exp(1j * (Theta_A - pi * arm_A + arm_A * lam * phi(r)))


def anchored_corners(x, arm_A, lam, phi, a=A_EXP):
    A = anchor_point(x[0], arm_A, lam, phi, a)
    w = complex(x[1], x[2])
    return [A] + [A * (1 + w * m) for m in M_CORNER]


def anchored_map(x, label, lam, phi, a=A_EXP):
    """The three residuals of the label map at x = (Theta_A, Re w, Im w); label = (arm_A, o_B, o_C, o_D)."""
    arm_A, offs = label[0], label[1:]
    Z = anchored_corners(x, arm_A, lam, phi, a)
    A = Z[0]
    phiA = phi(abs(A))
    out = []
    for X, o in zip(Z[1:], offs):
        arm_X = (arm_A + o) % 2
        out.append(abs(X) ** (-1 / a) - x[0] - np.angle(X / A) - pi * o + lam * (arm_X * phi(abs(X)) - arm_A * phiA))
    return np.array(out)


def anchored_jacobian_sign(x, label, lam, phi, a=A_EXP, h=1e-7):
    """Sign of det DF at a zero of the label map (local degree +-1; 0 means degenerate), and a conditioning number."""
    y = anchored_map(x, label, lam, phi, a)
    J = np.zeros((3, 3))
    steps = (h * x[0], h * abs(complex(x[1], x[2])), h * abs(complex(x[1], x[2])))
    for m in range(3):
        dx = np.zeros(3)
        dx[m] = steps[m]
        J[:, m] = (anchored_map(x + dx, label, lam, phi, a) - y) / steps[m]
    d = np.linalg.det(J)
    scale = np.prod(np.linalg.norm(J, axis=0))
    return int(np.sign(d)) if abs(d) > 1e-9 * scale else 0, float(abs(d) / scale)


def anchored_from_square(Z, ks, i=0, a=A_EXP):
    """Anchored coordinates and label of the square Z (cyclic corners) with corner i as A; ks are any consistent pass
    indices (only their differences and A's parity are used)."""
    o = [(i + q) % 4 for q in range(4)]
    A = Z[o[0]]
    x = np.array([abs(A) ** (-1 / a), 0.0, 0.0])
    w = (Z[o[1]] - A) / A
    x[1], x[2] = w.real, w.imag
    label = (int(ks[o[0]]) % 2,) + tuple(int(ks[o[q]] - ks[o[0]]) for q in (1, 2, 3))
    return x, label


def partner_anchored_label(label):
    """The half-turn swaps the arms and keeps (Theta_A, w) and the offsets."""
    return (1 - label[0],) + tuple(label[1:])


def phase_map(v, ks, lam, phi, a=A_EXP):
    """Chart-dependent absolute form (for comparison only): the four residuals of v = (Ax, Ay, angle, side) with
    corners on absolute passes ks, polar angles taken from v383's window (which has a branch cut)."""
    Z = twisted_corners(v)
    th, _, _ = window(Z, [0, 0, 0, 0])
    out = []
    for z, t, k in zip(Z, th, ks):
        val = abs(z) ** (-1 / a) - t - pi * k
        if k % 2:
            val += lam * phi(abs(z))
        out.append(val)
    return np.array(out)


def jacobian_sign(v, ks, lam, phi, a=A_EXP, h=1e-7):
    """Sign of det DF at a zero (local degree +-1; 0 means degenerate)."""
    y = phase_map(v, ks, lam, phi, a)
    J = np.zeros((4, 4))
    t = abs(v[3])
    for m in range(4):
        dv = np.zeros(4)
        dv[m] = h * t if m != 2 else h
        J[:, m] = (phase_map(v + dv, ks, lam, phi, a) - y) / dv[m]
    d = np.linalg.det(J)
    scale = np.prod(np.linalg.norm(J, axis=0))
    return int(np.sign(d)) if abs(d) > 1e-9 * scale else 0, float(abs(d) / scale)


def partner_label(Z, ks, a=A_EXP):
    """Absolute pass labels of the half-turn image (on the symmetric spiral)."""
    Zp = [-z for z in Z]
    th, _, _ = window(Zp, [0, 0, 0, 0])
    return [nearest_pass(z, t, 0.0, lambda r: 0.0, a) for z, t in zip(Zp, th)]


def labels_of_plain_squares(n_random=6000, seed=7, a=A_EXP):
    """Plain-spiral squares from the v380 search, each in all four anchorings, with residuals, local degrees and
    partner checks."""
    from general_comb_v380 import corners, search
    zero = lambda r: 0.0
    rows = []
    for v, Ths in search(0.0, 1.0, n_random=n_random, seed=seed):
        Z = corners(v)
        th, k, W = window(Z, Ths)
        ks = list(map(int, k))
        for i in range(4):
            x, lab = anchored_from_square(Z, ks, i, a)
            sgn, cond = anchored_jacobian_sign(x, lab, 0.0, zero, a)
            plab = partner_anchored_label(lab)
            rows.append(dict(anchor=i, label=list(lab), partner=list(plab),
                             residual=float(np.max(np.abs(anchored_map(x, lab, 0.0, zero, a)))),
                             partner_residual=float(np.max(np.abs(anchored_map(x, plab, 0.0, zero, a)))),
                             partner_corners_match=bool(np.allclose(sorted(np.round([-z for z in Z], 12), key=abs),
                                                                    sorted(np.round(anchored_corners(x, plab[0], 0.0, zero, a), 12), key=abs),
                                                                    atol=1e-9 * abs(Z[0]))),
                             degree_sign=sgn, conditioning=cond, sector=sector_two_arm(ks, Z, a),
                             t_over_A=float(abs(complex(x[1], x[2])))))
    return rows


def chart_check(a=A_EXP):
    """Move a square's A (and centre) across the negative real axis: the anchored residuals change continuously,
    the absolute (chart) residuals jump by 2 pi."""
    zero = lambda r: 0.0
    m = 40
    lab = (0, 2, 3, 1)
    q = 0.05
    w = q * (1 - 1j) / 2                      # (1 + i) w is real: the centre lies on A's ray
    out = []
    for eps in (-1e-9, 1e-9):
        x = np.array([pi * (2 * m + 1) + eps, w.real, w.imag])
        Z = anchored_corners(x, 0, 0.0, zero, a)
        ks = [2 * m] + [2 * m + o for o in lab[1:]]
        v = np.array([Z[0].real, Z[0].imag, np.angle(Z[1] - Z[0]), abs(Z[1] - Z[0])])
        out.append((anchored_map(x, lab, 0.0, zero, a), phase_map(v, ks, 0.0, zero, a)))
    return dict(anchored_jump=float(np.max(np.abs(out[0][0] - out[1][0]))),
                absolute_jump=float(np.max(np.abs(out[0][1] - out[1][1]))))


def run_checks(path=None):
    import json
    from collections import Counter
    from general_comb_v380 import corners, search
    from twisted_classification_v387 import on_branch
    rows = labels_of_plain_squares()
    res = dict(plain=dict(configurations=len(rows), squares=len(rows) // 4,
                          max_residual=max(r["residual"] for r in rows),
                          partner_is_zero_of_partner_label=sum(1 for r in rows if r["partner_residual"] < 1e-8
                                                               and r["partner_corners_match"]),
                          partner_same_label=sum(1 for r in rows if r["partner"] == r["label"]),
                          nondegenerate=sum(1 for r in rows if r["degree_sign"] != 0),
                          degree_signs=dict(Counter(r["degree_sign"] for r in rows)),
                          min_conditioning=min(r["conditioning"] for r in rows),
                          max_side_over_A=max(r["t_over_A"] for r in rows),
                          in_small_square_domain=sum(1 for r in rows if r["t_over_A"] < 1 / np.sqrt(2))))
    # one zero per anchored label (the uniqueness hypothesis U0), rerun with another seed: no label may be carried
    # by two different geometric squares, in any of the four anchorings
    sq = search(0.0, 1.0, n_random=6000, seed=19)
    owner = {}
    clash = 0
    for n, (v, Ths) in enumerate(sq):
        Z = corners(v)
        th, k, _ = window(Z, Ths)
        for i in range(4):
            _, lab = anchored_from_square(Z, list(map(int, k)), i)
            if owner.setdefault(lab, n) != n:
                clash += 1
    ok = [on_branch(corners(v))[1] for v, _ in sq]
    res["uniqueness"] = dict(squares=len(sq), anchored_labels=len(owner), labels_shared_by_two_squares=clash,
                             on_branch=int(sum(ok)))
    res["chart"] = chart_check()
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import json
    from pathlib import Path
    print(json.dumps(run_checks(Path(__file__).resolve().parents[1] / "results" / "v390" / "label_check.json"), indent=1))
