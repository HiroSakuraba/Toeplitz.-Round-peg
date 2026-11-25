"""v390: a square's label is a conserved quantity, and its degree survives any deformation that cannot push the
square out of its label region.

A label is an assignment of the four corners to passes (absolute pass indices). Along any continuous family of
curves gamma_lam, a square moves continuously with its label fixed: a corner cannot change pass without passing
through the wild point or through a fold of the pass structure. For a label region Omega (corners near their passes,
A in an annulus of radii, side in an interval), the mod-2 degree of the square map on Omega is constant in lam as long
as no square of any gamma_lam lies on the boundary of Omega (homotopy invariance). This replaces event counting:
  * parity bookkeeping alone is tautological (events at s = L(gamma_0, s) + L(gamma_1, s) mod 2);
  * a cascade through half-turn partners is impossible, because partners carry different labels and a branch keeps
    its label;
  * so a square of gamma_0 can disappear only by leaving its label region: collapsing into the wild point, or
    drifting out of its band of radii.
"""
import numpy as np
from math import pi

from smooth_cores_v386 import solve_twisted, sector_two_arm
from twisted_classification_v387 import nearest_pass, twisted_corners
from two_turn_squares_v383 import window

A_EXP = 0.5


def phase_map(v, ks, lam, phi, a=A_EXP):
    """The four comb-phase residuals of the square v = (Ax, Ay, angle, side) with corners on passes ks."""
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
    """Plain-spiral squares from the v380 search with their absolute labels, partner labels and local degrees."""
    from general_comb_v380 import corners, search
    rows = []
    for v, Ths in search(0.0, 1.0, n_random=n_random, seed=seed):
        Z = corners(v)
        th, k, W = window(Z, Ths)
        ks = list(map(int, k))
        w = np.array([Z[0].real, Z[0].imag, np.angle(Z[1] - Z[0]), abs(Z[1] - Z[0])])
        sgn, cond = jacobian_sign(w, ks, 0.0, lambda r: 0.0, a)
        pl = partner_label(Z, ks, a)
        rows.append(dict(ks=ks, partner=list(map(int, pl)), same_label=bool(pl == ks),
                         degree_sign=sgn, conditioning=cond, sector=sector_two_arm(ks, Z, a),
                         t_over_Z=float(abs(Z[1] - Z[0]) / abs(sum(Z) / 4))))
    return rows


def run_checks(path=None):
    import json
    from collections import Counter
    from general_comb_v380 import corners, search
    from twisted_classification_v387 import on_branch
    rows = labels_of_plain_squares()
    res = dict(plain=dict(squares=len(rows), partner_same_label=sum(r["same_label"] for r in rows),
                          nondegenerate=sum(1 for r in rows if r["degree_sign"] != 0),
                          degree_signs=dict(Counter(r["degree_sign"] for r in rows)),
                          min_conditioning=min(r["conditioning"] for r in rows)))
    # uniqueness per label on the plain spiral (v387 Proposition 387.2, rerun with another seed)
    sq = search(0.0, 1.0, n_random=6000, seed=19)
    labs = []
    for v, Ths in sq:
        Z = corners(v)
        th, k, _ = window(Z, Ths)
        labs.append(tuple(int(x) for x in k))
    ok = [on_branch(corners(v))[1] for v, _ in sq]
    res["uniqueness"] = dict(squares=len(sq), distinct_labels=len(set(labs)), on_branch=int(sum(ok)))
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import json
    from pathlib import Path
    print(json.dumps(run_checks(Path(__file__).resolve().parents[1] / "results" / "v390" / "label_check.json"), indent=1))
