"""v394: every twist below contact, at every slope, and only where each square lives.

1. In the two-arm power-spiral model the twist can't grow without limit. At radius r arm one sits at polar angle
   Theta - pi k (k even) and arm two at Theta - pi k (k odd) + phi(r), Theta = r^(-1/a): one point per arm per radius.
   They meet exactly when phi(r) = pi (mod 2 pi). A Jordan curve doesn't meet itself and phi is continuous, so once
   |phi| < pi somewhere, |phi| < pi everywhere. The open cases are arms approaching contact (sup |phi| = pi, not
   attained) and fast twisting (large |r phi'|).
2. Degree for every twist bound. v393's fine certificates never involved the twist. The coarse step is redone once,
   for the plain spiral: v392's tail certificate at zero twist excludes plain squares with scaled lift >= 1.5, and a
   branch-and-bound puts every plain zero with scaled lift in [1/4, 1.5] into v393's box B*. So every v392 region
   with s_lo >= 1/4, min offset >= 11 and B* inside it has degree +-1 (|j| = 1, M <= 2m).
3. A certified table of thresholds c(Phi, L) up to Phi = pi - 0.01 and L = 50: squares of every label with
   min offset m >= c (and >= 11), M <= 2m, |j| = 1, on every two-arm power spiral (a = 1/2) whose twist obeys
   |phi| <= Phi, |r phi'| <= L.
4. Local bounds: every certificate involves the twist only at the corners of squares whose corner A lies in the
   region, so the bounds are needed only on the label's annulus, not everywhere.
"""
from math import pi, sqrt

import numpy as np

from uniform_closure_v392 import (Q_MAX, certify_levels, resonance, scaled_side_bound, upper_face, plain_square)
from uniform_degree_v393 import (B_STAR, I, localize, patterns, star_inside_region)

A_EXP = 0.5
S_LO_MIN = 0.25
DELTAS_LO = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5)          # keeps s_lo = pi/(pi + Phi)(1 - delta) >= 1/4
N_DEG = 11                                           # v393's fine certificates cover eps <= 1/11, r in [1/2, 1]
RHO = 2.0


def arm_points(r, phi, a=A_EXP):
    """The two arms' points at radius r (twist phi there)."""
    T = r ** (-1 / a)
    return r * np.exp(1j * T), r * np.exp(1j * (T - pi + phi))


# ----------------------------------------------------------------------------------------------------------------
# Degree for every region: the plain-spiral tail and one coarse localization
# ----------------------------------------------------------------------------------------------------------------

def plain_tail(dhi=0.5):
    """v392's tail certificate with zero twist: no plain-spiral square of side <= Q_MAX |A| (any label with
    min offset >= 11, any signs) has scaled lift >= s_tail."""
    s_tail, ok = upper_face(0.0, 0.0, 1, 1 / N_DEG, 1 / N_DEG, dhi, uniform=True)
    return s_tail, ok


def plain_side_bound():
    """Side bound for plain squares at scaled lift >= 1/4 (uniform in eps <= 1/11): side <= eps Q0 |A|."""
    return scaled_side_bound(1 / N_DEG, 1 / N_DEG, S_LO_MIN, 0.0, 0.0, 1, uniform=True)


def coarse_plain(sig1, sig3, b_is_min, j, s_tail, Q0, depth=0, eps=None, r=None, max_depth=6, log=None):
    """Every plain zero with scaled lift in [1/4, s_tail] lies in B* (parameters split as needed)."""
    eps = eps if eps is not None else I(0, 1 / N_DEG)
    r = r if r is not None else I(1 / RHO, 1)
    par = (eps, r, sig1, sig3, b_is_min, j)
    t1 = float(r.a) if b_is_min else 1.0
    t3 = 1.0 if b_is_min else float(r.a)
    region = [I(S_LO_MIN, s_tail), I(-4 * Q0 / t1, 4 * Q0 / t1), I(-4 * Q0 / t3, 4 * Q0 / t3)]
    star = [I(1 - B_STAR, 1 + B_STAR)] * 3
    ok, n = localize(region, star, par, max_boxes=150000)
    if ok:
        return [dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], boxes=n)]
    if depth >= max_depth:
        if log is not None:
            log.append(dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)]))
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
        sub = coarse_plain(sig1, sig3, b_is_min, j, s_tail, Q0, depth + 1, e2, r2, max_depth, log)
        if sub is None:
            return None
        out.extend(sub)
    return out


def run_coarse(j):
    import time
    s_tail, ok = plain_tail()
    Q0 = plain_side_bound()
    out = dict(j=j, s_tail=s_tail, tail_certified=bool(ok), Q0=Q0, patterns={})
    for sig1, sig3, b_is_min in patterns(j):
        t0 = time.time()
        log = []
        co = coarse_plain(sig1, sig3, b_is_min, j, s_tail, Q0, log=log)
        out["patterns"][f"{sig1},{sig3},{b_is_min}"] = dict(certified=co is not None, pieces=len(co or []),
                                                           boxes=sum(c["boxes"] for c in co or []),
                                                           seconds=round(time.time() - t0, 1), failures=log[:5])
    return out


# ----------------------------------------------------------------------------------------------------------------
# Thresholds for every twist bound, and where they apply
# ----------------------------------------------------------------------------------------------------------------

def certify_entry(Phi, L, c):
    return certify_levels(Phi, L, c, 1, deltas_lo=DELTAS_LO)


def best_c(Phi, L, c0=N_DEG, c_cap=100000):
    """Smallest c (>= c0) with a certified all-levels region and s_lo >= 1/4 (doubling, then bisection)."""
    c = c0
    while not certify_entry(Phi, L, c)["certified"]:
        c *= 2
        if c > c_cap:
            return None
    lo, hi = max(c0 - 1, c // 2), c
    if hi == c0:
        return c0
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if certify_entry(Phi, L, mid)["certified"]:
            hi = mid
        else:
            lo = mid
    return hi


def entry_ok(cert):
    """Conditions for v393's degree to apply to this region: s_lo >= 1/4, c >= 11, B* inside the region."""
    st = star_inside_region(cert, RHO, 1, cert["c"])
    return bool(cert["certified"] and cert["s_lo"] >= S_LO_MIN and cert["c"] >= N_DEG and st["ok"]
                and cert["s_hi"] <= 1e9)


def annulus(cert, n1, n3, j, a=A_EXP):
    """Radii of every corner of every square of side <= Q_MAX |A| whose corner A lies in the label's region: the only
    place where the region's certificates look at the twist."""
    T = resonance(n1, n3, j, a)
    r_in = (cert["s_hi"] * T) ** (-a) * (1 - sqrt(2) * Q_MAX)
    r_out = (cert["s_lo"] * T) ** (-a) * (1 + sqrt(2) * Q_MAX)
    return r_in, r_out


DELTA_GRID = (1.5, 1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01)
L_GRID = (0, 2, 5, 10, 20, 30, 50)


def threshold_table():
    rows = []
    for d in DELTA_GRID:
        for L in L_GRID:
            Phi = pi - d
            c = best_c(Phi, L)
            cert = certify_entry(Phi, L, c) if c else dict(certified=False)
            if c:
                cert["c"] = c
            rows.append(dict(delta=d, Phi=Phi, L=L, c=c, cert=cert, degree_applies=bool(c and entry_ok(cert))))
    return rows


def local_scan(phi_fn, slope_fn, table, m_values, n_samples=400):
    """For the labels (m, m, 1), is there a table entry whose bounds hold on that entry's annulus for the label,
    with m >= c? Bounds are checked on a log-spaced sample of the annulus (a diagnostic; exact bounds come from the
    twist's formula)."""
    out = []
    for m in m_values:
        hit = None
        for row in sorted((r for r in table if r["degree_applies"]), key=lambda r: r["c"]):
            if m < row["c"]:
                continue
            r_in, r_out = annulus(row["cert"], m, m, 1)
            rs = np.geomspace(r_in, r_out, n_samples)
            if max(abs(phi_fn(x)) for x in rs) <= row["Phi"] and max(abs(slope_fn(x)) for x in rs) <= row["L"]:
                hit = row
                break
        out.append(dict(m=m, covered=hit is not None, delta=hit["delta"] if hit else None,
                        L=hit["L"] if hit else None, c=hit["c"] if hit else None))
    return out


# ----------------------------------------------------------------------------------------------------------------
# A strongly twisted example
# ----------------------------------------------------------------------------------------------------------------

def twisted_example(label, amp, freq, table):
    """Continue the plain square of label to the twist amp sin(freq log 1/r) (|phi| <= amp, |r phi'| <= amp freq)
    and check it stays inside the region of the table entry covering it."""
    from math import log, sin
    from label_degree_v390 import anchored_jacobian_sign
    from pinned_squares_v391 import continue_in_twist
    phi = lambda r: amp * sin(freq * log(1 / r))
    row = next(r for r in sorted(table, key=lambda r: r["c"] or 1e18)
               if r["degree_applies"] and r["Phi"] >= amp and r["L"] >= amp * freq and min(abs(label[0]), abs(label[1])) >= r["c"])
    cert = row["cert"]
    x0, lab = plain_square(*label)
    path = continue_in_twist(x0, lab, phi, 1, steps=400)
    T = resonance(*label)
    lifts = [x[0] / T for _, x in path]
    qs = [abs(complex(x[1], x[2])) for _, x in path]
    q_in = abs(label[2]) / min(abs(label[0]), abs(label[1])) * cert["side_floor"] / 2
    return dict(label=list(label), amp=amp, freq=freq, entry=dict(delta=row["delta"], L=row["L"], c=row["c"]),
                reached=path[-1][0] == 1.0, s_range=[min(lifts), max(lifts)], walls=[cert["s_lo"], cert["s_hi"]],
                q_range=[min(qs), max(qs)], q_in=q_in,
                inside=bool(cert["s_lo"] < min(lifts) and max(lifts) < cert["s_hi"] and q_in < min(qs)
                            and max(qs) < Q_MAX),
                degree_start=anchored_jacobian_sign(x0, lab, 0.0, lambda r: 0.0)[0],
                degree_end=anchored_jacobian_sign(path[-1][1], lab, 1.0, phi)[0] if path[-1][0] == 1.0 else None)


if __name__ == "__main__":
    import json
    import sys
    from math import log, sin
    from pathlib import Path
    out = Path(__file__).resolve().parents[1] / "results" / "v394"
    out.mkdir(parents=True, exist_ok=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode == "coarse":
        j = int(sys.argv[2])
        (out / f"coarse_{j}.json").write_text(json.dumps(run_coarse(j), indent=1))
    elif mode == "table":
        (out / "table.json").write_text(json.dumps(threshold_table(), indent=1))
    else:
        table = json.loads((out / "table.json").read_text())
        res = dict(coarse={str(j): json.loads((out / f"coarse_{j}.json").read_text()) for j in (1, -1)})
        res["degree_for_all_regions"] = all(c["tail_certified"] and all(p["certified"] for p in c["patterns"].values())
                                            for c in res["coarse"].values())
        res["table"] = [{k: v for k, v in r.items() if k != "cert"} |
                        dict(s_lo=r["cert"].get("s_lo"), s_hi=r["cert"].get("s_hi")) for r in table]
        ms = sorted(set(int(x) for x in np.geomspace(11, 2_000_000, 60)))
        res["scans"] = {
            "contact r^0.2": local_scan(lambda r: pi - 0.5 * r ** 0.2, lambda r: 0.1 * r ** 0.2, table, ms),
            "contact r^0.5": local_scan(lambda r: pi - 0.5 * r ** 0.5, lambda r: 0.25 * r ** 0.5, table, ms),
            "oscillating 2.9 sin(10 log 1/r)": local_scan(lambda r: 2.9 * sin(10 * log(1 / r)), lambda r: 29.0,
                                                           table, ms),
        }
        res["examples"] = [twisted_example((60, 61, 1), 2.0, 4.0, table),
                           twisted_example((-120, -130, 1), 2.6, 2.0, table)]
        (out / "contact_check.json").write_text(json.dumps(res, indent=1))
        print(json.dumps({k: v for k, v in res.items() if k not in ("table", "scans")}, indent=1)[:4000])
        for k, rows in res["scans"].items():
            cov = [r["m"] for r in rows if r["covered"]]
            print(k, len(cov), "of", len(rows), "covered; m range", cov[:1], cov[-1:])
