"""v396: any power-law winding rate.

The theorem so far assumed arms winding like r^(-2) (exponent a = 1/2), up to slowly varying angles. Here the exponent
is a parameter: every certificate is run with a in an interval, so one run covers every spiral r = Theta^(-a) with a in
that interval at once. (A constant factor in the winding, theta = C r^(-1/a), is a dilation of the picture by C^a, and
inscribed squares don't care about scale.)

  * Regions (v392, v394, v395): v392's estimates were written for general a and accept an interval a.
  * Degree (v393): the closed-form remainders were specific to a = 1/2. For general a, Taylor's formula with integral
    remainder gives (1 + rho)^(-p) = 1 - p rho + p(p+1)/2 rho^2 K(rho), p = 1/(2a), with K and G = (K - 1)/rho enclosed
    (values and derivatives) by mean values of positive-weight integrals (uniform_degree_v393.d_K, d_G). The limit
    system (sU - 1, sV - 1, +-(sUV - 1)) is the same for every a, because p(p+1) a^2 = (2a+1)/4.
  * Plain-spiral tail and coarse localization (v394): rerun with the interval a.
"""
from mpmath import iv

from uniform_degree_v393 import B_STAR, I, float_zero, krawczyk, localize, patterns

N_DEG = 11
RHO = 2.0


def _split(eps, r, aI):
    """Bisect the relatively widest parameter (eps weighted x8, as in v393; a relative to its own width)."""
    de, dr, da = float(eps.b) - float(eps.a), float(r.b) - float(r.a), (float(aI.b) - float(aI.a)) / float(aI.a)
    if de * 8 >= max(dr, da):
        m = (float(eps.a) + float(eps.b)) / 2
        return [(I(float(eps.a), m), r, aI), (I(m, float(eps.b)), r, aI)]
    if dr >= da:
        m = (float(r.a) + float(r.b)) / 2
        return [(eps, I(float(r.a), m), aI), (eps, I(m, float(r.b)), aI)]
    m = (float(aI.a) + float(aI.b)) / 2
    return [(eps, r, I(float(aI.a), m)), (eps, r, I(m, float(aI.b)))]


def certify_piece(eps, r, aI, sig1, sig3, b_is_min, j, outer, depth=0, max_depth=12,
                  halves=(0.06, 0.04, 0.025, 0.015, 0.01), log=None):
    """v393's fine step with the exponent as a parameter: Krawczyk on a small box around the piece's approximate zero,
    and localization of every zero in `outer` into it."""
    par = (eps, r, sig1, sig3, b_is_min, j, aI)
    mid = (I(float(eps.mid)), I(float(r.mid)), sig1, sig3, b_is_min, j, I(float(aI.mid)))
    z = float_zero(mid)
    for h in halves:
        box = [I(z[k] - h, z[k] + h) for k in range(3)]
        if not all(outer[k].a < box[k].a and box[k].b < outer[k].b for k in range(3)):
            continue
        kin, sgn, _ = krawczyk(box, par)
        if not kin:
            continue
        loc, n = localize(outer, box, par, max_boxes=20000)
        if loc:
            return [dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], a=[float(aI.a), float(aI.b)],
                         half=h, centre=[float(c) for c in z], det_sign=sgn, boxes=n)]
    if depth >= max_depth:
        if log is not None:
            log.append(dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], a=[float(aI.a), float(aI.b)]))
        return None
    out = []
    for e2, r2, a2 in _split(eps, r, aI):
        sub = certify_piece(e2, r2, a2, sig1, sig3, b_is_min, j, outer, depth + 1, max_depth, halves, log)
        if sub is None:
            return None
        out.extend(sub)
    return out


# ----------------------------------------------------------------------------------------------------------------
# The full run: five exponent pieces covering a in [1/4, 4]
# ----------------------------------------------------------------------------------------------------------------

import json
from math import pi, sqrt
from pathlib import Path

from uniform_closure_v392 import Q_MAX, certify_levels, scaled_side_bound, upper_face

A_PIECES = [(0.25, 0.4), (0.4, 0.7), (0.7, 1.2), (1.2, 2.0), (2.0, 4.0)]
# The coarse step and the table use a finer cut of the last piece: over all of [2, 4] the side bound Q0 is inflated by
# the interval dependency on a (6.4 instead of about 4.3 and 4.0), which makes the coarse search region far too big.
COARSE_PIECES = A_PIECES[:-1] + [(2.0, 3.0), (3.0, 4.0)]
N_A = 30                          # degree certificates cover min offset >= 30
S_LO_MIN = 0.25
DELTAS_LO = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
OUT = Path(__file__).resolve().parents[1] / "results" / "v396"


def run_fine(j):
    import time
    star = [I(1 - B_STAR, 1 + B_STAR)] * 3
    res = {}
    for a in A_PIECES:
        key_a = f"{a[0]}-{a[1]}"
        res[key_a] = {}
        for sig1, sig3, b_is_min in patterns(j):
            t0 = time.time()
            log = []
            leaves = certify_piece(I(0, 1 / N_A), I(1 / RHO, 1), I(*a), sig1, sig3, b_is_min, j, star, log=log)
            signs = sorted({lf["det_sign"] for lf in leaves}) if leaves else []
            res[key_a][f"{sig1},{sig3},{b_is_min}"] = dict(
                certified=leaves is not None and len(signs) == 1, det_sign=signs[0] if len(signs) == 1 else None,
                pieces=len(leaves or []), boxes=sum(lf["boxes"] for lf in leaves or []),
                seconds=round(time.time() - t0, 1), failures=log[:3], leaves=leaves or [])
    return res


def plain_tail_a(aI):
    """Zero-twist tail certificate (v394) for every exponent in aI: no plain square beyond scaled lift s_tail."""
    for dhi in (0.3, 0.5, 0.8, 1.2):
        s_tail, ok = upper_face(0.0, 0.0, 1, 1 / N_A, 1 / N_A, dhi, a=aI, uniform=True)
        if ok:
            return s_tail
    return None


def coarse_a(sig1, sig3, b_is_min, j, aI, s_tail, Q0, eps=None, r=None, depth=0, max_depth=6, log=None):
    """Every plain zero with scaled lift in [1/4, s_tail] lies in B*, for every exponent in aI.
    |omega| <= |j| Q0 gives |U| <= (2a+1) Q0 / (a t1), |V| <= (2a+1) Q0 / (a t3) (largest at the smallest a)."""
    eps = eps if eps is not None else I(0, 1 / N_A)
    r = r if r is not None else I(1 / RHO, 1)
    par = (eps, r, sig1, sig3, b_is_min, j, aI)
    t1 = float(r.a) if b_is_min else 1.0
    t3 = 1.0 if b_is_min else float(r.a)
    k = float(((2 * aI + 1) / aI).b) * Q0
    region = [I(S_LO_MIN, s_tail), I(-k / t1, k / t1), I(-k / t3, k / t3)]
    star = [I(1 - B_STAR, 1 + B_STAR)] * 3
    ok, n = localize(region, star, par, max_boxes=150000)
    if ok:
        return [dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], a=[float(aI.a), float(aI.b)],
                     boxes=n)]
    if depth >= max_depth:
        if log is not None:
            log.append(dict(eps=[float(eps.a), float(eps.b)], r=[float(r.a), float(r.b)], a=[float(aI.a), float(aI.b)]))
        return None
    out = []
    for e2, r2, a2 in _split(eps, r, aI):
        sub = coarse_a(sig1, sig3, b_is_min, j, a2, s_tail, Q0, e2, r2, depth + 1, max_depth, log)
        if sub is None:
            return None
        out.extend(sub)
    return out


def run_coarse(j):
    """Resumable: results are saved after every exponent piece and pattern (results/v396/coarse_{j}.json)."""
    import time
    path = OUT / f"coarse_{j}.json"
    res = json.loads(path.read_text()) if path.exists() else {}
    for a in COARSE_PIECES:
        key_a = f"{a[0]}-{a[1]}"
        aI = I(*a)
        d = res.get(key_a)
        if d is None:
            s_tail = plain_tail_a(aI)
            Q0 = scaled_side_bound(1 / N_A, 1 / N_A, S_LO_MIN, 0.0, 0.0, 1, a=aI, uniform=True)
            d = res[key_a] = dict(s_tail=s_tail, Q0=Q0, patterns={})
            path.write_text(json.dumps(res, indent=1))
        if d["s_tail"] is None or d["Q0"] is None:
            continue
        for sig1, sig3, b_is_min in patterns(j):
            key = f"{sig1},{sig3},{b_is_min}"
            if key in d["patterns"]:
                continue
            t0 = time.time()
            log = []
            co = coarse_a(sig1, sig3, b_is_min, j, aI, d["s_tail"], d["Q0"], log=log)
            d["patterns"][key] = dict(certified=co is not None, pieces=len(co or []),
                                      boxes=sum(c["boxes"] for c in co or []),
                                      seconds=round(time.time() - t0, 1), failures=log[:3])
            path.write_text(json.dumps(res, indent=1))
    return res


def star_inside_a(cert, aI, rho=RHO, j=1):
    """Excision containment (as v393's star_inside_region) for every exponent in aI: |w| / eps = c sqrt(t1^2 U^2 +
    t3^2 V^2) with c = a|j|/(2a+1) increasing in a."""
    lo, hi = I(1 - B_STAR), I(1 + B_STAR)
    c = I(abs(j)) * aI / (2 * aI + 1)
    w_lo = I(float(c.a)) * lo * iv.sqrt(I(1 / rho) ** 2 + 1)
    w_hi = I(float(c.b)) * hi * iv.sqrt(I(2))
    ok_s = cert["s_lo"] < float(lo.a) and float(hi.b) < cert["s_hi"]
    ok_in = (I(abs(j)) * cert["side_floor"] / 2).b < w_lo.a
    ok_out = (w_hi / cert["c"]).b < Q_MAX
    return bool(ok_s and ok_in and ok_out)


def best_c_a(Phi, L, aI, c0=N_A, c_cap=200000):
    def ok(c):
        return certify_levels(Phi, L, c, 1, a=aI, deltas_lo=DELTAS_LO)["certified"]
    c = c0
    while not ok(c):
        c *= 2
        if c > c_cap:
            return None
    lo, hi = max(c0 - 1, c // 2), c
    if hi == c0:
        return c0
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi


ENTRIES = [(0.4, 2.0), (1.2, 6.0), (pi - 0.3, 5.0), (pi - 0.1, (5.0, 2.0))]


def run_table():
    rows = []
    for a in COARSE_PIECES:
        aI = I(*a)
        for Phi, L in ENTRIES:
            c = best_c_a(Phi, L, aI)
            cert = certify_levels(Phi, L, c, 1, a=aI, deltas_lo=DELTAS_LO) if c else dict(certified=False)
            if c:
                cert["c"] = c
            ok = bool(c and cert["certified"] and cert["s_lo"] >= S_LO_MIN and star_inside_a(cert, aI))
            rows.append(dict(a=list(a), Phi=Phi, L=list(L) if isinstance(L, tuple) else L, c=c, degree_applies=ok,
                             s_lo=cert.get("s_lo"), s_hi=cert.get("s_hi"), side_floor=cert.get("side_floor")))
    return rows


def validate(fine_by_j, labels, a_values):
    """Plain squares for several exponents: scaled coordinates inside a certified fine box, local degree as certified."""
    from label_degree_v390 import anchored_jacobian_sign
    from uniform_closure_v392 import plain_square, resonance
    rows = []
    for a in a_values:
        key_a = next(f"{p[0]}-{p[1]}" for p in A_PIECES if p[0] <= a <= p[1])
        for n1, n3, j in labels:
            x, lab = plain_square(n1, n3, j, a)
            if x is None:
                rows.append(dict(a=a, label=[n1, n3, j], found=False))
                continue
            m, M = min(abs(n1), abs(n3)), max(abs(n1), abs(n3))
            eps, r = 1 / m, m / M
            b_is_min = abs(n1) == m
            s1, s3 = (1 if n1 > 0 else -1), (1 if n3 > 0 else -1)
            t1, t3 = (r, 1.0) if b_is_min else (1.0, r)
            c = a * abs(j) / (2 * a + 1)
            sc = (x[0] / resonance(n1, n3, j, a), x[1] / (-eps * s1 * c * t1), x[2] / (eps * s3 * c * t3))
            fine = fine_by_j[str(j)][key_a][f"{s1},{s3},{b_is_min}"]
            leaf = next((lf for lf in fine["leaves"] if lf["eps"][0] <= eps <= lf["eps"][1]
                         and lf["r"][0] <= r <= lf["r"][1] and lf["a"][0] <= a <= lf["a"][1]), None)
            inside = leaf is not None and all(abs(v - cc) <= leaf["half"] for v, cc in zip(sc, leaf["centre"]))
            rows.append(dict(a=a, label=[n1, n3, j], found=True, scaled=list(map(float, sc)), in_certified_box=bool(inside),
                             local_degree=anchored_jacobian_sign(x, lab, 0.0, lambda rr: 0.0, a)[0],
                             certified_sign=fine["det_sign"]))
    return rows


VAL_LABELS = [(30, 31, 1), (40, 70, 1), (-35, -50, 1), (45, -60, -1), (-60, 33, -1)]
VAL_A = [0.3, 0.55, 1.0, 1.7, 3.0]


if __name__ == "__main__":
    import sys
    OUT.mkdir(parents=True, exist_ok=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "merge"
    if mode == "fine":
        j = int(sys.argv[2])
        (OUT / f"fine_{j}.json").write_text(json.dumps(run_fine(j), indent=1))
    elif mode == "coarse":
        run_coarse(int(sys.argv[2]))
    elif mode == "table":
        (OUT / "table.json").write_text(json.dumps(run_table(), indent=1))
    else:
        fine = {str(j): json.loads((OUT / f"fine_{j}.json").read_text()) for j in (1, -1)}
        coarse = {str(j): json.loads((OUT / f"coarse_{j}.json").read_text()) for j in (1, -1)}
        table = json.loads((OUT / "table.json").read_text())
        summary = dict(
            fine={j: {ka: {k: {kk: vv for kk, vv in v.items() if kk != "leaves"} for k, v in d.items()}
                      for ka, d in fine[j].items()} for j in fine},
            coarse=coarse, table=table,
            validation=validate(fine, VAL_LABELS, VAL_A))
        def coarse_ok(j, kc):
            d = coarse[j].get(kc)
            return bool(d and d["s_tail"] is not None and len(d["patterns"]) == 4
                        and all(p["certified"] for p in d["patterns"].values()))

        # a fine piece is certified when its fine step holds and every coarse piece inside it is certified
        summary["degree_certified"] = {
            f"{fa[0]}-{fa[1]}": all(fine[j][f"{fa[0]}-{fa[1]}"][k]["certified"] for j in fine
                                    for k in fine[j][f"{fa[0]}-{fa[1]}"])
            and all(coarse_ok(j, f"{ca[0]}-{ca[1]}") for j in coarse for ca in COARSE_PIECES
                    if fa[0] <= ca[0] and ca[1] <= fa[1])
            for fa in A_PIECES}
        (OUT / "exponent_check.json").write_text(json.dumps(summary, indent=1))
        print(json.dumps(summary["degree_certified"]))
        for r in table:
            print(r)
        for v in summary["validation"]:
            print(v)
