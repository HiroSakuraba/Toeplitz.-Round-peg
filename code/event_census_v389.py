"""v389: graceful events along the symmetrizing path of a twisted smooth two-arm point.

Arms: Theta^(-a) e^{i Theta} and Theta^(-a) e^{i(Theta + pi + lam phi(r))}, phi = 0.4 sin(5 log 1/r) (v374's twist),
lam running from 1 (twisted) to 0 (half-turn symmetric). The large parity at cutoff s changes along the path exactly
when a graceful square's side passes s (v373). L(gamma_0, s) = 1 for the symmetric point (v370/v386), so
    L(gamma_1, s) = 1 + #(graceful events at cutoff s)   (mod 2).
Census: squares with side near s are found by random Newton seeds at several lam, each is continued in lam over [0, 1]
with its passes held fixed, and every crossing of side = s is refined and classified.
"""
import numpy as np
from math import log, pi, sin

from smooth_cores_v386 import sector_two_arm, solve_twisted
from twisted_classification_v387 import nearest_pass, twisted_corners
from two_turn_squares_v383 import split_type, window

A_EXP = 0.5


def phi0(r):
    return 0.4 * sin(5 * log(1 / r))


def seeds_at(lam, s, n_seeds, rng, r_range, a=A_EXP, max_ratio=1.0):
    found = []
    for _ in range(n_seeds):
        r = rng.uniform(*r_range)
        Th = r ** (-1 / a)
        arm2 = rng.random() < 0.5
        ang = Th + (pi + lam * phi0(r) if arm2 else 0.0)
        A = r * np.exp(1j * ang)
        v = np.array([A.real, A.imag, rng.uniform(0, 2 * pi), s * rng.uniform(0.6, 1.6)])
        Z = twisted_corners(v)
        th, _, _ = window(Z, [0, 0, 0, 0])
        ks = [nearest_pass(z, t, lam, phi0, a) for z, t in zip(Z, th)]
        res = solve_twisted(v, ks, lam, phi0, a)
        if res is None:
            continue
        w, Zs = res
        if not (0.5 * s < w[3] < 2 * s) or w[3] > max_ratio * abs(sum(Zs) / 4):
            continue
        found.append((lam, w, tuple(ks)))
    return found


def continue_branch(lam0, w0, ks, steps=60, a=A_EXP):
    """Follow a square in lam over [0, 1] from lam0 (passes fixed). Returns sorted list of (lam, w)."""
    out = {round(lam0, 12): w0}
    for direction in (1, -1):
        lam, w = lam0, w0.copy()
        h = 1.0 / steps
        while True:
            nl = lam + direction * h
            if nl < -1e-12 or nl > 1 + 1e-12:
                break
            nl = min(1.0, max(0.0, nl))
            res = solve_twisted(w.copy(), list(ks), nl, phi0, a, it=40)
            if res is None:
                if h > 1.0 / (steps * 64):
                    h /= 2
                    continue
                break
            w2, _ = res
            if abs(w2[3] - w[3]) > 0.2 * w[3]:      # jumped to another solution
                if h > 1.0 / (steps * 64):
                    h /= 2
                    continue
                break
            lam, w = nl, w2
            out[round(lam, 12)] = w
            if lam in (0.0, 1.0):
                break
    return sorted(out.items())


def refine_event(lam_lo, w_lo, lam_hi, w_hi, s, ks, a=A_EXP, it=50):
    """Bisection in lam for side = s along the branch."""
    for _ in range(it):
        lm = (lam_lo + lam_hi) / 2
        res = solve_twisted(w_lo.copy(), list(ks), lm, phi0, a, it=40)
        if res is None:
            break
        wm, _ = res
        if (wm[3] - s) * (w_lo[3] - s) > 0:
            lam_lo, w_lo = lm, wm
        else:
            lam_hi, w_hi = lm, wm
        if lam_hi - lam_lo < 1e-10:
            break
    return (lam_lo + lam_hi) / 2, w_lo, w_hi


def census(s, n_seeds=1500, lams=(0.0, 0.5, 1.0), seed=0, a=A_EXP, steps=60):
    rng = np.random.default_rng(seed)
    r_lo = max(1.5 * s, 0.02)
    cands = []
    for lam in lams:
        cands += seeds_at(lam, s, n_seeds, rng, (r_lo, 0.3), a)
    branches = []          # each: list of (lam, w), ks
    seen = set()
    for lam0, w0, ks in cands:
        key0 = (ks, tuple(np.round(w0, 7)))
        if key0 in seen:
            continue
        br = continue_branch(lam0, w0, ks, steps, a)
        # identify the branch by its state at its smallest lam
        lam_min, wm = br[0]
        key = (ks, round(lam_min, 6), tuple(np.round(wm[:3], 6)))
        if key in seen:
            continue
        seen.add(key)
        for lam, w in br:
            seen.add((ks, tuple(np.round(w, 7))))
        branches.append((br, ks))
    events = []
    for br, ks in branches:
        for (l1, w1), (l2, w2) in zip(br, br[1:]):
            if (w1[3] - s) * (w2[3] - s) < 0:
                le, wa, wb = refine_event(l1, w1, l2, w2, s, ks, a)
                Z = twisted_corners(wa)
                sec = sector_two_arm(list(ks), Z, a)
                Zc = abs(sum(Z) / 4)
                events.append(dict(lam=float(le), sector=sec, pattern=split_type(list(ks)), ks=list(map(int, ks)),
                                   direction=int(np.sign(w2[3] - w1[3]) * np.sign(l2 - l1)),
                                   centre=float(Zc), t_over_Z=float(s / Zc)))
    # dedupe events (same branch found twice through different seeds)
    uniq = []
    for e in sorted(events, key=lambda e: (e["lam"], e["centre"])):
        if not any(abs(e["lam"] - u["lam"]) < 1e-6 and abs(e["centre"] - u["centre"]) < 1e-6 for u in uniq):
            uniq.append(e)
    g = [e for e in uniq if e["sector"] == "graceful"]
    span = [br for br, _ in branches if br[0][0] == 0.0 and br[-1][0] == 1.0]
    return dict(s=s, candidates=len(cands), branches=len(branches), full_branches=len(span), events=uniq,
                graceful_events=len(g), crossed_events=len(uniq) - len(g),
                predicted_parity_at_lam1=int((1 + len(g)) % 2))


def band_census(s0, band=(0.9, 1.6), n_seeds=3000, lams=(0.0, 0.5, 1.0), seed=0, a=A_EXP, steps=40,
                r_hi_factor=1.3):
    """All branches whose side passes through the band [band0 s0, band1 s0] somewhere on the path."""
    from smooth_cores_v386 import r_smooth
    rng = np.random.default_rng(seed)
    lo, hi = band[0] * s0, band[1] * s0
    r_range = (max(1.2 * lo, 0.01), r_hi_factor * r_smooth(hi))
    branches, seen = [], set()
    per_lam = {}
    for lam in lams:
        cands = []
        for _ in range(n_seeds):
            cands += seeds_at(lam, s0 * rng.uniform(*band), 1, rng, r_range, a)
        per_lam[str(lam)] = len(cands)
        for lam0, w0, ks in cands:
            if not (lo <= w0[3] <= hi):
                continue
            if (ks, tuple(np.round(w0, 7))) in seen:
                continue
            br = continue_branch(lam0, w0, ks, steps, a)
            key = (ks, round(br[0][0], 6), tuple(np.round(br[0][1][:3], 6)))
            if key in seen:
                continue
            seen.add(key)
            for _, w in br:
                seen.add((ks, tuple(np.round(w, 7))))
            Z = twisted_corners(br[0][1])
            sides = np.array([w[3] for _, w in br])
            branches.append(dict(ks=list(map(int, ks)), lam_lo=float(br[0][0]), lam_hi=float(br[-1][0]),
                                 t_min=float(sides.min()), t_max=float(sides.max()),
                                 t_at_0=float(sides[0]) if br[0][0] == 0.0 else None,
                                 sector=sector_two_arm(list(ks), Z, a), pattern=split_type(list(ks)),
                                 centre=float(abs(sum(Z) / 4)), monotone=bool(np.all(np.diff(sides) >= -1e-12)
                                                                              or np.all(np.diff(sides) <= 1e-12)),
                                 sides=[float(x) for x in sides], lams=[float(l) for l, _ in br]))
    return dict(s0=s0, band=list(band), r_range=list(r_range), candidates=per_lam, branches=branches)


def events_at(branches, s):
    """Graceful and crossed crossings of side = s along the stored branch samples."""
    g = c = 0
    for b in branches:
        sd = np.array(b["sides"]) - s
        n = int(np.sum(sd[:-1] * sd[1:] < 0))
        if b["sector"] == "graceful":
            g += n
        else:
            c += n
    return g, c


def coverage(branches, s_lo, s_hi, n=4001):
    S = np.linspace(s_lo, s_hi, n)
    ev = [events_at(branches, s) for s in S]
    g = np.array([e[0] for e in ev])
    return dict(s_lo=s_lo, s_hi=s_hi, samples=n, zero_graceful=float(np.mean(g == 0)), odd_graceful=float(np.mean(g % 2 == 1)),
                max_graceful=int(g.max()), parity_odd_fraction=float(np.mean((1 + g) % 2 == 1)))


def half_turn_partner(w, ks, a=A_EXP):
    """Image of a lam = 0 square under the half-turn, polished on the symmetric spiral."""
    Z = [-z for z in twisted_corners(w)]
    th, _, _ = window(Z, [0, 0, 0, 0])
    ks2 = [nearest_pass(z, t, 0.0, phi0, a) for z, t in zip(Z, th)]
    v = np.array([Z[0].real, Z[0].imag, np.angle(Z[1] - Z[0]), abs(Z[1] - Z[0])])
    res = solve_twisted(v, ks2, 0.0, phi0, a)
    if res is None:
        return None
    return res[0], tuple(ks2)


def pair_census(s0, band=(0.9, 1.6), n_seeds=3000, seed=0, a=A_EXP, steps=40, r_hi_factor=1.3):
    """Squares at lam = 0 with side in the band, completed by their half-turn partners, each continued to lam = 1."""
    from smooth_cores_v386 import r_smooth
    rng = np.random.default_rng(seed)
    lo, hi = band[0] * s0, band[1] * s0
    r_range = (max(1.2 * lo, 0.01), r_hi_factor * r_smooth(hi))
    squares = {}
    for _ in range(n_seeds):
        for lam, w, ks in seeds_at(0.0, s0 * rng.uniform(*band), 1, rng, r_range, a):
            if not (lo <= w[3] <= hi):
                continue
            Z = twisted_corners(w)
            key = tuple(np.round(sorted((z.real, z.imag) for z in Z), 8).ravel())
            if key in squares:
                continue
            p = half_turn_partner(w, ks, a)
            if p is None:
                continue
            Zp = twisted_corners(p[0])
            keyp = tuple(np.round(sorted((z.real, z.imag) for z in Zp), 8).ravel())
            squares[key] = (w, ks, keyp)
            squares.setdefault(keyp, (p[0], p[1], key))
    pairs, done = [], set()
    for key, (w, ks, keyp) in squares.items():
        if key in done or keyp not in squares:
            continue
        done.update((key, keyp))
        wp, ksp, _ = squares[keyp]
        legs = []
        for ww, kk in ((w, ks), (wp, ksp)):
            br = continue_branch(0.0, ww, kk, steps, a)
            Z = twisted_corners(ww)
            legs.append(dict(ks=list(map(int, kk)), sector=sector_two_arm(list(kk), Z, a),
                             lam_hi=float(br[-1][0]), sides=[float(x[1][3]) for x in br],
                             lams=[float(x[0]) for x in br], centre=float(abs(sum(Z) / 4))))
        pairs.append(dict(t0=float(w[3]), legs=legs))
    return dict(s0=s0, band=list(band), r_range=list(r_range), squares=len(squares), pairs=pairs)


def pair_events(pairs, s):
    g = c = 0
    for p in pairs:
        for leg in p["legs"]:
            sd = np.array(leg["sides"]) - s
            n = int(np.sum(sd[:-1] * sd[1:] < 0))
            if leg["sector"] == "graceful":
                g += n
            else:
                c += n
    return g, c


def pair_summary(pc, n=4001):
    pairs = pc["pairs"]
    full = [p for p in pairs if all(l["lam_hi"] == 1.0 for l in p["legs"])]
    splits = [abs(p["legs"][0]["sides"][-1] - p["legs"][1]["sides"][-1]) / p["t0"] for p in full]
    moves = [abs(l["sides"][-1] - l["sides"][0]) / p["t0"] for p in full for l in p["legs"]]
    same_sector = sum(1 for p in pairs if p["legs"][0]["sector"] == p["legs"][1]["sector"])
    lo, hi = pc["band"][0] * pc["s0"], pc["band"][1] * pc["s0"]
    S = np.linspace(lo * 1.02, hi * 0.98, n)
    g = np.array([pair_events(full, s)[0] for s in S])
    return dict(pairs=len(pairs), full_pairs=len(full), partners_same_sector=same_sector,
                graceful_pairs=sum(1 for p in full if p["legs"][0]["sector"] == "graceful"),
                split_median=float(np.median(splits)) if splits else None,
                split_max=float(max(splits)) if splits else None,
                move_median=float(np.median(moves)) if moves else None, move_max=float(max(moves)) if moves else None,
                zero_graceful_fraction=float(np.mean(g == 0)), odd_fraction=float(np.mean(g % 2 == 1)),
                max_graceful_events=int(g.max()) if len(g) else 0)


def run_checks(path=None, n_seeds=1500):
    import json
    res = {}
    for s0 in (0.01, 0.02):
        A = pair_census(s0, n_seeds=n_seeds, seed=11)
        B = pair_census(s0, n_seeds=n_seeds, seed=12)
        ka = {round(p["t0"], 9) for p in A["pairs"]}
        kb = {round(p["t0"], 9) for p in B["pairs"]}
        ov = len(ka & kb)
        both = dict(A, pairs=A["pairs"] + [p for p in B["pairs"] if round(p["t0"], 9) not in ka])
        sA, sB, sAB = pair_summary(A), pair_summary(B), pair_summary(both)
        full = [p for p in both["pairs"] if all(l["lam_hi"] == 1.0 for l in p["legs"])]
        moves = [abs(l["sides"][-1] - l["sides"][0]) for p in full for l in p["legs"] if l["sector"] == "graceful"]
        width = (both["band"][1] - both["band"][0]) * s0
        N_est = len(ka) * len(kb) / ov if ov else None
        lower = len(ka) * len(kb) / max(1, ov)          # with no overlap: a lower bound in expectation
        mean_events = 2 * lower * (sAB["graceful_pairs"] / max(1, sAB["full_pairs"])) * float(np.mean(moves)) / width \
            if moves else None
        res[str(s0)] = dict(run_A=sA, run_B=sB, combined=sAB, overlap=ov, pairs_A=len(ka), pairs_B=len(kb),
                            lincoln_petersen=N_est, population_lower_bound=lower,
                            mean_graceful_move=float(np.mean(moves)) if moves else None,
                            estimated_events_per_cutoff=mean_events,
                            example_pairs=[dict(t0=p["t0"], legs=[dict(sector=l["sector"], t_end=l["sides"][-1],
                                                                       lam_hi=l["lam_hi"], centre=l["centre"])
                                                                  for l in p["legs"]]) for p in full[:20]])
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import json
    from pathlib import Path
    r = run_checks(Path(__file__).resolve().parents[1] / "results" / "v389" / "event_census.json")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "example_pairs"} for k, v in r.items()}, indent=1))
