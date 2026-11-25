"""v384: wild points accumulating at a limit point.

Test curves. The bottom edge of the rectangle [-1, 1] x [0, 1] runs along the x-axis through the limit point
p_inf = 0. Double spirals are spliced into it at p_k = x_k (k = k0, k0+1, ...), x_k -> 0 from the right. Each one
is a scaled copy of the v370 double spiral: the edge enters its disk at x_k - rho_k, arm one winds in from angle
pi (radius rho_k), a chord crosses the centre, and arm two winds out to x_k + rho_k. With arm two the exact
half-turn image of arm one, the whole local arc (edge pieces included) is invariant under the half-turn about
p_k. A twist turns arm two by pi + phi(r) instead (as in v374), which breaks the symmetry at every scale.

Models:
    sparse        x_k = 2^-k, rho_k = 4^-k            (each point isolated at its own scale)
    self-similar  x_k = 2^-k, rho_k = 2^-k / 8         (same picture at every scale)

A drawing keeps the spirals k0..K, each to a finite depth, and flattens the rest (the edge runs straight
through them and through p_inf). All vertices are rational.
"""
from fractions import Fraction as Q
from math import hypot, log, pi, sin, sqrt, tan

import numpy as np

from asymmetric_spiral_v374 import assert_jordan
from wild_point_parity_v370 import DEN, _unit, sector_counts

A_EXP = 0.5
THETA0 = 3 * pi          # arm one starts on the negative x-axis


def model(name, k):
    x = Q(1, 2 ** k)
    if name == "sparse":
        return x, Q(1, 4 ** k)
    if name == "self-similar":
        return x, x / 8
    raise ValueError(name)


def _round(v):
    return Q(round(v * DEN), DEN)


def spiral_piece(cx, rho, half_turns, per_turn=24, twist=None, a=A_EXP):
    """Vertices of one spliced double spiral, from (cx - rho0, 0) to (cx + rho0', 0), centre (cx, 0).

    rho0 = rounded rho; arm one: radius rho (theta/THETA0)^-a at angle theta, theta in [THETA0, THETA0 + n pi].
    twist: None (half-turn symmetric) or a function phi(r_rel) of the relative radius r/rho."""
    th1 = THETA0 + half_turns * pi
    steps = max(1, round((th1 - THETA0) * per_turn / (2 * pi)))
    arm1 = []
    for j in range(steps + 1):
        th = THETA0 + (th1 - THETA0) * j / steps
        r = _round(float(rho) * (th / THETA0) ** (-a))
        ux, uy = _unit(th)
        arm1.append((r * ux, r * uy))
    arm2 = []
    for x, y in reversed(arm1):
        if twist is None:
            arm2.append((-x, -y))
        else:
            rr = hypot(float(x), float(y)) / float(rho)
            u = Q(round(tan(twist(rr) / 2) * DEN), DEN)
            c, sn = (1 - u * u) / (1 + u * u), 2 * u / (1 + u * u)
            arm2.append((-(c * x - sn * y), -(sn * x + c * y)))
    pts = arm1 + arm2
    return [(cx + x, y) for x, y in pts]


def default_twist(r_rel):
    return 0.4 * sin(5 * log(1 / r_rel))


def drawing(name, k0, K, half_turns, per_turn=24, twist=None, kept=None):
    """Rectangle with spirals k0..K spliced into its bottom edge (or only those in `kept`), counterclockwise.

    half_turns: int or dict k -> depth."""
    ks = [k for k in range(K, k0 - 1, -1) if kept is None or k in kept]      # left to right: small x first
    V = [(Q(-1), Q(0))]
    for k in ks:
        cx, rho = model(name, k)
        n = half_turns[k] if isinstance(half_turns, dict) else half_turns
        V += spiral_piece(cx, rho, n, per_turn, twist)
    V += [(Q(1), Q(0)), (Q(1), Q(1)), (Q(-1), Q(1))]
    return V


def large_parity(V, s):
    assert_jordan(V)
    counts, degenerate, squares = sector_counts(V, s)
    if degenerate:
        raise ArithmeticError("degenerate small-square candidate")
    return (1 + counts["aligned"] + counts["reversed"]) % 2, counts, squares


def square_home(entry, name, ks):
    """Which spiral disks (index k, distance measured in units of rho_k) the square's corners are nearest to."""
    pts = entry[1]
    out = []
    for p in pts:
        x, y = float(p[0]), float(p[1])
        best = min(ks, key=lambda k: hypot(x - float(model(name, k)[0]), y) / float(model(name, k)[1]))
        out.append(best)
    return tuple(out)


# ---------------------------------------------------------------- the flattening lemma, made explicit

STRIP = 1 / (2 * sqrt(2))      # eta < s / (2 sqrt 2): no three corners of a side-s square fit in the strip


def strip_flattenable(name, ks, s, eta_frac=0.999 * STRIP):
    """Spirals (among ks) that Lemma 384.1 lets us flatten at cutoff s, all at once.

    eta = eta_frac * s. Candidates: disks inside the strip |y| <= eta (rho_k <= eta). Each candidate disk must be
    more than sqrt2 * s from every part of the curve outside the strip: the other spirals' disks with rho > eta,
    and the rectangle's vertical sides (distance 1 - x_k - rho_k). Returns (flattenable, kept)."""
    s = float(s)
    eta = eta_frac * s
    big = [k for k in ks if float(model(name, k)[1]) > eta]
    flat, kept = [], list(big)
    for k in ks:
        if k in big:
            continue
        x, rho = map(float, model(name, k))
        ok = 1 - x - rho > sqrt(2) * s and x + 1 - rho > sqrt(2) * s
        for b in big:
            xb, rb = map(float, model(name, b))
            if abs(x - xb) - rho - rb <= sqrt(2) * s:
                ok = False
        (flat if ok else kept).append(k)
    return sorted(flat), sorted(kept)


def separated(name, ks, s):
    """Additivity hypothesis: kept disks pairwise more than sqrt2 * s apart and from the rectangle's sides."""
    s = float(s)
    for i, k in enumerate(ks):
        x, rho = map(float, model(name, k))
        if 1 - x - rho <= sqrt(2) * s:
            return False
        for m in ks[i + 1:]:
            y, r2 = map(float, model(name, m))
            if abs(x - y) - rho - r2 <= sqrt(2) * s:
                return False
    return True


def domain_counts(squares, V, name, ks):
    """Sector counts of small squares split by the spiral disk nearest (in units of rho) to their corners.
    A square whose corners are attributed to two different spirals is reported under 'straddling'."""
    from wild_point_parity_v370 import orientation, sector
    o = orientation(V)
    out = {}
    for entry in squares.values():
        home = set(square_home(entry, name, ks))
        key = home.pop() if len(home) == 1 else "straddling"
        d = out.setdefault(key, {"aligned": 0, "reversed": 0, "crossed": 0})
        d[sector(entry, len(V), o)] += 1
    return out


# ---------------------------------------------------------------- experiments

def parity_record(name, ks, s, depth, per_turn, twist):
    V = drawing(name, min(ks), max(ks), depth, per_turn, twist, kept=set(ks))
    p, c, sq = large_parity(V, s)
    return dict(ks=list(ks), parity=p, counts=c, domains={str(k): v for k, v in domain_counts(sq, V, name, ks).items()})


def additivity_trial(name, s, k_all, depth=3, per_turn=12, twist=None, flat_depth=2):
    """Full drawing (all of k_all) versus the strip-kept drawing versus each kept spiral alone.

    Kept spirals are drawn to `depth` half-turns, flattenable ones to `flat_depth` (their content is irrelevant
    by Lemma 384.1, and shallow copies keep the exact count cheap)."""
    flat, kept = strip_flattenable(name, list(k_all), s)
    depth = {k: (depth if k in kept else flat_depth) for k in k_all}
    full = parity_record(name, k_all, s, depth, per_turn, twist)
    reduced = parity_record(name, kept, s, depth, per_turn, twist)
    singles = {str(k): parity_record(name, [k], s, depth, per_turn, twist) for k in kept}
    pred = (1 + sum(1 + v["parity"] for v in singles.values())) % 2
    return dict(model=name, s=str(s), twisted=twist is not None, flattenable=flat, kept=kept,
                separated=separated(name, kept, s), full=full, kept_only=reduced, singles=singles,
                additive_prediction=pred)


def cluster_split(name, ks, s, sep=sqrt(2)):
    """Split spirals into the 'cluster' at the limit point and 'big' isolated ones.

    Walking outward from the limit point, a spiral joins the cluster if it fits in the strip (rho <= s/(2 sqrt 2))
    or its disk comes within sep*s of a cluster disk; the first spiral that does neither, and every spiral beyond
    it, is big. Big disks are checked to be pairwise separated by `separated`."""
    s = float(s)
    cluster, big = [], []
    for k in sorted(ks, reverse=True):              # smallest x first
        x, rho = map(float, model(name, k))
        near = any(abs(x - float(model(name, m)[0])) - rho - float(model(name, m)[1]) <= sep * s for m in cluster)
        if not big and (not cluster or near or rho <= STRIP * s):
            cluster.append(k)
        else:
            big.append(k)
    return sorted(big), sorted(cluster)


def cluster_trial(name, s, k_all, depth=2, per_turn=12, twist=None):
    """Self-similar test: full drawing versus (big spirals alone) + (the cluster at the limit point alone)."""
    big, cluster = cluster_split(name, list(k_all), s)
    full = parity_record(name, k_all, s, depth, per_turn, twist)
    clus = parity_record(name, cluster, s, depth, per_turn, twist)
    singles = {str(k): parity_record(name, [k], s, depth, per_turn, twist) for k in big}
    pred = (1 + (1 + clus["parity"]) + sum(1 + v["parity"] for v in singles.values())) % 2
    # interaction probe: drop the largest spiral of the cluster that the strip lemma cannot flatten
    flat, _ = strip_flattenable(name, cluster, s)
    hard = [k for k in cluster if k not in flat]
    probe = None
    if len(hard) >= 2:
        drop = max(hard)                     # smallest of the unflattenable ones
        rest = [k for k in cluster if k != drop]
        probe = dict(dropped=drop, parity_without=parity_record(name, rest, s, depth, per_turn, twist)["parity"])
    return dict(model=name, s=str(s), twisted=twist is not None, big=big, cluster=cluster, flattenable=flat,
                full=full, cluster_only=clus, singles=singles, additive_prediction=pred, probe=probe)


def sparse_range_check(s_values, k0=2):
    """Corollary 384.6 hypotheses on the sparse model: at each cutoff, the spirals the strip lemma cannot
    flatten are pairwise separated by more than sqrt2*s (and from the rectangle's sides)."""
    bad = []
    for s in s_values:
        kmax = k0 + 4 + int(np.ceil(np.log2(1 / s)))      # far below the strip size
        flat, kept = strip_flattenable("sparse", list(range(k0, kmax + 1)), s)
        if not separated("sparse", kept, s) or max(flat, default=kmax) < kmax:
            bad.append(float(s))
    return bad


def interaction_trial(name, s, k_all, depth=2, per_turn=12, twist=None):
    """The cluster at the limit point with its strip-flattenable members removed (Lemma 384.2), against the sum of
    its members alone. A difference, or a square whose corners sit near two different spirals, is an interaction."""
    big, cluster = cluster_split(name, list(k_all), s)
    flat, hard = strip_flattenable(name, cluster, s)
    if not hard:      # everything in the cluster flattens: the plain rectangle, L = 1 (Lemma 370.1)
        joint = dict(ks=[], parity=1, counts={"aligned": 0, "reversed": 0, "crossed": 0}, domains={})
    else:
        joint = parity_record(name, hard, s, depth, per_turn, twist)
    singles = {str(k): parity_record(name, [k], s, depth, per_turn, twist) for k in hard}
    pred = (1 + sum(1 + v["parity"] for v in singles.values())) % 2
    return dict(model=name, s=str(s), twisted=twist is not None, kind="interaction", big=big, cluster=cluster,
                flattenable=flat, hard=hard, joint=joint, singles=singles, additive_prediction=pred,
                straddling=joint["domains"].get("straddling"))


# ---------------------------------------------------------------- resumable runner

SPARSE_CUTS = ["1/60", "1/120", "1/240", "1/480"]
SELF_CUTS = ["1/64", "1/76", "1/91", "1/108", "1/128"]


def run_checks(path, which=("sparse", "self-similar")):
    """Run every trial not yet stored in `path` (JSON list), saving after each one."""
    import json
    from pathlib import Path
    path = Path(path)
    done = json.loads(path.read_text()) if path.exists() else []
    have = {(r["model"], r["s"], r["twisted"], r.get("kind", "trial")) for r in done}
    jobs = []
    if "sparse" in which:
        jobs += [("sparse", c, tw, "trial") for tw in (True, False) for c in SPARSE_CUTS]
    if "self-similar" in which:
        jobs += [("self-similar", "1/64", False, "trial")]
        jobs += [("self-similar", c, tw, "interaction") for tw in (False, True) for c in SELF_CUTS]
    for job in jobs:
        name, c, tw = job[:3]
        kind = job[3] if len(job) > 3 else "trial"
        if (name, c, tw, kind) in have:
            continue
        twist = default_twist if tw else None
        if name == "sparse":
            r = additivity_trial(name, Q(c), range(2, 8), twist=twist)
        elif kind == "trial":
            r = cluster_trial(name, Q(c), range(2, 9), twist=twist)
        else:
            r = interaction_trial(name, Q(c), range(2, 9), twist=twist)
        done.append(r)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(done, indent=1))
        print(name, c, "twisted" if tw else "symmetric", "parity", (r.get("full") or r["joint"])["parity"],
              "prediction", r["additive_prediction"], flush=True)
    return done


if __name__ == "__main__":
    from pathlib import Path
    run_checks(Path(__file__).resolve().parents[1] / "results" / "v384" / "trials.json")
