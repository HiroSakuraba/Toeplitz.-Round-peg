"""v367: an adversarial refinement against the v366 certificates.

At each stage take the best certificate (largest side floor), find the outer locus endpoint whose
fourth corner D is closest to the boundary, and insert a small triangular bump on the nearest edge
that carries the boundary across D. That flips the side of D at that endpoint, so the certificate's
parity count changes unless something else compensates. New vertices inherit registration times
from the edge they subdivide, so the stretches keep their meaning on the original loop.
"""
from fractions import Fraction as Q

from sweep_search_v366 import RegisteredPolygon, certificate, locus, search


def _nearest_on_edge(a, b, p):
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)
    t = max(Q(0), min(Q(1), t))
    q = (a[0] + t * dx, a[1] + t * dy)
    return t, q, (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2


def weakest_endpoint(poly, cert, bins, offset=Q(0)):
    """The fourth corner D at an outer locus endpoint closest to the boundary."""
    arcs = [(offset + Q(x, bins), Q(cert['span'], bins)) for x in cert['arcs']]
    best = None
    for e0, e1 in locus(poly, *arcs, cert['orient']):
        for P, Qp, R, D, outer in (e0, e1):
            if outer and e0[:3] != e1[:3]:
                c = poly.clearance_sq(D)
                if best is None or c < best[0]:
                    best = (c, D)
    return best


def bump_across(poly, D):
    """Insert a triangular bump on the edge nearest to D whose apex lands beyond D, trying a few
    widths and overshoots. Returns a new RegisteredPolygon, or None if none is simple."""
    for width in (Q(3), Q(2), Q(1), Q(1, 2), Q(1, 4), Q(1, 8)):
        for overshoot in (Q(2), Q(3, 2), Q(5, 4)):
            new = _bump(poly, D, overshoot, width)
            if new is not None:
                return new
    return None


def _bump(poly, D, overshoot, half_width_factor):
    best = None
    for i, (a, b) in enumerate(poly.edges):
        t, q, d2 = _nearest_on_edge(a, b, D)
        if 0 < t < 1 and (best is None or d2 < best[2]):
            best = (i, t, d2, q)
    if best is None:
        return None
    i, t, d2, q = best
    a, b = poly.edges[i]
    apex = (q[0] + overshoot * (D[0] - q[0]), q[1] + overshoot * (D[1] - q[1]))
    # half width along the edge, as a fraction of the edge, proportional to the push
    edge_len2 = (b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2
    from math import isqrt
    ratio = Q(isqrt(int(d2 / edge_len2 * 10**12)) + 1, 10**6) * half_width_factor
    t0, t1 = t - ratio, t + ratio
    if not 0 < t0 < t < t1 < 1:
        return None
    p0 = (a[0] + t0 * (b[0] - a[0]), a[1] + t0 * (b[1] - a[1]))
    p1 = (a[0] + t1 * (b[0] - a[0]), a[1] + t1 * (b[1] - a[1]))
    ta, tb = poly.times[i], poly.times[i + 1]
    V = poly.V[:i + 1] + [p0, apex, p1] + poly.V[i + 1:]
    times = poly.times[:i + 1] + [ta + t0 * (tb - ta), ta + t * (tb - ta), ta + t1 * (tb - ta)] + poly.times[i + 1:]
    try:
        new = RegisteredPolygon(V, times)
    except ValueError:
        return None
    if new.side_of(D) == poly.side_of(D):
        return None
    return new


def run(poly, stages, bins=16, spans=(1, 2), gap=1, policy='best', log=None):
    """Run the adversary for a number of stages.

    policy 'best' attacks the certificate with the largest side floor; policy 'robust' attacks the
    certificate whose smallest endpoint clearance is largest, so it tries to drive the best margin
    over all certificates to zero. Returns one record per stage and the final polygon."""
    records = []
    for stage in range(stages):
        res = search(poly, bins=bins, spans=spans, gap=gap)
        if not res:
            records.append({'stage': stage, 'vertices': len(poly.V), 'certificates': 0})
            break
        margins = [c['clearance_sq'] for c in res]
        if policy == 'best':
            target = res[0]
        elif policy == 'robust':
            target = res[max(range(len(res)), key=lambda i: (margins[i], res[i]['side_sq_floor']))]
        else:
            raise ValueError("unknown policy")
        c2, D = weakest_endpoint(poly, target, bins)
        rec = {'stage': stage, 'vertices': len(poly.V), 'certificates': len(res),
               'best_floor_sq': res[0]['side_sq_floor'], 'best_margin_sq': max(margins),
               'target': (target['arcs'], target['span'], target['orient']),
               'target_floor_sq': target['side_sq_floor'], 'target_clearance_sq': c2,
               'triples': sorted({(c['arcs'], c['span'], c['orient']) for c in res})}
        records.append(rec)
        if log:
            log(rec)
        new = bump_across(poly, D)
        if new is None:
            rec['stopped'] = 'bump would not be simple'
            break
        poly = new
    return records, poly


def record(name, policy, stages, out_dir):
    """Run on the v247 decagon and write a JSON record (rationals as strings)."""
    import json
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    dec = [tuple(map(Q, v)) for v in json.loads((root / 'research/v247-refined-ten-gon.json').read_text())['vertices']]
    recs, poly = run(RegisteredPolygon(dec), stages, policy=policy)
    def enc(x):
        if isinstance(x, Q):
            return str(x)
        if isinstance(x, (list, tuple)):
            return [enc(y) for y in x]
        if isinstance(x, dict):
            return {k: enc(v) for k, v in x.items()}
        return x
    data = {'name': name, 'policy': policy, 'stages': enc(recs),
            'final_vertices': enc(poly.V), 'final_times': enc(poly.times)}
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / f'{name}.json').write_text(json.dumps(data, indent=1) + '\n')
    return data


def record_push(name, triple, orient, steps, out_dir):
    """Greedy parity push on the v247 decagon, with the union certificate before and after."""
    import json
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    dec = [tuple(map(Q, v)) for v in json.loads((root / 'research/v247-refined-ten-gon.json').read_text())['vertices']]
    start = RegisteredPolygon(dec)
    path, costs, poly = greedy_push(start, triple, orient, steps)
    bins = 16
    union = []
    for role in range(3):
        visited = sorted({t[role] for t in path})
        union.append((Q(visited[0], bins), Q(len(visited), bins)))
    before = certificate(start, *union, orient)
    after = certificate(poly, *union, orient)
    data = {'name': name, 'path': [list(t) for t in path], 'costs_sq': [str(c) for c in costs],
            'union': [[str(a), str(b)] for a, b in union],
            'union_before': {'degree': before['degree'], 'clearance_sq': str(before['clearance_sq']),
                             'squares': len(before['squares'])},
            'union_after': {'degree': after['degree'], 'clearance_sq': str(after['clearance_sq']),
                            'squares': len(after['squares'])},
            'final_vertices': [[str(x), str(y)] for x, y in poly.V],
            'final_times': [str(t) for t in poly.times]}
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / f'{name}.json').write_text(json.dumps(data, indent=1) + '\n')
    return data


def push_parity(poly, a_bin, b_bin, e_bin, orient, steps, bins=16, direction=-1):
    """Push the odd parity of the single-bin triple (a, b, e) step by step along the Q stretch by
    bumping across the fourth corner at the locus end on the face of b's bin facing `direction`.
    Records the clearance of the attacked end (the size of the bump needed) at each step."""
    out = []
    for _ in range(steps):
        arcs = [(Q(x, bins), Q(1, bins)) for x in (a_bin, b_bin, e_bin)]
        cert = certificate(poly, *arcs, orient)
        if not cert or cert['degree'] != 1:
            out.append({'b': b_bin, 'status': 'not odd'})
            break
        face = 'B0' if direction < 0 else 'B1'
        best = None
        for e0, e1 in locus(poly, *arcs, orient):
            for P, Qp, R, D, outer in (e0, e1):
                if face in outer and e0[:3] != e1[:3]:
                    c = poly.clearance_sq(D)
                    if best is None or c < best[0]:
                        best = (c, D)
        if best is None:
            out.append({'b': b_bin, 'status': 'no end on that face'})
            break
        new = bump_across(poly, best[1])
        if new is None:
            out.append({'b': b_bin, 'status': 'bump not simple', 'clearance_sq': best[0]})
            break
        nxt = (b_bin + direction) % bins
        moved = certificate(new, (Q(a_bin, bins), Q(1, bins)), (Q(nxt, bins), Q(1, bins)),
                            (Q(e_bin, bins), Q(1, bins)), orient)
        out.append({'b': b_bin, 'clearance_sq': best[0], 'moved': bool(moved and moved['degree'] == 1)})
        poly = new
        if not out[-1]['moved']:
            break
        b_bin = nxt
    return out, poly


def greedy_push(poly, triple, orient, steps, bins=16, gap=1):
    """Move the odd parity of a single-bin triple as far as possible: at each step attack the
    cheapest outer locus end whose face leads to an unvisited single-bin triple that still keeps the
    stretches separated by `gap` empty bins. Returns the path with the cost (clearance) of each move."""
    path = [tuple(triple)]
    costs = []
    for _ in range(steps):
        cur = path[-1]
        arcs = [(Q(x, bins), Q(1, bins)) for x in cur]
        options = []
        for e0, e1 in locus(poly, *arcs, orient):
            for P, Qp, R, D, outer in (e0, e1):
                if not outer or e0[:3] == e1[:3]:
                    continue
                for face in outer:
                    role = 'ABE'.index(face[0])
                    step = -1 if face[1] == '0' else 1
                    nxt = list(cur)
                    nxt[role] = (nxt[role] + step) % bins
                    nxt = tuple(nxt)
                    order = sorted(nxt)
                    sep = all(y - x >= 1 + gap for x, y in zip(order, order[1:] + [order[0] + bins]))
                    if sep and nxt not in path:
                        options.append((poly.clearance_sq(D), D, nxt))
        moved = False
        for c2, D, nxt in sorted(options, key=lambda o: o[0]):
            new = bump_across(poly, D)
            if new is None:
                continue
            cert = certificate(new, *[(Q(x, bins), Q(1, bins)) for x in nxt], orient)
            if cert and cert['degree'] == 1:
                poly = new
                path.append(nxt)
                costs.append(c2)
                moved = True
                break
        if not moved:
            break
    return path, costs, poly


if __name__ == '__main__':
    import sys
    if sys.argv[2] == 'push':
        record_push(sys.argv[1], (7, 3, 10), -1, int(sys.argv[3]), sys.argv[4])
    else:
        record(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])

