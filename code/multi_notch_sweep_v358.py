"""v358 sweeps: certify the robust disk bound on every integer notched rectangle in a box range and on
seeded random rational samples. Writes results/v358/<name>.json.

    python code/multi_notch_sweep_v358.py exhaustive K WMAX HMAX
    python code/multi_notch_sweep_v358.py random K N DENOMINATOR SEED
"""
import collections
import itertools
import json
import random
import sys
import time
from fractions import Fraction as Q
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from multi_notch_v358 import NotchedRectangle, certify  # noqa: E402


def integer_shapes(k, wmax, hmax):
    """Every integer rectangle with k roof cuts, 2 <= H < W <= wmax, H <= hmax, up to reflection."""
    for W in range(2 * k + 1, wmax + 1):
        for H in range(2, min(W - 1, hmax) + 1):
            for pos in itertools.combinations(range(1, W), 2 * k):
                for floors in itertools.product(range(1, H), repeat=k):
                    cuts = tuple((pos[2 * i], pos[2 * i + 1], floors[i]) for i in range(k))
                    mirror = tuple((W - b, W - a, h) for a, b, h in reversed(cuts))
                    if mirror >= cuts:
                        yield W, H, cuts


def rational_shapes(k, n, den, seed):
    rng = random.Random(seed)
    out = []
    while len(out) < n:
        H = Q(rng.randint(3 * den, 20 * den), den)
        W = Q(rng.randint(int(H * den) + 1, 60 * den), den)
        pts = sorted({Q(rng.randint(1, int(W * den) - 1), den) for _ in range(2 * k)})
        if len(pts) < 2 * k:
            continue
        out.append((W, H, tuple((pts[2 * i], pts[2 * i + 1], Q(rng.randint(1, int(H * den) - 1), den))
                                for i in range(k))))
    return out


def _work(shape):
    cert = certify(NotchedRectangle.make(*shape))
    if cert is None:
        return 'MISS', None, shape
    return cert.level, cert.kind, shape


def run(name, shapes):
    start = time.time()
    levels, kinds, misses = collections.Counter(), collections.Counter(), []
    with Pool() as pool:
        for level, kind, shape in pool.imap_unordered(_work, shapes, chunksize=16):
            levels[level] += 1
            kinds[kind] += 1
            if level == 'MISS':
                misses.append([str(x) for x in shape[:2]] + [[str(c) for c in cut] for cut in shape[2]])
    result = {'name': name, 'shapes': sum(levels.values()), 'levels': dict(sorted(levels.items())),
              'kinds': dict(sorted((str(k), v) for k, v in kinds.items())), 'misses': misses,
              'seconds': round(time.time() - start, 1)}
    out = Path(__file__).resolve().parents[1] / 'results' / 'v358'
    out.mkdir(parents=True, exist_ok=True)
    (out / f'{name}.json').write_text(json.dumps(result, indent=1) + '\n')
    print(json.dumps(result))
    return result


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'exhaustive':
        k, wmax, hmax = map(int, sys.argv[2:5])
        run(f'exhaustive_k{k}_W{wmax}_H{hmax}', integer_shapes(k, wmax, hmax))
    elif mode == 'random':
        k, n, den, seed = map(int, sys.argv[2:6])
        run(f'random_k{k}_n{n}_den{den}_seed{seed}', rational_shapes(k, n, den, seed))
    else:
        raise SystemExit(__doc__)
