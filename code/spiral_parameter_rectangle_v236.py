"""Verify a square-size bound over a rectangle of rough spiral parameters."""

import argparse
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor
import gzip
import json
from math import cos, sin
from multiprocessing import get_context
from pathlib import Path
import sys

from sage.all import QQ

from rough_spiral_family_v235 import _tent_terms, unit_phase_enclosure
from spiral_family_v234 import (
    FIELD, MANIFEST as SEED_MANIFEST, POSITION_TOLERANCE, SCALE, SIDE_FLOOR,
    derivatives, geometry, preconditioned_rows,
)
from spiral_pruned_v231 import _solve4_float, square_residual


MANIFEST = SEED_MANIFEST.with_name("v236-spiral-parameter-rectangle.json.gz")
MAX_AMPLITUDE = QQ(1)/100
DOMAIN = ((QQ(0), QQ(1)), (-MAX_AMPLITUDE, MAX_AMPLITUDE))
MAX_DEPTH = 19


def _sum_product(left, right):
    return sum(a*b for a, b in zip(left, right))


def _ray_acceleration(u, coefficient, amplitude, source_speed, deltas, chord):
    """Second derivative along the straight parameter ray from its midpoint."""
    lam_speed, amp_speed = deltas
    slope, intercept = chord
    inv = 1/u
    radius = u**2
    phase = inv**4 + coefficient*inv**2 + amplitude*(slope*u+intercept)
    phase_u = -4*inv**5 - 2*coefficient*inv**3 + amplitude*slope
    rate = phase_u*source_speed + inv**2*lam_speed + (slope*u+intercept)*amp_speed
    accel = ((20*inv**6+6*coefficient*inv**4)*source_speed**2
             + 2*(-2*inv**3*lam_speed+slope*amp_speed)*source_speed)
    real = 2*source_speed**2 - radius*rate**2
    imag = 4*u*source_speed*rate + radius*accel
    c, s = phase.cos(), phase.sin()
    return real*c-imag*s, real*s+imag*c


def verify_tile(bounds, orientation, center_nums, lambda_nums, amplitude_nums, width_num):
    """Prove an inscription above SIDE_FLOOR for every parameter pair in a tile."""
    if not (len(bounds) == 2 and orientation in (-1, 1)
            and all(len(data) == 4 for data in (center_nums, lambda_nums, amplitude_nums))
            and all(type(n) is int for n in (*center_nums, *lambda_nums, *amplitude_nums, width_num))
            and width_num > 0
            and all(c <= a < b <= d for (a, b), (c, d) in zip(bounds, DOMAIN))):
        raise ValueError("invalid tile or source box")

    midpoint = [(a+b)/2 for a, b in bounds]
    half = [(b-a)/2 for a, b in bounds]
    coefficient, amplitude = [FIELD(a, b) for a, b in bounds]
    deltas = [FIELD(-h, h) for h in half]
    center = [QQ(n)/SCALE for n in center_nums]
    speeds = [(QQ(l)/SCALE, QQ(a)/SCALE) for l, a in zip(lambda_nums, amplitude_nums)]
    width = QQ(width_num)/SCALE
    radius = [_sum_product(map(abs, velocity), half) for velocity in speeds]
    box = [(u-r-width, u+r+width) for u, r in zip(center, radius)]
    models = [unit_phase_enclosure(a, b) for a, b in box]
    if any(model is None for model in models):
        return "phase-outside-certified-domain"
    ordered = sorted(box)
    if any(b >= c for (_, b), (c, _) in zip(ordered, ordered[1:])):
        return "overlapping-source-intervals"

    max_amp = max(map(abs, bounds[1]))
    errors = [POSITION_TOLERANCE + max_amp*b*b*model[2]
              for (_, b), model in zip(box, models)]
    trig = (lambda x: x.cos(), lambda x: x.sin())
    mid_lambda, mid_amp = map(FIELD, midpoint)
    mid_data, full_data, along_rays = [], [], []
    amp_partials = []
    for u, r, velocity, (a, b), (slope, intercept, _) in zip(center, radius, speeds, box, models):
        data = geometry(FIELD(u), mid_lambda, *trig, mid_amp*slope, mid_amp*intercept)
        mid_data.append(data)
        amp_partials.append(tuple(x*(u*u*(slope*u+intercept)) for x in data[2]))
        full_data.append(geometry(FIELD(a, b), coefficient, *trig,
                                  amplitude*slope, amplitude*intercept))
        along_rays.append(_ray_acceleration(
            FIELD(u-r, u+r), coefficient, amplitude, _sum_product(velocity, deltas),
            deltas, (slope, intercept)))

    rows = preconditioned_rows(mid_data, orientation)
    if rows is None:
        return "singular-preconditioner"
    half_unit = FIELD(0, QQ(1)/2)
    for i, row in enumerate(rows):
        def dot(j, pair):
            return row[2*j]*pair[0] + row[2*j+1]*pair[1]

        value = sum(dot(j, d[0]) for j, d in enumerate(mid_data))
        slope_lam = sum(dot(j, d[2]) + speeds[j][0]*dot(j, d[1])
                        for j, d in enumerate(mid_data))
        slope_amp = sum(dot(j, amp_partials[j]) + speeds[j][1]*dot(j, d[1])
                        for j, d in enumerate(mid_data))
        offset = (value + slope_lam*deltas[0] + slope_amp*deltas[1]
                  + half_unit*sum(dot(j, a) for j, a in enumerate(along_rays)))
        bound = (abs(offset) + width*sum(abs(int(i == j)-dot(j, d[1]))
                                        for j, d in enumerate(full_data))
                 + sum(abs(weight)*errors[j//2] for j, weight in enumerate(row)))
        if not bound.upper().exact_rational() < width:
            return "inclusion-failed"
    a, b = full_data[0][0], full_data[1][0]
    side_error = FIELD(-errors[0]-errors[1], errors[0]+errors[1])
    squared = sum((y-x+side_error)**2 for x, y in zip(a, b))
    if not squared.lower().exact_rational() > SIDE_FLOOR**2:
        return "insufficient-side-bound"
    return "certified"


def _locate(bounds, orientation, seed):
    coefficient, amplitude = [float((a+b)/2) for a, b in bounds]
    root = seed
    for _ in range(16):
        if any(not .333 < u < .834 for u in root):
            return None
        terms = [_tent_terms(u) for u in root]
        data = [geometry(u, coefficient, cos, sin, amplitude*slope,
                         amplitude*(value-slope*u)) for u, (value, slope) in zip(root, terms)]
        residual = square_residual([d[0] for d in data], orientation)
        jacobian = derivatives(data, orientation, (0., 0.))
        if max(map(abs, residual)) < 2e-13:
            break
        correction = _solve4_float(jacobian, tuple(-x for x in residual))
        if correction is None:
            return None
        root = tuple(u+c for u, c in zip(root, correction))
    else:
        return None
    lambda_residual = square_residual([d[2] for d in data], orientation)
    amp_partials = [tuple(x*u*u*term[0] for x in d[2]) for u, term, d in zip(root, terms, data)]
    amp_residual = square_residual(amp_partials, orientation)
    lambda_speed = _solve4_float(jacobian, tuple(-x for x in lambda_residual))
    amp_speed = _solve4_float(jacobian, tuple(-x for x in amp_residual))
    if lambda_speed is None or amp_speed is None:
        return None
    return [[round(x*SCALE) for x in nums] for nums in (root, lambda_speed, amp_speed)]


def _split(bounds, axis):
    lo, hi = bounds[axis]
    mid = (lo+hi)/2
    return (bounds[:axis]+((lo, mid),)+bounds[axis+1:],
            bounds[:axis]+((mid, hi),)+bounds[axis+1:])


def _build_slab(bounds, start_depth):
    seeds = json.loads(SEED_MANIFEST.read_text(encoding="utf-8"))
    seed_starts = [QQ(row[0]) for row in seeds]

    def visit(bounds, depth):
        h_lam, h_amp = [float((b-a)/2) for a, b in bounds]
        if h_lam <= 1/256 and h_amp <= .005:
            center_lam = sum(bounds[0])/2
            seed_row = seeds[bisect_right(seed_starts, center_lam)-1]
            orientation = seed_row[2]
            data = _locate(bounds, orientation, tuple(x/SCALE for x in seed_row[3]))
            if data is not None:
                variation = h_lam*h_lam + (h_amp/4)**2
                for padding in (1, 4, 16, 64, 256, 1024):
                    width = max(64, round(padding*variation*SCALE))
                    if verify_tile(bounds, orientation, *data, width) == "certified":
                        return [orientation, *data, width]
        if depth == MAX_DEPTH:
            raise ValueError(f"search inconclusive for parameter tile {bounds}")
        axis = 0 if h_lam >= h_amp/4 else 1
        left, right = _split(bounds, axis)
        return [axis, visit(left, depth+1), visit(right, depth+1)]

    return visit(bounds, start_depth)


def _build_subtree(task):
    return _build_slab(*task)


def generate(jobs=1):
    depth = 5
    slabs = [DOMAIN]
    for _ in range(depth):
        slabs = [child for slab in slabs for child in _split(slab, 0)]
    tasks = [(slab, depth) for slab in slabs]

    def collect(results):
        trees = []
        for i, result in enumerate(results, 1):
            trees.append(result)
            print(f"Certified coefficient slab {i}/{len(slabs)}", file=sys.stderr, flush=True)
        return trees

    if jobs == 1:
        trees = collect(map(_build_subtree, tasks))
    else:
        with ProcessPoolExecutor(max_workers=jobs, mp_context=get_context("spawn")) as pool:
            trees = collect(pool.map(_build_subtree, tasks))
    while len(trees) > 1:
        trees = [[0, a, b] for a, b in zip(trees[::2], trees[1::2])]
    return trees[0]


def replay(tree, bounds=DOMAIN, depth=0):
    if depth > MAX_DEPTH or not isinstance(tree, list):
        raise ValueError("invalid certificate tree")
    if len(tree) == 3 and type(tree[0]) is int and tree[0] in (0, 1):
        axis, left, right = tree
        a, b = _split(bounds, axis)
        return replay(left, a, depth+1) + replay(right, b, depth+1)
    if len(tree) != 5:
        raise ValueError("invalid certificate leaf")
    outcome = verify_tile(bounds, *tree)
    if outcome != "certified":
        raise ValueError(f"tile {bounds}: {outcome}")
    return 1


def _format_tree(tree, depth=0):
    indent = " " * depth
    if len(tree) == 3:
        axis, left, right = tree
        return (f"{indent}[{axis},\n{_format_tree(left, depth+1)},\n"
                f"{_format_tree(right, depth+1)}\n{indent}]")
    return indent + json.dumps(tree, separators=(",", ":"))


def load_tree(path):
    data = path.read_bytes()
    if path.suffix == ".gz":
        data = gzip.decompress(data)
    return json.loads(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    tree = generate(args.jobs) if args.generate else load_tree(args.manifest)
    count = replay(tree)
    if args.generate:
        data = (_format_tree(tree)+"\n").encode("utf-8")
        if args.manifest.suffix == ".gz":
            data = gzip.compress(data, mtime=0)
        args.manifest.write_bytes(data)
    print(f"Certified [0,1] x [{-MAX_AMPLITUDE},{MAX_AMPLITUDE}] in {count} tiles."
          f" Square side > {SIDE_FLOOR} with coordinate error <= {POSITION_TOLERANCE}.")


if __name__ == "__main__":
    main()
