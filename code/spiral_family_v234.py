"""Interval certificates for squares on the spiral family r^-2 + lambda/r."""

import argparse
import json
import sys
from fractions import Fraction
from math import cos, sin
from pathlib import Path

from sage.all import QQ, RDF, RealIntervalField, matrix, vector

from mixed_power_spiral_v233 import GRID, WITNESSES, PhaseLaw, _sample
from spiral_multicandidate_v232 import scan_sampled_candidates
from spiral_pruned_v231 import _solve4_float, refine_square, square_residual


FIELD = RealIntervalField(100)
SCALE = 2**40
SIDE_FLOOR = QQ(1) / 2
POSITION_TOLERANCE = QQ(1) / 10**9
MANIFEST = Path(__file__).resolve().parents[1] / "research" / "v234-spiral-family.json"


def geometry(u, coefficient, cosine, sine):
    inverse = 1 / u
    theta = inverse**4 + coefficient * inverse**2
    radial = u*u
    c, s = cosine(theta), sine(theta)
    angular_speed = -4 * inverse**3 - 2 * coefficient * inverse
    return ((radial*c, radial*s),
            (2*u*c-angular_speed*s, 2*u*s+angular_speed*c),
            (-s, c))


def derivatives(data, orientation, zero):
    cols = []
    for j, (_, tangent, _) in enumerate(data):
        points = [zero] * 4
        points[j] = tangent
        cols.append(square_residual(points, orientation))
    return tuple(zip(*cols))


def path_acceleration(u, coefficient, speed):
    inv = 1/u
    radial = u*u
    theta = inv**4 + coefficient*inv**2
    rate = inv**2 - speed*(4*inv**5 + 2*coefficient*inv**3)
    acceleration = speed**2*(20*inv**6 + 6*coefficient*inv**4) - 4*speed*inv**3
    real = 2*speed**2 - radial*rate**2
    imag = 4*u*speed*rate + radial*acceleration
    c, s = theta.cos(), theta.sin()
    return real*c-imag*s, real*s+imag*c


def verify_box(lo, hi, orientation, numerators, speed_numerators, width_numerators,
               position_tolerance=POSITION_TOLERANCE):
    """Prove a distinct square above the side floor for all lambda in [lo, hi]."""
    if not (0 <= lo < hi <= 1 and orientation in (-1, 1)
            and len(numerators) == len(speed_numerators) == len(width_numerators) == 4
            and all(type(n) is int for n in (*numerators, *speed_numerators, *width_numerators))
            and all(n > 0 for n in width_numerators) and position_tolerance >= 0):
        raise ValueError("invalid coefficient interval or source box")
    mid = (lo+hi) / 2
    half = (hi-lo) / 2
    center = [QQ(n) / SCALE for n in numerators]
    speed = [QQ(n) / SCALE for n in speed_numerators]
    widths = [QQ(n) / SCALE for n in width_numerators]
    reach = [abs(v)*half+w for v, w in zip(speed, widths)]
    bounds = [(v-w, v+w) for v, w in zip(center, reach)]
    if not all(0 < a < b < 1 for a, b in bounds):
        return "source-box-exits-arm"
    ordered = sorted(bounds)
    if any(b >= c for (_, b), (c, _) in zip(ordered, ordered[1:])):
        return "overlapping-source-intervals"
    coefficient = FIELD(lo, hi)
    delta = FIELD(-half, half)
    remainder = FIELD(0, half**2/2)
    point = [FIELD(u) for u in center]
    path = [u+v*delta for u, v in zip(point, speed)]
    box = [FIELD(a, b) for a, b in bounds]
    trig = (lambda v: v.cos(), lambda v: v.sin())
    zero = (FIELD(0), FIELD(0))
    mid_data = [geometry(u, FIELD(mid), *trig) for u in point]
    acceleration = [path_acceleration(u, coefficient, v) for u, v in zip(path, speed)]
    full = [geometry(u, coefficient, *trig) for u in box]
    nominal = matrix(RDF, [[x.center() for x in row]
                           for row in derivatives(mid_data, orientation, zero)])
    try:
        numerical_inverse = nominal.inverse()
    except (ArithmeticError, ValueError):
        return "singular-nominal-jacobian"
    inverse = matrix(QQ, 4, 4, [QQ(round(float(x)*SCALE)) / SCALE
                                 for x in numerical_inverse.list()])
    if not inverse.det():
        return "singular-preconditioner"
    # Compose the exact linear square residual with the preconditioner before
    # evaluating intervals. Each column depends on just one vertex tangent.
    basis = []
    for j in range(4):
        for axis in ((1, 0), (0, 1)):
            unit = [(0, 0)] * 4
            unit[j] = axis
            basis.append(inverse * vector(QQ, square_residual(unit, orientation)))
    rows = matrix(QQ, basis).transpose().rows()
    for i, width in enumerate(widths):
        row = rows[i]
        weights = tuple((row[2*j], row[2*j+1]) for j in range(4))
        def dot(j, pair):
            x, y = weights[j]
            return x*pair[0] + y*pair[1]
        offset = (sum(dot(j, d[0]) for j, d in enumerate(mid_data))
                  + sum(dot(j, d[2]) + speed[j]*dot(j, d[1])
                        for j, d in enumerate(mid_data))*delta
                  + sum(dot(j, a) for j, a in enumerate(acceleration))*remainder)
        bound = (abs(offset)
                 + sum(abs(int(i == j)-dot(j, d[1]))*widths[j]
                       for j, d in enumerate(full))
                 + position_tolerance*sum(map(abs, row)))
        if not bound.upper().exact_rational() < width:
            return "inclusion-failed"
    a, b = full[0][0], full[1][0]
    side_error = FIELD(-2*position_tolerance, 2*position_tolerance)
    length_sq = sum((y-x+side_error)**2 for x, y in zip(a, b))
    if not length_sq.lower().exact_rational() > SIDE_FLOOR**2:
        return "insufficient-side-bound"
    return "certified"


def _correct(seed, orientation, coefficient):
    root = tuple(seed)
    for _ in range(12):
        if not all(.05 < u < 1 for u in root):
            return None
        data = [geometry(u, coefficient, cos, sin) for u in root]
        residual = square_residual([d[0] for d in data], orientation)
        if max(map(abs, residual)) < 2e-13:
            return root
        step = _solve4_float(derivatives(data, orientation, (0., 0.)),
                             tuple(-v for v in residual))
        if step is None:
            return None
        root = tuple(u+du for u, du in zip(root, step))
    return None


def _candidates(coefficient):
    terms = ((4, 1), (2, Fraction(coefficient))) if coefficient > 0 else ((4, 1),)
    law = PhaseLaw(terms)
    for inner in (.55, .45, .36, .30, .25, .22):
        samples = _sample(law, inner, .22, 400)
        if not samples:
            continue
        scan = scan_sampled_candidates(samples, coefficient, inner, .22, max_pairs=30_000)
        for candidate in scan.distinct_candidates:
            try:
                root, _, _ = refine_square(
                    candidate, lambda u: geometry(u, coefficient, cos, sin)[0],
                    lambda u: geometry(u, coefficient, cos, sin)[1])
            except ValueError:
                continue
            yield candidate["sigma"], root


def _propose(lo, hi, orientation, root):
    mid = float((lo+hi)/2)
    root = _correct(root, orientation, mid)
    if root is None:
        return None, None
    data = [geometry(u, mid, cos, sin) for u in root]
    slope = square_residual([d[2] for d in data], orientation)
    speed = _solve4_float(derivatives(data, orientation, (0., 0.)),
                          tuple(-v for v in slope))
    if speed is None:
        return None, None
    nums = [round(u*SCALE) for u in root]
    speeds = [round(v*SCALE) for v in speed]
    half = float((hi-lo)/2)
    for padding in (1, 4, 16, 64, 256):
        width = max(64, round(padding*half**2*SCALE))
        widths = [width]*4
        if verify_box(lo, hi, orientation, nums, speeds, widths) == "certified":
            return [str(lo), str(hi), orientation, nums, speeds, widths], root
    return None, root


def generate():
    orientation, endpoints = WITNESSES["quadratic-plus-inverse"]
    seed = tuple((a+b)/(2*GRID) for a, b in zip(endpoints[::2], endpoints[1::2]))
    right = QQ(1)
    step = QQ(1)/32
    rows = []
    while right > 0:
        left = max(QQ(0), right-step)
        row, root = _propose(left, right, orientation, seed)
        if row:
            rows.append(row)
            seed = root
            right = left
            if len(rows) % 50 == 0:
                print(f"Certified [{right}, 1] in {len(rows)} boxes", file=sys.stderr, flush=True)
            step = min(QQ(1)/32, step*2)
            continue
        step /= 2
        if step >= QQ(1)/4096:
            continue
        step = QQ(1)/4096
        left = max(QQ(0), right-step)
        print(f"Searching for another square near coefficient {right}", file=sys.stderr, flush=True)
        for sign, candidate in _candidates(float((left+right)/2)):
            row, root = _propose(left, right, sign, candidate)
            if row:
                rows.append(row)
                orientation, seed = sign, root
                right = left
                step *= 2
                break
        else:
            raise ValueError(f"search inconclusive below coefficient {right}")
    return rows[::-1]


def replay(rows):
    through = QQ(0)
    for i, (a, b, sign, center, speed, widths) in enumerate(rows):
        lo, hi = QQ(a), QQ(b)
        if lo != through:
            raise ValueError(f"interval {i} does not start at {through}")
        outcome = verify_box(lo, hi, sign, center, speed, widths)
        if outcome != "certified":
            raise ValueError(f"interval {i} [{lo}, {hi}]: {outcome}")
        through = hi
    if through != 1:
        raise ValueError(f"last certified endpoint is {through}, expected 1")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true", help="find and certify a new cover")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    args = parser.parse_args()
    if args.generate:
        rows = generate()
        replay(rows)
        args.manifest.write_text("[\n" + ",\n".join(json.dumps(row) for row in rows) + "\n]\n",
                                 encoding="utf-8")
    else:
        rows = json.loads(args.manifest.read_text(encoding="utf-8"))
        replay(rows)
    print(f"Certified [0, 1] with {len(rows)} boxes and square side > {SIDE_FLOOR}"
          f" under coordinate perturbations <= {POSITION_TOLERANCE}.")
