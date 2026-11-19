"""Certify square sizes with an infinite Takagi phase on the witness arcs."""

import argparse
import json
from math import ceil, cos, floor, sin
from pathlib import Path

from sage.all import QQ

from spiral_family_v234 import (
    MANIFEST as SEED_MANIFEST, POSITION_TOLERANCE, SCALE, SIDE_FLOOR,
    derivatives, geometry, replay, verify_box,
)
from spiral_pruned_v231 import _solve4_float, square_residual


MANIFEST = SEED_MANIFEST.with_name("v235-rough-spiral-family.json")
BASE = 32
TERMS = 6
AMPLITUDE = QQ(1)/100
UNIT_TAIL_BOUND = 1/(2*BASE**TERMS*(1-QQ(1)/BASE))
MIN_U, MAX_U = QQ(1)/3, QQ(5)/6


def _tent(x):
    fraction = x-floor(x)
    return min(fraction, 1-fraction)


def _scaled_tent(u, frequency):
    return _tent(u*frequency)/frequency


def tent_chord_error(a, b, frequency):
    """Exact maximum error against a tent term's rational endpoint chord."""
    start, end = _scaled_tent(a, frequency), _scaled_tent(b, frequency)
    slope = (end-start)/(b-a)
    first, last = ceil(2*frequency*a), floor(2*frequency*b)
    errors = [QQ(0)]
    for parity in (0, 1):
        # Within each parity the tent's value is constant at the knots.
        # The signed error is affine there, so only the extremes are needed.
        low = first+(parity-first)%2
        high = last-(last-parity)%2
        if low <= high:
            knot_value = QQ(parity)/(2*frequency)
            for index in (low, high):
                knot = QQ(index)/(2*frequency)
                errors.append(abs(knot_value-start-slope*(knot-a)))
    return start, end, max(errors)


def unit_phase_enclosure(a, b):
    """Affine chord and phase error for the unit Takagi series on [a,b]."""
    if not MIN_U <= a < b <= MAX_U:
        return None
    rows = [tent_chord_error(a, b, BASE**k) for k in range(TERMS)]
    start = sum(row[0] for row in rows)
    end = sum(row[1] for row in rows)
    slope = (end-start)/(b-a)
    angle_error = sum(row[2] for row in rows) + UNIT_TAIL_BOUND
    return slope, start-slope*a, angle_error


def phase_enclosure(a, b):
    """Affine reference angle and uniform Cartesian error for the full series."""
    model = unit_phase_enclosure(a, b)
    if model is None:
        return None
    slope, intercept, angle_error = model
    return AMPLITUDE*slope, AMPLITUDE*intercept, b*b*AMPLITUDE*angle_error


def _tent_terms(u):
    value, slope = 0., 0.
    for k in range(TERMS):
        frequency = BASE**k
        frac = u*frequency-floor(u*frequency)
        value += min(frac, 1-frac)/frequency
        slope += 1 if frac < .5 else -1
    return value, slope


def _rough_geometry(u, coefficient):
    """A floating-point locator using the truncated phase and one-sided slopes."""
    value, slope = _tent_terms(u)
    amplitude = float(AMPLITUDE)
    return geometry(u, coefficient, cos, sin, amplitude*slope,
                    amplitude*(value-slope*u))


def verify_rough_box(lo, hi, orientation, center, speed, width):
    return verify_box(lo, hi, orientation, center, speed, width,
                      phase_enclosure=phase_enclosure)


def _propose(lo, hi, orientation, seed):
    coefficient = float((lo+hi)/2)
    root = seed
    for _ in range(16):
        if any(not .333 < u < .834 for u in root):
            return None
        data = [_rough_geometry(u, coefficient) for u in root]
        discrepancy = square_residual([d[0] for d in data], orientation)
        jacobian = derivatives(data, orientation, (0., 0.))
        if max(map(abs, discrepancy)) < 2e-13:
            break
        step = _solve4_float(jacobian, tuple(-x for x in discrepancy))
        if step is None:
            return None
        root = tuple(u+du for u, du in zip(root, step))
    else:
        return None
    slope = square_residual([d[2] for d in data], orientation)
    speed = _solve4_float(jacobian, tuple(-x for x in slope))
    if speed is None:
        return None
    center = [round(u*SCALE) for u in root]
    speed_numerators = [round(v*SCALE) for v in speed]
    half = float((hi-lo)/2)
    for width in (max(64, round(padding*half**2*SCALE))
                  for padding in (1, 4, 16, 64, 256, 1024, 4096)):
        widths = [width]*4
        if verify_rough_box(lo, hi, orientation, center, speed_numerators, widths) == "certified":
            return [str(lo), str(hi), orientation, center, speed_numerators, widths]
    return None


def _subdivide(lo, hi, orientation, seed, remaining):
    row = _propose(lo, hi, orientation, seed)
    if row is not None:
        yield row
        return
    if not remaining:
        raise ValueError(f"search inconclusive on [{lo}, {hi}]")
    mid = (lo+hi)/2
    yield from _subdivide(lo, mid, orientation, seed, remaining-1)
    yield from _subdivide(mid, hi, orientation, seed, remaining-1)


def generate():
    for a, b, orientation, center, _, _ in json.loads(SEED_MANIFEST.read_text(encoding="utf-8")):
        seed = tuple(u/SCALE for u in center)
        yield from _subdivide(QQ(a), QQ(b), orientation, seed, 5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", action="store_true")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    args = parser.parse_args()
    if args.generate:
        rows = list(generate())
        replay(rows, verify_rough_box)
        args.manifest.write_text("[\n" + ",\n".join(json.dumps(row) for row in rows) + "\n]\n",
                                 encoding="utf-8")
    else:
        rows = json.loads(args.manifest.read_text(encoding="utf-8"))
        replay(rows, verify_rough_box)
    print(f"Certified [0, 1] in {len(rows)} boxes for the infinite base-{BASE} Takagi phase."
          f" Square side > {SIDE_FLOOR} with extra coordinate error <= {POSITION_TOLERANCE}.")


if __name__ == "__main__":
    main()
