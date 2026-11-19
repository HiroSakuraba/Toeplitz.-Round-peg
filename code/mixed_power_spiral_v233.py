"""Blind square discovery and exact replay for mixed inverse-power phases."""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as Q
from math import ceil, cos, isfinite, sin, sqrt

from spiral_multicandidate_v232 import scan_sampled_candidates
from spiral_pruned_v231 import (
    Sample, _R_matrix, _inverse_and_root, _round_fraction,
    _trig_rational_enclosure, _side_floor_holds, refine_square,
)


GRID = 2**30


@dataclass(frozen=True)
class PhaseLaw:
    """Phase in u = sqrt(radius): sum of c*u**(-n), with positive c."""

    terms: tuple[tuple[int, Q], ...]

    def __post_init__(self):
        terms = tuple((n, Q(c)) for n, c in self.terms)
        if (not terms or any(type(n) is not int or n < 1 or c <= 0 for n, c in terms)
                or len({n for n, _ in terms}) != len(terms)):
            raise ValueError("expected distinct positive integer powers and positive coefficients")
        object.__setattr__(self, "terms", terms)

    def phase(self, u):
        return sum(c * u**(-n) for n, c in self.terms)

    def phase_slope(self, u):
        return -sum(n * c * u**(-n-1) for n, c in self.terms)

    def point(self, u):
        theta = self.phase(u)
        return u*u*cos(theta), u*u*sin(theta)

    def tangent(self, u):
        theta = self.phase(u)
        speed = u*u*self.phase_slope(u)
        return 2*u*cos(theta)-speed*sin(theta), 2*u*sin(theta)+speed*cos(theta)

    def chord_error(self, a, b):
        """Coordinate interpolation error in u, for rational 0 < a < b <= 1."""
        a, b = Q(a), Q(b)
        if not 0 < a < b <= 1:
            raise ValueError("expected an increasing source interval in (0,1]")
        angular = sum(abs(n*(n-3))*c*a**(-n) for n, c in self.terms)
        curvature = 2 + (a*self.phase_slope(a))**2 + angular
        return (b-a)**2 * curvature / 8


@dataclass(frozen=True)
class Certificate:
    status: str
    boundary_margin: Q = Q(0)
    coordinate_error: Q = Q(0)

    @property
    def certified(self):
        return self.status == "certified"


def _endpoint(law, u):
    sine, cosine, error = _trig_rational_enclosure(law.phase(u))
    actual_approximation = (u*u*cosine, u*u*sine)
    center = tuple(_round_fraction(x) for x in actual_approximation)
    rounding_error = max(abs(v-c) for v, c in zip(actual_approximation, center))
    return center, u*u*error + rounding_error


def certify(law, orientation, endpoints, side_floor=Q(1, 4), radius_floor=Q(1, 5)):
    """Replay four rational intervals encoded as u numerators over GRID."""
    side_floor, radius_floor = Q(side_floor), Q(radius_floor)
    if (orientation not in (-1, 1) or len(endpoints) != 8
            or any(type(x) is not int for x in endpoints)
            or side_floor <= 0 or radius_floor < 0):
        raise ValueError("expected an orientation, eight integer endpoints and valid floors")
    intervals = tuple((Q(a, GRID), Q(b, GRID))
                      for a, b in zip(endpoints[::2], endpoints[1::2]))
    ordered = sorted(intervals)
    if (any(not 0 < a < b <= 1 for a, b in ordered)
            or any(b >= c for (_, b), (c, _) in zip(ordered, ordered[1:]))):
        raise ValueError("source intervals must increase, be disjoint and lie in (0,1]")

    starts, ends = [], []
    error = Q(0)
    for a, b in intervals:
        start, start_error = _endpoint(law, a)
        end, end_error = _endpoint(law, b)
        starts.append(start)
        ends.append(end)
        error = max(error, max(start_error, end_error) + law.chord_error(a, b))

    matrix, rhs = _R_matrix(starts, ends, orientation)
    try:
        root, inverse_norm = _inverse_and_root(matrix, rhs)
    except ValueError:
        return Certificate("singular-chords")
    margin = min(min(v, 1-v) for v in root)
    if margin <= 0:
        return Certificate("affine-root-outside")
    boundary = margin / inverse_norm
    if 4*error >= boundary:
        return Certificate("insufficient-degree-margin", boundary, error)
    if not _side_floor_holds(starts, ends, root, error, side_floor):
        return Certificate("insufficient-side-bound", boundary, error)
    if any(a*a <= radius_floor for a, _ in intervals):
        return Certificate("insufficient-radius-bound", boundary, error)
    return Certificate("certified", boundary, error)


def _sample(law, inner_radius, step, max_segments):
    inner = sqrt(inner_radius)
    start, end = law.phase(1.), law.phase(inner)
    requested = (end-start)/step
    if requested > max_segments or max_segments < 2:
        return ()
    count = max(2, ceil(requested))
    samples = [Sample(1., *law.point(1.))]
    for i in range(1, count):
        target = start + (end-start)*i/count
        lo, hi = inner, 1.
        for _ in range(42):
            middle = (lo+hi)/2
            if law.phase(middle) > target:
                lo = middle
            else:
                hi = middle
        u = (lo+hi)/2
        samples.append(Sample(u, *law.point(u)))
    samples.append(Sample(inner, *law.point(inner)))
    return tuple(samples)


def discover(law, radii=(.55, .45, .36, .30, .22), step=.22,
             max_segments=400, max_pairs=12_000, max_checks=16,
             side_floor=Q(1, 4), radius_floor=Q(1, 5)):
    """Bounded numerical search. Only the independent replay can certify."""
    radii = tuple(radii)
    if (not isfinite(step) or step <= 0 or min(max_segments, max_pairs, max_checks) < 1
            or any(not 0 < r < 1 for r in radii)):
        raise ValueError("expected positive search budgets and radii in (0,1)")
    ledger = []
    for radius in radii:
        samples = _sample(law, radius, step, max_segments)
        if not samples:
            ledger.append({"radius": radius, "status": "sample-budget"})
            continue
        scan = scan_sampled_candidates(samples, law, radius, step, max_pairs=max_pairs)
        failures = Counter()
        row = {"radius": radius, "pairs": scan.pairs_examined,
               "pair_budget_exhausted": scan.pair_budget_exhausted,
               "checks": 0, "rejections": failures}
        ledger.append(row)
        for candidate in scan.distinct_candidates:
            if row["checks"] >= max_checks:
                break
            try:
                root, _, _ = refine_square(candidate, law.point, law.tangent)
            except ValueError:
                failures["newton-failure"] += 1
                continue
            for half_angle in (.015, .0075, .003):
                if row["checks"] >= max_checks:
                    break
                row["checks"] += 1
                endpoints = []
                for u in root:
                    width = half_angle / abs(law.phase_slope(u))
                    endpoints.extend((round((u-width)*GRID), round((u+width)*GRID)))
                orientation = candidate["sigma"]
                try:
                    proof = certify(law, orientation, endpoints, side_floor, radius_floor)
                except ValueError:
                    failures["invalid-source-box"] += 1
                    continue
                if proof.certified:
                    return {"status": "certified", "orientation": orientation,
                            "endpoints": tuple(endpoints), "ledger": ledger}
                failures[proof.status] += 1
    return {"status": "inconclusive", "ledger": ledger}


LAWS = {
    "quadratic-plus-inverse": PhaseLaw(((4, Q(1)), (2, Q(1)))),
    "five-halves-plus-square": PhaseLaw(((5, Q(1)), (4, Q(1, 2)))),
    "cubic-plus-three-halves": PhaseLaw(((6, Q(1)), (3, Q(3, 2)))),
}

# Boxes found by discover() on each law. No earlier power-spiral endpoints
# or numerical root are used during replay.
WITNESSES = {
    "quadratic-plus-inverse": (1, (771308297, 771921610, 537868722, 537981623,
                                  660074117, 660371755, 823475927, 824303326)),
    "five-halves-plus-square": (1, (964463260, 965712473, 928184417, 929186262,
                                   653424741, 653556394, 743638559, 743917212)),
    "cubic-plus-three-halves": (-1, (1048970763, 1050317323, 813306561, 813596486,
                                    720802410, 720936907, 1010775551, 1011860630)),
}

FLOORS = {
    "quadratic-plus-inverse": (Q(1, 2), Q(1, 4)),
    "five-halves-plus-square": (Q(1, 4), Q(9, 25)),
    "cubic-plus-three-halves": (Q(2, 5), Q(9, 20)),
}


def report(rediscover=False):
    rows = {}
    for name, law in LAWS.items():
        side, radius = FLOORS[name]
        proof = certify(law, *WITNESSES[name], side, radius)
        row = {"status": proof.status, "side_strict_floor": str(side),
               "radius_strict_floor": str(radius),
               "boundary_margin_approx": float(proof.boundary_margin),
               "residual_error_bound_approx": float(4*proof.coordinate_error)}
        if rediscover:
            row["discovery"] = discover(law, side_floor=side, radius_floor=radius)
        rows[name] = row
    return rows


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rediscover", action="store_true", help="rerun bounded numerical search")
    args = parser.parse_args()
    result = report(args.rediscover)
    print(json.dumps(result, indent=2))
    if any(row["status"] != "certified" or
           (args.rediscover and row["discovery"]["status"] != "certified")
           for row in result.values()):
        raise SystemExit("At least one positive control was inconclusive")
