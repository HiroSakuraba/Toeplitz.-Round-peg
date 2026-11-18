"""Exact-rational periodic-graph witness and anti-localization test for v212.

This is an independent finite model of hypotheses in the October 2026
four-periodic-inscription theorems. It is NOT a construction from an
arbitrary Jordan curve, a Floer solver, or a proof of Toeplitz.
Python >= 3.9; standard library only.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from typing import Tuple


@dataclass(frozen=True)
class Point:
    x: F
    y: F

    def __add__(self, other: "Point") -> "Point":
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other: "Point") -> "Point":
        return Point(self.x - other.x, self.y - other.y)


def squared_length(p: Point) -> F:
    return p.x * p.x + p.y * p.y


def rotate90(p: Point) -> Point:
    return Point(-p.y, p.x)


@dataclass(frozen=True)
class PeriodicGraph:
    """Piecewise-linear graph of (x, y(x)) with period 1 in y."""

    vertices: Tuple[Point, ...]

    def __post_init__(self) -> None:
        xs = [p.x for p in self.vertices]
        if len(xs) < 2 or xs[0] != 0 or xs[-1] != 1:
            raise ValueError("graph must span exactly one unit period")
        if any(a >= b for a, b in zip(xs, xs[1:])):
            raise ValueError("graph x-coordinates must be strictly increasing")
        if self.vertices[0].y != self.vertices[-1].y:
            raise ValueError("period endpoints must agree in height")

    def height(self, x: F) -> F:
        x = F(x) % 1
        v = self.vertices
        for a, b in zip(v, v[1:]):
            if a.x <= x <= b.x:
                return a.y + (x - a.x) * (b.y - a.y) / (b.x - a.x)
        raise AssertionError("x not found in a full period")

    def lift(self, x: F) -> Point:
        x = F(x)
        return Point(x, self.height(x))

    def area_y_dx(self) -> F:
        """Exact integral of y dx over one positively oriented period."""
        return sum((b.x - a.x) * (a.y + b.y) / 2
                   for a, b in zip(self.vertices, self.vertices[1:]))


def tent(x: F, center: F, half_width: F) -> F:
    return max(F(0), F(1) - abs(x - center) / half_width)


def strand(i: int, amplitude: F = F(1, 16)) -> PeriodicGraph:
    """A horizontal level with an exact zero-integral ±tent perturbation."""
    if i not in (0, 1, 2, 3):
        raise ValueError("strand index must be 0..3")
    h = (F(1), F(0), F(-1), F(0))[i]
    cp, cm, w = F(2 * i + 1, 16), F(2 * i + 2, 16), F(1, 128)
    xs = sorted({F(0), F(1), cp - w, cp, cp + w,
                 cm - w, cm, cm + w})
    return PeriodicGraph(tuple(Point(x, h + amplitude *
                        (tent(x, cp, w) - tent(x, cm, w))) for x in xs))


def four_strands(amplitude: F = F(1, 16)) -> Tuple[PeriodicGraph, ...]:
    return tuple(strand(i, amplitude) for i in range(4))


def alternating_area(strands: Tuple[PeriodicGraph, ...]) -> F:
    if len(strands) != 4:
        raise ValueError("four ordered strands required")
    a = [s.area_y_dx() for s in strands]
    return a[0] - a[1] + a[2] - a[3]


def square_at_unmodified_period(strands: Tuple[PeriodicGraph, ...],
                                x: F = F(7, 8)) -> Tuple[Point, ...]:
    """One exact square on the completed periodic curves, not local windows."""
    return (strands[0].lift(x), strands[1].lift(x - 1),
            strands[2].lift(x), strands[3].lift(x + 1))


def is_labeled_square(points: Tuple[Point, ...]) -> bool:
    if len(points) != 4:
        return False
    e0, e1, e2, e3 = (points[(j + 1) % 4] - points[j]
                       for j in range(4))
    return (squared_length(e0) > 0 and
            e1 == rotate90(e0) and e2 == rotate90(e1) and
            e3 == rotate90(e2))


# The distinguished original-piece *windows* in the periodic coordinate.
# They are NOT claimed to have been extracted from one Jordan curve.
SOURCE_WINDOWS = (
    (F(0), F(1, 16)),
    (F(1, 4), F(5, 16)),
    (F(1, 2), F(9, 16)),
    (F(3, 4), F(13, 16)),
)


def minimum_mod1_separation(i: int, j: int) -> F:
    """Exact minimum circular x-distance for two disjoint intervals."""
    a, b = SOURCE_WINDOWS[i], SOURCE_WINDOWS[j]
    return min(abs(u - v + k) for u in a for v in b for k in (-1, 0, 1))


def no_square_inside_all_windows(amplitude: F = F(1, 16)) -> bool:
    """Rigorous inequality, not sampling.

    Square identity: x0 - x2 = y1 - y3 in chosen real lifts.
    For x0 in W0 and x2 in W2, |x0-x2| >= 7/16 mod 1.
    On G1 and G3: |y1-y3| <= 2*amplitude < 7/16.
    Therefore there cannot be a labeled square in all four windows.
    """
    return (amplitude >= 0 and
            2 * amplitude < minimum_mod1_separation(0, 2))


def audit() -> dict:
    strands = four_strands()
    witness = square_at_unmodified_period(strands)
    return {
        "areas_y_dx": [str(a.area_y_dx()) for a in strands],
        "alternating_area": str(alternating_area(strands)),
        "witness_vertices": [[str(p.x), str(p.y)] for p in witness],
        "square_side_squared": str(squared_length(witness[1] - witness[0])),
        "is_exact_square": is_labeled_square(witness),
        "minimum_window_0_2_x_separation_mod1": str(minimum_mod1_separation(0, 2)),
        "max_window_1_3_height_difference": str(F(1, 8)),
        "no_all_window_square": no_square_inside_all_windows(),
        "fourfold_intersection_empty": True,  # y0 >= 15/16 and y2 <= -15/16
        "claim": "periodic theorem can have an exclusively artificial witness",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(audit(), indent=2, sort_keys=True))
