"""Exact moving charts for the v278 concave quadrilateral."""
import argparse
from copy import copy
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path

from sage.all import AA, QQ, PolynomialRing, matrix, vector

from adversarial_rectangles_v278 import concave_polygon
from polygon_orbit_census_v223 import census_polygon, pt
from rectangle_atlas_v274 import _candidates, _inner_sample

TIME = PolynomialRing(QQ, "t")
T = TIME.gen()
RING = PolynomialRing(TIME, "x")
X = RING.gen()
ANGLE = PolynomialRing(QQ, "x")
END = QQ(1)/100
POINTS = ((0, 0), (4, T), (1, 1), (0, 3))
DIRECTIONS = tuple(tuple(b-a for a, b in zip(POINTS[i], POINTS[(i+1)%4]))
                   for i in range(4))
ARTIFACT = Path(__file__).resolve().parents[1] / "results/v280/moving-concave-atlas.json"
PERSISTENT = (0, 1, 2, 3)
FIRST_LOWER = (0, 0, 2, 3)
FIRST_UPPER = (0, 0, 1, 3)
SECOND_LOWER = (0, 1, 3, 3)
SECOND_UPPER = (0, 2, 3, 3)
HALF = QQ(1)/2


def rotations(values):
    return tuple(values[i:]+values[:i] for i in range(4))


def residual(points):
    a, b, c, d = points
    mid = a+c-b-d
    side = b-a
    turn = vector(RING, [-side[1], side[0]])
    return vector(RING, tuple(mid) + tuple(X*(d-a)-(1-X)*turn))


@dataclass(frozen=True)
class MovingChart:
    ids: tuple
    det: object
    numerators: tuple

    @property
    def walls(self):
        return self.numerators + tuple(self.det-n for n in self.numerators)


def moving_chart(ids):
    points = [vector(RING, POINTS[i]) for i in ids]
    constant = residual(points)
    columns = []
    for j, i in enumerate(ids):
        displaced = points.copy()
        displaced[j] = displaced[j]+vector(RING, DIRECTIONS[i])
        columns.append(residual(displaced)-constant)
    system = matrix(RING, columns).transpose()
    rhs = -constant
    det = system.det()
    if not det:
        raise ArithmeticError(f"entirely singular chart {ids}")
    numerators = []
    for j in range(4):
        changed = copy(system)
        changed.set_column(j, rhs)
        numerators.append(changed.det())
    return MovingChart(ids, det, tuple(numerators))


def representative_charts():
    return tuple(moving_chart(ids) for ids in _candidates(4)
                 if ids == min(rotations(ids)))


def specialized(poly, time):
    return ANGLE([c(time) for c in poly.list()])


def projection_factors(chart):
    """Nonzero projection polynomials for roots, signs, and wall coincidence."""
    field = TIME.fraction_field()
    field_ring = PolynomialRing(field, "x")
    factors = set()
    polys = set()
    for p in (chart.det, *chart.walls, X, 2*X-1, X-1):
        if not p:
            continue
        polys.add(p.leading_coefficient())
        factors.update(f.monic() for f, _ in field_ring(p).factor())
    for f in factors:
        polys.update(fi.denominator() for fi in f.list())
        polys.add(f.discriminant().numerator())
    for f, g in combinations(factors, 2):
        polys.add(f.resultant(g).numerator())
    return {p.monic() for p in polys if p.degree() > 0}


def projection_times(chart):
    return sorted({AA(-END), AA(END)} |
                  {root for p in projection_factors(chart) for root in p.roots(AA, multiplicities=False)
                   if -END < root < END})


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def contact_walls(parameters):
    return tuple((i, int(p == 1)) for i, p in enumerate(parameters) if p in (0, 1))


def reduced_columns(chart):
    # All four charts with a repeated edge have a removable flat factor x.
    factor = RING(1) if chart.ids == PERSISTENT else X
    columns = (chart.det, *chart.numerators)
    reduced = tuple(p.quo_rem(factor) for p in columns)
    require(all(not remainder for _, remainder in reduced),
            f"nonremovable denominator in {chart.ids}")
    return tuple(p for p, _ in reduced)


def slice_chart(chart, time):
    """Isolate every wall and determinant root in this rational-time fiber."""
    det = specialized(chart.det, time)
    ns = tuple(specialized(n, time) for n in chart.numerators)
    require(det, f"entirely singular slice {chart.ids} at {time}")
    polys = (det, *ns, *(det-n for n in ns))
    events = sorted({AA(0), AA(HALF), AA(1)} |
                    {root for p in polys if p.degree() > 0
                     for root in p.roots(AA, multiplicities=False) if 0 < root < 1})
    singular = 0

    def feasible(a):
        determinant = det(a)
        if not determinant:
            return False
        return all(0 <= n(a)/determinant <= 1 for n in ns)

    bands = []
    for a, b in zip(events, events[1:]):
        x = _inner_sample(a, b)
        if feasible(x):
            walls = contact_walls(tuple(n(x)/det(x) for n in ns))
            if bands and bands[-1][1] == a:
                require(feasible(a), "feasible bands meet at an unaudited singularity")
                p, _, old_walls, inside = bands[-1]
                require(walls == old_walls, "persistent contacts changed inside a band")
                at_event = contact_walls(tuple(n(a)/det(a) for n in ns))
                bands[-1] = (p, b, walls, inside | set(at_event))
            else:
                bands.append((a, b, walls, set()))
    for a in events[1:-1]:
        if not det(a):
            singular += 1
            # Cramer's augmented determinant gives an exact inconsistency witness.
            # This is invariant within a projection chamber by simultaneous root order.
            require(any(n(a) for n in ns), f"unresolved singular fiber {chart.ids} at {(time, a)}")
        elif feasible(a):
            require(any(p <= a <= q for p, q, _, _ in bands),
                    f"unresolved isolated fiber {chart.ids} at {(time, a)}")

    if not bands:
        return (), singular
    den, *numerators = (specialized(p, time) for p in reduced_columns(chart))
    require(not any(0 <= x <= 1 for x in den.roots(AA, multiplicities=False)),
            f"reduced denominator vanishes in {chart.ids} at {time}")
    details = []
    for a, b, walls, inside in bands:
        ends = tuple(tuple(n(x)/den(x) for n in numerators) for x in (a, b))
        require(all(0 <= y <= 1 for end in ends for y in end),
                f"limit exits cube in {chart.ids} at {time}")
        require(inside <= set(walls), f"extra interior edge contact in {chart.ids} at {time}")
        details.append((a, b, walls, tuple(map(contact_walls, ends))))
    return tuple(details), singular


def first_seam(t):
    return (t*t+16)/(t*t-4*t+32)


def second_seam(t):
    return (4-t)/(8-t)


def expected_bands(time):
    result = {PERSISTENT: ((0, 1, (), (((0, 1), (1, 0), (2, 1), (3, 0)),
                                      ((0, 0), (1, 1), (2, 0), (3, 1)))),)}
    if time < 0:
        return result
    a, b = first_seam(time), second_seam(time)
    if time:
        result.update({
            FIRST_LOWER: ((0, a, (), (((2, 1), (3, 0)), ((2, 0),))),),
            FIRST_UPPER: ((a, 1, (), (((2, 1),), ((0, 0), (1, 1), (2, 0), (3, 1)))),),
            SECOND_LOWER: ((0, b, (), (((0, 1), (1, 0)), ((1, 1),))),),
            SECOND_UPPER: ((b, 1, (), (((1, 0),), ((0, 0), (1, 1), (2, 0), (3, 1)))),),
        })
    else:
        result.update({
            FIRST_LOWER: ((0, a, ((0, 0),), (((0, 0), (1, 0), (2, 1), (3, 0)),
                                             ((0, 0), (2, 0)))),),
            FIRST_UPPER: ((a, 1, ((0, 0),), (((0, 0), (2, 1)),
                                             ((0, 0), (1, 1), (2, 0), (3, 1)))),),
            SECOND_LOWER: ((0, b, ((3, 1),), (((0, 1), (1, 0), (2, 1), (3, 1)),
                                              ((1, 1), (3, 1)))),),
            SECOND_UPPER: ((b, 1, ((3, 1),), (((1, 0), (3, 1)),
                                              ((0, 0), (1, 1), (2, 0), (3, 1)))),),
        })
    return result


def substitute(poly, t, x):
    return sum(coefficient(t)*x**i for i, coefficient in enumerate(poly.list()))


def sources(chart, time, angle):
    """Lifted source coordinates; 4 and 0 name the same polygon vertex."""
    d, *ns = reduced_columns(chart)
    det = substitute(d, time, angle)
    return tuple(i+substitute(n, time, angle)/det for i, n in zip(chart.ids, ns))


def mod_four(positions):
    return tuple(0 if p == 4 else p for p in positions)


def verify_gluing(charts):
    c = {chart.ids: chart for chart in charts}
    time = TIME.fraction_field()(T)
    for lower, upper, event in ((FIRST_LOWER, FIRST_UPPER, first_seam(time)),
                                (SECOND_LOWER, SECOND_UPPER, second_seam(time))):
        require(sources(c[lower], time, event) == sources(c[upper], time, event),
                f"corner transition fails at {lower} to {upper}")

    angle = ANGLE.fraction_field().gen()
    for first, second in ((FIRST_LOWER, SECOND_UPPER), (FIRST_UPPER, SECOND_LOWER)):
        a = sources(c[first], QQ(0), angle)
        b = sources(c[second], QQ(0), 1-angle)
        require(a[1:]+a[:1] == mod_four(b),
                f"cyclic identification fails at t=0 for {first} and {second}")

    endpoint_sources = {
        "persistent_at_zero": sources(c[PERSISTENT], time, 0),
        "persistent_at_one": sources(c[PERSISTENT], time, 1),
        "first_at_zero": sources(c[FIRST_LOWER], time, 0),
        "first_at_one": sources(c[FIRST_UPPER], time, 1),
        "second_at_zero": sources(c[SECOND_LOWER], time, 0),
        "second_at_one": sources(c[SECOND_UPPER], time, 1),
    }
    u = 3*time/(time*time+16)
    v = 4-time/3
    require({key: mod_four(p) for key, p in endpoint_sources.items()} == {
        "persistent_at_zero": (1, 1, 3, 3), "persistent_at_one": (0, 2, 2, 0),
        "first_at_zero": (u, u, 3, 3), "first_at_one": (0, 1, 1, 0),
        "second_at_zero": (1, 1, v, v), "second_at_one": (0, 3, 3, 0),
    }, "flat contacts changed")
    # u is in (0,1), and v is in (3,4) for positive t in the window.
    return {key: list(map(str, mod_four(p))) for key, p in endpoint_sources.items()}


def points_at(chart, time, angle):
    return tuple(tuple(a(time)+d(time)*(u-i) for a, d in
                       zip(map(TIME, POINTS[i]), map(TIME, DIRECTIONS[i])))
                 for u, i in zip(sources(chart, time, angle), chart.ids))


def verify_input_geometry():
    # Nonadjacent pairs have the coordinate gaps proved in v280's note.
    # Check no adjacent overlap and exclude all parallel edge pairs.
    for u, v in combinations(DIRECTIONS, 2):
        det = TIME(u[0]*v[1]-u[1]*v[0])
        require(det and not any(-END <= r <= END for r in det.roots(AA, multiplicities=False)),
                "edge directions become parallel in the window")


def certify_atlas():
    verify_input_geometry()
    charts = representative_charts()
    require(len(charts) == 42 and all(len(set(rotations(c.ids))) == 4 for c in charts),
            "chart enumeration lost a cyclic orbit")
    events = set()
    for chart in charts:
        events.update(projection_times(chart))
    require(sorted(events) == [AA(-END), AA(0), AA(END)], "new projection event in the time window")
    flat_data = verify_gluing(charts)
    audits = []
    by_id = {chart.ids: chart for chart in charts}
    for t in (-END, -END/2, QQ(0), END/2, END):
        found = {}
        singular = 0
        for chart in charts:
            bands, count = slice_chart(chart, t)
            singular += count
            if bands:
                found[chart.ids] = bands
        require(found == expected_bands(t), f"unexpected feasible chart or contact at {t}")
        active_squares = (PERSISTENT,) + ((FIRST_LOWER,) if t == 0 else
                                          (FIRST_LOWER, SECOND_UPPER) if t > 0 else ())
        squares = {tuple(sorted(points_at(by_id[ids], t, HALF))) for ids in active_squares}
        independent = {tuple(sorted((QQ(p.x), QQ(p.y)) for p in square.vertices))
                       for square in census_polygon(tuple(pt(*p) for p in concave_polygon(str(t)))).squares}
        require(squares == independent and len(squares) == len(active_squares),
                f"independent square orbit census differs at {t}")
        audits.append({"time": str(t), "active_representative_charts": len(found),
                       "inconsistent_singular_representative_fibers": singular,
                       "geometric_square_count": len(squares)})
    return {
        "time_window": [str(-END), str(END)],
        "ordered_contact_charts": 4*len(charts),
        "projection_event_times": list(map(str, sorted(events))),
        "regular_and_event_fiber_audits": audits,
        "circumradius_lower": str((1-END)/2),
        "positive_time_seams": [str(first_seam(TIME.fraction_field()(T))),
                                str(second_seam(TIME.fraction_field()(T)))],
        "lifted_flat_contacts": flat_data,
        "geometric_squares_negative_zero_positive": [1, 2, 3],
        "framed_relative_intervals_negative_zero_positive": [[1], [1, 1], [1, 1, 1]],
        "moving_components": [
            {"time_map": "t = s", "interval_coordinate": "x", "relative_index": 1},
            {"time_map": "t = abs(s)", "interval_coordinate": "x for first strip, 1-x for second",
             "relative_index": 1, "regular_time_total_parity": 0},
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    audit = certify_atlas()
    data = json.dumps(audit, indent=2, sort_keys=True) + "\n"
    if args.write:
        ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
        ARTIFACT.write_bytes(data.encode("utf-8"))
    else:
        require(ARTIFACT.read_text(encoding="utf-8") == data, "moving atlas artifact differs from replay")
    print("Verified 168 ordered charts, the complete moving atlas and its fold.")


if __name__ == "__main__":
    main()
