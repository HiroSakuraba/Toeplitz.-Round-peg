"""Find tip chords certifying the v302 continuum for an unmarked rational Jordan polygon."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from polygon_orbit_census_v223 import cross, dot, orientation_area, pt, sub, validate_jordan_polygon
from rough_crossing_v302 import DISP, ERROR, height_band

RECORD = Path(__file__).resolve().parents[1] / "results/v303/polygon-directions.json"


def arc_reserve(arc, upper):
    """Certify each entire segment by splitting the tent at its sole breakpoint."""
    lo, hi = height_band(upper)
    slope = Q(3, 2) if upper else Q(-3, 4)
    margin = ERROR

    def height_margin(x, y):
        delta = y-slope*min(x, 1-x)
        return min(delta-lo, hi-delta, DISP-abs(delta))

    for a, b in zip(arc, arc[1:]):
        if a.x >= b.x:
            return None
        margin = min(margin, height_margin(a.x, a.y), height_margin(b.x, b.y))
        if a.x < Q(1, 2) < b.x:
            y = a.y + (b.y-a.y)*(Q(1, 2)-a.x)/(b.x-a.x)
            margin = min(margin, height_margin(Q(1, 2), y))
        if margin < 0:
            return None
    return margin


def direction_certificates(polygon):
    polygon = validate_jordan_polygon(polygon)
    if orientation_area(polygon) < 0:
        polygon = polygon[::-1]
    found = []
    for i, start in enumerate(polygon):
        loop = polygon[i:] + polygon[:i] + (start,)
        for j in range(1, len(polygon)):
            axis = sub(loop[j], start)
            length_sq = dot(axis, axis)
            coords = [sub(p, start) for p in loop]
            frame = [pt(dot(p, axis)/length_sq, cross(axis, p)/length_sq) for p in coords]
            lower = arc_reserve(frame[:j+1], False)
            if lower is None:
                continue
            upper = arc_reserve(frame[j:][::-1], True)
            if upper is not None:
                found.append({"start": start, "end": loop[j], "length_squared": length_sq,
                              "normalized_tube": min(lower, upper)})
    return found


def example():
    # Unequal upper and lower deviations on independently chosen intervals.
    return tuple(pt(x, y) for x, y in (
        (0, 0), (Q(1, 4), Q(-19, 100)), (Q(1, 2), Q(-3, 8)),
        (1, 0), (Q(3, 4), Q(39, 100)), (Q(1, 2), Q(77, 100)),
        (Q(1, 4), Q(3, 8))))


def moved(p):
    return pt(7+3*p.x-4*p.y, -2+4*p.x+3*p.y)


def serial(rows):
    return [{"start": [str(r["start"].x), str(r["start"].y)],
             "end": [str(r["end"].x), str(r["end"].y)],
             "length_squared": str(r["length_squared"]), "normalized_tube": str(r["normalized_tube"])}
            for r in rows]


def example_report():
    polygon = example()
    square = (pt(0, 0), pt(1, 0), pt(1, 1), pt(0, 1))
    return {"independent_arcs": serial(direction_certificates(polygon)),
            "rotated_scaled_translated": serial(direction_certificates(tuple(map(moved, polygon)))),
            "square_boundary": serial(direction_certificates(square))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--polygon", type=Path, help="JSON array of cyclic [x,y] pairs")
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.polygon:
        points = json.loads(args.polygon.read_text(encoding="utf-8"), parse_float=Q)
        rows = direction_certificates(tuple(pt(x, y) for x, y in points))
        print(json.dumps(serial(rows), indent=2))
    else:
        data = json.dumps(example_report(), indent=2, sort_keys=True) + "\n"
        if args.write:
            RECORD.write_text(data, encoding="utf-8")
        elif RECORD.read_text(encoding="utf-8") != data:
            raise SystemExit("v303 record differs from exact replay")
        print("Verified the polygon direction search examples.")
