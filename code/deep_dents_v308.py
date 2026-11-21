"""Exact strip-cover and large-square parity certificate for rational Jordan polygons."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from polygon_orbit_census_v223 import pt, validate_jordan_polygon
from convex_shell_v307 import inside, pairs


def segment_interval(a, b, box):
    """Closed parameter interval cut out by a closed axis-aligned box."""
    lo, hi = Q(0), Q(1)
    for start, end, lower, upper in ((a.x,b.x,box[0],box[1]), (a.y,b.y,box[2],box[3])):
        if end == start:
            if not lower <= start <= upper:
                return None
            continue
        u, v = (lower-start)/(end-start), (upper-start)/(end-start)
        lo, hi = max(lo,min(u,v)), min(hi,max(u,v))
        if lo > hi:
            return None
    return lo, hi


def covers_segment(a, b, boxes):
    intervals = [segment_interval(a,b,box) for box in boxes]
    end = Q(0)
    for lo, hi in sorted(i for i in intervals if i is not None):
        if lo > end:
            return False
        end = max(end,hi)
    return end == 1


def certify(polygon, lanes, ymin, ymax, cap):
    """Outer endpoints of the first and last lanes are the left and right walls."""
    polygon = validate_jordan_polygon(polygon)
    lanes = [tuple(map(Q,lane)) for lane in lanes]
    ymin, ymax, cap = map(Q,(ymin,ymax,cap))
    if len(lanes)<2 or any(len(lane)!=2 for lane in lanes):
        raise ValueError("need two or more ordered vertical strips")
    gaps = [(b,c) for (_,b),(c,_) in zip(lanes,lanes[1:])]
    if any(a>=b for a,b in lanes+gaps) or not 0<2*cap<ymax-ymin:
        raise ValueError("invalid strip arrangement")
    width = max(b-a for a,b in lanes)
    cutoff = min(ymax-ymin-2*cap, min(b-a for a,b in gaps))/2
    if cutoff<=width or cutoff**2<=2*cap**2:
        raise ValueError("no open interval of excluded sizes")
    left, right = lanes[0][0], lanes[-1][1]
    boxes = [(a,b,ymin,ymax) for a,b in lanes]
    boxes += [(left,right,ymin,ymin+cap), (left,right,ymax-cap,ymax)]
    if not all(covers_segment(a,b,boxes) for a,b in pairs(polygon)):
        raise ValueError("uncovered polygon edge")
    occupied = [inside(pt((a+b)/2,(ymin+ymax)/2),polygon) for a,b in gaps]
    runs = sum(flag and (i==0 or not occupied[i-1]) for i,flag in enumerate(occupied))
    return {"occupied_cells": list(map(int,occupied)), "runs":runs, "parity":runs%2,
            "forbidden_half_diagonal_squared_lower_open":str(max(width**2,2*cap**2)),
            "forbidden_half_diagonal_upper_open":str(cutoff),
            "forced_side_squared_at_least":str(2*cutoff**2) if runs%2 else None}


def example(kind):
    cap = Q(1,32)
    left, right = (-10,10) if kind=="three_runs" else (-6,6)
    mids = (-6,-2,2,6) if kind=="three_runs" else (-2,2)
    lanes = [(left,left+cap), *((x-Q(1,8),x+Q(1,8)) for x in mids), (right-cap,right)]
    top, bottom = 2-cap/2, -2+cap/2
    if kind=="opposite_dents":
        a,b = lanes[1]
        c,d = lanes[2]
        vertices = [(left,-2),(c,-2),(c,top),(d,top),(d,-2),(right,-2),
                    (right,2),(b,2),(b,bottom),(a,bottom),(a,2),(left,2)]
    elif kind in ("two_runs","three_runs"):
        vertices = [(left,-2)]
        for a,b in zip(mids[::2],mids[1::2]):
            vertices.extend([(a,-2),(a,top),(b,top),(b,-2)])
        vertices.extend([(right,-2),(right,2),(left,2)])
    else:
        raise ValueError("unknown example")
    return tuple(pt(*v) for v in vertices), lanes, -2, 2, cap


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="JSON keys polygon, lanes, ymin, ymax, cap")
    args = parser.parse_args()
    if args.input:
        data = json.loads(args.input.read_text(),parse_float=Q)
        result = certify(tuple(pt(*v) for v in data["polygon"]),data["lanes"],
                         data["ymin"],data["ymax"],data["cap"])
    else:
        result = {name:certify(*example(name)) for name in ("opposite_dents","two_runs","three_runs")}
    print(json.dumps(result,indent=2))
