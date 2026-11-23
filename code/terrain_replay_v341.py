"""Reproduce v341 component witnesses for the v336 single and glued samples."""
import csv
from fractions import Fraction as Q
from pathlib import Path
import random
import sys

from polygon_component_v341 import PolygonComponents, UnsupportedComponent, degree, segment_floor
from polygon_orbit_census_v223 import _side_sq
from robust_squares_v336 import terrain_polygon
from terrain_gluing_v334 import column_maps, interior_square, polygon, terrain_i
from terrain_ledge_v341 import ledge_squares

PREVIOUS = Path(__file__).resolve().parents[1] / "results/v340/end-wall.csv"


def fixtures():
    for kind, count in (("terrain",300), ("union",200)):
        for seed in range(count):
            rng = random.Random(seed)
            upper = [rng.randint(1,8) for _ in range(rng.randint(2,8 if kind == "terrain" else 7))]
            lower, start = [], 0
            if kind == "terrain":
                i, boundary = terrain_i(upper), terrain_polygon(upper)
            else:
                lower = [rng.randint(1,8) for _ in range(rng.randint(1,6))]
                start = rng.randint(-len(lower)+1,len(upper)-1)
                columns = column_maps(upper, lower, 0, start)
                i, _ = interior_square(*columns)
                boundary = polygon(*columns)
            if i >= 2:
                yield kind, seed, upper, lower, start, i, boundary


def witness(boundary, cutoff_squared):
    engine = PolygonComponents(boundary)
    for seed in engine.lattice_seeds():
        if _side_sq(seed) < cutoff_squared:
            break
        try:
            cells = engine.trace(seed)
            floor = min(segment_floor(c.squares) for c in cells.values())
        except UnsupportedComponent:
            continue
        if floor >= cutoff_squared:
            k, hits = degree(cells)
            index = sum(sign for _,sign in hits)
            if index:
                return seed, cells, floor, k, hits, index
    return None


def write_report(stream):
    with PREVIOUS.open() as f:
        previous = {int(r["seed"]):r["method"] for r in csv.DictReader(f)}
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("kind","seed","upper","lower","lower_start","i","v340","ledge_floor2",
                     "marked_square","cells","floor2","ray_k","signed_charts","degree"))
    for kind, n, upper, lower, start, i, boundary in fixtures():
        old = previous.get(n,"interior") if kind == "terrain" else ""
        threshold = Q(i*i,2)
        ledge = ""
        if old == "unresolved":
            sides = [(s[1][0]-s[0][0])**2+(s[1][1]-s[0][1])**2
                     for s,*_ in ledge_squares(boundary)]
            ledge = max((x for x in sides if x >= threshold), default="")
        prefix = (kind,n," ".join(map(str,upper))," ".join(map(str,lower)),
                  start if lower else "",i,old,ledge)
        result = witness(boundary,threshold)
        if result is None:
            writer.writerow((*prefix,"","","","","",""))
            continue
        square,cells,floor,k,hits,index = result
        charts = " ".join(f"{sign:+d}@"+"/".join(map(str,ids)) for ids,sign in hits)
        writer.writerow((*prefix," ".join(f"{p.x}:{p.y}" for p in square),
                         len(cells),floor,k,charts,index))


if __name__ == "__main__":
    write_report(sys.stdout)
