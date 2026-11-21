"""Independent odd square degrees on the two sides of a bridged, deeply cut polygon."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from deep_dents_v308 import certify, covers_segment, example, pt
from convex_shell_v307 import pairs
from polygon_orbit_census_v223 import _affine_system, _rank, census_polygon, validate_jordan_polygon


def certify_runs(polygon, lanes, low, high, cap, displacement=0, inset=0):
    """Localize one degree per occupied run using a single cap across each empty gap."""
    poly = validate_jordan_polygon(polygon)
    lanes = [tuple(map(Q,lane)) for lane in lanes]
    low, high, cap, eta, inset = map(Q,(low,high,cap,displacement,inset))
    base = certify(poly,lanes,low,high,cap)
    marks = base["occupied_cells"]
    if not marks[0] or not marks[-1] or any(a==b==0 for a,b in zip(marks,marks[1:])):
        raise ValueError("need occupied outer cells and isolated empty cells")
    if min(eta,inset)<0:
        raise ValueError("negative allowance")

    # Each empty gap can use its own bridging cap. Test whole segments in either cover.
    bridges = []
    left_edge, right_edge = lanes[0][0], lanes[-1][1]
    vertical = [(l,r,low,high) for l,r in lanes]
    for i,mark in enumerate(marks):
        if mark:
            continue
        lo,hi = lanes[i][1],lanes[i+1][0]
        for name,near,far in (("top",(high-cap,high),(low,low+cap)),
                              ("bottom",(low,low+cap),(high-cap,high))):
            boxes = vertical + [(left_edge,right_edge,*near), (left_edge,lo,*far),
                                (hi,right_edge,*far)]
            if all(covers_segment(p,q,boxes) for p,q in pairs(poly)):
                bridges.append((i,name))
                break
        else:
            raise ValueError("empty gap uses both caps")

    lanes = [(l-eta,r+eta) for l,r in lanes]
    height, cap = high-low+2*eta, cap+2*eta
    gaps = [r[0]-l[1] for l,r in zip(lanes,lanes[1:])]
    b = min(*gaps,height-2*cap)/2
    width = max(r-l for l,r in lanes)
    if b<=width or b*b<=2*cap*cap:
        raise ValueError("size gap unavailable after displacement")
    mids = [(l+r)/2 for l,r in lanes]
    walls = []
    run_start = 0
    for i,bridge in bridges:
        if mids[i+1]-mids[i] <= height/2:
            raise ValueError("nonadjacent runs could interact across short gaps")
        left, right = mids[i]-inset, mids[i+1]+inset
        if not (left+height/2 < lanes[i+1][0] and right-height/2 > lanes[i][1]):
            raise ValueError("wall can reach the far lane")
        strip = max(lanes[i][1]-left,right-lanes[i+1][0])
        wall_bound = 2*(strip+cap)**2
        if wall_bound >= b*b:
            raise ValueError("one-cap wall clearance exhausted")
        if 8*inset**2 >= min(height,mids[i]-mids[run_start])**2:
            raise ValueError("wall could cut the pure reference square set")
        walls.append(dict(left_wall=str(left),right_wall=str(right),bridge=bridge,
                          wall_radius_squared_upper=str(wall_bound)))
        run_start = i+1
    if bridges and 8*inset**2 >= min(height,mids[-1]-mids[run_start])**2:
        raise ValueError("wall could cut the final pure reference square set")
    return dict(occupied_cells=marks,walls=walls,large_squares_at_least=1+len(walls),
                gap_upper=str(b),side_squared_floor=str(2*b*b),
                registered_budget=str(eta))


def perturbed_fixture(name="two_runs"):
    poly,*_ = example(name)
    return tuple(pt(p.x+Q((i*i+3*i+1)%19-9,32768),
                    p.y+Q((2*i*i+i+4)%23-11,32768)) for i,p in enumerate(poly))


def census_record(name="two_runs"):
    points = perturbed_fixture(name)
    squares = census_polygon(points,cutoff_squared=Q(2)).squares
    walls = certify_runs(*example(name),displacement=Q(1,1024),inset=1)["walls"]
    rows = []
    for square in squares:
        ids = [int(t) for t in square.curve_positions]
        edges = [(points[i],points[(i+1)%len(points)]) for i in ids]
        if square.has_polygon_vertex or _rank(_affine_system(edges,1)[0])!=4:
            raise ArithmeticError("singular census zero")
        cx = sum(p.x for p in square.vertices)/4
        cy = sum(p.y for p in square.vertices)/4
        rows.append(dict(center_x=str(cx),center_y=str(cy),side_squared=str(square.side_squared),
                         region=sum(cx>Q(w[key]) for w in walls for key in ("left_wall","right_wall"))))
    rows.sort(key=lambda s:Q(s["center_x"]))
    for region in range(2*len(walls)+1):
        count = sum(s["region"]==region for s in rows)
        if count%2 != (1-region%2):
            raise ArithmeticError("census disagrees with localized degree",region)
    return rows


def report():
    data = example("two_runs")
    return {"base":certify_runs(*data,inset=1),
            "quarter_unit_tube":certify_runs(*data,displacement=Q(1,4)),
            "v308_perturbation_tube":certify_runs(*data,displacement=Q(1,1024),inset=1),
            "three_runs":certify_runs(*example("three_runs"),inset=1),
            "census_above_unit_half_diagonal":census_record(),
            "three_run_census":census_record("three_runs")}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,help="v308 polygon input; optional displacement and inset")
    args = parser.parse_args()
    if args.input:
        j = json.loads(args.input.read_text(),parse_float=Q)
        result = certify_runs(tuple(pt(*p) for p in j["polygon"]),j["lanes"],j["ymin"],
                             j["ymax"],j["cap"],j.get("displacement",0),j.get("inset",0))
    else:
        result = report()
    print(json.dumps(result,indent=2))
