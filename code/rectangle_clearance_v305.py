"""Physical rectangle boxes found from unmarked rational Jordan polygons."""
import argparse
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json
from pathlib import Path

from polygon_orbit_census_v223 import (pt, sub, dot, cross, _affine_system, _solve,
                                       _evaluate, validate_jordan_polygon, orientation_area)

T = Q(99, 100)
BOX = Q(1, 64)
SPAN = Q(1, 6)
HEIGHT = Q(1, 4)
UNIT = (pt(0, 0), pt(1, 0), pt(1, 1), pt(0, 1))
RECORD = Path(__file__).resolve().parents[1] / "results/v304-v306/crossings.json"


def fixture():
    raw = [(-4, '-1/8'), ('-1/4', 0), (0, 0), ('1/4', '1/128'),
           ('3/4', '-1/16'), ('5/4', '1/16'), (5, 0), (5, 1),
           ('5/4', '7/8'), ('3/4', '9/8'), ('1/4', '9/8'),
           ('-1/4', '7/8'), (-4, 1)]
    return tuple(pt(Q(x)-Q(y), Q(x)+Q(y)) for x, y in raw)


def local_graph(edges, p, slope):
    pieces = []
    for a, b in edges:
        if a.x > b.x:
            a, b = b, a
        l, r = max(a.x, p.x-SPAN), min(b.x, p.x+SPAN)
        if l > r:
            continue
        if a.x == b.x:
            if max(a.y, b.y) >= p.y-HEIGHT and min(a.y, b.y) <= p.y+HEIGHT:
                return None
            continue
        m = (b.y-a.y)/(b.x-a.x)
        yl, yr = a.y+m*(l-a.x), a.y+m*(r-a.x)
        if min(yl, yr) > p.y+HEIGHT or max(yl, yr) < p.y-HEIGHT:
            continue
        if max(abs(yl-p.y), abs(yr-p.y)) >= HEIGHT:
            return None
        if l < r:
            pieces.append((l, r, yl, yr, m))
        elif l not in (p.x-SPAN, p.x+SPAN):
            return None
    pieces.sort()
    if not pieces or pieces[0][0] != p.x-SPAN or pieces[-1][1] != p.x+SPAN:
        return None
    if any(a[1] != b[0] or a[3] != b[2] for a, b in zip(pieces, pieces[1:])):
        return None
    defect = max(abs(y-p.y-slope*(x-p.x)) for l,r,yl,yr,m in pieces for x,y in ((l,yl),(r,yr)))
    return defect, max(abs(row[4]) for row in pieces)


def matrix(slopes, k):
    a,b,c,d = slopes
    return ((-a,1,0,0),(-b,1,-b,1),(-c,1,k-c,1+c*k),(-d,1,k,d*k))


def certify(poly, corners, chart):
    axis = sub(corners[1], corners[0])
    size = dot(axis, axis)
    if not size:
        return None
    points = [sub(v, corners[0]) for v in poly]
    points = [pt(dot(p,axis)/size, cross(axis,p)/size) for p in points]
    edges = list(zip(points, points[1:]+points[:1]))
    slopes = []
    for i in chart:
        a,b = edges[i]
        if a.x == b.x:
            return None
        slopes.append((b.y-a.y)/(b.x-a.x))
    A = matrix(slopes, 1)
    cols = [_solve(A, [int(i==j) for i in range(4)]) for j in range(4)]
    if None in cols:
        return None
    motion = BOX*max(sum(abs(z) for z in row) for row in zip(*cols))
    # Worst displacement of a corner from its unit-square location in either coordinate.
    offset = (2+1/T)*motion + 1/T-1
    clearance = min(SPAN,HEIGHT)-offset
    for i,(p,m) in enumerate(zip(UNIT, slopes)):
        graph = local_graph(edges,p,m)
        if graph is None:
            return None
        defect, lip = graph
        ratio_error = (1/T-1)*(1+motion+abs(m)*motion) if i>=2 else 0
        clearance = min(clearance,(BOX-defect-ratio_error)/(1+lip))
    if clearance <= 0 or motion >= 1:
        return None
    return dict(chart=chart, origin=corners[0], axis=axis, slopes=slopes, motion=motion,
                clearance=clearance, side_floor=T*(1-motion), anchor_side_squared=size,
                physical_side_floor_squared=size*(T*(1-motion))**2,
                physical_clearance_squared=size*clearance**2)


def search(polygon):
    p = validate_jordan_polygon(polygon)
    if orientation_area(p)<0:
        p = p[::-1]
    edges = list(zip(p,p[1:]+p[:1]))
    # Gracefully ordered isolated edge-chart squares are candidate anchors.
    # Singular charts can still be covered by a neighboring regular candidate.
    for ids in combinations_with_replacement(range(len(p)),4):
        if len(set(ids))<2 or any(ids.count(i)>2 for i in ids):
            continue
        ed = [edges[i] for i in ids]
        tags = _solve(*_affine_system(ed,1))
        if tags is None or min(tags)<0 or max(tags)>1:
            continue
        q = _evaluate(ed,tags)
        if len(set(q))!=4:
            continue
        for j in range(4):
            cert = certify(p,q[j:]+q[:j],ids[j:]+ids[:j])
            if cert:
                yield cert


def endpoint(k, cert, kink):
    slopes = cert['slopes'].copy()
    slopes[0] = kink
    return _solve(matrix(slopes,k), (0,0,1-k,1-k))


def report():
    cert = next(c for c in search(fixture()) if c['origin']==pt(0,0) and c['axis']==pt(1,1))
    if cert['clearance'] <= Q(1,400) or cert['side_floor'] <= Q(94,100):
        raise ArithmeticError("insufficient physical clearance")
    ends = []
    for k in (T,1/T):
        solutions = []
        for slope in (Q(0),Q(1,32)):
            x = endpoint(k,cert,slope)
            if (x[0] >= 0) == bool(slope):
                A = matrix(cert['slopes'],1)
                if any(abs(sum(a*b for a,b in zip(row,x))) >= BOX for row in A):
                    raise ArithmeticError("endpoint outside box")
                solutions.append(x)
        if len(solutions)!=1:
            raise ArithmeticError("endpoint multiplicity")
        ends.append(solutions[0])
    if not (ends[0][0]<0<ends[1][0] and ends[1][3]<0<ends[0][3]):
        raise ArithmeticError("missing support or orientation exchange")
    return dict(chart=cert['chart'], normalized_clearance=str(cert['clearance']),
                side_floor=str(cert['side_floor']), endpoints=[list(map(str,x)) for x in ends],
                tail_budget="1/1024", radial_budget="1/100000")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--polygon', type=Path, help='cyclic JSON [x,y] pairs, fractions accepted')
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--write', action='store_true')
    args = parser.parse_args()
    if args.polygon:
        points = json.loads(args.polygon.read_text(),parse_float=Q)
        for certificate in search(tuple(pt(x,y) for x,y in points)):
            for key in ("origin", "axis"):
                point = certificate[key]
                certificate[key] = [str(point.x), str(point.y)]
            print(json.dumps(certificate,default=str))
    else:
        data = json.dumps(report(),indent=2)+'\n'
        if args.write:
            RECORD.write_text(data)
        elif RECORD.read_text()!=data:
            raise SystemExit('certificate differs')
