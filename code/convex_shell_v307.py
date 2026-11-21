"""A sufficient rational convex-hull certificate for a macroscopic inscribed square."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from polygon_orbit_census_v223 import cross, dot, pt, sub, validate_jordan_polygon


def convex_hull(points):
    def half(vertices):
        chain = []
        for p in vertices:
            while len(chain) >= 2 and cross(sub(chain[-1], chain[-2]), sub(p, chain[-1])) <= 0:
                chain.pop()
            chain.append(p)
        return chain[:-1]

    ordered = sorted(points)
    return half(ordered) + half(reversed(ordered))


def pairs(loop):
    return zip(loop, loop[1:] + loop[:1])


def inside(p, loop):
    winding = 0
    for a, b in pairs(loop):
        turn = cross(sub(b, a), sub(p, a))
        if turn == 0 and dot(sub(p, a), sub(p, b)) <= 0:
            return False
        if a.y <= p.y < b.y and turn > 0:
            winding += 1
        if b.y <= p.y < a.y and turn < 0:
            winding -= 1
    return winding != 0


def covered(a, b, hull, depth):
    # Cover the WHOLE input segment by a union of support strips.
    # max(|dx|,|dy|) <= Euclidean edge length makes each strip conservative.
    intervals = []
    for p, end in pairs(hull):
        edge = sub(end, p)
        bound = depth * max(abs(edge.x), abs(edge.y))
        start = cross(edge, sub(a, p)) - bound
        finish = cross(edge, sub(b, p)) - bound
        if start <= 0 and finish <= 0:
            return True
        if start <= 0 < finish:
            intervals.append((Q(0), -start / (finish - start)))
        if finish <= 0 < start:
            intervals.append((start / (start - finish), Q(1)))
    reach = Q(0)
    for lo, hi in sorted(intervals):
        if lo > reach:
            return False
        reach = max(reach, hi)
    return reach == 1


def certify(polygon, center, r, R, depth, q):
    polygon = validate_jordan_polygon(polygon)
    hull = convex_hull(polygon)
    gap = q*q - (r-depth)*q + depth*R
    if not (0 < q < r <= R and depth >= 0 and gap < 0):
        raise ValueError("no strict size gap at the requested cutoff")
    if not inside(center, polygon):
        raise ValueError("the loop must enclose the disk center")
    for a, b in pairs(hull):
        edge = sub(b, a)
        slack = cross(edge, sub(center, a))
        if slack <= 0 or slack**2 < r*r*dot(edge, edge):
            raise ValueError("inner disk exceeds hull")
        if dot(sub(a, center), sub(a, center)) > R*R:
            raise ValueError("outer disk misses hull")
    if not all(covered(a, b, hull, depth) for a, b in pairs(polygon)):
        raise ValueError("support strips do not certify the whole boundary")
    return {"cutoff_polynomial": str(gap), "half_diagonal_above": str(q),
            "large_root_linear_term": str(r-depth),
            "large_root_discriminant": str((r-depth)**2 - 4*depth*R)}


def fixture(notch=False):
    d = Q(1,32)
    top = [(1,1),(1,1-d),(1-d,1-d),(1-d,1)] if notch else [(1,1),(Q(1,2),1-d),(-1,1)]
    return tuple(pt(x,y) for x,y in [(-4,-1),(4,-1),(4,1),*top,(-4,1)])


def replay():
    return certify(fixture(),pt(0,0),Q(1),Q(5),Q(1,32),Q(3,4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="JSON keys polygon, center, r, R, depth, q")
    args = parser.parse_args()
    if args.input:
        data = json.loads(args.input.read_text(), parse_float=Q)
        result = certify(tuple(pt(*xy) for xy in data["polygon"]), pt(*data["center"]),
                         *(Q(data[k]) for k in ("r","R","depth","q")))
    else:
        result = replay()
    print(json.dumps(result,indent=2))
