"""Certified marked ledge contacts in the four quarter-turned terrain frames."""
from itertools import product


def rotate(point, turns):
    x, y = point
    for _ in range(turns % 4):
        x, y = -y, x
    return x, y


def wall_status(p, j, y, slope):
    """Does the marked wall contain the point or cross the required line at its end?"""
    start, end = p[j], p[(j+1) % len(p)]
    if min(start[1], end[1]) < y < max(start[1], end[1]):
        return "open"
    if y == start[1]:
        before, after = p[j-1], end
    elif y == end[1]:
        before, after = start, p[(j+2) % len(p)]
    else:
        return None
    level = start[0] + slope*y
    signs = [x+slope*h-level for x,h in (before, after)]
    return "corner" if signs[0]*signs[1] < 0 else None


def ledge_squares(polygon):
    """Yield (labeled square, original base edge, frame, wall-contact type).

    Input is a simplified orthogonal Jordan polygon without collinear subdivisions.
    """
    for turns in range(4):
        p = [rotate(v, turns) for v in polygon]
        edges = list(zip(p, p[1:]+p[:1]))
        walls = [j for j,(a,b) in enumerate(edges) if a[0] == b[0]]
        for j, ((x1,y),(x2,y2)) in enumerate(edges):
            if y != y2:
                continue
            for left_edge, top_edge, right_edge in product(walls, repeat=3):
                a, k, b = (p[edge][0] for edge in (left_edge, top_edge, right_edge))
                base_x = a+b-k
                if a == k == b or not min(x1,x2) < base_x < max(x1,x2):
                    continue
                if wall_status(p,top_edge,y+b-a,1) != "open":
                    continue
                left = wall_status(p,left_edge,y+k-a,1)
                right = wall_status(p,right_edge,y+b-k,-1)
                if not left or not right or "open" not in (left, right):
                    continue
                square = ((base_x,y),(b,y+b-k),(k,y+b-a),(a,y+k-a))
                method = "corner" if "corner" in (left,right) else "strict"
                yield tuple(rotate(v,-turns) for v in square), j, turns, method
