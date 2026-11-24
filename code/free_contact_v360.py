"""Three traversals and a free fourth contact. All geometric tests are rational."""
from fractions import Fraction as Q

from polygon_orbit_census_v223 import pt, validate_jordan_polygon


def _inside_clear(point, edges, clearance):
    x, y = point
    inside = False
    for a, b in edges:
        dx, dy = b.x-a.x, b.y-a.y
        t = ((x-a.x)*dx+(y-a.y)*dy)/(dx*dx+dy*dy)
        t = max(0, min(1, t))
        if (x-a.x-t*dx)**2+(y-a.y-t*dy)**2 <= clearance**2:
            return None
        if (a.y > y) != (b.y > y) and x < a.x+(y-a.y)*dx/dy:
            inside = not inside
    return inside


def free_contact_certificate(vertices, windows, L, error):
    """Scan one marked frame. Windows are cyclic vertex-index walks (first,last),
    following polygon order, for bottom, right/red, and left/blue. Return opposite
    probes and squared side floor, or None when this sufficient scan misses.
    """
    L, e = Q(L), Q(error)
    poly = validate_jordan_polygon([pt(*v) for v in vertices])
    if not 0 <= 2*e < L or len(windows) != 3:
        raise ValueError("invalid strip width or window count")
    edges = list(zip(poly, poly[1:]+poly[:1]))
    n = len(poly)
    ends = []
    used = set()
    for (i, j), axis, level in zip(windows, (1,0,0), (0,L,0)):
        if not (0 <= i < n and 0 <= j < n and i != j):
            raise ValueError("invalid window indices")
        indices = [(i+k)%n for k in range((j-i)%n+1)]
        if used.intersection(indices):
            raise ValueError("marked windows must be disjoint")
        used.update(indices)
        coords = [(poly[k].x, poly[k].y) for k in indices]
        if any(abs(v[axis]-level)>e for v in coords):
            return None
        ends.append(sorted((coords[0][1-axis], coords[-1][1-axis])))
    bottom, right, blue = ends
    lo = max(blue[0], L-bottom[1]+4*e, L-right[1]+8*e)
    hi = min(blue[1], L-bottom[0]-4*e, L-right[0]-8*e)
    if lo >= hi:
        return None
    breaks = {lo, hi}
    for a, b in edges:
        if a.y == b.y == L:
            breaks.update(x for x in (a.x,b.x) if lo < x < hi)
        elif a.y != b.y and min(a.y,b.y) <= L <= max(a.y,b.y):
            x = a.x+(b.x-a.x)*(L-a.y)/(b.y-a.y)
            if lo < x < hi:
                breaks.add(x)
    ordered = sorted(breaks)
    probes = {}
    for a, b in zip(ordered, ordered[1:]):
        x = (a+b)/2
        side = _inside_clear((x,L), edges, 8*e)
        if side is not None:
            probes[side] = x
        if len(probes) == 2:
            return (*sorted(probes.values()), (L-2*e)**2/2)
    return None
