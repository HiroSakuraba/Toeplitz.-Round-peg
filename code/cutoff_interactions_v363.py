"""Finite differences of filtered square parity and an exact four-edit example."""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import comb

from polygon_orbit_census_v223 import _solve, pt, validate_jordan_polygon


def reconstruction_masks(count, order=4):
    """Masks with nonzero coefficient in the full-mask parity formula over F2."""
    if count < 0 or order < 0:
        raise ValueError("negative count or order")
    for size in range(min(count, order)+1):
        coefficient = sum(comb(count-size, j) for j in range(min(count, order)-size+1))
        if coefficient % 2:
            for subset in combinations(range(count), size):
                yield sum(1 << i for i in subset)


# The triangle is (0,0), (4,0), (1,3). Contact lines in order around the square:
# y=0, y=0, x+y=4, 3x-y=0.
SIDE = Q(12,7)
CENTERS = ((Q(4,7),0),(Q(16,7),0),(Q(16,7),SIDE),(Q(4,7),SIDE))
WIDTH = Q(1,64)
EPS = Q(1,2**18)
CUTOFF = (12+Q(19,2)*EPS)/7
NORMALS = ((0,1),(1,1),(3,-1))
LEVELS = (0,4,0)
# Each corner's x and y rows in variables (origin x,y; first side u,v).
COORDS = (((1,0,0,0),(0,1,0,0)), ((1,0,1,0),(0,1,0,1)),
          ((1,0,1,-1),(0,1,1,1)), ((1,0,0,-1),(0,1,1,0)))


def _size_interval(constraints):
    # At least two distinct assigned corners bound the side vector in both
    # coordinates. Test vertices and perpendicular feet for the norm extrema.
    def feasible(u,v):
        return all(a*u+b*v <= c for a,b,c in constraints)

    sizes = []
    for (a,b,c),(d,e,f) in combinations(constraints,2):
        det = a*e-b*d
        if det:
            u,v = (c*e-b*f)/det, (a*f-c*d)/det
            if feasible(u,v):
                sizes.append(u*u+v*v)
    if not sizes:
        return None
    largest = max(sizes)
    if feasible(0,0):
        return Q(0),largest
    for a,b,c in constraints:
        norm = a*a+b*b
        if norm:
            u,v = a*c/norm,b*c/norm
            if feasible(u,v):
                sizes.append(u*u+v*v)
    return min(sizes),largest


def box_size_intervals(boxes):
    """Exact squared sides of squares hitting each box at a distinct corner.

    Each closed box is (xmin, xmax, ymin, ymax). Zero in the result permits
    the collapsed configuration in overlapping boxes.
    """
    boxes = [tuple(map(Q,box)) for box in boxes]
    if not 2 <= len(boxes) <= 4 or any(
            len(b)!=4 or b[0]>b[1] or b[2]>b[3] for b in boxes):
        raise ValueError("expected two to four ordered coordinate boxes")
    intervals = []
    for rest in permutations(range(1,4),len(boxes)-1):
        corners = (0,)+rest
        bounds = {}
        for axis in range(2):
            for i,ci in enumerate(corners):
                for j,cj in enumerate(corners):
                    x,y = COORDS[ci][axis][2:]
                    xx,yy = COORDS[cj][axis][2:]
                    normal = (xx-x,yy-y)
                    cap = boxes[j][2*axis+1]-boxes[i][2*axis]
                    bounds[normal] = min(bounds.get(normal,cap),cap)
        constraints = [(a,b,c) for (a,b),c in bounds.items()]
        interval = _size_interval(constraints)
        if interval is not None:
            intervals.append(interval)
    merged = []
    for low,high in sorted(intervals):
        if merged and low <= merged[-1][1]:
            merged[-1] = (merged[-1][0],max(high,merged[-1][1]))
        else:
            merged.append((low,high))
    return merged


def clear_squared_cutoff(band, patch_groups):
    """Choose a squared cutoff in the widest open gap, or None if covered."""
    low,high = map(Q,band)
    if not 0 < low < high:
        raise ValueError("expected an increasing positive squared-size band")
    forbidden = sorted(interval for group in patch_groups for interval in box_size_intervals(group))
    gaps = []
    cursor = low
    for a,b in forbidden:
        if cursor < min(a,high):
            gaps.append((cursor,min(a,high)))
        cursor = max(cursor,b)
        if cursor >= high:
            break
    if cursor < high:
        gaps.append((cursor,high))
    if gaps:
        a,b = max(gaps,key=lambda g:g[1]-g[0])
        return (a+b)/2
    return None


def triangle(mask):
    """Four disjoint trapezoidal graph edits. Bit order follows the square."""
    if not 0 <= mask < 16:
        raise ValueError("invalid edit mask")
    vertices = [pt(0,0)]
    for i, tangent, shift, endpoint in (
            (0,(1,0),(0,-1),None), (1,(1,0),(0,-1),(4,0)),
            (2,(-1,1),(0,1),(1,3)), (3,(-1,-3),(-1,0),None)):
        if mask & (1 << i):
            cx,cy = CENTERS[i]
            for k, active in ((-2,0),(-1,1),(1,1),(2,0)):
                vertices.append(pt(cx+k*WIDTH*tangent[0]+active*EPS*shift[0],
                                   cy+k*WIDTH*tangent[1]+active*EPS*shift[1]))
        if endpoint:
            vertices.append(pt(*endpoint))
    return validate_jordan_polygon(vertices)


def moving_square(mask, time=Q(1)):
    """The four central parallel pieces determine the square exactly."""
    t = Q(time)
    if not 0 <= t <= 1 or not 0 <= mask < 16:
        raise ValueError("invalid motion")
    a,b,c,d = [int(bool(mask & (1 << i))) for i in range(4)]
    u = (12+t*EPS*(a+3*b+3*c+3*d))/7
    v = t*EPS*(a-b)
    y = -a*t*EPS
    x = v+(y+u)/3-d*t*EPS
    return ((x,y),(x+u,y+v),(x+u-v,y+u+v),(x-v,y+u)), u*u+v*v


def cutoff_chart_cover():
    """Enclose every cutoff square on the swept triangle, including its connectors."""
    survivors = []
    for lines in product(range(3),repeat=4):
        if len(set(lines)) == 1:
            # A strip |normal.p-level| <= 3 EPS has projection width 6 EPS.
            # A square projects to at least its side on every unit direction.
            if 6*EPS >= CUTOFF:
                raise ArithmeticError("single strip could contain a cutoff square")
            continue
        system = []
        for i,(rx,ry) in zip(lines,COORDS):
            a,b = NORMALS[i]
            system.append([a*x+b*y for x,y in zip(rx,ry)])
        base = _solve(system,[LEVELS[i] for i in lines])
        columns = [_solve(system,[int(i==j) for i in range(4)]) for j in range(4)]
        if base is None:
            raise ArithmeticError("unexpected singular mixed chart")
        errors = [3*EPS*sum(abs(x) for x in row) for row in zip(*columns)]
        lower = sum(max(0,abs(base[j])-errors[j])**2 for j in (2,3))
        upper = sum((abs(base[j])+errors[j])**2 for j in (2,3))
        if not lower <= CUTOFF**2 <= upper:
            continue
        # Even allowing independent errors, all corners stay on the central
        # portions: the connectors have at least WIDTH-EPS x displacement.
        for rx,ry in COORDS:
            center = tuple(sum(q*t for q,t in zip(row,base)) for row in (rx,ry))
            if center not in CENTERS:
                raise ArithmeticError("unexpected cutoff location")
            if sum(abs(q)*e for q,e in zip(rx,errors)) >= WIDTH-EPS:
                raise ArithmeticError("cutoff might touch a connector")
        survivors.append(lines)
    return tuple(survivors)


def distant_bays():
    """Exact fixed-cutoff accounting on four unequal bays moved to new sites."""
    from half_turn_degree_v313 import ear_census, profile
    from polygon_orbit_census_v223 import add, mul

    height = Q(2)
    cutoff = height/5
    sites = (-9,-4,1,6)
    scales = (Q(19,20),Q(24,25),Q(49,50),Q(1))
    roots = ear_census()
    events = sorted((cutoff**2/(size*scale**2), i) for size,_,_ in roots
                    for i,scale in enumerate(scales) if size*scale**2 > cutoff**2)
    if any(size*scale**2 == cutoff**2 for size,_,_ in roots for scale in scales):
        raise ArithmeticError("endpoint lies on the cutoff")
    shift = 2
    shape = [pt(12,1)]
    for x,scale in reversed(list(zip(sites,scales))):
        shape.extend(add(pt(x+shift,1),mul(p,scale)) for p in reversed(profile()[1:-1]))
    shape.extend([pt(-12,1),pt(-12,-1),pt(-1,-1),pt(0,Q(-9,8)),
                  pt(Q(-1,4),Q(-19,16)),pt(1,-1),pt(12,-1)])
    validate_jordan_polygon(shape)
    # The whole square at an old site lies within this horizontal aperture.
    aperture = Q(7,16)+Q(17,32)
    clear = min(abs(a-(b+shift)) for a in sites for b in sites)-aperture-Q(7,16)
    if clear < cutoff or min(min(x+12,12-x-shift) for x in sites) <= aperture+cutoff:
        raise ArithmeticError("old event may meet the new bays or the end walls")
    return events, min(clear,cutoff)
