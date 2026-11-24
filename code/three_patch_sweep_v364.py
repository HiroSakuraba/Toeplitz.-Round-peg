"""Fourth-corner enclosures and rational replays for the three-patch sweep."""
from fractions import Fraction as Q

from cutoff_interactions_v363 import CENTERS, EPS, SIDE, WIDTH
from free_contact_v360 import _inside_clear
from polygon_orbit_census_v223 import pt, validate_jordan_polygon


def completion_boxes(boxes, side):
    """Outer enclosures for the unmarked corner. Every feasible assignment is covered."""
    boxes = [tuple(map(Q, b)) for b in boxes]
    side = Q(side)
    if side <= 0 or len(boxes) != 3 or any(
            len(b) != 4 or b[0]>b[1] or b[2]>b[3] for b in boxes):
        raise ValueError("expected positive side and three closed boxes")

    def possible(i,j,squared):
        low = high = 0
        for a in (0,2):
            left = boxes[i][a]-boxes[j][a+1]
            right = boxes[i][a+1]-boxes[j][a]
            low += max(0,left,-right)**2
            high += max(abs(left),abs(right))**2
        return low <= squared <= high

    enclosures = []
    for pivot in range(3):
        i,j = [k for k in range(3) if k != pivot]
        if not (possible(i,j,2*side**2) and possible(i,pivot,side**2)
                and possible(j,pivot,side**2)):
            continue
        box = []
        for a in (0,2):
            box.extend((boxes[i][a]+boxes[j][a]-boxes[pivot][a+1],
                        boxes[i][a+1]+boxes[j][a+1]-boxes[pivot][a]))
        enclosures.append(tuple(box))
    return tuple(enclosures)


def _edge_hits_box(p,q,box):
    lo,hi = Q(0),Q(1)
    for x,y,a,b in ((p.x,q.x,box[0],box[1]),(p.y,q.y,box[2],box[3])):
        if x == y:
            if not a <= x <= b:
                return False
        else:
            start,end = sorted(((a-x)/(y-x),(b-x)/(y-x)))
            lo,hi = max(lo,start),min(hi,end)
            if lo > hi:
                return False
    return True


def completion_side(vertices, boxes, side):
    """Return inside/outside only if every surviving enclosure has this strict sign."""
    poly = validate_jordan_polygon([pt(*p) for p in vertices])
    edges = tuple(zip(poly,poly[1:]+poly[:1]))
    signs = set()
    for box in completion_boxes(boxes,side):
        if any(_edge_hits_box(p,q,box) for p,q in edges):
            return None
        signs.add(_inside_clear((box[0],box[2]),edges,0))
    return signs.pop() if len(signs)==1 else None


def target_roof():
    mid = Q(23,28)
    return ((0,0),(4,0),(4,SIDE),(2,SIDE),
            (Q(3,2),SIDE+mid-Q(3,2)),(0,SIDE+mid))


def fixture_boxes(depth, lam):
    lam = Q(lam)
    if depth < 0 or not 0 <= lam <= 1:
        raise ValueError("invalid sweep parameter")
    r = 3*WIDTH/2**depth
    return [(x+lam/2-r,x+lam/2+r,y-r,y+r) for x,y in CENTERS[:3]]


def fixture_cutoff(depth):
    if depth < 0:
        raise ValueError("negative depth")
    return (12+Q(13,2)*EPS/2**depth)/7


def crossing_times(lam):
    """All distinct roots of the full-mask schedule; a turning-point root is counted once."""
    lam = Q(lam)
    if not 0 <= lam <= 1:
        raise ValueError("invalid sweep parameter")
    levels = (Q(0),Q(9,10)+lam/10,Q(17,20)+lam/10,Q(1))
    root = Q(13,14)
    times = set()
    for i,(a,b) in enumerate(zip(levels,levels[1:])):
        if min(a,b) <= root <= max(a,b):
            times.add((i+(root-a)/(b-a))/3)
    return tuple(sorted(times))
