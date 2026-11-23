"""Exact fixed-axis square-annulus selection for a column polygon, using Matschke 2.1."""
from dataclasses import dataclass
from fractions import Fraction as Q
import csv
import sys

from terrain_gluing_v334 import column_maps
from terrain_replay_v341 import fixtures


@dataclass(frozen=True)
class Frame:
    radius: Q
    center: tuple
    outer_radius: Q

    @property
    def registration_margin(self):
        # Epsilon below this keeps the enlarged outer / reduced inner ratio < 3.
        return (3*self.radius-self.outer_radius)/4


def _interval_frame(left, right, bottom, top, bounds):
    radius = Q(min(right-left,top-bottom),2)
    if radius <= 0:
        return None
    xmin,xmax,ymin,ymax = bounds
    # Closest point of the admissible center rectangle to the outer box midpoint
    # minimizes the largest of the four extreme coordinate distances.
    def clamp(v,lo,hi):
        return min(hi,max(lo,v))
    x = clamp(Q(xmin+xmax,2),left+radius,right-radius)
    y = clamp(Q(ymin+ymax,2),bottom+radius,top-radius)
    outer = max(x-xmin,xmax-x,y-ymin,ymax-y)
    return Frame(radius,(x,y),outer) if outer <= 3*radius else None


def best_frame(left, right, upper, lower):
    bounds = left,right,min(lower.values()),max(upper.values())
    best = None
    for a in range(left,right):
        top,bottom = upper[a],lower[a]
        for b in range(a+1,right+1):
            top = min(top,upper[b-1])
            bottom = max(bottom,lower[b-1])
            candidate = _interval_frame(a,b,bottom,top,bounds)
            if candidate and (best is None or (candidate.radius,-candidate.outer_radius)
                               > (best.radius,-best.outer_radius)):
                best = candidate
    return best


def write_report(stream):
    writer = csv.writer(stream,lineterminator="\n")
    writer.writerow(("kind","seed","i","radius","center_x","center_y","outer_radius",
                     "registration_margin","at_target"))
    for kind,seed,upper,lower,start,i,_ in fixtures():
        if lower:
            columns = column_maps(upper,lower,0,start)
        else:
            columns = 0,len(upper),dict(enumerate(upper)),dict.fromkeys(range(len(upper)),0)
        frame = best_frame(*columns)
        if frame is None:
            writer.writerow((kind,seed,i,"","","","","","none"))
            continue
        status = "below"
        if 2*frame.radius == i:
            status = "strict" if frame.registration_margin > 0 else "equality"
        writer.writerow((kind,seed,i,frame.radius,*frame.center,frame.outer_radius,
                         frame.registration_margin,status))


if __name__ == "__main__":
    write_report(sys.stdout)
