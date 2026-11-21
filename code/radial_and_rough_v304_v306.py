"""Rational bounds for the radial four-contact chart and thin rough-strip target."""
from fractions import Fraction as Q

from polygon_orbit_census_v223 import _solve
from rectangle_clearance_v305 import T, matrix

ETA = Q(1,100000)
WEIGHTS = (Q(1), Q(-3,4), Q(1,2), Q(-1,4))
XS = (0,1,1,0)
YS = (0,0,1,1)


def radial_bounds():
    radius = Q(989,2000)
    rmax = Q(1011,2000)
    speed_r = Q(21,1000)
    acceleration = Q(64,100)
    bend = Q(1,100)
    cos_phi = (1-T*T)/(1+T*T)
    leakage = rmax*(cos_phi+bend)+speed_r
    midpoint_pivot = radius-leakage-acceleration*bend
    diagonal_error = 2*ETA+acceleration*(5*ETA)**2
    radial_error = diagonal_error+speed_r*5*ETA
    angular_wall = (1-cos_phi)*(Q(1,2000)-Q(1,50)*cos_phi)
    short_sq = 4*(radius-ETA-acceleration*(5*ETA)**2/2)**2*T*T/(1+T*T)
    if not (rmax+2*speed_r+Q(82,1000)<acceleration and midpoint_pivot>Q(2,5)
            and diagonal_error < radius*Q(99,100)*5*ETA
            and angular_wall > radial_error and short_sq > Q(69,100)**2):
        raise ArithmeticError("radial box lost its clearance")
    return dict(midpoint_pivot=midpoint_pivot, diagonal_error=diagonal_error,
                angular_wall=angular_wall, radial_error=radial_error, side_squared=short_sq)


def wave(x, depth):
    value = Q(0)
    for j in range(1,depth+1):
        phase = (4**j*x) % 1
        value += 2*min(phase,1-phase)/2**j
    return value


def rough_shift(i,x,depth):
    distance = abs(x-XS[i])
    cutoff = min(Q(1),max(Q(0),12*(Q(1,4)-distance)))
    plateau_x = min(XS[i]+Q(1,6),max(XS[i]-Q(1,6),x))
    return cutoff*WEIGHTS[i]*(1+wave(plateau_x,depth))/2048


def local_height(i,x,depth):
    dx = x-XS[i]
    slope = (Q(1,32) if dx>0 else Q(0), Q(1,4), Q(-1,2), Q(1,2))[i]
    return YS[i]+slope*dx+rough_shift(i,x,depth)


def rough_witness(k,depth):
    """Solve candidate local affine charts and verify all four actual graph contacts."""
    xvals = list(map(Q,XS))
    step = Q(1,2*4**depth)
    for _ in range(24):
        slopes, offsets = [], []
        for i,x in enumerate(xvals):
            l = (x//step)*step
            # On the inner plateau the wave and the base are affine on each grid cell.
            l,r = max(l,XS[i]-Q(1,6)), min(l+step,XS[i]+Q(1,6))
            if not l<=x<=r or r<=l:
                raise ArithmeticError("chart left plateau")
            m = (local_height(i,r,depth)-local_height(i,l,depth))/(r-l)
            slopes.append(m)
            offsets.append(local_height(i,l,depth)+m*(XS[i]-l)-YS[i])
        offsets[2] += 1-k
        offsets[3] += 1-k
        x,y,h,v = _solve(matrix(slopes,k),offsets)
        xvals = (x,x+1+h,x+1+h-k*v,x-k*v)
        heights = (y,y+v,y+v+k*(1+h),y+k*(1+h))
        if all(abs(z-u)<Q(1,6) and local_height(i,z,depth)==b
               for i,(z,u,b) in enumerate(zip(xvals,XS,heights))):
            return (x,y,h,v)
    raise ArithmeticError("sample chart did not settle")


def verification_record():
    return {
        "radial_margins": {name: str(value) for name, value in radial_bounds().items()},
        "rough_contacts": [
            {"depth": n, "ratio": str(k), "x_y_h_v": list(map(str, rough_witness(k,n)))}
            for n in range(1,5) for k in (T,Q(1),1/T)
        ],
    }


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = Path(__file__).resolve().parents[1] / "results/v304-v306/radial-and-rough.json"
    payload = json.dumps(verification_record(), indent=2) + "\n"
    if args.write:
        output.write_text(payload, encoding="utf-8")
    elif output.read_text(encoding="utf-8") != payload:
        raise SystemExit("radial bounds or rough contact records differ")
