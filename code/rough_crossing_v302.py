"""Exact margins for the protected near-square continuum on the v301 rough family."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from balanced_sweep_v300 import profile_value
from polygon_orbit_census_v223 import _affine_system, _evaluate, _solve, add, dot, mul, pt, sub, validate_jordan_polygon
from tilted_sweep_v301 import BERNSTEIN, SAMPLE_TIMES, WINNER, audit_charts, determinant, seed, winning_tags

R_IN, R_OUT = Q(29, 100), Q(69, 100)
T_MIN = Q(99, 100)
ROUGH, SHEAR, ERROR = Q(1, 32), Q(1, 100), Q(1, 1000)
DISP = 4 * ROUGH / 3 + SHEAR + ERROR
SIDE = Q(407, 1000)
PACKET = Q(1, 12)
ORIGIN = pt(Q(1, 2), Q(1, 8))
RECORD = Path(__file__).resolve().parents[1] / "results/v302/rough-crossing.json"


def require(ok, reason):
    if not ok:
        raise ArithmeticError(reason)


def height_band(upper):
    """Allowed y minus seed height at the source coordinate u."""
    free = SHEAR + ERROR
    if upper:
        return -free, 4 * ROUGH / 3 + free
    return -2 * ROUGH / 3 - free, free


def annulus():
    r, R, t = R_IN, R_OUT, T_MIN
    require((r + DISP)**2 < min(Q(4, 25), Q(25, 208)) and Q(5, 8) + DISP < R,
            "seed containment")
    clearance = r*r + 2*t*r*R - R*R
    require(clearance > 0 and 4*t*t > 1 + t*t, "annular separation")
    side_sq = 4*r*r*t*t / (1 + t*t)
    width = 2*r*t / (1 + t*t)
    require(side_sq > SIDE**2 and width > 2*DISP, "size and triple orientation")
    require(SIDE > 2*DISP and (SIDE - 2*DISP)**2 > 13*PACKET**2, "source spacing")
    return dict(wall_clearance=clearance, side_squared=side_sq, triple_width=width,
                source_spacing_squared_margin=(SIDE-2*DISP)**2-13*PACKET**2)


def orientation_gate():
    t = T_MIN
    slope = (1-t)/(1+t)
    a_min, b_max = Q(289, 1000), Q(1, 250)
    require(a_min**2*(1+slope**2) < R_IN**2 and R_OUT*slope < b_max, "tip positions")
    up, down = height_band(True), height_band(False)
    m, n = Q(3, 2), Q(3, 4)
    # Allow source-to-position error and then the horizontal offset of a vertical tip.
    band = lambda raw, speed, shift: (raw[0]-speed*shift, raw[1]+speed*shift)
    at_x = band(up, m, ERROR)
    top = band(at_x, m, b_max)
    bottom = band(band(down, n, ERROR), n, b_max)
    U, L = (Q(-19, 1000), Q(61, 1000)), (Q(-36, 1000), Q(15, 1000))
    require(U[0] <= top[0] <= top[1] <= U[1] and L[0] <= bottom[0] <= bottom[1] <= L[1],
            "vertical tip height envelopes")
    require(2*a_min > max(U[1]-U[0], L[1]-L[0], down[1]-up[0]), "vertical tip assignment")
    a_lo, a_hi = (U[0]-L[1])/2, (U[1]-L[0])/2
    mid_lo = (U[0]+L[0])/2
    tent_lo = (a_min-a_hi)/((m+n)/2)
    require(tent_lo > Q(21, 100) and (m-n)/2*tent_lo+mid_lo-b_max > down[1],
            "horizontal tips lie on upper branch")
    diff = (2*b_max+at_x[1]-at_x[0])/m
    require(diff < min(2*a_min, Q(1, 20)), "horizontal tips straddle the crest")
    a_forced = (m+n)/2*Q(19, 40)+a_lo
    a_cap = Q(1, 2)+ERROR
    require(a_forced > a_cap, "orientation wall")
    return dict(diagonal_offset_ratio=slope, tent_difference_upper=diff,
                forced_half_span=a_forced, available_half_span=a_cap, excess=a_forced-a_cap)


def anchor():
    d = SHEAR
    polygon = seed(d)
    edges = [(polygon[i], polygon[(i+1) % 4]) for i in WINNER]
    tags = winning_tags(d)
    corners = _evaluate(edges, tags)
    edges = [edges[-1], *edges[:-1]]
    a, b, c, e = (corners[-1], *corners[:-1])
    v = sub(b, a)
    center = mul(add(a, c), Q(1, 2))
    offset, diagonal = sub(center, ORIGIN), sub(c, center)
    D = determinant(_affine_system(edges, 1)[0])
    require(D == 9*(4*d*d+12*d+17)/256 and D != 0, "seed transversality")
    require(all(0 < s < 1 for s in tags) and tags[2] < tags[3], "interior anchor")
    require(v.x > abs(v.y) and dot(offset, offset) < dot(diagonal, diagonal), "anchor sector")
    require(sub(c, b) == sub(e, a) == pt(-v.y, v.x), "anchor equations")
    return dict(sorted_charts=len(audit_charts()), determinant=D,
                sector_margin=dot(diagonal, diagonal)-dot(offset, offset),
                horizontal_margin=v.x-abs(v.y))


def rough_chart():
    """A single rational contact branch on a perturbed sixteen-edge rough target."""
    polygon = []
    for i in range(16):
        x = Q(min(i, 16-i), 8)
        height = profile_value(1, x)
        y = (Q(-2, 3) if i <= 8 else Q(4, 3))*height + SHEAR*x
        dx, dy = Q((3*i % 7)-3, 10000), Q((2*i % 5)-2, 10000)
        require(dx*dx + dy*dy <= ERROR**2, "rough polygon perturbation")
        polygon.append(pt(x+dx, y+dy))
    validate_jordan_polygon(polygon)
    chart = (15, 4, 8, 13)
    edges = [(polygon[i], polygon[(i+1) % 16]) for i in chart]
    clearance = Q(1, 100)
    samples = []
    for u in SAMPLE_TIMES:
        k = T_MIN + u*(1/T_MIN-T_MIN)
        matrix, rhs = _affine_system(edges, k)
        det = determinant(matrix)
        tags = _solve(matrix, rhs)
        require(tags is not None, "singular sampled chart")
        # These are the determinant and eight cleared inequalities for Cramer tags.
        samples.append([det, *(det*(v-clearance) for v in tags),
                        *(det*(1-clearance-v) for v in tags)])
    orientation = 1 if samples[0][0] > 0 else -1
    for values in zip(*samples):
        coeffs = _solve(BERNSTEIN, values)
        require(min(orientation*c for c in coeffs) > 0, "rough chart Bernstein inequality")
    square = _evaluate(edges, _solve(*_affine_system(edges, 1)))
    v = sub(square[1], square[0])
    center = mul(add(square[0], square[2]), Q(1, 2))
    offset, diagonal = sub(center, ORIGIN), sub(center, square[0])
    require(v.x > abs(v.y) and dot(offset, offset) < dot(diagonal, diagonal), "rough chart sector")
    return dict(ordered_edges=",".join(map(str, chart)), tag_clearance=clearance,
                positive_bernstein_coefficients=5*len(samples[0]))


def packet_start(s1, red):
    """Circle phase of a packet ending red overlap at the first red/blue junction."""
    return (s1-PACKET*(1-red)) % 1


def report():
    return {name: {key: str(value) for key, value in section().items()}
            for name, section in (("annulus", annulus), ("orientation", orientation_gate),
                                  ("anchor", anchor), ("rough_example", rough_chart))}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = json.dumps(report(), indent=2, sort_keys=True) + "\n"
    if args.write:
        RECORD.write_text(payload, encoding="utf-8")
    elif RECORD.read_text(encoding="utf-8") != payload:
        raise SystemExit("v302 record differs from exact replay")
    print("Verified the annulus, orientation exclusion, seed degree inputs, rough polygon branch and source spacing.")
