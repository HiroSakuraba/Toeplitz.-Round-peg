"""Exact finite certificates and analytic-model audits for research v269-v273.

This module contains *restricted* ellipse, combinatorial sector, and polygon
certificates. It does not compute a global rectangle cobordism, a Floer complex,
or an all-Jordan obstruction. Mathematical assertions that cannot be decided by
finite arithmetic (Brouwer degree, Floer realization, limit compactness) are
proved or carefully delimited in the accompanying research notes.
"""
from fractions import Fraction as Q
from itertools import permutations
from math import cos, pi, sin


def _q(x):
    return Q(x)


def _det(m):
    """Exact determinant by Leibniz formula for small square rational matrices."""
    n = len(m)
    if n == 0 or any(len(r) != n for r in m):
        raise ValueError("nonempty square matrix required")
    total = Q(0)
    for p in permutations(range(n)):
        inv = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = Q(-1 if inv % 2 else 1)
        for i, j in enumerate(p):
            term *= m[i][j]
        total += term
    return total


def _vadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _vsub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _dot(a, b):
    return a[0]*b[0] + a[1]*b[1]


def _J(a):
    return (-a[1], a[0])


def square_residual(points):
    """Two vector equations for a CCW square, in four scalar coordinates."""
    p0, p1, p2, p3 = points
    edge = _vsub(p1, p0)
    first = _vsub(_vsub(p2, p1), _J(edge))
    second = _vsub(_vsub(p3, p0), _J(edge))
    return (*first, *second)


def ellipse_square_jacobian(a, b):
    """Exact determinant of dF/d(phi0..phi3) at ellipse's positive square.

    gamma(phi)=(a cos(phi), b sin(phi)); s=ab/sqrt(a²+b²).
    Each column has a common factor s. Removing these yields the
    rational normalized tangent vectors used below; s^4 is rational.
    """
    a, b = _q(a), _q(b)
    if a <= 0 or b <= 0 or a == b:
        raise ValueError("positive noncircular ellipse axes required")
    t = [(-a/b, b/a), (-a/b, -b/a), (a/b, -b/a), (a/b, b/a)]
    # F1=p2-p1-J(p1-p0); F2=p3-p0-J(p1-p0).
    cols = []
    for j, v in enumerate(t):
        if j == 0:
            x, y = _J(v), _vsub(_J(v), v)
        elif j == 1:
            x, y = _vsub((-v[0], -v[1]), _J(v)), (-_J(v)[0], -_J(v)[1])
        elif j == 2:
            x, y = v, (Q(0), Q(0))
        else:
            x, y = (Q(0), Q(0)), v
        cols.append((*x, *y))
    jac = [[cols[j][i] for j in range(4)] for i in range(4)]
    s2 = a*a*b*b/(a*a+b*b)
    exact = _det(jac) * s2*s2
    simplified = -4*a*a*b*b*(a*a-b*b)/(a*a+b*b)
    if exact != simplified:
        raise AssertionError("independent determinant and formula disagree")
    return exact


def ellipse_relative_certificate(a, b):
    a, b = _q(a), _q(b)
    determinant = ellipse_square_jacobian(a, b)
    square_half_side_sq = a*a*b*b/(a*a+b*b)
    return {
        "flat_end_dot_signs": (-4*a*a, 4*b*b),
        "square_side_squared": 4*square_half_side_sq,
        "square_circumradius_squared": 2*square_half_side_sq,
        "local_jacobian_det": determinant,
        "local_degree_sign": 1 if determinant > 0 else -1,
        "global_ellipse_rectangle_locus": "open interval (proved analytically)",
        "scope": "relative endpoint and local C0 stability near one noncircular ellipse",
    }


def relative_interval_parity(f0, f1):
    """Mod-2 simple-zero parity *if* endpoint values are nonzero and zeros transverse."""
    f0, f1 = _q(f0), _q(f1)
    if f0 == 0 or f1 == 0:
        raise ValueError("boundary section is zero: relative parity is undefined")
    return int((f0 > 0) != (f1 > 0))


def cyclic_dot_character(k):
    if not isinstance(k, int):
        raise TypeError("integer cyclic relabeling required")
    return -1 if k & 1 else 1


def ellipse_angle_certificate(a, b, c):
    """Full noncircular-ellipse angle locus: c=cos(angle of oriented diagonals).

    All results describe squared distances, so no irrational coordinate
    or angle approximation is used. The angle itself lies in (0, pi).
    """
    a, b, c = _q(a), _q(b), _q(c)
    if a <= 0 or b <= 0 or a == b or not (-1 < c < 1):
        raise ValueError("positive unequal axes and -1 < cos(theta) < 1 required")
    den = b*b*(1-c) + a*a*(1+c)
    x2 = a*a*b*b*(1-c)/den
    y2 = a*a*b*b*(1+c)/den
    r2 = x2+y2
    assert x2/a**2 + y2/b**2 == 1
    assert (y2-x2)/r2 == c
    return {"cos_diagonal_angle":c,"x_squared":x2,"y_squared":y2,
            "circumradius_squared":r2,"side_squared":(4*x2,4*y2)}


def ellipse_angle_band_floors(a, b, eta):
    """Exact attained lower bounds for |cos theta|<=eta<1 on ellipse."""
    a, b, eta = _q(a), _q(b), _q(eta)
    if not 0 <= eta < 1 or a <= 0 or b <= 0 or a == b:
        raise ValueError("noncircular ellipse and 0 <= eta < 1 needed")
    dmax = a*a+b*b+abs(a*a-b*b)*eta
    r_floor = 2*a*a*b*b/dmax
    x_floor = ellipse_angle_certificate(a,b,eta)["x_squared"]
    y_floor = ellipse_angle_certificate(a,b,-eta)["y_squared"]
    return {"circumradius_squared_floor":r_floor,
            "minimum_side_squared_floor":4*min(x_floor,y_floor)}


def _all_pairings():
    return (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def _rectangle_pairing_certificate(p, pairing):
    (i,j),(k,l) = pairing
    d0, d1 = _vsub(p[j],p[i]), _vsub(p[l],p[k])
    return {
        "shared_midpoint":_vadd(p[i],p[j]) == _vadd(p[k],p[l]),
        "equal_diagonal_squared":_dot(d0,d0) == _dot(d1,d1),
        "perpendicular":_dot(d0,d1) == 0,
        "dot":_dot(d0,d1),
    }


def permutation_sector_ledger(a=Q(3), b=Q(4)):
    """All 24 labellings / C4 = 6 sectors; exactly one pairing per ordering.

    Does *not* enumerate the geometric rectangle locus of a Jordan curve.
    """
    a,b = _q(a),_q(b)
    if a <= 0 or b <= 0:
        raise ValueError("positive side half-lengths required")
    corners = ((a,b),(-a,b),(-a,-b),(a,-b))
    reps = {}
    rows = []
    for p in permutations(range(4)):
        orbit = tuple(p[i:]+p[:i] for i in range(4))
        representative = min(orbit)
        reps.setdefault(representative, []).append(p)
        q = tuple(corners[i] for i in p)
        matches = [j for j, pair in enumerate(_all_pairings()) if
                   (lambda z: z["shared_midpoint"] and z["equal_diagonal_squared"])(
                     _rectangle_pairing_certificate(q,pair))]
        if len(matches)!=1:
            raise AssertionError("missing or duplicated diagonal pairing")
        rows.append({"ordering":p,"cyclic_representative":representative,
                     "opposite_pairing_index":matches[0]})
    if len(reps)!=6 or any(len(v)!=4 for v in reps.values()):
        raise AssertionError("not six free C4 orbits")
    counts={k:sum(r["opposite_pairing_index"]==k for r in rows) for k in range(3)}
    if sorted(counts.values()) != [8,8,8]:
        raise AssertionError("pairing coverage not balanced")
    return {"six_sectors":len(reps),"labellings":len(rows),
            "orbit_sizes":tuple(sorted(len(v) for v in reps.values())),
            "pairing_counts":counts,"rows":rows}





def dihedral_sector_ledger():
    """Unrestricted 24 orderings modulo all eight D4 word symmetries.

    This is a *combinatorial* three-class quotient, not a claim that
    all six curve-order sectors are geometrically interchangeable.
    """
    d={}
    for p in permutations(range(4)):
        reflect=tuple(reversed(p))
        orbit={q[i:]+q[:i] for q in (p,reflect) for i in range(4)}
        d.setdefault(min(orbit),set()).add(p)
    if len(d)!=3 or any(len(z)!=8 for z in d.values()):
        raise AssertionError("D4 must give three free eightfold orbits")
    return {"dihedral_classes":3,"eightfold_orbit_sizes":(8,8,8),
            "representatives":tuple(sorted(d))}


def cyclic_square_residual_equivariance():
    """Exact universal linear D4 test on all 8 Cartesian coordinate axes.

    The four-cycle acts by a source matrix of determinant -1; the
    displayed target R has determinant +1. Therefore a transverse
    *labeled* local square index changes sign under T. This is no
    global degree theorem for an arbitrary curve.
    """
    R = ((-1,1,1,0),(-1,-1,0,1),(0,1,0,0),(-1,0,0,0))
    if _det(R)!=1:
        raise AssertionError("target action has wrong determinant")
    for k in range(8):
        v=[Q(0)]*8
        v[k]=Q(1)
        points=tuple((v[2*i],v[2*i+1]) for i in range(4))
        f=square_residual(points)
        shifted=square_residual(points[1:]+points[:1])
        if tuple(sum(R[i][j]*f[j] for j in range(4)) for i in range(4))!=shifted:
            raise AssertionError("target action fails on Cartesian basis")
    def matmul(A,B):
        return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4))
                           for j in range(4)) for i in range(4))
    power=tuple(tuple(Q(i==j) for j in range(4)) for i in range(4))
    for i in range(4):
        power=matmul(power,R)
    if power!=tuple(tuple(Q(i==j) for j in range(4)) for i in range(4)):
        raise AssertionError("R^4 != I")
    source4cycle=tuple(tuple(Q(j==(i+1)%4) for j in range(4)) for i in range(4))
    if _det(source4cycle)!=-1:
        raise AssertionError("source four-cycle must reverse orientation")
    return {"basis_checks":8,"target_determinant":1,
            "source_determinant":-1,"target_order_four":True,
            "labeled_index_character":"(-1)^k"}


def rectangle_boundary_stratum(side1_squared, side2_squared):
    """Classification of any Euclidean rectangle limit with finite side lengths."""
    x, y = _q(side1_squared), _q(side2_squared)
    if x < 0 or y < 0:
        raise ValueError("squared lengths must be nonnegative")
    if x == y == 0:
        return "total_four_corner_collision"
    if x == 0 or y == 0:
        return "two_pair_flattening_nonzero_diagonal"
    return "genuine_noncollapsed_rectangle"


def _orientation(p,q,r):
    v, w = _vsub(q,p),_vsub(r,p)
    return v[0]*w[1]-v[1]*w[0]


def _on_segment(p,q,r):
    return (_orientation(p,q,r)==0 and
            min(p[0],q[0])<=r[0]<=max(p[0],q[0]) and
            min(p[1],q[1])<=r[1]<=max(p[1],q[1]))


def _segments_intersect(p,q,r,s):
    a,b,c,d=(_orientation(p,q,r),_orientation(p,q,s),
               _orientation(r,s,p),_orientation(r,s,q))
    if (a>0 and b<0 or a<0 and b>0) and (c>0 and d<0 or c<0 and d>0):
        return True
    return ((a==0 and _on_segment(p,q,r)) or
            (b==0 and _on_segment(p,q,s)) or
            (c==0 and _on_segment(r,s,p)) or
            (d==0 and _on_segment(r,s,q)))


def is_simple_polygon(vertices):
    n=len(vertices)
    if n<3 or len(set(vertices))!=n:
        return False
    for i in range(n):
        if vertices[i]==vertices[(i+1)%n]:
            return False
        for j in range(i+1,n):
            if j==i+1 or (i==0 and j==n-1):
                continue
            if _segments_intersect(vertices[i],vertices[(i+1)%n],
                                   vertices[j],vertices[(j+1)%n]):
                return False
    return True


def notch_square_certificate(n, anchor=Q(1,3)):
    """Fixed-size Jordan boundaries with a selected square shrinking to zero.

    The upper boundary of the unit square takes an outward square bump of
    height=width=eps. This demonstrates *selected branch shrinkout*, not
    absence of other large squares or a counterexample to Toeplitz.
    """
    if not isinstance(n,int) or n<4:
        raise ValueError("integer n>=4 required")
    anchor=_q(anchor); eps=Q(1,n)
    if anchor<=0 or anchor+eps>=1:
        raise ValueError("outward bump must stay in upper edge interior")
    p=((Q(0),Q(0)),(Q(1),Q(0)),(Q(1),Q(1)),
       (anchor+eps,Q(1)),(anchor+eps,Q(1)+eps),
       (anchor,Q(1)+eps),(anchor,Q(1)),(Q(0),Q(1)))
    if not is_simple_polygon(p):
        raise AssertionError("outward bump is not a Jordan polygon")
    sq=(p[4],p[5],p[6],p[3])
    if square_residual(sq)!=(Q(0),)*4:
        raise AssertionError("selected bump corners not a square")
    blowup = tuple(((x-anchor)/eps,(y-1)/eps) for x,y in sq)
    if blowup!=((1,1),(0,1),(0,0),(1,0)):
        raise AssertionError("nonunit blowup")
    return {"polygon":p,"selected_square":sq,"side_squared":eps*eps,
            "circumradius_squared":eps*eps/2,
            "hausdorff_distance_to_unit_square_upper_bound":eps,
            "normalized_square":blowup,"simple_Jordan_polygon":True}


def spectral_axiom_countermodel(area_over_radius_squared):
    """Analytic countermodel, NOT an actual Floer spectral invariant.

    ell(theta)=A/pi * (theta + sin(2 theta)/2), theta in [0,pi].
    A/R²=B with 0<B<=pi/2 gives the needed R² Lipschitz bound.
    The triangle gap is (2A/pi) sin(x+y) sin(x) sin(y)>=0.
    ell'(pi/2)=0 although ell is strictly increasing.

    Input B is numerical for the simple dimensionless bound only; the
    symbolic trigonometric proof is in research/v270.
    """
    B=float(area_over_radius_squared)
    if not (0<B<=pi/2):
        raise ValueError("countermodel takes 0 < A/R² <= pi/2")
    return {"area_over_radius_squared":B,"maximum_derivative_over_radius_squared":2*B/pi,
            "midpoint_derivative_over_radius_squared":0.0,
            "derivative_normalized":lambda theta: 2*cos(theta)**2,
            "normalized_action":lambda theta: (theta+sin(2*theta)/2)/pi,
            "normalized_triangle_gap":lambda x,y:2*sin(x+y)*sin(x)*sin(y)/pi}


if __name__ == "__main__":
    print("v269 ellipse:",ellipse_relative_certificate(3,4))
    print("v271 six-sector counts:", {k:v for k,v in permutation_sector_ledger().items() if k!='rows'})
    print("v272 collapse:",notch_square_certificate(100)['side_squared'])
    print("v273 angle floors:",ellipse_angle_band_floors(3,4,Q(1,2)))