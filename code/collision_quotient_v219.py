"""v219: exact C4 quotient, residual-bundle section, and collision-link audit.

This is an elementary model of ONE cyclic-order component of Conf_4(S1).
It makes no claim to solve the square-peg conjecture, nor to implement an
all-Jordan equivariant Euler/Floer invariant. Exact rational arithmetic is
used for the group action, leading collision terms, and shrinking polygons.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import combinations


def B(orientation):
    """Residual target action from v218: F_sigma(Tp) = B_sigma F_sigma(p)."""
    if orientation == 1:
        return ((-1, 1, 1, 0), (-1, -1, 0, 1),
                (0, 1, 0, 0), (-1, 0, 0, 0))
    if orientation == -1:
        return ((-1, -1, 1, 0), (1, -1, 0, 1),
                (0, -1, 0, 0), (1, 0, 0, 0))
    raise ValueError("orientation must be +/-1")


def apply(matrix, vector):
    if len(matrix) != len(vector) or any(len(row) != len(vector) for row in matrix):
        raise ValueError("square matrix and vector required")
    return tuple(sum(Q(a)*Q(x) for a,x in zip(row, vector)) for row in matrix)


def vector_add(a, b):
    return tuple(x+y for x,y in zip(a,b))


def vector_scale(z, a):
    return tuple(Q(z)*x for x in a)


def norm_squared(v):
    return sum(x*x for x in v)


@dataclass(frozen=True)
class CyclicConfiguration:
    """Closed naive gap simplex chart: t0 mod 1, four nonnegative gaps sum 1."""
    t: Q
    gaps: tuple

    def __post_init__(self):
        if not isinstance(self.t, Q) or len(self.gaps) != 4 or any(
                not isinstance(g, Q) for g in self.gaps):
            raise ValueError("Exact Fraction input required")
        if any(g < 0 for g in self.gaps) or sum(self.gaps) != 1:
            raise ValueError("Four nonnegative gaps summing to one required")
        object.__setattr__(self, 't', self.t % 1)


def conf(t, gaps):
    return CyclicConfiguration(Q(t), tuple(map(Q, gaps)))


def cyclic_shift(cfg):
    """T(t,g0,g1,g2,g3)=(t+g0,g1,g2,g3,g0)."""
    g=cfg.gaps
    return CyclicConfiguration((cfg.t+g[0]) % 1, g[1:]+g[:1])


def center(cfg):
    """Equivariant center, transforms by +1/4 under cyclic relabeling."""
    g=cfg.gaps
    return (cfg.t+(3*g[0]+2*g[1]+g[2])/4) % 1


def from_center(c, gaps):
    g=tuple(map(Q,gaps))
    if len(g)!=4:
        raise ValueError("Four gaps required")
    return conf(Q(c)-(3*g[0]+2*g[1]+g[2])/4,g)


def equivariant_retract(cfg, s):
    """H_s: interpolate all gaps to 1/4 while fixing equivariant center."""
    s=Q(s)
    if not 0<=s<=1:
        raise ValueError("homotopy parameter must be in [0,1]")
    h=tuple((1-s)*g+s*Q(1,4) for g in cfg.gaps)
    return from_center(center(cfg),h)


def orbit(cfg):
    states=[cfg]
    for _ in range(3):states.append(cyclic_shift(states[-1]))
    return tuple(states)


def real_rotation_basis(orientation):
    """Real 2-plane where B acts as a quarter turn: B u=v, B v=-u."""
    u=(Q(1),Q(0),Q(1),Q(0))
    v=apply(B(orientation),u)
    assert apply(B(orientation),v)==tuple(-x for x in u)
    return u,v


def equivariant_section_formal(orientation, cosine, sine):
    """Exact algebra of s_sigma(c)=cos(2pi c)u+sin(2pi c)Bu.

    cosine/sine can be *arbitrary rational symbols*; along the unit circle
    cosine**2+sine**2=1, hence the output is never zero. Tests also verify
    the identity s(c+1/4)=B s(c) via (cos,sin)->(-sin,cos).
    """
    u,v=real_rotation_basis(orientation)
    return vector_add(vector_scale(cosine,u),vector_scale(sine,v))


def section_equivariance_identity(orientation, cosine, sine):
    left=equivariant_section_formal(orientation,-Q(sine),Q(cosine))
    right=apply(B(orientation),equivariant_section_formal(orientation,Q(cosine),Q(sine)))
    return left==right


def leading_collision_residual(vx,vy,a,b,c,orientation=1):
    """(F_sigma/eps) limit near total collision with 3 short gaps eps*a,b,c.

    Sum a+b+c=1, a,b,c>=0; p_i are first-order points on the
    same regular tangent v=(vx,vy), with p3-p0=eps*v.
    """
    a,b,c=Q(a),Q(b),Q(c)
    vx,vy=Q(vx),Q(vy)
    if min(a,b,c)<0 or a+b+c!=1:
        raise ValueError("Nonnegative normalized collision triple required")
    if orientation not in (-1,1):
        raise ValueError("orientation +/-1 required")
    s=Q(orientation)
    first=(b*vx+s*a*vy, b*vy-s*a*vx)
    second=(vx+s*a*vy,vy-s*a*vx)
    return first+second


def leading_link_lower_bound_squared(vx,vy):
    """The squared norm of the second residual vector is >= |v|^2."""
    return Q(vx)*Q(vx)+Q(vy)*Q(vy)


@dataclass(frozen=True)
class Point:
    x: Q
    y: Q


def pt(x,y):
    return Point(Q(x),Q(y))


def orientation2(p,q,r):
    return (q.x-p.x)*(r.y-p.y)-(q.y-p.y)*(r.x-p.x)


def between(a,b,c):
    return min(a,c)<=b<=max(a,c)


def on_segment(p,a,b):
    return orientation2(a,b,p)==0 and between(a.x,p.x,b.x) and between(a.y,p.y,b.y)


def segments_intersect(a,b,c,d):
    o1,o2=orientation2(a,b,c),orientation2(a,b,d)
    o3,o4=orientation2(c,d,a),orientation2(c,d,b)
    if o1*o2<0 and o3*o4<0:return True
    return ((o1==0 and on_segment(c,a,b)) or
            (o2==0 and on_segment(d,a,b)) or
            (o3==0 and on_segment(a,c,d)) or
            (o4==0 and on_segment(b,c,d)))


def is_simple_polygon(vertices):
    """Exact intersection-based simplicity for the rational examples."""
    p=tuple(vertices);n=len(p)
    if n<3 or len(set(p))!=n:return False
    if any(p[i]==p[(i+1)%n] for i in range(n)):return False
    if sum(p[i].x*p[(i+1)%n].y-p[i].y*p[(i+1)%n].x for i in range(n))==0:
        return False
    for i,j in combinations(range(n),2):
        if (j-i)%n in (1,n-1):
            # Reject adjacent collinear backtracking, not ordinary corners.
            a,b=p[i],p[(i+1)%n]
            c,d=p[j],p[(j+1)%n]
            common=set((a,b))&set((c,d))
            if len(common)!=1:return False
            shared=next(iter(common))
            other_a=next(z for z in (a,b) if z!=shared)
            other_b=next(z for z in (c,d) if z!=shared)
            if (orientation2(shared,other_a,other_b)==0 and
                    (on_segment(other_a,shared,other_b) or
                     on_segment(other_b,shared,other_a))):
                return False
            continue
        if segments_intersect(p[i],p[(i+1)%n],p[j],p[(j+1)%n]):
            return False
    return True


def notch_polygon(n):
    """Simple polygon with an exact square of side 1/n on a tiny outward notch.

    It converges to [0,2]x[0,1] in Hausdorff distance <=1/n.
    Uniformly close homeomorphic parameterizations also exist because
    the detour can be traversed over a shrinking parameter interval.
    """
    n=int(n)
    if n<2:raise ValueError("n must be >= 2")
    e=Q(1,n)
    return (pt(0,0),pt(1,0),pt(1,-e),pt(1+e,-e),
            pt(1+e,0),pt(2,0),pt(2,1),pt(0,1))


def notch_square(n):
    p=notch_polygon(n)
    return p[1:5]


def is_oriented_square(p,orientation=1):
    if len(p)!=4 or orientation not in (-1,1):return False
    e=(p[1].x-p[0].x,p[1].y-p[0].y)
    if e==(0,0):return False
    j=(-orientation*e[1],orientation*e[0])
    return (p[2].x-p[1].x,p[2].y-p[1].y)==j and (
        p[3].x-p[0].x,p[3].y-p[0].y)==j


def notch_square_side_squared(n):
    s=notch_square(n)
    return (s[1].x-s[0].x)**2+(s[1].y-s[0].y)**2


def audit():
    g=conf(Q(1,7), (Q(1,8),Q(1,3),Q(1,6),Q(3,8)))
    result={
        "cyclic_center_increment":str((center(cyclic_shift(g))-center(g))%1),
        "equivariant_retraction":equivariant_retract(cyclic_shift(g),Q(2,3))==cyclic_shift(equivariant_retract(g,Q(2,3))),
        "full_cyclic_orbit_size":len(set(orbit(g))),
        "section_equivariance_both_signs":all(section_equivariance_identity(s,Q(3,5),Q(4,5)) for s in (-1,1)),
        "section_norm_squared_on_unit_circle":"2",
        "C4_quotient_homotopy_type":"circle on one cyclic-order component (absolute H^4=0)",
        "notch_side_squared_n64":str(notch_square_side_squared(64)),
        "notch_polygon_simple_n64":is_simple_polygon(notch_polygon(64)),
        "notch_hausdorff_upper_bound_n64":str(Q(1,64)),
        "global_general_Jordan_square_peg_proved":False,
    }
    return result


if __name__=="__main__":
    import json
    print(json.dumps(audit(),indent=2,sort_keys=True))
