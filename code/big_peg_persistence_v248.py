"""v248: an exact quantitative NONCOLLAPSE theorem for inscribed squares.

Let four distinct selected segments e_i=A_i+t E_i (0<=t<=1) contain
the vertices of a nondegenerate square, in cyclic order, with parameters
strictly inside their segments. If the 4x4 square-incidence line matrix
M is invertible, its exact inverse yields an explicit epsilon>0 such
that *every* continuous perturbation of the four segments as normal
graphs, with normal displacement bounded by epsilon, still has an
inscribed square of an explicitly positive minimum side length.

The proof uses Brouwer's fixed-point theorem on a rational 4-cube,
not differentiation of the perturbed arcs. The certificate uses only
Fraction; no interval or floating-point decisions. This is a local
theorem for a controlled family, not a uniform Big Peg constant for all
Jordan curves. In particular it does not solve Toeplitz.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from pathlib import Path
import json

from polygon_orbit_census_v223 import pt, add, sub, dot, validate_jordan_polygon
from big_peg_exact_v243 import read_polygon

ROOT=Path(__file__).resolve().parents[1]
OCTAGON_GRID=((39,21),(18,40),(-17,39),(-38,20),
              (-43,-19),(-20,-38),(21,-38),(39,-23))
OCTAGON=tuple(pt(Q(x,20),Q(y,20)) for x,y in OCTAGON_GRID)
OCTAGON_ASSIGNMENT=(0,2,4,6)
ORIGINAL_ASSIGNMENT=(2,4,5,6)
V247_ASSIGNMENT=(1,4,5,6)

# C_k(z)=C_x[k] . z, C_y[k] . z, z=(p_x,p_y,s_x,s_y).
# Vertices run around the square P, P+S, P+S+JS, P+JS.
CORNER_COEFFS=(
    ((1,0,0,0),(0,1,0,0)),
    ((1,0,1,0),(0,1,0,1)),
    ((1,0,1,-1),(0,1,1,1)),
    ((1,0,0,-1),(0,1,1,0)),
)


def invert_and_solve(matrix, vector):
    """Exact Gauss-Jordan inverse and root, fail on singular incidence."""
    n=4
    rows=[[Q(x) for x in matrix[i]]+
          [Q(int(i==j)) for j in range(n)]+[Q(vector[i])]
          for i in range(n)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if rows[i][j]),None)
        if pivot is None:
            raise ValueError("singular four-edge square-incidence matrix")
        rows[j],rows[pivot]=rows[pivot],rows[j]
        scale=rows[j][j]
        rows[j]=[x/scale for x in rows[j]]
        for i in range(n):
            if i==j:
                continue
            factor=rows[i][j]
            rows[i]=[u-factor*v for u,v in zip(rows[i],rows[j])]
    inverse=tuple(tuple(rows[i][n+j] for j in range(n)) for i in range(n))
    root=tuple(rows[i][-1] for i in range(n))
    for i in range(n):
        assert sum(matrix[i][j]*root[j] for j in range(n)) == vector[i]
        for k in range(n):
            assert sum(matrix[i][j]*inverse[j][k] for j in range(n))==Q(i==k)
    return inverse,root


def square_corner(root,k):
    """Exact corner of a (counterclockwise) Euclidean square."""
    cx,cy=CORNER_COEFFS[k]
    return pt(sum(Q(cx[j])*root[j] for j in range(4)),
              sum(Q(cy[j])*root[j] for j in range(4)))


@dataclass(frozen=True)
class PersistenceCertificate:
    assignment: tuple
    epsilon: Q
    root: tuple
    parameters: tuple
    parameter_sensitivities: tuple
    normal_response_bound: Q
    cube_radius: Q
    minimum_square_side_squared: Q
    original_square_side_squared: Q

    @property
    def nondegenerate(self):
        return self.minimum_square_side_squared>0


def certify_persistence(vertices, assignment, epsilon):
    """Prove square survives ALL continuous normal-graph perturbations.

    A selected segment A_i+t E_i may be replaced by
        A_i+t E_i + g_i(t) J(E_i)/||E_i||
    with continuous real g_i:[0,1]->R and |g_i(t)|<=epsilon.
    The four perturbed arcs need not be differentiable; if they lie on
    a Jordan curve, the resulting square is boundary-inscribed there.

    This routine checks the exact sufficient conditions for a *uniform*
    square-side bound over every such choice of four g_i. The other
    portions of the Jordan curve are unrestricted.
    """
    P=validate_jordan_polygon(tuple(vertices))
    if type(epsilon) is not Q or epsilon<=0:
        raise ValueError("epsilon must be an exact positive Fraction")
    if len(assignment)!=4 or any(type(i) is not int or i<0 or i>=len(P)
                                 for i in assignment):
        raise ValueError("choose exactly four valid ordered edge indices")
    if len(set(assignment))!=4:
        raise ValueError("this certificate requires four distinct edges")
    M=[];b=[];starts=[];directions=[]
    for k,i in enumerate(assignment):
        a=P[i]
        e=sub(P[(i+1)%len(P)],a)
        n=pt(-e.y,e.x)
        cx,cy=CORNER_COEFFS[k]
        M.append(tuple(n.x*cx[j]+n.y*cy[j] for j in range(4)))
        b.append(dot(n,a))
        starts.append(a)
        directions.append(e)

    inverse,root=invert_and_solve(M,b)
    parameters=[];sensitivities=[]
    for k,(a,e) in enumerate(zip(starts,directions)):
        corner=square_corner(root,k)
        ee=dot(e,e)
        t=dot(e,sub(corner,a))/ee
        if not (Q(0)<t<Q(1)):
            raise ValueError("square corner touches/misses an edge endpoint")
        parameters.append(t)
        cx,cy=CORNER_COEFFS[k]
        sensitivity=sum(abs(e.x*cx[j]+e.y*cy[j])/ee
                        for j in range(4))
        sensitivities.append(sensitivity)

    # For normal displacement g_i, the incidence equation shifts
    # by h_i(t)=g_i(t)*||E_i||. We bound ||E_i|| by |E_ix|+|E_iy|,
    # keeping every decision rational and allowing arbitrary continuity.
    norms=tuple(abs(e.x)+abs(e.y) for e in directions)
    response=max(sum(abs(inverse[j][i])*norms[i] for i in range(4))
                 for j in range(4))
    cube_radius=2*epsilon*response
    if not all(cube_radius*l<min(t,Q(1)-t)
               for t,l in zip(parameters,sensitivities)):
        raise ValueError("normal-graph square chart is not stable at this epsilon")

    # Each square side-vector component changes by at most cube_radius.
    side_lower2=sum(max(abs(root[j])-cube_radius,Q(0))**2
                    for j in (2,3))
    if side_lower2<=0:
        raise ValueError("no strictly positive square size certified")

    return PersistenceCertificate(
        tuple(assignment),epsilon,root,tuple(parameters),
        tuple(sensitivities),response,cube_radius,side_lower2,
        root[2]**2+root[3]**2
    )


def report(c):
    """Exact output. Approximate displays are not used in the theorem."""
    return {
        "assignment":list(c.assignment),
        "epsilon_normal_distance":str(c.epsilon),
        "root_square_z":[str(v) for v in c.root],
        "edge_parameters":[str(v) for v in c.parameters],
        "edge_parameter_sensitivities":[str(v) for v in c.parameter_sensitivities],
        "normal_response_bound":str(c.normal_response_bound),
        "fixed_point_cube_radius":str(c.cube_radius),
        "original_square_side_squared":str(c.original_square_side_squared),
        "guaranteed_perturbed_square_side_squared":str(c.minimum_square_side_squared),
        "guarantee":"For EVERY four continuous normal-graph perturbations with sup norm <= epsilon, a square persists with side squared >= the exact bound.",
    }


def audit():
    """Reproduce a robust rough-curve family and two fragile Big Peg examples."""
    octagon=certify_persistence(OCTAGON,OCTAGON_ASSIGNMENT,Q(1,512))
    assert octagon.root == (
        Q(21561,15136), Q(808131,529760),
        Q(-44719,15136), Q(-20205,105952))
    assert octagon.parameters == (
        Q(13257,26488), Q(8577,13244),
        Q(18673,26488), Q(8317,13244))
    assert octagon.cube_radius == Q(367709,3390464)
    assert octagon.minimum_square_side_squared == Q(
        46593827703305,5747623067648)
    # Any Jordan curve lying inside the radius-3 ball and containing
    # these four perturbed arcs has r_in<=3. Thus J>9/10:
    assert octagon.minimum_square_side_squared>Q(27,10)**2

    original,_,_=read_polygon(
        ROOT/"research"/"v244-pochuev-original-ten-gon.json")
    v247,_,_=read_polygon(
        ROOT/"research"/"v247-refined-ten-gon.json")
    poc=certify_persistence(original,ORIGINAL_ASSIGNMENT,Q(1,10**6))
    fine=certify_persistence(v247,V247_ASSIGNMENT,Q(1,10**9))
    # The maximizer on the optimized decagon nearly meets an endpoint:
    assert min(min(t,1-t) for t in fine.parameters)<Q(1,10**7)
    try:
        certify_persistence(v247,V247_ASSIGNMENT,Q(1,10**8))
    except ValueError:
        pass
    else:
        raise AssertionError("fragile near-endpoint chart unexpectedly certified")

    return {
        "status":"PASS: exact quantitative persistence for continuous rough normal graphs",
        "robust_octagon":report(octagon),
        "restricted_rough_curve_class_big_peg_bound":"J > 9/10 if the Jordan curve stays inside the radius-3 ball and contains the four epsilon-perturbed octagon arcs",
        "pochuev_original":report(poc),
        "v247_near_endpoint":report(fine),
        "scope":"LOCAL theorem only: does not provide a positive universal c for arbitrary Jordan curves or resolve Toeplitz.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True))
