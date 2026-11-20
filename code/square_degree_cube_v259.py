"""v259: exact Poincare-Miranda / local Brouwer-degree square certificate.

Four continuous boundary arcs p_k(s_k), s_k in [-rho,rho], give the
four-variable square-defect map F(s):
 F_0:1 = p2-p1-J(p1-p0)
 F_2:3 = p3-p0-J(p1-p0).
A zero means the four actual boundary points are a counterclockwise square.

The reference arcs are continuous piecewise affine and MAY have genuine
polygon corners. No differentiability at the kink, no individual Cramer
square root, and no individual cell-Jacobian invertibility is assumed.
Choose an invertible rational PRECONDITIONER B, not necessarily a square
incidence matrix. Prove strict signs of G=B F on all 8 cube faces by
checking EXACT rational vertices of every affine cell on each face.

If all face signs point OUTWARD, the homotopy (1-t)G+t*identity has no
boundary zeros, so deg(G,cube,0)=+1. This is the Poincare-Miranda theorem
and guarantees at least one square, without selecting a 4-edge chart.
The same survives ANY continuous coordinatewise arc displacements
bounded by epsilon if B's exact global perturbation response times
epsilon is smaller than the verified face margin.

The example uses the ACTUAL noncollinear vertex from v250's convex
9-gon and additionally a genuine nearby ASYMMETRIC Jordan polygon
whose reference arc centers do NOT form an exact square. Its degree
certificate crosses the real corner and proves an exact, nonzero
inscribed square for all allowed continuous arc perturbations, subject
to those four arcs belonging to one Jordan boundary.

Limit: LOCAL parameter cube only; not a full global degree invariant
under arbitrary polygon deformations or a proof of Toeplitz.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
import json

from polygon_orbit_census_v223 import (
    pt, add, sub, mul, dot, cross, _on_seg, validate_jordan_polygon,
    census_polygon,
)
from big_peg_persistence_v248 import (
    square_corner, invert_and_solve,
)
from big_peg_corner_v250 import (
    family_vertices as moving_corner_family, prove_jordan_family,
    verify_three_edge_kernel,
)
from big_peg_multichart_v249 import HALF, ppoint_at
from central_symmetry_square_v256 import J

RADIUS=Q(1,100)
EPSILON=Q(1,100000)
ASYMMETRIC_EDGE_SHIFT=pt(Q(1,20000),-Q(1,25000))
ACTIVE_EDGES=(3,5,7)


@dataclass(frozen=True)
class PiecewiseArc:
    """Continuous PL arc, with a genuine kink allowed at parameter 0."""
    center: object
    negative_slope: object
    positive_slope: object

    def point(self,s):
        if type(s) is not Q:
            raise ValueError("exact rational arc parameter required")
        return add(self.center,mul(
            self.negative_slope if s<=0 else self.positive_slope,s))


def square_defect(points):
    """Zero iff the 4 ordered corners form a CCW square (possibly zero)."""
    if len(points)!=4:
        raise ValueError("four ordered square candidates required")
    p0,p1,p2,p3=points
    u=sub(p1,p0)
    x=sub(sub(p2,p1),J(u))
    y=sub(sub(p3,p0),J(u))
    return (x.x,x.y,y.x,y.y)


def column_for_arc(k,e):
    zero=pt(0,0)
    v=[zero]*4
    v[k]=e
    return square_defect(v)


def _invert4(matrix):
    """Exact 4x4 preconditioner inverse; handles no near-singular floats."""
    inverse,root=invert_and_solve(matrix,(Q(0),)*4)
    if any(root):
        raise AssertionError("homogeneous exact matrix must have zero root")
    return inverse


def _precondition(inverse,defect):
    return tuple(sum(inverse[i][j]*defect[j] for j in range(4))
                 for i in range(4))


def _average_preconditioner(arcs):
    """One convenient invertible preconditioner. NO cell Jacobian gate."""
    minus=tuple(column_for_arc(k,arcs[k].negative_slope)
                for k in range(4))
    plus=tuple(column_for_arc(k,arcs[k].positive_slope)
               for k in range(4))
    average=tuple(
        tuple((minus[j][i]+plus[j][i])/2 for j in range(4))
        for i in range(4))
    return average


def _face_vertices(axis,sign,radius):
    """Every vertex of the (possibly kink-split) 4D cube face.

    The only potential internal break is s_k=0; including 0 for every
    unpinned parameter makes the enumeration valid for any of 4 kinked
    input arcs, not just the first. Each component of F is affine on
    each corresponding subbox; extrema are therefore at its vertices.
    """
    choices=[(sign*radius,) if k==axis else
             (-radius,Q(0),radius) for k in range(4)]
    return product(*choices)


def _perturbation_response(inverse):
    """Exact L_infty norm for 8 independent coordinate perturbations."""
    zero=pt(0,0)
    effects=[]
    for k in range(4):
        for unit in (pt(1,0),pt(0,1)):
            points=[zero]*4
            points[k]=unit
            effects.append(_precondition(inverse,square_defect(points)))
    return max(sum(abs(effects[j][i]) for j in range(8))
               for i in range(4))


@dataclass(frozen=True)
class CubeDegreeCertificate:
    radius:Q
    continuous_arc_error:Q
    exact_face_margins:tuple
    minimum_face_margin:Q
    perturbation_response:Q
    maximum_transformed_perturbation:Q
    certified_square_side_squared_lower:Q
    reference_center_square_defect:tuple
    preconditioner:tuple
    local_degree:int
    checked_face_vertices:int


def certify_cube_degree(arcs,radius,epsilon,preconditioner=None):
    """A *strict* Poincare-Miranda and Brouwer-degree certificate.

    Exact reference F(s) is continuous and piecewise-affine on 4-cube.
    An arbitrary perturbed arc is reference_arc(s)+h_k(s), with h_k
    continuous and |h_k(s).x|,|h_k(s).y|<=epsilon. Both arbitrary
    nonsmooth perturbation and an actual polygon-vertex kink are allowed.

    For every face s_i= +/- radius, certify +/- (B F)_i>0 on every
    reference subbox vertex. Piecewise affinity extends these signs
    to the ENTIRE face. The perturbed sign remains strict if the
    exact 8-coordinate operator infinity norm times epsilon is smaller
    than every face margin. Thus the degree of the perturbed map is +1.

    This proof NEVER tests determinant of either PL cell Jacobian; it
    uses only ONE optional nonsingular coordinate preconditioner B.
    For the nonzero-side guarantee, compare p1-p0 to its central value
    using all coordinate slope bounds plus 2*epsilon.
    """
    if len(arcs)!=4 or any(not isinstance(a,PiecewiseArc) for a in arcs):
        raise ValueError("four piecewise-affine boundary arcs required")
    if type(radius) is not Q or not Q(0)<radius<Q(1):
        raise ValueError("rational radius must lie strictly between 0 and 1")
    if type(epsilon) is not Q or epsilon<0:
        raise ValueError("nonnegative rational continuous arc error required")
    P=tuple(arcs)
    matrix=(_average_preconditioner(P) if preconditioner is None
            else tuple(tuple(Q(x) for x in row) for row in preconditioner))
    if len(matrix)!=4 or any(len(row)!=4 for row in matrix):
        raise ValueError("preconditioner must be an exact four-by-four matrix")
    inverse=_invert4(matrix)

    def transformed(s):
        return _precondition(
            inverse,square_defect(tuple(a.point(s[k])
                                       for k,a in enumerate(P))))
    faces=[]
    checked=0
    for axis in range(4):
        for sign in (-1,1):
            values=[]
            for s in _face_vertices(axis,sign,radius):
                values.append(sign*transformed(s)[axis])
                checked+=1
            margin=min(values)
            if margin<=0:
                raise ValueError(
                    "strict Poincare-Miranda face sign not certified")
            faces.append(margin)
    global_margin=min(faces)
    response=_perturbation_response(inverse)
    if epsilon*response>=global_margin:
        raise ValueError(
            "continuous boundary perturbation exceeds face-sign margin")

    s0=sub(P[1].center,P[0].center)
    bounds=tuple(
        radius*(max(abs(getattr(P[0].negative_slope,axis)),
                    abs(getattr(P[0].positive_slope,axis)))+
                max(abs(getattr(P[1].negative_slope,axis)),
                    abs(getattr(P[1].positive_slope,axis))))+2*epsilon
        for axis in ("x","y"))
    lower=sum(max(abs(getattr(s0,axis))-bounds[i],Q(0))**2
              for i,axis in enumerate(("x","y")))
    if lower<=0:
        raise ValueError("no positive-size square guaranteed")
    return CubeDegreeCertificate(
        radius,epsilon,tuple(faces),global_margin,response,
        epsilon*response,lower,
        square_defect(tuple(a.center for a in P)),
        matrix,1,checked)


def _parameter_on_segment(q,a,b):
    e=sub(b,a)
    if dot(e,e)<=0:
        raise ValueError("degenerate edge")
    if not _on_seg(q,a,b):
        raise ValueError("central arc contact not on actual polygon edge")
    return dot(e,sub(q,a))/dot(e,e)


def _seam_polygon_and_reference_root():
    family=moving_corner_family()
    prove_jordan_family(family)
    P=tuple(ppoint_at(p,HALF) for p in family)
    validate_jordan_polygon(P)
    _,z=verify_three_edge_kernel()
    if square_corner(z,0)!=P[1]:
        raise AssertionError("square must hit genuine v250 polygon kink")
    if cross(sub(P[1],P[0]),sub(P[2],P[1]))<=0:
        raise AssertionError("v250 kink flattened or reversed")
    return P,z


def build_v250_kink_arcs(shifted=False,radius=RADIUS):
    """Construct four ACTUAL polygon subarcs, crossing a non-collinear kink.

    Optional tiny rational shift moves BOTH endpoints of edge #3 of the
    nine-gon, producing a genuinely asymmetric polygon; the original
    four central points then are only a NEAR-square, not an exact one.
    """
    P,z=_seam_polygon_and_reference_root()
    P=list(P)
    if shifted:
        P[3]=add(P[3],ASYMMETRIC_EDGE_SHIFT)
        P[4]=add(P[4],ASYMMETRIC_EDGE_SHIFT)
    P=validate_jordan_polygon(tuple(P))
    centers=[square_corner(z,k) for k in range(4)]
    if shifted:
        centers[1]=add(centers[1],ASYMMETRIC_EDGE_SHIFT)
    e0minus=sub(P[1],P[0])
    e0plus=sub(P[2],P[1])
    if centers[0]!=P[1]:
        raise AssertionError("first arc is not centered on kink")
    arcs=[PiecewiseArc(P[1],e0minus,e0plus)]
    if type(radius) is not Q or not 0<radius<1:
        raise ValueError("invalid local radius")
    for k,i in enumerate(ACTIVE_EDGES,start=1):
        a,b=P[i],P[(i+1)%len(P)]
        t=_parameter_on_segment(centers[k],a,b)
        if not (radius<t<Q(1)-radius):
            raise ValueError("local square parameter box exits a selected edge")
        e=sub(b,a)
        arcs.append(PiecewiseArc(centers[k],e,e))
    if not (0<radius<1):
        raise ValueError("local kink parameter leaves source segments")
    return P,tuple(arcs)


def _linear_square_witness(arcs,cert):
    """Independent exact piecewise-affine root on one of two kink cells.

    Not used by the DEGREE existence proof. This verifies that a
    concrete exact square can be recovered in our sample geometry.
    """
    if len(arcs)!=4:
        raise ValueError("requires four arcs")
    centers=tuple(a.center for a in arcs)
    const=square_defect(centers)
    roots=[]
    for positive in (False,True):
        slopes=tuple((a.positive_slope if positive else a.negative_slope)
                     for a in arcs)
        columns=tuple(column_for_arc(k,slopes[k]) for k in range(4))
        M=tuple(tuple(columns[k][j] for k in range(4))
                for j in range(4))
        try:
            _,s=invert_and_solve(M,tuple(-x for x in const))
        except ValueError:
            continue  # square may sit in a singular PL cell; degree still valid
        if not all(-cert.radius<s[k]<cert.radius for k in range(4)):
            continue
        if s[0]>0 and not positive or s[0]<0 and positive:
            continue
        points=tuple(arcs[k].point(s[k]) for k in range(4))
        if any(square_defect(points)):
            raise AssertionError("independent square root fails original constraints")
        roots.append((s,points))
    if not roots:
        raise AssertionError("could not extract simple exact witness from degree")
    return roots[0]


def audit():
    base,arcs=build_v250_kink_arcs()
    base_cert=certify_cube_degree(arcs,RADIUS,EPSILON)
    if base_cert.minimum_face_margin!=Q(17387,2696700):
        raise AssertionError("unexpected v250 kink degree margin")
    if base_cert.perturbation_response!=Q(89890,3311):
        raise AssertionError("unexpected reference response")
    if base_cert.certified_square_side_squared_lower<=Q(9):
        raise AssertionError("large square size bound failed")

    target,arcs2=build_v250_kink_arcs(shifted=True)
    target_cert=certify_cube_degree(arcs2,RADIUS,EPSILON)
    if target_cert.minimum_face_margin!=Q(11853176941,1785754740000):
        raise AssertionError("unexpected near-square degree margin")
    if not any(target_cert.reference_center_square_defect):
        raise AssertionError("asymmetric boundary contacts were already an exact square")
    if target_cert.certified_square_side_squared_lower<=Q(9):
        raise AssertionError("asymmetric target large square bound failed")

    # Independent exact root extraction and independent COMPLETE v223
    # all-geometric-square census on actual perturbed polygon.
    parameters,corners=_linear_square_witness(arcs2,target_cert)
    census=census_polygon(target)
    if not any(set(corners)==set(s.vertices) for s in census.squares):
        raise AssertionError("degree-extracted square absent in v223 exhaustive census")
    if len(set(corners))!=4:
        raise AssertionError("collapsed square")
    return {
        "status":"PASS: exact 4D Poincare-Miranda degree +1 across genuine kink, robust under arbitrary continuous small arc perturbations",
        "polygon_time":"v250 moving convex 9-gon at t=1/2",
        "reference_corner_noncollinear":True,
        "local_parameter_cube_radius":str(RADIUS),
        "arbitrary_continuous_coordinatewise_arc_displacement":str(EPSILON),
        "reference_exact_face_sign_margin":str(base_cert.minimum_face_margin),
        "asymmetric_target_exact_face_sign_margin":str(target_cert.minimum_face_margin),
        "preconditioned_perturbation_operator_norm":str(base_cert.perturbation_response),
        "maximum_continuous_residual_perturbation":str(base_cert.maximum_transformed_perturbation),
        "strict_square_side_squared_lower_reference":str(base_cert.certified_square_side_squared_lower),
        "strict_square_side_squared_lower_asymmetric":str(target_cert.certified_square_side_squared_lower),
        "exact_face_vertices_checked_each":base_cert.checked_face_vertices,
        "local_degree":base_cert.local_degree,
        "asymmetric_target_center_points_are_not_a_square":True,
        "independent_actual_asymmetric_square_parameters":[str(x) for x in parameters],
        "independent_actual_asymmetric_square_points":[[str(p.x),str(p.y)] for p in corners],
        "independent_v223_all_square_orbit_count":census.count,
        "scope":"Local 4-arc 4D degree, not a global square-orbit invariant or general singular-chart continuation. Four perturbed arcs must lie on a Jordan boundary for a boundary-inscribed square claim. No universal positive Big Peg constant or Toeplitz proof.",
    }


if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
