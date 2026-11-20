"""v267: exact secant-spin / diagonal-rectangle square-peg audit.

For ANY Jordan embedding gamma, the antipodal secant
v(t)=gamma(t+1/2)-gamma(t) is nonzero and satisfies v(t+1/2)=-v(t).
Its degree is odd; it is a robust twist, NOT an unconditional
square-forcing invariant. The unordered-pair chord-orientation line
bundle has w1=1 over the core circle of the open Moebius band.

This program EXACTLY checks the PL rational cases, the full D4
transformation of the diagonal-square residual, and both a protected
rectangle-path sign crossing and a collapsing-square counterstress
test. Not a universal Toeplitz proof.
"""
from fractions import Fraction as Q
from itertools import combinations
import json

from polygon_orbit_census_v223 import (
    pt,add,sub,mul,dot,cross,validate_jordan_polygon,
    orientation_area,
)
from collision_quotient_v219 import (
    notch_polygon,notch_square,notch_square_side_squared,
)
from configuration_symmetry_v218 import determinant


class SecantAuditInconclusive(ValueError):
    """A malformed or unsupported source is not a zero winding result."""


def _gamma(P,t):
    """Equally timed edges of rational Jordan polygon, exact for rational t."""
    t=Q(t)%1
    n=len(P);u=t*n
    j=int(u);f=u-j
    return add(P[j],mul(sub(P[(j+1)%n],P[j]),f))


def secant(P,t):
    return sub(_gamma(P,Q(t)+Q(1,2)),_gamma(P,t))


def _breakpoints(P):
    n=len(P)
    return tuple(sorted(
        {Q(0),Q(1)}
        |{Q(i,n) for i in range(n)}
        |{(Q(i,n)-Q(1,2))%1 for i in range(n)}
    ))


def _origin_on_segment(v,w):
    """Zero on closed segment [v,w], using exact rational geometry."""
    if v==pt(0,0) or w==pt(0,0):
        return True
    return cross(v,w)==0 and dot(v,w)<0


def _min_segment_norm_sq(v,w):
    d=sub(w,v);den=dot(d,d)
    if den==0:return dot(v,v)
    s=max(Q(0),min(Q(1),-dot(v,d)/den))
    x=add(v,mul(d,s))
    return dot(x,x)


def exact_secant_winding(vertices):
    """Exact integer winding of t -> gamma(t+1/2)-gamma(t), rational PL.

    The breakpoints include all polygon vertex and shifted-vertex times.
    Hence secant is EXACTLY LINEAR between consecutive breakpoints.
    A half-open positive-x-axis crossing convention yields winding
    without atan2/float approximations or sampling.
    """
    P=validate_jordan_polygon(vertices)
    times=_breakpoints(P)
    vectors=tuple(secant(P,t) for t in times)
    if vectors[0]!=vectors[-1]:
        raise AssertionError("secant winding path not closed")
    crossings=0;min_sq=None;positive_ray_hits=0
    for a,b in zip(vectors,vectors[1:]):
        if _origin_on_segment(a,b):
            raise SecantAuditInconclusive("antipodal chord hit origin")
        floor=_min_segment_norm_sq(a,b)
        min_sq=floor if min_sq is None else min(min_sq,floor)
        if a.y<=0<b.y or b.y<=0<a.y:
            x=(a.x*b.y-b.x*a.y)/(b.y-a.y)
            if x>0:
                crossings+=1 if a.y<b.y else -1
                positive_ray_hits+=1
    if min_sq is None or min_sq<=0:
        raise AssertionError("missing strict chord floor")
    if crossings%2==0:
        raise AssertionError("nonzero antipodal secant must have odd degree")
    if crossings!=(1 if orientation_area(P)>0 else -1):
        raise AssertionError(
            "Jordan polygon antipodal chord degree differs from orientation")
    return {
        "secant_winding":crossings,
        "odd_degree":crossings%2!=0,
        "antiperiodic_reversal_exact":all(
            secant(P,t+Q(1,2))==mul(secant(P,t),Q(-1))
            for t in times
        ),
        "positive_ray_crossings_checked":positive_ray_hits,
        "secant_piecewise_affine_segments":len(times)-1,
        "minimum_antipodal_secant_length_squared":min_sq,
        # If each curve coordinate is uniformly displaced by epsilon,
        # secants move by at most 2*sqrt(2)*epsilon. Thus the norm^2
        # bound implies a positive C0 stability budget for *this* curve.
        "uniform_coordinate_displacement_epsilon_squared_budget":
            min_sq/Q(8),
        "w1_of_unoriented_chord_orientation_line_on_core":1,
    }


def _square_pairing_residual(corners,pairing=((0,2),(1,3))):
    """(mA-mB in R2, |dA|^2-|dB|^2, dA.dot(dB)).

    Four zeros with distinct square corners force same-center,
    equal-length and perpendicular diagonals, hence a square.
    """
    if len(corners)!=4 or len(pairing)!=2:
        raise ValueError("four points and two chord pairs required")
    inds=tuple(j for pair in pairing for j in pair)
    if sorted(inds)!=[0,1,2,3]:
        raise ValueError("the two chords must partition four endpoints")
    (a,b),(c,d)=pairing
    A,B,C,D=(corners[j] for j in (a,b,c,d))
    da,db=sub(B,A),sub(D,C)
    m= sub(add(A,B),add(C,D))
    return (m.x/Q(2),m.y/Q(2),
            dot(da,da)-dot(db,db),dot(da,db))


def all_pairings(corners):
    """Three unordered partitions; do not assume graceful curve order."""
    return {((0,1),(2,3)):_square_pairing_residual(corners,((0,1),(2,3))),
            ((0,2),(1,3)):_square_pairing_residual(corners,((0,2),(1,3))),
            ((0,3),(1,2)):_square_pairing_residual(corners,((0,3),(1,2)))}


def is_positive_square_diagonals(corners):
    if len(corners)!=4 or len(set(corners))!=4:
        return False
    return sum(v==(Q(0),)*4 for v in all_pairings(corners).values())==1


ROT=(1,2,3,0)
REFLECT=(0,3,2,1)
ROT_TARGET=(-1,-1,-1,-1)
REFLECT_TARGET=(1,1,1,-1)


def _source_sign(perm):
    if sorted(perm)!=list(range(4)):
        raise ValueError("invalid corner permutation")
    return (-1)**sum(perm[i]>perm[j] for i in range(4)
                     for j in range(i+1,4))


def _eval_perm(coords,perm):
    points=tuple(pt(coords[2*i],coords[2*i+1]) for i in range(4))
    return _square_pairing_residual(tuple(points[j] for j in perm))


def _quadratic_unisolvent_points():
    """45 exact nodes determining any degree<=2 polynomial on Q^8."""
    basis=[tuple(Q(1) if i==j else Q(0) for i in range(8))
           for j in range(8)]
    zero=(Q(0),)*8
    return (zero,
            *(tuple(v) for v in basis),
            *(tuple(2*x for x in v) for v in basis),
            *(tuple(a+b for a,b in zip(basis[i],basis[j]))
              for i,j in combinations(range(8),2)))


def certify_diagonal_D4_equivariance():
    """Exhaust coefficient basis for every quadratic residual identity.

    F is polynomial degree <=2 in the eight Cartesian corner coordinates.
    Values at 0, e_i, 2e_i and e_i+e_j reconstruct *every* constant,
    linear, square and mixed coefficient. Thus verification on these
    45 unisolvent nodes is a UNIVERSAL polynomial-identity check,
    rather than a finite random numerical test.
    """
    nodes=tuple(_quadratic_unisolvent_points())
    if len(nodes)!=45 or len(set(nodes))!=45:
        raise AssertionError("quadratic identity basis incomplete")
    for perm,target in ((ROT,ROT_TARGET),(REFLECT,REFLECT_TARGET)):
        for c in nodes:
            f=_eval_perm(c,(0,1,2,3))
            permuted=_eval_perm(c,perm)
            expected=tuple(t*x for t,x in zip(target,f))
            if permuted!=expected:
                raise AssertionError("D4 target action identity failed")
    if not (_source_sign(ROT)==-1
            and _source_sign(REFLECT)==-1
            and _product(ROT_TARGET)==1
            and _product(REFLECT_TARGET)==-1):
        raise AssertionError("wrong source/target orientation characters")
    return {
        "quadratic_coefficients_probed_per_generator":45,
        "total_exact_polynomial_identity_probes":90,
        "cyclic_rotation":dict(
            source_orientation=-1,target_orientation=1,
            local_index_transport=-1,
        ),
        "diagonal_endpoint_reversal":dict(
            source_orientation=-1,target_orientation=-1,
            local_index_transport=1,
        ),
        "simple_absolute_integer_index_from_chord_spin_forces_square":False,
        "comparison":"same alternating cyclic sign obstruction as v218; w1 chord twist alone does not remove v219 absolute Euler no-go",
    }


def _product(v):
    x=1
    for z in v:x*=z
    return x


def diagonal_nowherezero_section_formal(cosine,sine):
    """Values s(c)=(cos(4*pi*c),sin(4*pi*c),0,0) at a unit phase pair.

    On the one graceful cyclic configuration component, v219's
    center c obeys c(Tx)=c(x)+1/4. The diagonal square residual's
    target representation under T is -I4, so s(Tx)=-s(x)
    exactly, and |s|^2=1. Therefore its rank-four associated
    bundle admits an EXPLICIT nowhere-zero equivariant section.
    This is a direct absolute-obstruction no-go for the *new*
    diagonal target representation, not an inference from v219's
    different target representation.
    """
    a,b=Q(cosine),Q(sine)
    if a*a+b*b!=Q(1):
        raise ValueError("exact rational phase must lie on unit circle")
    return (a,b,Q(0),Q(0))


def certify_diagonal_absolute_no_go():
    samples=((Q(1),Q(0)),(Q(0),Q(1)),
             (Q(3,5),Q(4,5)),(Q(-5,13),Q(12,13)))
    for a,b in samples:
        v=diagonal_nowherezero_section_formal(a,b)
        rotated=diagonal_nowherezero_section_formal(-a,-b)
        if rotated!=tuple(z*x for z,x in zip(ROT_TARGET,v)):
            raise AssertionError("diagonal section fails C4 equivariance")
        if sum(z*z for z in v)!=Q(1):
            raise AssertionError("diagonal section vanished")
    return {
        "global_section_formula":"s(c)=(cos(4*pi*c), sin(4*pi*c), 0, 0)",
        "center_quarter_shift":"c(Tx)=c(x)+1/4",
        "target_R":"-I4",
        "nonzero_section_squared_norm":"1",
        "absolute_top_w4_obstruction_on_one_cyclic_quotient":0,
        "no_unrestricted_Jordan_square_conclusion":True,
    }


def tensor_perpendicularity_character(flip_first,flip_second,swap_chords=False):
    """Monodromy of signed diagonal dot: (-1)^(flip_first+flip_second).

    The scalar is a section of L_first tensor L_second; simultaneous
    reversal is TRIVIAL, even if both individual chord lines twist.
    Interchanging the two diagonals also leaves the dot unchanged.
    """
    if type(flip_first) is not bool or type(flip_second) is not bool or type(swap_chords) is not bool:
        raise ValueError("exact Z2 monodromy flags required")
    return -1 if flip_first!=flip_second else 1


def tensor_w1_on_loop(first_w1,second_w1):
    """The *relative* obstruction is w1(L1 tensor L2)=w1(L1)+w1(L2)."""
    if first_w1 not in (0,1) or second_w1 not in (0,1):
        raise ValueError("w1 inputs must be 0 or 1")
    return (first_w1+second_w1)%2


def rectangle_family(x,height=Q(1)):
    """Four distinct corners (±x,±height) on fixed rectangle boundary."""
    x,h=Q(x),Q(height)
    if x<=0 or h<=0 or x>2:
        raise ValueError("need 0<x<=2, positive height")
    return (
        pt(-x,-h),pt(x,-h),pt(x,h),pt(-x,h),
    )


def certify_rectangle_orthogonality_crossing(left,right):
    """A *conditional* IVT square certificate on a known rectangle family.

    The four points stay distinct on a common Jordan rectangle for
    positive rational x <=2 and fixed h=1. Midpoint and diagonal-length
    constraints vanish identically as polynomials in x. A change in the
    signed diagonal dot product forces an intermediate square.
    """
    a,b=Q(left),Q(right)
    if not 0<a<b<=2:
        raise ValueError("require 0<a<b<=2")
    f0=_square_pairing_residual(rectangle_family(a))
    f1=_square_pairing_residual(rectangle_family(b))
    if f0[:3]!=(Q(0),)*3 or f1[:3]!=(Q(0),)*3:
        raise AssertionError("rectangle path failed its three constraints")
    if not f0[3]*f1[3]<0:
        raise SecantAuditInconclusive(
            "orthogonality sign is not forced to change on this family")
    if not a<Q(1)<b:
        raise AssertionError("axis-aligned example must cross x=1")
    witness=rectangle_family(Q(1))
    if not is_positive_square_diagonals(witness):
        raise AssertionError("intermediate square identity failed")
    return {
        "rectangle_path_interval":(a,b),
        "orthogonality_at_ends":(f0[3],f1[3]),
        "exact_square_at_parameter":Q(1),
        "strict_chord_half_length_lower_bound_squared":
            min(a*a+Q(1),b*b+Q(1)),
        "positive_size_protected_on_chosen_family":True,
        "this_is_a_conditional_1d_winding_or_sign_crossing":True,
    }


def notch_shrinkout_comparison(n):
    """Same odd universal secant winding, but a *chosen* square shrinks."""
    # v219 and v223 define distinct Point dataclasses; normalize into
    # v223's exact polygon vertex type, without changing coordinates.
    P=tuple(pt(p.x,p.y) for p in notch_polygon(n))
    v=exact_secant_winding(P)
    sq=tuple(pt(p.x,p.y) for p in notch_square(n))
    if not is_positive_square_diagonals(sq):
        raise AssertionError("notch's four corners not a square")
    size=notch_square_side_squared(n)
    if size!=Q(1,int(n)**2):
        raise AssertionError("wrong shrinkout scale")
    # Strong exact non-implication: for these eight equally timed
    # edges every antipodal secant has |v| >=1 (one coordinate has
    # magnitude >=1), and equality occurs at t=1/8. The *selected*
    # tiny square nevertheless has side 1/n.
    if v["minimum_antipodal_secant_length_squared"]!=Q(1):
        raise AssertionError("notch should have UNIFORM unit antipodal chord floor")
    return {
        "n":int(n),"secant_winding":v["secant_winding"],
        "minimum_antipodal_chord_squared":
            v["minimum_antipodal_secant_length_squared"],
        "chosen_square_side_squared":size,
        "limit_rectangle_still_has_macroscopic_squares":True,
        "no_uniform_square_size_floor_deduced_from_winding":True,
    }


def audit():
    ccw=(pt(0,0),pt(2,0),pt(2,2),pt(0,2))
    clockwise=tuple(reversed(ccw))
    diag=certify_diagonal_D4_equivariance()
    rotated=tuple(ccw[j] for j in ROT)
    if not is_positive_square_diagonals(ccw):
        raise AssertionError("square missing")
    if not is_positive_square_diagonals(rotated):
        raise AssertionError("relabelled square missing")
    return {
        "ccw_square_secant":_format(exact_secant_winding(ccw)),
        "clockwise_square_secant":_format(exact_secant_winding(clockwise)),
        "chord_swap_monodromy_w1":1,
        "D4_target_character":diag,
        "diagonal_target_absolute_no_go":certify_diagonal_absolute_no_go(),
        "tensor_diagonal_twist":{
            "one_chord_endpoint_reversal":tensor_perpendicularity_character(True,False),
            "both_chord_endpoint_reversals":tensor_perpendicularity_character(True,True),
            "w1_of_perpendicularity_line_if_two_chord_w1_are_1_0":
                tensor_w1_on_loop(1,0),
            "w1_of_perpendicularity_line_if_two_chord_w1_are_1_1":
                tensor_w1_on_loop(1,1),
        },
        "rectangle_sign_crossing":_format(
            certify_rectangle_orthogonality_crossing(Q(1,2),Q(3,2))),
        "notch_shrinkout_n_64":_format(notch_shrinkout_comparison(64)),
        "unconditional_all_Jordan_square_existence_inferred":False,
        "actionable_conditional_theorem":"a compact 1-cycle of nondegenerate inscribed rectangles with w1(L_first tensor L_second)=1 has a zero of its diagonal-dot-product section, hence a square",
        "remaining_obstruction":"construct such a protected nondegenerate rectangle cycle for every Jordan curve, or a boundary-relative class that cannot disappear into zero chord length",
    }


def _format(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {k:_format(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [_format(v) for v in x]
    return x


if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
