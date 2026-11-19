"""v236: rational checks for a geometric NO-SQUARE theorem on open logarithmic spirals.

The theorem is analytical; executing rational inequalities is a regression audit,
not a substitute for the geometric width / graph arguments in the research note.

Gamma_a(t) = exp((-a + i)*t), t >= 0, plus limiting origin.
When a >= 4, Gamma_a contains no nondegenerate inscribed Euclidean square.
The theorem does NOT assert anything about a closed Jordan curve obtained by
joining the spiral endpoints, and does NOT classify 0<a<4.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from math import cos, exp, hypot, log, sin

EPSILON=F(1,4)
A_THRESHOLD=F(4)


def exp_positive_partial_sum(x, last_degree):
    """Exact rational strict lower bound for e^x when x>0."""
    x=F(x)
    if x<=0 or last_degree<0:raise ValueError('positive x and nonnegative degree required')
    total=F(1)
    term=F(1)
    for n in range(1,last_degree+1):
        term*=x/n
        total+=term
    return total


@dataclass(frozen=True)
class ExclusionCertificate:
    openness:F
    epsilon:F
    strict_delta_upper:F
    triangle_min_width_lower:F
    tangent_angle_upper:F
    two_inner_points_impossible:bool
    three_outer_points_impossible:bool
    four_outer_points_impossible:bool
    status:str

    @property
    def excluded(self):
        return self.status=='square-free-open-spiral'


def exclusion_certificate(a):
    """Certify a>=4 using exact rational inequalities only.

    Geometry assumptions (proved in the note):
    * at most one square vertex in radius <=epsilon R, because pairwise
      square vertex separation >=s >=(1-epsilon)*R/sqrt(2);
    * any 3 square corners form a right-isosceles triangle of minimum
      strip width s/sqrt(2);
    * the remaining spiral annulus is a strict 1-Lipschitz graph in an axis
      through the outermost vertex, so cannot contain all four corners.
    """
    a=F(a)
    if a<=0:raise ValueError('positive openness required')
    if a<A_THRESHOLD:
        return ExclusionCertificate(a,EPSILON,F(3,2)/a,(1-EPSILON)/2,
                                    F(3,2)/a+1/a,False,False,False,
                                    'not-established')
    # Exponential partial sum: exp(3/2) > 67/16 > 4, so log 4 < 3/2.
    assert exp_positive_partial_sum(F(3,2),3)==F(67,16)>4
    # Delta = (log 4)/a < (3/2)/a.  Strictness matters at a=4.
    delta_upper=F(3,2)/a
    # A single inner vertex and an outermost vertex imply side
    # s >= (1-eps)R/sqrt(2), so the 3 other vertices have strip
    # width >= s/sqrt(2) >= (1-eps)R/2.
    width_lower=(1-EPSILON)/2
    assert delta_upper<=width_lower # strict real Delta < upper
    # If two vertices are in the inner disk, their separation <=2 eps R,
    # whereas square-vertex separation >=s.
    two_small=2*(2*EPSILON)**2<(1-EPSILON)**2
    # Angular deviation of the remaining outer spiral <=delta, and its
    # derivative is rotated by -atan(1/a), whose absolute angle <1/a.
    # Delta + atan(1/a) < delta_upper + 1/a <=5/8 < pi/4.
    tangent_upper=delta_upper+1/a
    assert tangent_upper<=F(5,8)<F(3,4) # pi>3
    assert two_small
    return ExclusionCertificate(a,EPSILON,delta_upper,width_lower,
                                tangent_upper,True,True,True,
                                'square-free-open-spiral')


def gamma(a,t):
    """Numerical plotting/sanity ONLY. Never used for certification."""
    a=float(a);t=float(t)
    if a<=0 or t<0:raise ValueError('a>0, t>=0')
    r=exp(-a*t)
    return (r*cos(t),r*sin(t))


def v99_special_pitch_locator():
    """Historical numerical locator, NOT a proof of its pitch value.

    The exact v99 Miranda face checker supplies a special a in (1/20,1/4).
    Its actual isolated parameter is located numerically near 1/9.9726096.
    """
    return {'k_approx':9.972609617563232,
            'a_approx':1/9.972609617563232,
            'certified_outer_bracket':'1/20 < a < 1/4',
            'existence_source':'restored v99 exact Poincare-Miranda theorem',
            'number_of_squares':'infinite self-similar descendants if one exists',
            'scope':'not a positive-measure interval of openness values'}


def report():
    certs=[exclusion_certificate(a) for a in (F(4),F(5),F(10),F(101,4))]
    return {'theorem':'For all a>=4, the compactified open logarithmic spiral has no nondegenerate inscribed square',
            'curve':'Gamma_a(t)=exp((-a+i)t), t>=0, with endpoint at 0',
            'rational_checks':{'exp_3over2_lower':str(exp_positive_partial_sum(F(3,2),3)),
                               'eps':str(EPSILON),
                               'delta_upper_at_a4':str(certs[0].strict_delta_upper),
                               'three_corner_min_strip_width_lower_at_a4':str(certs[0].triangle_min_width_lower),
                               'tangent_angle_upper_at_a4':str(certs[0].tangent_angle_upper)},
            'certified_a':list(map(lambda x:str(x.openness),certs)),
            'v99_special_pitch':v99_special_pitch_locator(),
            'limits':'Does not rule out squares on Jordan closures, or determine all 0<a<4 pitches.'}

if __name__=='__main__':
    import json
    print(json.dumps(report(),indent=2,sort_keys=True))
