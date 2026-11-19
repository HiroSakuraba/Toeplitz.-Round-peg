# v240: Square-free compactified open logarithmic spirals for a >= 17/20

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** A new elementary geometric strengthening of v239, with fully rational numerical inequality audits. The geometric argument, not test execution, proves the theorem. External mathematical review remains appropriate.

## Theorem and scope

For real a>0, define the compactified open logarithmic spiral

$$
S_a=\{e^{(-a+i)t}:t\geq0\}\cup\{0\}.
$$

**Theorem.** For every real a>=17/20, no four distinct points of S_a are the corners of a nondegenerate Euclidean square.

Unlike a Jordan curve, S_a is an arc with the origin as its accumulation endpoint. That endpoint is included in the argument. No claim is made for 0<a<17/20; in particular the older v99 Poincare-Miranda theorem supplies a certified square-bearing log-spiral pitch near a=0.1002747. This is *not* a solution or counterexample to Toeplitz's square-peg problem.

## Exact rational prerequisites

Let a0=17/20, q=4/15, and x=108/125. The certificate program verifies all strict comparisons below using only Python Fraction values.

**(1) Stronger two-inner observation.** Since 3q=4/5<1, we have 2q<1-q. In a square, among any two vertices other than a fixed vertex, at least one is *adjacent* to that fixed vertex (there is only one opposite corner). This improves v239's general diagonal lower bound when at least two corners are inner.

**(2) Rational lower bound pi>25/8.** Inscribe a regular 24-gon in a unit circle. Comparing its perimeter with circumference gives

$$
\pi>24\sin(\pi/24)=12\sqrt{2-\sqrt{2+\sqrt3}}.
$$

The rational inequalities

$$
(1733/1000)^2>3,\qquad
(17807/9216)^2>2+1733/1000
$$

imply that sqrt(2+sqrt 3)<17807/9216. Since 2-17807/9216=625/9216, it follows that

$$
\pi>12\sqrt{625/9216}=25/8.
$$

This bound is strict and uses only a standard perimeter-versus-circle inequality and nested-radical algebra.

**(3) Logarithmic tangent-span bound.** Because the exact positive Taylor sum through degree five gives

$$
e^{53/40}>\sum_{k=0}^{5}\frac{(53/40)^k}{k!}>\frac{15}{4},
$$

we know log(15/4)<53/40. Therefore, for any a>=a0, the angular span of the radial annulus [qR,R] is

$$
\delta=\frac{\log(15/4)}{a}
<\frac{53/40}{17/20}=\frac{53}{34}
<\frac{25}{16}<\frac\pi2.
$$

All comparisons between rational fractions are strict.

**(4) Arctangent bound.** For x=108/125 in (0,1), the alternating Taylor inequalities give

$$
\sin x < S=x-x^3/6+x^5/120,\quad
\cos x > C=1-x^2/2+x^4/24-x^6/720>0.
$$

Exact multiplication verifies 17S<20C; hence tan x<20/17 and

$$
\arctan(20/17)>108/125.
$$

**(5) Maximum transverse displacement.** At a0=17/20, the function f(theta)=e^{-a0 theta}sin theta on [0,pi/2] is maximized at theta*=arctan(20/17), with

$$
f(\theta_*)=\frac{e^{-a0\theta_*}}{\sqrt{1+a0^2}}.
$$

The exact squared inequality (131/100)^2 < 1+(17/20)^2 and (4) imply

$$
f(\theta_*) < \frac{100}{131}\,e^{-459/625}.
$$

The positive Taylor sum through degree four proves

$$
e^{459/625}>\sum_{k=0}^{4}\frac{(459/625)^k}{k!}>\frac{3000}{1441}.
$$

Since (100/131)/(3000/1441)=11/30, it follows strictly that

$$
\max_{0\leq\theta\leq\pi/2} f(\theta)<\frac{11}{30}.
$$

For *every* a>=a0 and theta>=0 with theta<=delta<pi/2, pointwise

$$
0\leq e^{-a\theta}\sin\theta
\leq e^{-a0\theta}\sin\theta<\frac{11}{30}.
$$

## Exclusion proof

Assume four points of S_a form a square with side s>0, where a>=17/20. Let R>0 be the *largest* of their distances from the origin, and rotate the plane so that a farthest corner A lies at (R,0). Set q=4/15.

Call a square corner *inner* if r<qR, and *outer* otherwise. The corner A is outer. Every outer corner lies on the spiral's **single radial annulus** qR<=r<=R, with representation

$$
z=R e^{-a\theta}(\cos\theta+i\sin\theta),
\quad 0\leq\theta\leq\delta=\log(15/4)/a<\pi/2.
$$

### Case A: At least two inner vertices

Among two inner vertices at least one must be **adjacent to A**, since A has only one opposite vertex in the square. Their distance from A is exactly s and, by the reverse triangle inequality, strictly exceeds (1-q)R. But their mutual distance is less than 2qR. With 2q<1-q, this pair of inner corners would be separated by less than s, impossible for a square.

### Case B: Exactly one inner vertex

The farthest-to-inner distance exceeds (1-q)R and is at most sqrt(2)s, so

$$
s>\frac{(1-q)R}{\sqrt2}.
$$

The other three vertices form a right-isosceles triangle with legs s. Its width in *any* fixed direction is at least s/sqrt(2), the altitude to its hypotenuse. Its vertical width is therefore strictly greater than

$$
\frac{(1-q)R}{2}=\frac{11R}{30}.
$$

But by the radial representation, the three outer corners have vertical coordinates satisfying

$$
0\leq y/R=e^{-a\theta}\sin\theta<11/30
$$

by prerequisite (5). They lie inside a strip of vertical width less than 11R/30. Contradiction.

### Case C: No inner vertex

All four corners lie in the same radial subarc qR<=r<=R. Write theta(r)=log(R/r)/a. Then

$$
z(r)=r e^{i\theta(r)},\qquad
\frac{dz}{dr}=e^{i\theta(r)}(1-i/a).
$$

The tangent angle theta(r)-arctan(1/a) varies through an interval of width delta<pi/2. Rotate a single fixed coordinate axis (X,Y) to the midpoint of that tangent interval. Every tangent makes an angle strictly less than pi/4 with its positive X-axis. Hence dX/dr>0 and |dY/dr|<dX/dr throughout the interval. Integration over any pair of radii shows the whole subarc obeys the strict *global* 1-Lipschitz property

$$
|\Delta Y|<|\Delta X|.
$$

A square has two perpendicular adjacent edges. In any fixed coordinates, at least one of two perpendicular equal-length edge vectors satisfies |Delta Y|>=|Delta X|. Therefore no complete square can be contained in this subarc.

The three cases exhaust all possibilities, including 0 as a possible inner vertex. QED.

## Research contribution and limits

v236 supplied a>=4, v238 supplied a>=1, and v239 improved the sufficient condition to a>=9/10. This v240 proof strengthens it to **a>=17/20**, chiefly by using a sharper inner-pair argument and narrow exact rational margins for the transverse maximum and pi. The proof does not classify a<17/20 or establish optimality.

The exact Python auditor validates the numerical inequalities in prerequisites (1)–(5). The geometric triangle-width and tangent-graph lemmas are written here explicitly rather than encoded as false claims of automatic formal verification. Numerical tests sample diagrams only and are not proof inputs.

Run:

~~~sh
python code/log_spiral_exclusion_v240.py
python -m unittest discover -s tests -p 'test_log_spiral_exclusion_v240.py' -v
~~~
