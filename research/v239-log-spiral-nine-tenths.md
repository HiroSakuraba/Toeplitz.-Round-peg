# v239: Square-free compactified open logarithmic spirals for a >= 9/10

**Research program:** Benjamin John Schulz  
**19 November 2026**  
**Status:** Elementary analytic strengthening of v238, with exact-rational audits. Independent mathematical review is welcome. This proves no result for arbitrary closed Jordan curves.

## Statement

Let S_a = {exp((-a+i)t): t >= 0} union {0}, with a>0. This is a compactified **open arc**. Its four square vertices may include the limiting origin, and need not occur consecutively along the arc.

**Theorem:** For every real a >= 9/10, there do not exist four distinct vertices of S_a forming a nondegenerate Euclidean square.

The threshold is only sufficient. No result is asserted for 0<a<9/10. In fact, the separately certified v99 special-pitch case near a=0.1002747 **does** contain an inscribed square.

## Exact rational lemmas

All inequalities below are strict where necessary.

1. Put q=13/50. Then 2(2q)^2 < (1-q)^2: in exact fractions,
   338/625 < 1369/2500. This ensures that two vertices in the inner disk of radius qR are impossible once the farthest-to-inner chord bound is applied.
2. The positive Taylor sum of exp(27/20) through degree 6 exceeds 50/13. Thus log(50/13)<27/20. For a>=9/10 the angular span delta=log(50/13)/a is strictly less than (27/20)/(9/10)=3/2<pi/2. Here pi>3 follows from the regular hexagon inscribed in a circle.
3. At x=4/5, the standard alternating Taylor estimates give sin x < x-x^3/6+x^5/120 and cos x > 1-x^2/2=17/25 > 0. An exact rational comparison proves sin(x) < (10/9)(17/25), so tan(4/5)<10/9. Since 0<4/5<pi/2 and tan is increasing, arctan(10/9)>4/5.
4. sqrt(1+(9/10)^2)=sqrt(181)/10>4/3, because 181/100>16/9.
5. The positive Taylor sum of exp(18/25) through degree 3 exceeds 75/37. Hence (3/4)exp(-18/25)<37/100.

The Python checker uses exact Fraction arithmetic for all five comparisons. This is not a substitution for the geometric proof.

## Geometric proof

Assume a square with side s>0 has all four vertices on S_a. Choose R>0, the greatest of their distances from the origin, and rotate so that a farthest vertex lies at (R,0). Any other spiral point with radius between qR and R can be written as

$$
z=R e^{-a\theta}(\cos\theta+i\sin\theta),\qquad
0\le\theta\le\delta=\frac{\log(1/q)}a=\frac{\log(50/13)}a<\frac32<\frac\pi2.
$$

Call any square vertex with radius less than qR *inner*, and the others *outer*. The farthest vertex is outer.

### Case A: two or more inner vertices

Distance from the farthest vertex to an inner vertex is greater than (1-q)R. Any pair of vertices in a square has distance at most sqrt(2)s. Therefore

$$
s\ge \frac{(1-q)R}{\sqrt2}.
$$

But any two inner vertices are separated by less than 2qR. The exact inequality 2(2q)^2 < (1-q)^2 implies 2qR < (1-q)R/sqrt(2) <= s. This contradicts their being two distinct corners of the square.

### Case B: exactly one inner vertex

The other three vertices form a right-isosceles triangle with legs s. In any strip direction, its width is at least its altitude to the hypotenuse, s/sqrt(2). So its **vertical** width is at least

$$
\frac s{\sqrt2}\ge\frac{(1-q)R}{2}=\frac{37R}{100}.
$$

All three outer vertices satisfy 0<=theta<pi/2, so their y coordinates are nonnegative:

$$
0\le\frac yR=e^{-a\theta}\sin\theta.
$$

Since a>=9/10 and 0<=theta<=delta<pi/2,

$$
e^{-a\theta}\sin\theta
\le e^{-(9/10)\theta}\sin\theta
\le \frac{e^{-(9/10)\arctan(10/9)}}{\sqrt{1+(9/10)^2}}.
$$

The last inequality is an elementary one-variable maximization: the derivative of e^{-b theta}sin theta vanishes at theta=arctan(1/b), which is its maximum on [0,pi/2].

By Lemmas 3–5, the right side is strictly less than

$$
\frac{e^{-18/25}}{4/3}
=\frac34e^{-18/25}
<\frac{37}{100}.
$$

Thus every outer vertex lies in the horizontal strip 0<=y<37R/100, contradicting the required minimum vertical width of the triangle. This works for *every* a>=9/10, not just rational values.

### Case C: no inner vertices

All four lie on one radial segment qR<=r<=R. With theta(r)=log(R/r)/a, differentiation gives

$$
\frac{dz}{dr}=e^{i\theta(r)}(1-i/a).
$$

The tangent argument ranges from -arctan(1/a) to delta-arctan(1/a). Choose a fixed axis whose direction is the midpoint of this interval. Every tangent then differs from that axis by at most delta/2 < 3/4 < pi/4. The rotated arc is a strict graph with |dY/dr|<dX/dr throughout and hence |Delta Y|<|Delta X| for any two different arc points. A square has perpendicular adjacent edges, and at least one of them violates this inequality. Contradiction.

The three cases are exhaustive, including whenever one of the corners is the limiting point 0. QED.

## Relation to v238 and research limitations

v238 proved a>=1 using q=1/4 and the generic bound y/R<=1/(ae). This proof optimizes q to 13/50 and uses the sharper actual maximum of e^{-a theta}sin theta. It strengthens only the sufficient **local arc exclusion threshold**, not the unrestricted Toeplitz square-peg conjecture.

The exact-rational checker validates the sufficient constants at a=9/10 and fails closed below that threshold. Numerical tests of transverse coordinates and tangent angles are sanity checks, not proof inputs. No novelty claim against the full mathematical literature is made.

Run:

~~~sh
python code/log_spiral_exclusion_v239.py
python -m unittest discover -s tests -p 'test_log_spiral_exclusion_v239.py' -v
~~~
