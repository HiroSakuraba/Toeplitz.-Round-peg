# v238: Open logarithmic spirals have no inscribed square for a >= 1

**Research program:** Benjamin John Schulz  
**Status:** New elementary geometric proof strengthening the v236 sufficient condition from a >= 4 to a >= 1. Exact-rational regression checks accompany the proof. This is a result about compactified open arcs, not arbitrary Jordan curves.

## Theorem

For real a > 0, put

$$
S_a=\{e^{(-a+i)t}:t\geq0\}\cup\{0\}.
$$

**Theorem (v238).** For every real a >= 1, no four distinct points of S_a form a nondegenerate Euclidean square.

An inscribed square may use the origin as a vertex; its vertices need not occur consecutively on the arc. This theorem does **not** imply anything about the existence of squares in a closed Jordan curve formed by adding a second arm or a connector. It makes **no claim** for 0 < a < 1. In particular, the existing v99 Poincare-Miranda certificate proves that a special logarithmic spiral with a approximately 0.1002747 *does* contain a square.

## Proof

Suppose a square of side length s > 0 has all four corners on S_a. Let R > 0 be the largest distance of any corner from the origin, and rotate coordinates so that one corner of radius R lies on the positive horizontal axis. On the spiral segment with radii r in [R/4,R], write

$$
r=R e^{-a\vartheta},\qquad
\vartheta=\frac{1}{a}\log\frac{R}{r},\qquad
0\leq\vartheta\leq\delta=\frac{\log4}{a}.
$$

For a >= 1, the strict bound log 4 < 3/2 gives delta < 3/2 < pi. Here log 4 < 3/2 follows from the strictly positive exponential tail and the exact rational lower sum

$$
e^{3/2}>1+\frac32+\frac98+\frac9{16}
=\frac{67}{16}>4.
$$

Call a vertex **inner** if its radius is < R/4, and **outer** otherwise. The farthest vertex is outer.

### Case 1: At least two inner vertices

The distance between the farthest vertex and any inner one is at least 3R/4; every pair of square vertices is separated by at most sqrt(2) s. Hence

$$
s\geq\frac{3R}{4\sqrt2}>\frac R2
$$

(the strict inequality is equivalent to 8 < 9). But two inner points are less than R/2 apart, since each has radius less than R/4. Every pair of distinct square vertices is at least s apart. Contradiction.

### Case 2: Exactly one inner vertex

The other three square vertices form a right-isosceles triangle, with legs s. Its minimum width in any direction is its altitude to the hypotenuse, s/sqrt(2). Therefore its width in the vertical direction, perpendicular to the farthest vertex's radius, is at least

$$
\frac{s}{\sqrt2}\geq\frac{3R}{8}.
$$

But all three outer vertices lie on the angular span 0 <= vartheta <= delta < pi. Their vertical coordinates are nonnegative and satisfy

$$
0\leq\frac yR
=e^{-a\vartheta}\sin\vartheta
\leq\vartheta e^{-a\vartheta}
\leq\frac1{ae}
<\frac{3}{8a}
\leq\frac38.
$$

The third inequality follows by maximizing x exp(-ax) over x >= 0: its global maximum is 1/(ae). The **strict** fourth inequality is certified without floating point because

$$
e>1+1+\frac12+\frac16=\frac83.
$$

Thus all three outer vertices fit in a horizontal strip of width strictly less than 3R/8, contradicting the triangle-width floor. This remains strict when a = 1.

### Case 3: No inner vertices

All four vertices lie on the single radial interval [R/4,R]. As a function of radius,

$$
z(r)=r e^{i\vartheta(r)},\qquad
\frac{dz}{dr}=e^{i\vartheta(r)}(1-i/a).
$$

The tangent direction is therefore

$$
\phi(r)=\vartheta(r)-\arctan(1/a).
$$

These tangent directions fill an interval of angular width delta = log(4)/a. Choose the **fixed** coordinate axis with direction at the midpoint of this interval, namely

$$
\phi_*=\frac{\delta}{2}-\arctan(1/a).
$$

Every tangent lies strictly within pi/4 of this axis, because

$$
|\phi(r)-\phi_*|
\leq\frac{\delta}{2}
=\frac{\log4}{2a}
<\frac{3}{4a}
\leq\frac34
<\frac\pi4.
$$

For the last strict inequality we use pi > 3, which follows geometrically by comparing the circumference of a circle with its inscribed regular hexagon. Thus in the rotated (X,Y) coordinates dX/dr > 0 and |dY/dr| < dX/dr on the **whole** radial interval. Integrating between any two radii shows

$$
|\Delta Y|<|\Delta X|
$$

for any two distinct points on that subarc. But if one edge of a square has |Delta Y| < |Delta X|, its perpendicular adjacent edge has |Delta Y| > |Delta X| in the same coordinates. It is impossible for all four square vertices to lie on the subarc.

The three exhaustive cases prove the theorem. The origin is included in Case 1 or Case 2 whenever it is a vertex. QED.

## What changed from v236

The earlier v236 theorem used the looser transverse estimate R*delta and kept the original radial direction as the graph axis. Those choices demanded a >= 4. v238 uses the exact radial-decay envelope R/(ae) for three outer vertices and bisects the tangent-angle range for four outer vertices. Both improvements have strict elementary margins at a = 1.

The threshold 1 is **sufficient, not asserted optimal**. Neither this proof nor the v99 special-pitch existence proof classifies the intermediate interval (0,1); neither yields a counterexample to the unrestricted Toeplitz conjecture.

## Reproducibility and verification scope

Run:

~~~sh
python code/log_spiral_exclusion_v238.py
python -m unittest discover -s tests -p 'test_log_spiral_exclusion_v238.py' -v
~~~

Every theorem-critical numerical comparison is rational, using Python's Fraction arithmetic. The code checks the exponential Taylor sums and the strict inequalities at the boundary a = 1, and deliberately reports **not established** for 0 < a < 1. The facts about minimum triangle width, monotone projection of a tangent cone, and pi > 3 are proven mathematically above; no claim is made that executing a program alone establishes those geometric facts. The floating-point tests are diagnostic only. The v236 proof and tests remain intact as a historical, weaker result.
