# v258: Quantitative transversality lifts an approximate square to an exact inscribed square

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact-rational, conditional implicit-correction theorem for four distinct **polygonal boundary edges**. An explicit nonconvex centrally symmetric polygon is perturbed into a genuinely asymmetric one; v257 supplies a near-square, and the new v258 verifier supplies a **fully exact inscribed square** and positive lower side bound. A separate v248 Brouwer certificate protects a square under an infinite-dimensional class of rough normal-graph boundary perturbations.

**Not established:** Near-symmetry by itself implies an exact large square on every Jordan curve; a universal Big Peg constant; or unrestricted Toeplitz.

## 1. Why the extra transversality condition is necessary

v256 proves an exact centered square for any centrally symmetric Jordan curve. v257 proves that **an arbitrary** Jordan curve with a disk centered at \(0\) has two boundary points forming a right angle and, if the curve's Hausdorff half-turn defect is \(\delta\), four boundary points within \(\delta\) of an ideal centered square.

There is a logical gap between *four boundary points approximately forming a square* and *four boundary points exactly forming a square*. A naive implicit-function argument fails if the square constraints are singular or if a required boundary point lies at the endpoint of a supporting edge.

A particularly instructive counterexample to a **proposed proof mechanism**, not to square existence, is v256's notched twelve-gon. Its centered diamond has vertices \((4,0),(0,4),(-4,0),(0,-4)\) on four axis-aligned edges. The four-by-four supporting-line square-incidence matrix is **singular**. It does possess an exact square, but this particular four-line chart does not select an isolated solution. No inverse-based transversality estimate can be applied to it.

v258 imposes precisely the missing conditions: an invertible square-incidence matrix, four distinct interior edge contacts, and a quantitative correction smaller than each contact's endpoint margin.

## 2. A general exact rational near-square correction theorem

Let \(P\) be a Jordan polygon (not necessarily convex). Let

\[
z_0=(p_x,p_y,u,v)\in\mathbf R^4
\]

describe an **ideal positive-side Euclidean square** with four counterclockwise corners

\[
C_0(z)=p,\quad
C_1(z)=p+(u,v),\quad
C_2(z)=p+(u-v,u+v),\quad
C_3(z)=p+(-v,u).
\]

Choose four distinct polygon edges \(e_k=\{A_k+tE_k:0\le t\le1\}\), \(k=0,1,2,3\). Assume actual boundary points

\[
q_k=A_k+t_kE_k,\qquad 0<t_k<1,
\qquad |q_k-C_k(z_0)|\le\delta.
\]

Define the exact rational \(4\times4\) matrix \(M\) and vector \(b\) from the supporting-line incidence equations

\[
F_k(z)=JE_k\cdot(C_k(z)-A_k)=0,\qquad
Mz=b.
\]

Suppose \(\det M\ne0\).

Write \(m_{kj}=\partial C_k/\partial z_j\), understood as a two-component constant vector from the displayed linear square parametrization. Set

\[
L=\max_{j=0,\dots,3}\sum_{k=0}^3
  |(M^{-1})_{jk}|\,\|E_k\|_1,\qquad
R=\delta L,
\]

and define the edge-parameter sensitivity

\[
T_k=\frac{\displaystyle\sum_{j=0}^3|E_k\cdot m_{kj}|}
              {\|E_k\|_2^2}.
\]

### Exact lifting criterion

Assume the four **strict rational inequalities**

\[
\boxed{
\delta\,\frac{\|E_k\|_1}{\|E_k\|_2^2}
+R\,T_k < \min(t_k,1-t_k),
\qquad k=0,1,2,3.}
\]

Then **there exists an exact, positive, boundary-inscribed square**, provided the explicit side lower bound below is positive.

**Proof.** Since \(q_k\) lies on its edge,

\[
F_k(z_0)=JE_k\cdot(C_k(z_0)-q_k).
\]

The right-hand side has magnitude at most \(\|E_k\|_1\delta\). Since \(M\) is invertible, the unique exact supporting-line square is

\[
z_*=M^{-1}b,\qquad
z_*-z_0=-M^{-1}F(z_0),
\]

and consequently

\[
\boxed{\|z_*-z_0\|_\infty\le R.}
\]

The edge parameter of its corrected corner \(C_k(z_*)\) differs from \(t_k\) by at most

\[
\left|\frac{E_k\cdot(C_k(z_*)-A_k)}{\|E_k\|_2^2}-t_k\right|
\le
\delta\frac{\|E_k\|_1}{\|E_k\|_2^2}+RT_k.
\]

By the four strict margin inequalities this corrected parameter lies in \((0,1)\). Moreover \(F_k(z_*)=0\), so each corrected square corner lies **on the corresponding actual polygon segment** rather than merely its supporting line.

Finally, the squared side of the corrected square satisfies the completely explicit lower bound

\[
\boxed{\operatorname{side}(z_*)^2
\ge
\max(|u|-R,0)^2+\max(|v|-R,0)^2.}
\]

If this expression is strictly positive, the corrected square is nondegenerate. QED.

**Nature of the argument:** This is a finite-dimensional linear inverse/implicit-correction theorem, not a new general topological degree theorem. It is substantially more targeted than the v248 Brouwer result: it turns the **actual four near-square points** provided by v257 into a checkable exact polygonal square, without needing to parametrize all nearby polygon edges as normal graphs of a common reference curve.

The code independently checks the determinant by a 24-term permutation formula, recomputes \(M^{-1}\) by the v248 exact elimination routine, verifies all four corrected edge memberships and confirms the stated side bound using only \`fractions.Fraction\`.

## 3. A nonconvex, centrally symmetric reference with transverse corners

Define the counterclockwise ten-gon \(P_0\) with vertices

\[
\begin{aligned}
&(4,-1),(4,1),(3/2,3/2),(2,5),(-2,3),\\
&(-4,1),(-4,-1),(-3/2,-3/2),(-2,-5),(2,-3).
\end{aligned}
\]

It has exact antipodal vertex pairing \(P_{i+5}=-P_i\). At the vertices \((3/2,3/2)\) and its opposite, the signed successive turn is negative, so **the polygon is nonconvex**. It nevertheless contains the centered disk of radius 1, as certified by its exact edge-to-origin distances and Jordan winding.

The centered diamond with vertices

\[
(4,0),(0,4),(-4,0),(0,-4)
\]

meets the four edges in cyclic order \((0,3,5,8)\), with each contact **strictly interior** to its edge. These edge normals are not in the singular alignment of the axis-aligned notched twelve-gon; their four-by-four incidence determinant is exactly \(+64\).

The pre-existing **v248 Brouwer theorem** now supplies an additional, genuinely infinite-dimensional statement: replace each of these four chosen segments by any **continuous** normal graph of amplitude at most \(\varepsilon=1/100\). Provided all four perturbed graph arcs lie on one Jordan boundary, at least one exact nondegenerate square persists, with a separately computed explicit positive squared-side bound. Differentiability is not required.

This is a conditional local robustness result, not a theorem for arbitrary perturbations of the entire curve.

## 4. Break symmetry and solve the actual v257 near-square

Produce an asymmetric polygon \(P_1\) by applying these **rational vertex displacements** to \(P_0\):

| Vertex indices | Displacement |
|---|---|
| 0 and 1 | \((1/160,0)\) |
| 2 | \((1/200,-1/240)\) |
| 3 and 4 | \((1/140,-1/70)\) |
| 5 and 6 | \((-1/180,0)\) |
| 7 | \((0,0)\) |
| 8 and 9 | \((1/190,-1/95)\) |

The exact code checks that \(P_1\) is still a simple nonconvex Jordan polygon, contains the closed centered unit disk, and **is no longer centrally symmetric**.

Run the actual **v257 quarter-turn intersection and antipodal indexed-edge construction** on this target. It returns an ideal centered square and four actual boundary points, with the four contacts lying strictly inside the selected edges

\[
\boxed{(3,5,8,0)}.
\]

Their indexed antipodal mismatch and maximum actual corner error are proved **strictly below** \(\delta^2=(1/30)^2\).

The v258 near-square-lifting checker then verifies:

| Exact quantity | Certificate |
|---|---|
| Four-edge incidence determinant | \(-64\) |
| Rational upper corner displacement | \(\delta=1/30\) |
| Inverse response bound \(L\) | \(6\) |
| Four-parameter correction bound \(R=\delta L\) | \(1/5\) |
| Edge endpoint margins | Each strictly exceeds its independently calculated upper parameter correction |
| Corrected exact squared side | \(94179510917/2934392832\) |
| Guaranteed squared side lower bound | Strictly greater than **28** |

The corrected four boundary vertices are rational and lie on **four actual segments of the asymmetric polygon**, not on an idealized circle, a nearby reference polygon, or extrapolated infinite supporting lines.

An independent v223 all-orbit exact polygon census checks that the resulting geometric square occurs in the full list of inscribed square configurations. The 22 new tests additionally cover determinant degeneracy, wrong edge assignment, insufficient or overly loose perturbation bounds, invalid disks, nonsquare ideals, false geometry and floating-point inputs.

**Interpretation:** This is a genuine example in which the approximate-square information of v257, supplemented by four-edge transversality and conditioning, yields an exact large square on a **nonconvex asymmetric Jordan polygon**.

## 5. Why this still does not solve Toeplitz

The essential additional hypotheses are not automatic:

1. The four near-square boundary points must each lie in a **strictly interior** contact position on four **distinct** edges.
2. The four supporting-line square-incidence constraints must be **invertible**; symmetry may create a singular chart.
3. The inverse's condition number and the actual boundary proximity must together permit a correction smaller than all four endpoint margins.

An arbitrary Jordan curve has no canonical finite edge chart, and a polygonal approximation may have increasingly ill-conditioned incidence matrices and vanishing endpoint margins as its complexity grows. Furthermore, small central-symmetry defect by itself only controls geometric proximity, **not** these three transversality margins.

**Therefore v258 proves a conditional exact-square lifting principle and an explicit nonsymmetric nonconvex instance. It does not produce a complexity-independent positive Big Peg constant, nor a proof of Toeplitz for arbitrary Jordan curves.**

A worthwhile v259 would develop a **chart-independent topological degree** or a finite-atlas index robust under singular four-edge contacts, so the persistence of at least one large square could be proved without selecting a well-conditioned chart in advance.

## Reproduction

~~~sh
python code/near_square_transversality_v258.py
python -m unittest -v tests.test_near_square_transversality_v258
~~~

The dedicated workflow also replays v256-v257 and v248 predecessor examples. All formal proof decisions use Python \`Fraction\`, with no floating-point or sampled-time assumptions.
