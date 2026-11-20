# v261: A nonzero oriented multi-region atlas, and an obstruction on unlabeled square space

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Result:** The v254 moving convex pentagon admits **two disjoint, uniformly isolated rational four-dimensional regions** containing all three *tracked* square branches. The glued two-square collision region has degree **0** as rigorously established in v260. The third, persistent square has a separately certified local Brouwer degree **+1** in exactly the same oriented boundary-parameter convention. Consequently the degree on the **disconnected union of these two particular regions** is **+1 throughout the exact interval that contains the irrational collision**, including the collision instant. A second result demonstrates why this is **not** yet a canonical global integer index on the **unlabeled** square-configuration quotient: cyclic relabeling reverses the local signed index.

## 1. The correct claim and the part still missing

v260 proved one true positive-size square pair death:

- The two disappearing square configurations occupy a continuously glued 4D local source box; their indices are **-1** and **+1** before the collision.
- They meet on an **internal** interface at the unique irrational algebraic collision time and disappear afterward. The whole pair box remains uniformly boundary-zero-free, with degree **0**.
- A distinct third square on edge assignment **C=(0,1,2,4)** survives, but v260 did not calculate its degree or place it in the same global orientation.

The present milestone makes this last calculation *without pretending that the selected neighborhoods are the complete configuration space*. The exact, reproducible result is

\[
\deg(F_t,U_A\cup U_C,0)
=\underbrace{0}_{\text{glued collision region }U_A}
 +\underbrace{1}_{\text{persistent region }U_C}=1
\]

for **every real** \(t\in[a,b]\), with

\[
a=64722/100000,\qquad b=64723/100000.
\]

Here \(U_A\) means v260's **one whole glued pair region**, not an isolated branch, and the union is disconnected.

**This is an oriented relative degree in a fixed labeled and unwrapped sector.** It is not a universal Toeplitz theorem, not an all-time exhaustive census of every square on the pentagon, and not a well-defined *integer* degree on the unlabeled quotient.

## 2. A single square-defect convention

Let the five counterclockwise moving pentagon vertices \(P_0(t),\ldots,P_4(t)\) be exactly those of [v254](v254-genuine-square-pair-death.md). The family is strictly convex and Jordan throughout the larger original interval \([3/5,13/20]\), by v254's exact certificate.

Let each of four selected boundary points be \(p_0,p_1,p_2,p_3\), in counterclockwise order, and define

\[
F(p_0,p_1,p_2,p_3)=
\begin{pmatrix}
p_2-p_1-J(p_1-p_0)\\
p_3-p_0-J(p_1-p_0)
\end{pmatrix}\in\mathbb R^4,
\qquad J(x,y)=(-y,x).
\]

Its zeros are counterclockwise squares, with positive side whenever \(p_1\ne p_0\). Use four **increasing** boundary parameters, with a global unwrapped coordinate \(q_i=e_i+u_i\) on edge \(e_i\), so that

\[
0<q_0<q_1<q_2<q_3<5.
\]

The persistent chart uses

\[
q_C=(u_0,1+u_1,2+u_2,4+u_3).
\]

The colliding pair uses exactly the v260 hinge coordinates

\[
q_A=(s_0,2+h,3+s_2,4+s_3),
\quad h\in[-1/1000,1/1000].
\]

For \(h<0\), the second corner lies on edge 1; for \(h>0\), it lies on edge 2; at \(h=0\), both expressions are the same genuine polygon vertex \(P_2(t)\). On each side, **the derivative from local parameters to these globally increasing parameters is the identity**. Therefore the signed determinant orientations in v260 and this paper can be compared without silently flipping a chart.

## 3. A fixed rational box for the surviving third square

The persistent-square parameters \(u=(u_0,u_1,u_2,u_3)\), along edges \((0,1,2,4)\), lie in

\[
\boxed{
U_C =
(14/100,17/100)\times
(30/100,34/100)\times
(78/100,84/100)\times
(27/100,31/100).
}
\]

The new [v261 verifier](../code/square_oriented_atlas_v261.py) constructs this exact affine four-arc square-defect map and calls the inherited v224 rational-polynomial Cramer chart engine. Write its affine equation as

\[
F_t(u)=M_C(t)u+c_C(t),\qquad
D_C(t)=\det M_C(t),\quad
r_{C,j}(t)=N_{C,j}(t)/D_C(t).
\]

All real-time inequalities below are proved by **strict rational Bernstein coefficient positivity over the entire normalized time interval** \([0,1]\), not by floating-point root samples:

1. \(D_C(t)>1/5\) throughout \([a,b]\), so the unique candidate root is regular and has orientation **+1**.
2. The \(4\times4\) matrix entries are bounded in absolute value by **2** at the two rational endpoints. They are affine in time, so this bounds the entries at *every real time* in the interval. The determinant polynomial is checked against an independent permutation determinant at both endpoints.
3. For each of the four parameters, two strict Bernstein-certified polynomial inequalities prove the **extended Cramer root** lies between the box faces, with a **minimum distance of 1/1000** from *each* face, for every real time. That is eight root-coordinate sign inequalities in addition to the determinant lower-bound inequality.
4. Direct Fraction substitution at both rational endpoints independently checks the exact square equations and the globally ordered boundary parameters.

Because the only root of this affine map lies uniformly inside the box, its boundary is zero-free for *all* real times. Thus

\[
\boxed{\deg(F_t,U_C,0)=+1\quad(t\in[a,b]).}
\]

The approximate endpoint parameters are given only for readability:

| Parameter | \(t=a\) | \(t=b\) |
|---|---:|---:|
| \(u_0\) | 0.152553796 | 0.152554470 |
| \(u_1\) | 0.316701884 | 0.316704018 |
| \(u_2\) | 0.811744300 | 0.811745849 |
| \(u_3\) | 0.289138950 | 0.289139832 |

All assertions and decisions in the program use exact Fraction expressions, not these decimals.

## 4. External boundary, noncollapse and rough perturbation budgets

The exact determinant and entry bounds permit the same explicit inverse-matrix estimate as v260. Each 3-by-3 cofactor has magnitude at most \(6\cdot 2^3\); summing four cofactors per row, with \(D_C>1/5\), yields

\[
\|M_C^{-1}\|_\infty<24\cdot2^3/(1/5)=960.
\]

The root is more than \(1/1000\) away from every exterior box face, so

\[
\boxed{\|F_t\|_{\infty,\partial U_C}>1/960000}.
\]

v260's separate pair region has the smaller certified residual bound

\[
\|F_t\|_{\infty,\partial U_A}>1/19200000.
\]

Consequently **the union of both outer boundaries** has the uniform conservative lower bound

\[
\boxed{\|F_t\|_{\infty,\partial(U_A\cup U_C)}>1/19200000}.
\]

For an arbitrary **continuous**, coordinatewise displacement of each used boundary arc by at most \(\varepsilon=10^{-9}\), the square residual changes by at most \(4\varepsilon\), which is strictly smaller than this certified boundary margin. Therefore the degree remains **+1**, forcing *at least one exact square* somewhere in the selected union for every such perturbation.

There is also a certified **size floor** independent of which point in the persistent box is chosen. The vector from first to second corner satisfies

\[
(p_1-p_0)_x<-1,\qquad(p_1-p_0)_y>9/100.
\]

These inequalities follow by checking the eight rational time/parameter vertices because both components are multi-affine in \((t,u_0,u_1)\). After the allowed coordinatewise perturbation, any root in the persistent box has

\[
\operatorname{side}^2>
(1-2\varepsilon)^2+(9/100-2\varepsilon)^2>1.
\]

The colliding box has its own certified positive side lower bound in v260. Thus **every zero of the selected perturbed atlas corresponds to a nondegenerate square**, and the more direct persistent-region degree already guarantees such a square in that box.

**Caveat:** This is a small fixed-parameter and small \(C^0\) neighborhood theorem, not a guarantee for an arbitrary homotopy of all Jordan curves. Without a boundary-zero-free homotopy, the relative degree may change as roots move out of the selected regions.

## 5. The complete local ledger around the event

The two boxes are **disjoint in the same global configuration coordinates**. In particular their first corner parameters are separated:

\[
q_{C,0}<17/100<9/10<q_{A,0}.
\]

They are separated in their third corner parameters as well:

\[
q_{C,2}<2+84/100<3+43/100<q_{A,2}.
\]

So additivity across the disconnected union is unambiguous.

| Time | Persistent region \(U_C\) | Pair region \(U_A\) | Total |
|---|---:|---:|---:|
| Before \(\alpha\) | +1 | -1 + 1 = 0 | **+1** |
| At \(\alpha\) | +1 | 0 (root at inner hinge) | **+1** |
| After \(\alpha\) | +1 | 0 (no pair roots) | **+1** |

This is a mathematical certificate for all real times in the interval. It is consistent with the independent v223 complete whole-pentagon geometric census **only at the two rational times**: three geometric squares just before, one just after. It does **not** turn two complete endpoint censuses into an all-real-time exhaustive classification.

## 6. Fundamental obstruction: an unlabeled square has no canonical first corner

An inscribed geometric square admits four cyclic choices of which corner is called \(p_0\). The rotation

\[
(p_0,p_1,p_2,p_3)\mapsto(p_1,p_2,p_3,p_0)
\]

acts on the ordered *unwrapped* source parameters by

\[
(q_0,q_1,q_2,q_3)\mapsto(q_1,q_2,q_3,q_0+5).
\]

The corresponding \(4\times4\) permutation has determinant **-1**. But the square-residual vector transforms by the constant exact matrix

\[
R=
\begin{pmatrix}
-1&1&1&0\\
-1&-1&0&1\\
0&1&0&0\\
-1&0&0&0
\end{pmatrix},
\qquad \det R=+1.
\]

In fact \(F(p_1,p_2,p_3,p_0)=R F(p_0,p_1,p_2,p_3)\) **for arbitrary four points**, not just squares. The verifier proves this symbolic linear identity by testing all eight standard coordinate basis vectors, which suffices by linearity.

For an isolated regular square, the local degree consequently changes by

\[
\operatorname{sgn}(\det R)\,
\operatorname{sgn}(\det P)^{-1}
=(+1)(-1)=\boxed{-1}
\]

under a cyclic change of the first corner.

**Therefore a naively summed integer orientation index on the quotient by cyclic relabeling is not well-defined.** Our +1 is valid for the *specified labeled lift*. It cannot be advertised as a canonical integer invariant of unlabeled geometric squares. A mod-2 index does survive this orientation reversal, and a twisted-orientation local system may permit a more refined quotient-level construction. Neither candidate alone supplies uniform protection against small-square collapse along arbitrary Jordan limits.

This is not a failure of the v260/v261 certificates. It is a **structural restriction on how their signed local indices may be globally assembled**.

## 7. Exact checks and replay

The program and regression tests depend only on the Python standard library and the pre-existing exact Fraction/Bernstein modules.

```sh
python code/square_oriented_atlas_v261.py
python -m unittest -v tests.test_square_oriented_atlas_v261
```

The [new tests](../tests/test_square_oriented_atlas_v261.py) verify the time-uniform strict root enclosure, exact endpoint squares, the common global ordered-coordinate convention, the genuine vertex-glued pair, disjointness of source regions, strict determinant and residual bounds, noncollapse, signed index addition, the cyclic source/target orientation discrepancy, and fail-closed controls for altered geometry, boxes and perturbation budgets.

### Next stage (v262)

1. **Choose the correct quotient invariant.** Build the fourfold cyclic action on a complete positive-side-size configuration space. Use mod-2 degree or the appropriate orientation local system. Prove transition compatibility across cyclic relabeling, repeated-edge strata and genuine polygon-vertex events.
2. **Close the exhaustive pentagon event ledger.** Extend v223/v224 beyond rational snapshots. Certify all relevant singular/rank-deficient charts and repeated-edge assignments over the real parameter interval. Only then decide whether the selected three-region atlas represents every square at every real time.
3. **Control the collapse boundary.** Derive a quantitative positive side-length floor or an alternative compactness invariant that survives families with arbitrarily many vertices. This is the key missing step to apply an existence invariant to a limiting arbitrary continuous Jordan curve.
4. **Cross-check against known parity topology.** Any new mod-2 theorem on polygons must be reconciled with classical odd-orbit parity; merely recovering that known generic parity is not by itself a solution of unrestricted Toeplitz.

**Bottom line:** v261 certifies a real **nonzero relative degree** for the three tracked configurations through a genuine pair annihilation, and identifies a concrete **orientation-twisting obstruction** to promoting that result naively to all geometric square orbits. The unrestricted Toeplitz conjecture remains open.

### References

- [v254: genuine positive-square pair death](v254-genuine-square-pair-death.md).
- [v259: four-dimensional local square degree](v259-piecewise-arc-local-degree.md).
- [v260: exact glued degree zero for the annihilating pair](v260-relative-degree-pair-death.md).
- B. Matschke, *A Survey on the Square Peg Problem*, *Notices of the AMS* 61 (2014).
- D. Salamon, *Lectures on Differential Topology*, chapter on Brouwer degree, https://people.math.ethz.ch/~salamon/PREPRINTS/difftop.pdf.
