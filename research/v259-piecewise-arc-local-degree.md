# v259: A four-dimensional local degree across a genuine polygon corner

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** A rigorous, **exact rational**, local Poincaré–Miranda/Brouwer-degree certificate for square existence on four continuous boundary arcs that may have **genuine nonsmooth corners**. It is robust under a quantified class of arbitrary continuous coordinatewise perturbations, and it also lifts a **nonexact near-square** on a slightly asymmetric Jordan polygon to a genuine square without tracking any chosen edge-chart root.

**Scope:** A compact local four-arc parameter cube, not an index for the whole square-configuration space. The method certifies face inequalities directly and does **not** require invertibility of any of the individual piecewise-affine cell Jacobians. It does require that one chosen rational *coordinate preconditioner* is invertible, and that the face inequalities hold. This is **not** an unconditional Toeplitz proof.

## 1. From an edge-chart root to a four-parameter square map

v248–v258 mostly certify an exact square using a selected **four-edge** incidence matrix. That approach can fail at square corners which pass between adjacent polygon edges, and at configurations where the incidence matrix is singular.

Instead, start with four **continuous** arcs

\[
p_k:[-\rho,\rho]\longrightarrow \mathbb R^2,\qquad k=0,1,2,3.
\]

They may be nonsmooth or piecewise affine. We do not assume four supporting lines or solve their Cramer equations.

Let \(J(x,y)=(-y,x)\). Define the continuous *square defect*

\[
F(s_0,s_1,s_2,s_3)=
\begin{pmatrix}
p_2(s_2)-p_1(s_1)-J\big(p_1(s_1)-p_0(s_0)\big)\\
p_3(s_3)-p_0(s_0)-J\big(p_1(s_1)-p_0(s_0)\big)
\end{pmatrix}
\in \mathbb R^4.
\]

If \(F(s)=0\), set \(u=p_1(s_1)-p_0(s_0)\). Then

\[
p_2=p_1+Ju,\qquad p_3=p_0+Ju.
\]

Thus the four points form a counterclockwise Euclidean square, possibly collapsed if \(u=0\). A separate quantitative side bound below excludes collapse.

**Key advantage:** \(F\) is continuous even when one of the boundary parameters traverses a genuine polygon vertex. No piecewise derivative need agree across the vertex. A degree argument for the entire four-cube can see through the chart change.

## 2. The general face-sign/degree theorem

Let \(K=[-\rho,\rho]^4\), and let \(B\in GL_4(\mathbb R)\) be **any fixed invertible coordinate preconditioner**. Write

\[
G(s)=B F(s).
\]

Suppose, for each coordinate \(i=0,1,2,3\), that

\[
\boxed{
G_i(s)>0\ \text{on the entire face }s_i=+\rho,
\qquad
G_i(s)<0\ \text{on the entire face }s_i=-\rho.}
\]

Then **there is some \(s\in\operatorname{int}K\) with \(F(s)=0\)**. Indeed, the strict face conditions are the Poincaré–Miranda theorem. More precisely, for the homotopy

\[
H_\tau(s)=(1-\tau)G(s)+\tau s,\quad 0\le \tau\le1,
\]

if \(s\in\partial K\), at least one coordinate \(s_i=\pm\rho\) and the same coordinate of \(H_\tau(s)\) has strict sign \(\pm\). Thus the entire homotopy is boundary-zero-free and

\[
\boxed{\deg(G,\operatorname{int}K,0)
=\deg(\mathrm{id},\operatorname{int}K,0)=+1.}
\]

This is a **local, coordinate-dependent degree certificate**. The integer degree is independent of *how many* isolated square roots there are in the cube, and remains meaningful if individual square-root charts merge or cease to be differentiable. It is not a global degree for all four-tuples on an arbitrary Jordan curve.

For the general Poincaré–Miranda theorem and its Brouwer-degree interpretation, see *Variations on the Brouwer Fixed Point Theorem: A Survey*, *Mathematics* 8(4), 501 (2020), https://www.mdpi.com/2227-7390/8/4/501 .

### Arbitrary continuous arc perturbations

Suppose each reference arc \(p_k^0\) is replaced by an **arbitrary continuous** arc \(p_k\) satisfying the **coordinatewise** uniform bound

\[
|p_{k,x}(t)-p_{k,x}^0(t)|\le\varepsilon,\qquad
|p_{k,y}(t)-p_{k,y}^0(t)|\le\varepsilon
\]

for all \(k,t\). No differentiability or frequency bound is assumed.

The square defect \(F\) is a linear expression in the eight Cartesian input coordinates. Let \(D:\mathbb R^8\to\mathbb R^4\) denote this linear map, and set

\[
\Lambda=\max_{i=0,\dots,3}
\sum_{j=0}^7 |(BD)_{ij}|.
\]

If the reference face signs have a common **strict margin** \(m>0\), then

\[
\|BF-BF^0\|_\infty\le\Lambda\varepsilon.
\]

Consequently, the exact inequality

\[
\boxed{\Lambda\varepsilon<m}
\]

protects all eight strict face signs for **every such continuous perturbation**, giving **degree +1 and at least one square** on the perturbed four arcs.

If all four perturbed arcs belong to one Jordan boundary, that square is boundary-inscribed there. Without the Jordan-boundary condition, the unconditional statement is only that a square exists on the *union of those four arcs*.

### Positive square-size floor

Let \(u_0=p_1^0(0)-p_0^0(0)\). Suppose each reference \(p_k^0\) is affine on the two halves of its parameter interval, with negative and positive slope vectors \(E_{k,-},E_{k,+}\). Then a solution \(s\in K\) satisfies the coordinatewise estimates

\[
|(p_1(s_1)-p_0(s_0))_\alpha-(u_0)_\alpha|
\le d_\alpha
\]

for \(\alpha=x,y\), where

\[
d_\alpha=\rho\left(
\max(|E_{0,-,\alpha}|,|E_{0,+,\alpha}|)
+\max(|E_{1,-,\alpha}|,|E_{1,+,\alpha}|)
\right)+2\varepsilon.
\]

Therefore any square guaranteed by the degree theorem has

\[
\boxed{\operatorname{side}^2
\ge\max(|u_{0,x}|-d_x,0)^2
+\max(|u_{0,y}|-d_y,0)^2.}
\]

A strictly positive right-hand side proves **noncollapse**, independently of the number, smoothness, or identities of the square roots in \(K\).

## 3. How the computer proves all face inequalities, not just samples

For the actual polygon data, the reference arcs are piecewise affine. The square-defect map and its linear transform \(G\) are then affine on every **subcube** cut by the hyperplanes \(s_k=0\), because each \(p_k^0(s_k)\) depends only on its *own* coordinate.

On any of the eight faces of \(K\), the extrema of each component of \(G\) on a subcube occur at its vertices. It is consequently **rigorous** to test the exact values at the Cartesian product of parameter values \(-\rho,0,\rho\) for each free coordinate.

The v259 verifier checks **216 exact rational face vertices**: on each of the eight faces there are \(3^3=27\) enumerated vertices. Redundancy is deliberate: some input arcs are straight, but the same code handles kinks in *any of the four arcs*.

A strict common positive rational margin proves the signs at **all real points of all faces**, with no numerical sampling and no differential or Cramer-root assumptions. The coefficient matrices and \(BD\) perturbation norm use exact rational Gauss–Jordan elimination and direct coordinate-basis evaluation.

If a face-sign proof fails, the verifier **fails closed**, not assuming that another chart must save the square.

## 4. Exact genuine-corner reference: v250's moving convex nine-gon

Use v250's exact moving-kink family at time \(\lambda=1/2\). Its first square vertex is the **noncollinear polygon vertex** indexed 1, while its remaining square vertices lie strictly inside polygon edges \((3,5,7)\).

Let \(V=P_1\) be the actual kink. Define the first reference arc for \(-1/100\le s\le1/100\) by

\[
p_0^0(s)=
\begin{cases}
V+s(P_1-P_0),&s\le0,\\
V+s(P_2-P_1),&s\ge0.
\end{cases}
\]

The two directions are **not collinear**, and this arc genuinely crosses a polygon corner. The other three reference arcs are the actual straight subarcs on edges \((3,5,7)\) centered at the remaining three inscribed square corners.

The entire nine-gon is Jordan, strictly convex and contains all these subarcs. This is reverified using v250's exact all-time Jordan certificate. The cube radius \(\rho=1/100\) is small enough that all four selected subarcs remain inside their actual assigned edges.

Choose the preconditioner \(B\) as the exact inverse of the average of the two reference affine defect matrices, one for \(s_0<0\) and one for \(s_0>0\). **No determinant of either side's supporting-edge incidence chart is tested or required by the degree theorem.**

The exact verifier obtains

\[
\boxed{
m=\frac{17387}{2696700}>0,\qquad
\Lambda=\frac{89890}{3311}.
}
\]

For coordinatewise continuous arc perturbation of at most

\[
\boxed{\varepsilon=\frac1{100000}},
\]

the maximum change of the preconditioned square defect is

\[
\boxed{
\Lambda\varepsilon=\frac{8989}{33110000}
<m.}
\]

Thus every such continuous perturbation has a square on the four arcs, and if those arcs are part of a Jordan curve it is an **inscribed** square. The exact side bound computed from the slopes and cube radius is **strictly greater than \(9\)**. This is stronger than a mere claim of nonzero side.

## 5. A genuinely nonexact near-square on an asymmetric Jordan polygon

To demonstrate that the result is **not** merely recovering an exact square already present at the center of the parameter cube, modify the v250 polygon by translating **both endpoints of edge 3** by

\[
\Delta=(1/20000,-1/25000).
\]

The modified polygon is independently checked to remain Jordan, and the selected subarcs still lie on its actual boundary. The first corner stays at the genuine vertex kink; the second center is moved by \(\Delta\). Their four central boundary points therefore have a **nonzero** four-component square defect.

Nevertheless, the exact 4D face inequalities remain valid, with the **new** strict margin

\[
\boxed{m_{\rm shifted}=
\frac{11853176941}{1785754740000}>0.}
\]

The preconditioned response \(\Lambda\) is unchanged, and again

\[
\Lambda/100000<m_{\rm shifted}.
\]

Hence this actual **asymmetric, nonsquare-center** polygon has an exact inscribed square of side squared \(>9\), even after any additional continuous coordinatewise displacement of the four chosen arcs within \(\varepsilon=10^{-5}\) (provided they remain on a Jordan boundary).

As an independent sanity check, the verifier solves the piecewise affine equations on the correct side of the kink and compares the resulting exact rational square with the **v223 complete geometric-square orbit census** on the modified polygon. This independent witness extraction is **not used in the topological existence proof**, which relies only on strict boundary signs.

## 6. What this achieves, and what it does not

**New v259 capabilities:**
- A local degree certificate across a true polygon-vertex kink without selecting a single four-edge Cramer branch.
- Existence survives a **near-square** at the cube center, not only a square at that center.
- Any additional continuous, not necessarily differentiable, perturbations of the four selected arcs are covered by a quantitative \(\varepsilon\) bound.
- Explicit positive square-side floor, guaranteeing that the found root cannot collapse.
- Exact rational face-verification algorithm that never requires invertibility of each individual piecewise-affine Jacobian.

**Essential limitations:**
- The preconditioning matrix must be invertible, and all eight face signs must be verified. The code refuses singular choices and any failed sign or perturbation margin. It does **not** prove that such a box exists around every square on every Jordan curve.
- This is a **local parameter-cube degree**, not a chart-independent global configuration-space degree class for arbitrary curves. Although the local result does not rely on four-edge incidence charts, the parameterization still chooses four boundary arcs.
- The arcs must be contained in the same Jordan boundary for the usual Toeplitz notion of a boundary-inscribed square.
- The side bound and perturbation budget depend on the *specific polygon geometry*. There is no complexity-independent universal Big Peg constant.
- The current example is a **vertex crossing**, not a complete size-filtered parity theory handling arbitrary unfiltered orbit births/deaths such as v254 at once.

### v260 recommendation

Define a genuine **finite union of local degree regions** and a rigorous gluing law under parameter changes. Then test it on v254's positive-square pair-death example. At an actual two-to-zero bifurcation the local signed degree must remain consistent, even when the individual root count changes by two. Prove the boundary nonvanishing hypotheses or report that the proposed invariant is not yet available. Do **not** infer a global Toeplitz theorem from local degrees.

## Reproduce

~~~sh
python code/square_degree_cube_v259.py
python -m unittest -v tests.test_square_degree_cube_v259
~~~

The new 25 regression tests include the reference 216-vertex face proof, the asymmetric nonexact-center version, exact noncollapsing square bounds, arbitrary continuous perturbation budgets, a separate v223 geometric orbit census, independently recovered exact root, and negative controls for singular preconditioners, overly large \(\varepsilon\), excessively small cubes and false boundary contacts.

**Sources:** Classical Poincaré–Miranda theorem and Brouwer degree, as surveyed by *Variations on the Brouwer Fixed Point Theorem: A Survey*, *Mathematics* **8** (2020), 501 (https://www.mdpi.com/2227-7390/8/4/501); Benjamin Matschke, *A Survey on the Square Peg Problem*, *Notices of the AMS* **61** (2014), https://www.ams.org/notices/201404/rnoti-p346.pdf. The specialized rational example and verifier are original to this repository; no historical priority is claimed for the general topological theorem.
