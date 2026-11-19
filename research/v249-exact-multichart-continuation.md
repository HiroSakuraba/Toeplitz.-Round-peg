# v249: Exact multi-chart continuation across a moving square-vertex event

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Constructive, exact-rational *finite-atlas* theorem and a nontrivially changing 1-parameter family of Jordan polygons. It proves a noncollapsing square branch across an edge-chart change, including the exact event when a square corner coincides with a moving polygon vertex. **Not a general polygonal cobordism or unrestricted Toeplitz theorem.**

## Why this is different from v248

v248's local Brouwer certificate assumes all four square corners remain strictly inside four selected edges. That sufficient condition naturally fails when a corner hits a polygon vertex: a single fixed edge assignment cannot describe the square on both sides of the event.

v249 replaces this one-chart assumption with a *finite cover of the entire deformation parameter interval*. Different intervals can be certified by different four-edge incidence charts, with no gaps. A square branch can cross an edge endpoint because both adjacent **closed** edges contain the shared vertex at the transition.

Unlike a numerical scan, this work establishes a theorem **for every real parameter value**, using exact polynomial-sign certificates, including endpoints.

## Structural finite-atlas lemma

Let \(P(\lambda)\), \(0\le\lambda\le1\), be a family of simple polygonal Jordan curves whose vertices depend **affinely with rational coefficients** on \(\lambda\). For each ordered choice of four boundary edges \((i_0,i_1,i_2,i_3)\), parametrize an oriented Euclidean square by

\[
z=(p_x,p_y,u,v),\quad
C_0=p,\quad C_1=p+(u,v),\quad
C_2=p+(u-v,u+v),\quad C_3=p+(-v,u).
\]

The condition that each corner \(C_k(z)\) lies on the **supporting line** of edge \(i_k\) gives

\[
M(\lambda)z=b(\lambda),
\]

with polynomial coefficients in \(\lambda\). Define the Cramer determinants

\[
D(\lambda)=\det M(\lambda),\qquad
N_j(\lambda)=\det M^{(j)}_b(\lambda).
\]

On any closed time interval \(I=[a,b]\) where \(D\) is nonzero, the unique square candidate is the rational function \(z_j=N_j/D\), continuous on \(I\). To certify that it is a *real boundary-inscribed square* and has squared side greater than a prescribed rational \(s_*^2>0\), we must prove for **every** \(\lambda\in I\):

1. \(D(\lambda)\) has a fixed nonzero sign.
2. For each edge, the tangential parameter \(t_k(\lambda)\) satisfies **inclusive** bounds \(0\le t_k\le1\). These become signs of polynomials after multiplying by \(D\) with its certified sign.
3. The polynomial \(N_u^2+N_v^2-s_*^2D^2\) is **strictly positive**.

For a rational polynomial \(f(\lambda)\), transform it to the Bernstein basis on \([a,b]\). Every Bernstein basis function is nonnegative and their sum is one. Therefore if every Bernstein coefficient is nonnegative, \(f\ge0\) on that entire interval; if all are strictly positive, \(f>0\). Optional rational bisection and re-evaluation of the Bernstein coefficients handles narrower positivity regions. The verifier **fails closed** if the signs cannot be certified within its limit.

Now suppose a finite collection of such certified chart intervals **covers the whole** \([0,1]\). Then \(P(\lambda)\) has a boundary-inscribed square of side \(>s_*\) for every \(\lambda\). If the square-parameter functions of consecutive charts agree exactly at each shared endpoint, they also produce a **continuous, noncollapsing square branch** throughout the deformation.

This is a rigorous, computationally checkable *finite-atlas continuation lemma*. The explicit coverage and seam checks are essential: the existence of a square at a sampled set of times is **not** sufficient.

## An explicit moving octagon with a chart handoff

Begin with the rational convex octagon used by v248, its vertices \(A_0,\ldots,A_7\) given by integer pairs divided by 20:

| i | X | Y |
|---:|---:|---:|
| 0 | 39 | 21 |
| 1 | 18 | 40 |
| 2 | -17 | 39 |
| 3 | -38 | 20 |
| 4 | -43 | -19 |
| 5 | -20 | -38 |
| 6 | 21 | -38 |
| 7 | 39 | -23 |

This time **change its geometry**: move vertex \(A_3\) as

\[
A_3(\lambda)=A_3+\lambda\left(\frac1{20},-\frac1{30}\right).
\]

Introduce a new subdivision vertex on the unchanged original edge \(A_0A_1\):

\[
V(\lambda)=A_0+q(\lambda)(A_1-A_0),\qquad
q(\lambda)=\frac{39906}{73777}+\frac{\lambda-\frac12}{2}.
\]

The resulting 9-vertex polygon has cyclic vertices

\[
(A_0,V(\lambda),A_1,A_2,A_3(\lambda),
  A_4,A_5,A_6,A_7).
\]

The inserted vertex remains strictly between \(A_0\) and \(A_1\) for all \(\lambda\in[0,1]\). It does not change the geometric boundary, but it changes which *edge* supports the moving square corner. Moving \(A_3\) does change the geometric polygon and the square's position.

### Every polygon in the family is Jordan

This statement is proven for **all real times**, not checked on a finite grid. The original octagon is exactly simple, and each non-adjacent pair of original edges has an axis-aligned coordinate bounding-box separation of at least \(9/10\). Only vertex \(A_3\) moves, by at most \(1/20\) in the infinity norm. Thus no pair of originally non-adjacent edges can meet, because the separation exceeds twice this maximum displacement.

Every successive turn remains strictly counterclockwise through the entire homotopy, established by exact Bernstein signs of the cross products of consecutive affine edges. Every nonzero edge length is similarly certified. The inserted point simply splits the first edge into two nonbacktracking collinear pieces. All vertices stay in the radius-3 Euclidean disk for all \(\lambda\), as shown by checking the affine vertex endpoints and convexity of the disk.

Consequently every \(P(\lambda)\) is a simple closed Jordan polygon lying in the closed radius-3 disk, so its true inradius is **at most 3**.

### The two certified square charts

The underlying octagon's square uses old edges \((0,2,4,6)\). After inserting \(V(\lambda)\), these become chart assignments:

| Time interval | Four ordered supporting edges | First square corner |
|---|---|---|
| \(0\le\lambda\le1/2\) | \((1,3,5,7)\) | Second part of split edge |
| \(1/2\le\lambda\le1\) | \((0,3,5,7)\) | First part of split edge |

The square root is a rational function of \(\lambda\). Exact Bernstein checks prove that *all four* square corners lie in their **closed** assigned segments across each complete half-interval. The determinant is nonzero throughout both intervals, and

\[
u(\lambda)^2+v(\lambda)^2 > \frac{729}{100}
=\left(\frac{27}{10}\right)^2
\]

for **every real \(\lambda\in[0,1]\)**.

At \(\lambda=1/2\), the first corner lies exactly at the moving subdivision vertex \(V(1/2)\). It is the endpoint of both split edges, so both chart descriptions are valid. The Cramer roots agree at that instant. Stronger still, the implementation verifies the **polynomial identity**

\[
N^{\rm left}_j D^{\rm right}
= N^{\rm right}_j D^{\rm left}
\]

for each of the four square parameters \(j\), proving the two rational descriptions actually represent the same continuous branch wherever both are algebraically defined.

A chart that uses the *wrong* side of the subdivision fails its segment-membership constraints. Thus the handoff is a real **edge-chart event**, not a numerical artifact. The full time cover has no unexamined gap.

### Quantitative corollary

For every polygon in this family,

\[
s_{\max}(P(\lambda))>\frac{27}{10},
\qquad
r_{\mathrm{in}}(P(\lambda))\le3.
\]

It follows that

\[
\boxed{J(P(\lambda))>9/10
\quad\text{for every}\ \lambda\in[0,1].}
\]

This is a uniform Big Peg bound on an explicit **nonconstant one-parameter family**, surviving a supporting-edge handoff at a polygon vertex.

## What this does, and does NOT, establish

**Established:**
- A reusable, exact rational **finite-chart certificate** for affine polygon homotopies.
- Proof of square existence **at every real deformation time**, including the edge-transition instant.
- An explicitly continuous nondegenerate square branch, with a positive uniform minimum size.
- A globally valid Jordan-domain check and a restricted uniform Big Peg ratio, without assuming the square stays interior to one four-edge chart.
- Negative controls for invalid interval covers, wrong charts, negative polynomial values, and oversized square cutoffs.

**Not established:**
- A generic square-orbit birth/death transition, a singular determinant crossing, or independent geometric square branches exchanging existence. This example is a **moving subdivision-vertex** event on an evolving boundary, not an arbitrary singular cobordism.
- A finite atlas valid for *every* Jordan curve or every polygon homotopy.
- A global positive Big Peg constant independent of polygon complexity.
- Toeplitz for arbitrary Jordan curves.

The next structural task is to certify genuine **non-collinear vertex crossings** and potentially singular multiple-orbit transitions, with a size cutoff, and to bound the possibility that *all* chart branches simultaneously escape into small squares.

## Reproduce

The certificate needs only Python's standard library, no floating-point decisions and no third-party interval package:

~~~sh
python code/big_peg_multichart_v249.py
python -m unittest -v tests.test_big_peg_multichart_v249
~~~

The exact output includes each chart's closed time interval, ordered edge assignment, number of Bernstein sign checks, and rational endpoint parameters, plus the verified handoff and full-family Jordan property.

## Background and context

- v223 exact complete single-polygon square orbit census (this repository) provides geometric square census machinery but not a general homotopy theorem.
- v248 Brouwer normal-graph persistence (this repository) provides robust nonsmooth local persistence only while all four contacts remain inside selected edges.
- B. Matschke, *A Survey on the Square Peg Problem*, Notices AMS (2014), discusses polygonal approximations and the quantitative Big Peg route.
- D. Bernier's 2026 Big Peg challenge and the Pochuev counterexample: https://dbernier.ca/2026/06/16/the-sharp-big-peg-conjecture/ .
