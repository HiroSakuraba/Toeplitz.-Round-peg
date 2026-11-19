# v223 — Complete rational Jordan-polygon square-orbit census and the vertex-event obstruction

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Scope:** An exact decision/census theorem for one **fixed rational simple polygon**, plus a complete geometric-orbit size-cutoff ledger for **uniform similarity deformations** of that polygon. An explicit strictly convex polygonal deformation demonstrates nongeneric vertex-event behavior. **Neither a general moving-polygon cobordism nor a proof of Toeplitz for arbitrary Jordan curves is claimed.**

## 1. Why v222 was not yet a one-curve square census

[v222](v222-multibranch-cutoff-and-orbit-audit.md) glued square branches on *four independently selected translating polygonal arcs*. Its sources need not arise from one Jordan curve, and its generator counts retained source labels, a fixed orientation and deck lifts. It did not enumerate every square-orbit of a whole closed polygon: a square may use two vertices from the *same polygon edge* or land exactly on a polygon vertex. Furthermore, a geometric square has eight cyclic/reversal ordered descriptions. Counting those eight separately modulo two would erase it.

Here the four points are required to lie on **one verified embedded polygon**. Enumerating every oriented four-edge product, including repeated indices and endpoints, gives a complete finite rational decision on that input polygon. The result is stronger than v218's intentionally limited four-*distinct-open-edge* diagnostic. The procedure also detects **positive-dimensional square families**, which make a naive finite orbit parity undefined.

## 2. Theorem A: complete exact finite-model decision, with continuum detection

Let \(P=(v_0,\ldots,v_{n-1})\) be an oriented simple closed polygon whose vertex coordinates are rational. Let \(\kappa\ge0\) be a rational squared-side cutoff. Its positively sized inscribed squares are the sets of four distinct polygon points forming a Euclidean square with side-square \(q>\kappa\). They need not be gracefully ordered along the curve.

**Theorem A.** There is a finite rational-arithmetic decision procedure that returns precisely one of the following:

1. A complete list of **all geometric square orbits** with \(q>\kappa\), when there are finitely many;
2. A certificate of **no** such squares; or
3. A rational, nondegenerate witness on an explicitly identified **positive-dimensional family** of such squares, correctly refusing to assign finite orbit parity.

A geometric orbit is a *set of four planar vertices*, not an ordering of them. The checker returns for finite orbits their exact four coordinates, normalized polygon parameters, canonical six-sector curve order, strict size, whether any vertex lies on a polygon corner, and the number of redundant edge-chart labelings.

### Proof: exhaustive edge products

Every polygon point belongs to one or two closed original edges. Fix an ordered tuple of four edge indices \((i_0,i_1,i_2,i_3)\), **allowing repetition**, and a square orientation \(\sigma\in\{+1,-1\}\). Each chosen edge has an exact affine parameterization

\[
p_j(u_j)=v_{i_j}+u_j(v_{i_j+1}-v_{i_j}),\qquad0\le u_j\le1.
\]

The four ordered-square equations are exactly the four scalar affine equations

\[
F_\sigma(u)=\bigl(p_2-p_1-\sigma J(p_1-p_0),\ p_3-p_0-\sigma J(p_1-p_0)\bigr)
=Au-b=0,
\]

where \(A,b\) are rational. There are only \(2n^4\) ordered edge/orientation combinations. Triples on one and the same *straight* edge are safely skipped, since three distinct vertices of a nondegenerate square are never collinear.

If \(\operatorname{rank}A=4\), Gaussian elimination gives the only possible solution. Check all four parameters in the closed box, the exact residual and \(q>\kappa\). Thus edge endpoints are included, not thrown away.

If \(\operatorname{rank}A<4\), the feasible set

\[
K=\{u\in[0,1]^4:Au=b\}
\]

is a compact rational polytope. Its vertices are found by selecting the independent rows of \(A\), adding \(4-\operatorname{rank}A\) active unit-box constraints, solving rational square systems, and checking *all* original equations. A feasible polytope is the convex hull of its enumerated vertices. The squared side length \(q(u)=\|p_1(u_1)-p_0(u_0)\|^2\) is a **convex quadratic** on \(K\), hence its maximum is attained at a vertex (since any convex combination has \(q\) no greater than the largest value at its vertices).

- If every vertex satisfies \(q\le\kappa\), the entire chart contributes **no** square above cutoff.
- If \(K\) is a singleton, its rational vertex is one isolated square candidate, even though \(A\) is singular.
- If \(K\) has more than one vertex and one has \(q>\kappa\), the open subset of \(K\) near that vertex also satisfies \(q>\kappa\): there are **infinitely many different square configurations**. Since a simple Jordan polygon has a one-to-one parameterization, varying a parameter genuinely varies the polygon point. The implementation raises `PositiveDimensionalSquares` and supplies a rational witness rather than inventing a finite count.

Every inscribed square belongs to at least one tuple of closed edges, so this finite procedure is **complete** for rational polygonal inputs. The inference is independent of numerical sampling or a global square-peg theorem. Because affine constraints have rational coefficients and a rational feasible polytope has rational vertices, the existence of *any* square above cutoff entails an exact rational-coordinate witness in this finite model.

### Geometric-orbit normalization and six sectors

Every nondegenerate geometric square has **eight ordered representations** (four cyclic rotations times two orientations). At a polygon corner, several *edge charts* may describe the very same point and the same ordered square. The checker normalizes each vertex's curve parameter modulo the total number of edges and keys each orbit by its lexicographically sorted set of four planar points. On an embedded simple polygon this identifies one *geometric* square exactly once. It reconstructs the positive-orientation vertex order starting at the lexicographically least geometric vertex, then sorts the other labels by cyclic polygon parameter. This gives one of the **six** cyclic order sectors; it does **not** assume every square is graceful.

**Caveat:** The number of *edge-chart presentations* can exceed eight when a square vertex hits a polygon corner. That number is diagnostic metadata, **not** the physical orbit multiplicity.

## 3. Exact positive and adversarial examples

| One rational Jordan polygon | Certified outcome | Scientific point |
|---|---|---|
| \((0,0),(2,0),(3,3),(0,2)\) | Exactly **one** geometric square, \(q=72/25\), 8 edge-chart labelings | Agreement with independent v218 four-distinct-open-edge case |
| \((0,0),(2,0),(2,3/2),(0,1/2)\) | Exactly **two** squares, \(q=45/49\) and \(q=1\), with 8 and 16 edge-chart labelings | Repeated-edge and polygon-vertex cases can alter finite counts |
| \((0,0),(2,0),(3,1),(2,2),(0,1)\) | Exactly **one**, \(q=20/9\), 16 chart presentations | A five-edge input, square touching a vertex |
| Unit square \((0,0),(1,0),(1,1),(0,1)\) | **Continuum** of inscribed squares with positive side | Finite raw orbit parity not defined without regularization |
| \(2\times1\) rectangle | **Continuum**, including products with repeated edges | Repeated-edge singular families must not be ignored |

The checker also proves there are **no** squares exceeding a sufficiently large positive size cutoff in these examples, even when unfiltered configurations form a continuum. A strict filter at the exact side-square \(q=\kappa\) excludes that orbit, by definition.

## 4. Theorem B: exact geometric-orbit cutoff incidence for a genuine one-polygon similarity deformation

Let a *fixed rational Jordan polygon* \(P\) have a finite set of geometric square orbits of sides squared \(q_1,\ldots,q_m>0\) (certified by Theorem A). Define, for \(t\in[a,b]\),

\[
P(t)=s(t)P+z(t),\qquad s(t)=s_0+v t>0,
\]

where \(z(t)\) is any rational affine translation. Every \(P(t)\) is an **embedded Jordan polygon**. The similarity gives a *bijection of geometric square orbits* and preserves their curve-order sectors. In particular each orbit has squared side \(q_j(t)=q_j s(t)^2\).

Fix \(\kappa>0\), with no orbit exactly at cutoff at \(a\) or \(b\). If \(v\ne0\), the positive scale is strictly monotone, and each orbit crosses the cutoff at most once. The exact crossing condition is

\[
q_j(s_0+vt)^2=\kappa,
\qquad t_j=\frac{\sqrt{\kappa/q_j}-s_0}{v}.
\]

A time is returned *exactly* if rational; otherwise the implementation gives a **certified rational sign-changing isolating interval** narrower than \(2^{-40}\). This uses only `Fraction` arithmetic and integer square roots. With the strict cutoff count \(N_\kappa(t)=\#\{j:q_js(t)^2>\kappa\}\),

\[
\boxed{N_\kappa(a)+N_\kappa(b)\equiv
\#\{\hbox{cutoff-crossing geometric orbits}\}\pmod2.}
\]

This is an exact **geometric-orbit**, rather than eightfold labeled, incidence identity on a family of **one embedded polygon**. It includes no untracked source-arc boundary because the entire polygon is transformed together.

Examples:

- For the one-square quadrilateral above, \(q=72/25\), \(s(t)=1/2+t\), \(t\in[0,1]\), \(\kappa=1\): the count changes **0 to 1** at the irrational time \(t=\sqrt{25/72}-1/2\), with one certified crossing.
- For the two-square trapezoid, \(s(t)=1/4+t\), \(\kappa=1/4\), \(t\in[0,1]\): the count changes **0 to 2**, with one rational crossing \(t=1/4\) for \(q=1\), and one irrational crossing for \(q=45/49\). The total event parity is zero.
- Under pure translation, all sizes and orbits are unchanged and there are **no** cutoff events.

This family is deliberately restricted: similarities preserve the entire orbit set, so **no birth/death or internal corner event can occur**. It is not a proof of a general polygon-homotopy incidence theorem.

## 5. A nongeneric vertex-event obstruction inside a strictly convex homotopy

Consider the **same** 4-vertex polygon moving affinely:

\[
P(t)=((0,0),(2,0),(3-t,3-3t/2),(0,2-3t/2)),
\quad 0\le t<4/3.
\]

The four consecutive oriented corner cross-products are

\[
6-3t,\qquad \frac{3t^2-13t+16}{2},\qquad
\frac{(t-3)(3t-4)}{2},\qquad4-3t.
\]

All are **strictly positive** throughout this time interval: the quadratic has negative discriminant, and the other expressions have the stated signs. Thus this is a homotopy through **strictly convex embedded Jordan quadrilaterals**, without an embedding failure.

Exact complete census snapshots show:

| Time | Geometric square orbits | Polygon-vertex square at that time? |
|---|---:|---|
| \(0\) | 1 | No |
| \(63/64\) | 1 | No |
| \(1\) | **2** | **Yes**, one of the two |
| \(33/32\) | **3** | No |

One square of side-square \(q=1\) passes through an exact polygon vertex at the special time \(t=1\), and the integer orbit count there has parity **even**, whereas at both nearby sampled *regular* times the counts are odd. This is a **nongeneric polygonal corner event**, not evidence against the published generic smooth parity theorem and not a complete event classification on this time interval. The endpoint census is exact; the current module does **not** trace every solution branch or isolate every intervening algebraic event time.

This example demonstrates why a whole-polygon geometric orbit **census** still does not automatically deliver a **homotopy invariant**. Time-varying edge directions cause the square equations to be bilinear in the time and edge parameters; determinant singularities, vertex contacts, repeated-edge families and size-filter crossings all require rigorous event processing.

## 6. What v223 genuinely establishes versus what remains open

| Statement | Status |
|---|---|
| Complete exact decision/census of inscribed square orbits on a fixed rational Jordan polygon, including repeated edges and corners | **Proved and implemented**; continuum cases reported explicitly |
| Every rational-polygon square existence claim has a rational-coordinate witness | **Proved in this PL finite model** |
| Symmetry-free geometric orbit count, with six-sector classification | **Implemented**; no forced graceful assumption |
| Complete geometric-orbit cutoff ledger for one-polygon positive similarity homotopies | **Proved and implemented** (explicit restricted family) |
| Strict convexity of the above vertex-event family and the four exact snapshot counts | **Proved/tested**; no generic parity claim at the singular time |
| Complete event ledger for a general vertexwise-affine polygon deformation | **Not implemented** |
| Well-defined cutoff-independent nonzero mod-2 class for any Jordan curve | **Not established** |
| Uniform positive-size witness under arbitrary continuous approximation | **Not established** |
| Unrestricted Toeplitz square-peg conjecture | **Open** |

## 7. Recommended v224: exact algebraic event isolator for a moving Jordan polygon

For \(v_i(t)=v_i^{(0)}+t w_i\), each product-edge square residual is affine in the four edge parameters but has coefficients **linear in \(t\)**. Accordingly the unique regular root is a rational function of \(t\), with denominator \(\det A(t)\). A complete homotopy ledger would need certified root/event isolation for:

1. \(\det A(t)=0\) and feasible singular square families;
2. parameters \(u_i(t)=0,1\), including simultaneous polygon-vertex crossings;
3. positive size-cutoff crossings, obtained by clearing denominator in \(q(t)-\kappa\);
4. joint events between different charts, repeated-edge assignments, and all eight geometric labelings;
5. polygon simplicity/convexity (or general embedding safety) **for the entire time interval**, not only at selected snapshots.

Rational Sturm sequences or other exact algebraic-number sign methods are promising, but the implementation must be genuinely executable in CI. Any introduction of a computer algebra dependency must be explicit. A local exact census at finitely many sample times **cannot** silently substitute for a certified exhaustive cobordism.

Even completing this exact polygonal homotopy ledger would not solve the **rough-curve compactness** barrier exposed by [v220](v220-fixed-rectifiable-collision-link-obstruction.md). The v218–v222 symmetry/size-boundary results remain mandatory negative controls for any proposed universal invariant.

## 8. Reproduction and source context

```bash
python -m unittest discover -s tests -v
python code/polygon_orbit_census_v223.py
```

The solver is Python standard-library-only, uses exact `fractions.Fraction` and `math.isqrt`, and deliberately runs no floating-point or random search to *certify* a result. The tests exercise valid/invalid Jordan polygons, isolated and continuum square families, strictly positive size cutoffs, corner and repeated-edge duplication, cyclic starting-point invariance, reversal, geometric reflection, all-orbit similarity incidence, rational/irrational cutoff times, and the convex vertex-event snapshots. Independent comparisons to the narrower v218 distinct-open-edge checker and v216 polygon validator are included in the **full-repository** suite. These tests supplement, but do not replace, the mathematical proof above.

**Prior art/context, not claims of originality or proof of Toeplitz:**

- Cantarella, Denne and McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, [arXiv:2103.07506](https://arxiv.org/abs/2103.07506), on generic smooth-curve square counts and transversality.
- Matschke, *On the Square Peg Problem and Its Relatives*, [arXiv:1001.0186](https://arxiv.org/abs/1001.0186), on square-peg classes and parity/background.
- Internal [v218 global symmetry audit](v218-global-configuration-symmetry-and-collision-audit.md), [v220 fixed-curve collision obstruction](v220-fixed-rectifiable-collision-link-obstruction.md), [v221 positive-cutoff incidence](v221-positive-size-cutoff-incidence.md), [v222 multi-branch PL incidence](v222-multibranch-cutoff-and-orbit-audit.md).
