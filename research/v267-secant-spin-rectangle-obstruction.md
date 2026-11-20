# v267: Secant–spin, the Möbius chord twist, and a conditional rectangle-to-square obstruction

**Research program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** An exact topological lemma for **every continuous Jordan curve**, complete rational verification on polygonal examples, a precise **conditional** square-forcing statement using a twisted line bundle over the **nondegenerate rectangle locus**, and an exact shrinkout stress test. **Not a proof of the unrestricted Toeplitz conjecture.**

This is deliberately **v267**, not v266: [v266](v266-two-corner-stability-and-convex-cutoff.md) already contains the convex two-near-corner, stability, and geometric size-cutoff lemmas.

## 1. A true spin-like invariant that needs no tangent

Let \(\gamma:\mathbb R/\mathbb Z\to\mathbb R^2\) be an arbitrary **continuous injection**. Set

\[
v(t)=\gamma(t+\tfrac12)-\gamma(t).
\]

A Jordan embedding makes **all** these secants nonzero. Moreover,

\[
\boxed{v(t+\tfrac12)=-v(t).}
\]

Normalize \(u(t)=v(t)/|v(t)|\in S^1\) and choose a continuous lifted angle \(\theta:\mathbb R\to\mathbb R\) satisfying \(u(t)=\exp(i\theta(t))\). The equality \(u(t+1/2)=-u(t)\) implies

\[
\theta(t+\tfrac12)-\theta(t)=(2k+1)\pi
\]

for some integer \(k\); because the left-hand side is continuous in \(t\) but belongs to a discrete set, \(k\) is constant. Applying the identity twice gives

\[
\theta(t+1)-\theta(t)=2(2k+1)\pi.
\]

Hence

\[
\boxed{\deg u=2k+1\equiv1\pmod2.}
\]

This is an **actual, universal, nonzero topological winding invariant** of every Jordan embedding, including fractal and nowhere-differentiable examples. No derivative, curvature, or ambient polygon approximation is assumed. In the examples of Jordan polygons checked by v267, its integer winding agrees with the curve's orientation, \(+1\) for CCW and \(-1\) for CW.

Defining a half-angle lift

\[
\psi(t)=\exp(i\theta(t)/2)
\]

gives \(\psi(t+1)=-\psi(t)\). This is **spinor-like monodromy**, not a claim that the curve's electron spin is \(1/2\), nor a canonical physical spin structure or chirality of its Euclidean shape.

In particular a perfectly symmetric circle satisfies the same odd-degree property. We are detecting **ordered-pair orientation**, not visual mirror asymmetry.

### The precise line bundle

The space of **unordered distinct pairs** of points on \(S^1\), \(\mathrm{UConf}_2(S^1)\), is an **open Möbius strip**. The particular path

\[
t\in[0,\tfrac12]\mapsto\{\gamma(t),\gamma(t+\tfrac12)\}
\]

is a closed loop in that space, because its endpoints are the *same unordered pair*. An oriented chord vector along the path starts at \(v(0)\) and ends at \(-v(0)\). Thus the double cover choosing a **direction** for the unordered chord is nontrivial over this loop.

If \(L\) is its associated real sign line,

\[
\boxed{\langle w_1(L),[\text{antipodal chord loop}]\rangle=1.}
\]

This is the mathematical Möbius twist under chord-endpoint exchange. It exists regardless of how irregular the Jordan embedding is.

**Prior art:** This unordered-chord Möbius geometry is **not new**. It is central to Vaughan's theorem that every Jordan curve contains an inscribed rectangle, as described in B. Matschke, *A survey on the Square Peg Problem*, Notices AMS **61** (2014), 346–352, DOI [10.1090/noti1100](https://doi.org/10.1090/noti1100). The new use here is to audit the *different* twist attached to perpendicularity in the full four-point square equations.

## 2. Which chord equations characterize a square?

Take four distinct points \(p_0,p_1,p_2,p_3\in\gamma(S^1)\), grouped into **two candidate diagonals**. Write

\[
d_1=p_2-p_0,\qquad d_2=p_3-p_1,
\]
\[
m_1=\tfrac12(p_0+p_2),\qquad m_2=\tfrac12(p_1+p_3).
\]

Define the four-component *diagonal square residual*

\[
\boxed{\mathcal F(p_0,p_1,p_2,p_3)=
\left(
m_{1x}-m_{2x},\
m_{1y}-m_{2y},\
|d_1|^2-|d_2|^2,\
d_1\cdot d_2
\right).}
\]

For four pairwise distinct points, \(\mathcal F=0\) if and only if the two chords are the **diagonals of a nondegenerate square**. The first two conditions make a parallelogram, the third makes it a rectangle, and the fourth makes the rectangle a square.

If all four points collide, this residual also vanishes. Therefore the **nondegeneracy requirement is essential** and must not be removed by compactifying the source.

For **arbitrary nonconvex Jordan curves**, the curve's cyclic order need not agree with the square's geometric order. There are precisely **three** ways of partitioning four curve-labeled points into two unordered diagonal pairs:

\[
(01|23),\quad(02|13),\quad(03|12).
\]

The exact checker tries all three, and verifies that an actual positive square has precisely one square-diagonal pairing. We do **not** assume a graceful sector.

## 3. The exact dihedral symmetry audit shows the old sign obstruction survives

Let \(R\) be cyclic relabeling \(p\mapsto(p_1,p_2,p_3,p_0)\) in a fixed geometric-square labeling. Then

\[
d_1\mapsto d_2,\qquad d_2\mapsto-d_1,\qquad
m_1-m_2\mapsto-(m_1-m_2),
\]

so

\[
\boxed{\mathcal F(Rp)=-\mathcal F(p).}
\]

The target matrix is \(-I_4\), with determinant **+1**. However, the four-cycle on the **source** parameters has determinant **−1**. Thus, for a transverse zero, cyclic relabeling changes its ordinary integer Jacobian index's sign:

\[
\boxed{\operatorname{ind}(Rp)=-\operatorname{ind}(p).}
\]

The four cyclic lift signs alternate \((+,-,+,-)\). This reproduces the exact cancellation already established for a different square residual in [v218](v218-global-configuration-symmetry-and-collision-audit.md) and [v262](v262-cyclic-quotient-relative-parity.md).

Under reflection \(S:p\mapsto(p_0,p_3,p_2,p_1)\), the *second diagonal reverses*. The first three residual components stay fixed, the last changes sign:

\[
\boxed{\mathcal F(Sp)=
\operatorname{diag}(1,1,1,-1)\mathcal F(p).}
\]

Both source and target transformations now have determinant **−1**, so the local Jacobian index is unchanged by this reflection.

The checker proves both transformations as **universal rational polynomial identities**, not by sampling shapes. Every component of \(\mathcal F\) is a polynomial of total degree at most **2** in eight Cartesian coordinates. The set

\[
0,\quad e_i,\quad2e_i,\quad e_i+e_j\ (i<j)
\]

has **45** elements and uniquely determines every constant, linear, square and mixed quadratic coefficient. Evaluating the two proposed identities at these 45 exact rational nodes verifies all coefficients. **90 exact polynomial-identity probes** cover both generators.

### An explicit nowhere-zero section for this *new* target representation

We should not simply transfer v219's section to the diagonal residual, because their target representations are not identical. Instead, construct a **new section directly**.

On the positive cyclic-order four-point configuration component, v219 established a center coordinate \(c\in\mathbb R/\mathbb Z\) satisfying

\[
c(Tx)=c(x)+\tfrac14.
\]

Our new target generator is precisely \(R_{\mathrm{target}}=-I_4\). Define

\[
\boxed{s(x)=\bigl(\cos(4\pi c(x)),\ \sin(4\pi c(x)),\ 0,\ 0\bigr).}
\]

Then

\[
s(Tx)=\bigl(-\cos(4\pi c(x)),-\sin(4\pi c(x)),0,0\bigr)
      =R_{\mathrm{target}}s(x),
\]
\[
\boxed{\|s(x)\|^2=1\quad\text{for every configuration }x.}
\]

Thus the rank-four associated bundle for the **diagonal** square residual on this cyclic quotient also admits a **globally nowhere-zero equivariant section**, entirely independently of the specific section constructed for v219's other residual. Its absolute top \(w_4\) obstruction is zero.

The code verifies this quarter-turn equivariance and unit norm with exact rational unit-phase pairs; the full all-phase formula and nonvanishing follow analytically from trigonometry. As in v219, this says nothing about a *relative* section constrained by physically meaningful collision-boundary data and makes no assertion about an all-sector compactification.

**Conclusion:** A different square-defect formula does **not** evade the v218 signed-degree cancellation or the v219-style absolute-bundle no-go. The Möbius sign twist is genuine, but the naive inference **“odd secant winding therefore a square” is unproved**.

## 4. The relevant twist for perpendicularity is a tensor product

Suppose we are **already on the locus of nondegenerate inscribed rectangles**:

\[
\mathcal R_\gamma=
\left\{
\{c_1,c_2\}:
\begin{array}{l}
c_i \text{ are two disjoint, nonzero, unordered chords on }\gamma,\\
m_1=m_2,\ |d_1|=|d_2|
\end{array}
\right\},
\]

with suitable identification under exchanging the two chord labels.

On this locus, the three midpoint/equal-length residual equations are already zero. The only missing square equation is

\[
f=d_1\cdot d_2.
\]

Let \(L_1,L_2\) be the two **orientation-sign lines** of the unordered diagonals. Reversing only the first chord changes \(f\mapsto-f\); reversing only the second does the same. **Reversing both gives \(f\mapsto f\).**

Consequently, \(f\) is a continuous section of the **tensor product**

\[
\boxed{L_\perp=L_1\otimes L_2,}
\]

not of either chord line by itself. Its Stiefel–Whitney class obeys

\[
\boxed{w_1(L_\perp)=w_1(L_1)+w_1(L_2)\pmod2.}
\]

Exchanging the diagonal labels also preserves their dot product and canonically swaps the tensor factors. The exact verifier checks the full elementary monodromy table:

| Chord orientation flips along a loop | \(w_1(L_1)\) | \(w_1(L_2)\) | \(w_1(L_\perp)\) |
|---|---:|---:|---:|
| Neither | 0 | 0 | 0 |
| First only | 1 | 0 | **1** |
| Second only | 0 | 1 | **1** |
| Both | 1 | 1 | **0** |

That last row is crucial. **The fact that each individual chord is Möbius-twisted does not imply that perpendicularity is twisted**: their twists may cancel.

### A rigorous conditional square-forcing theorem

**Proposition.** Let \(C\subset\mathcal R_\gamma\) be a compact loop of **strictly nondegenerate** inscribed rectangles, with chord pairings tracked continuously and the corresponding tensor sign line \(L_\perp|_C\) well-defined. If

\[
\boxed{\langle w_1(L_\perp),[C]\rangle=1,}
\]

then \(\gamma\) inscribes a positive-size square.

**Proof.** The perpendicularity \(d_1\cdot d_2\) is a continuous section of \(L_\perp|_C\). If it never vanished, division by its norm would define a nowhere-zero section and **trivialize the real line bundle** along \(C\). But a real line bundle with nonzero \(w_1\) on \(S^1\) cannot be trivial. Therefore the section has a zero; by the rectangle equations, the two nonzero disjoint chords are equal-length, same-midpoint and perpendicular, giving a nondegenerate inscribed square. QED.

The hypothesis is stronger than odd antipodal winding. **Nothing in v267 proves that such a compact, protected, nontrivially tensor-twisted rectangle loop exists on every Jordan curve.** That is precisely the new explicit obligation. A loop allowed to enter collapsed-chord strata could satisfy the same algebraic zero with *no positive square*.

This proposition is a standard line-bundle obstruction applied to square diagonals, **not** a newly solved unrestricted case.

## 5. Exact example of the conditional square mechanism

Consider the **fixed** Jordan rectangle \(\partial([-2,2]\times[-1,1])\). For each rational \(0<x\le2\), choose four boundary points

\[
p_0=(-x,-1),\ p_1=(x,-1),\
p_2=(x,1),\ p_3=(-x,1).
\]

They are four *distinct* corners of an inscribed rectangle. Their diagonal midpoints agree, and the diagonals have equal squared length. Direct substitution gives

\[
\boxed{\mathcal F(p_0,p_1,p_2,p_3)=
(0,0,0,4(1-x^2)).}
\]

At \(x=1/2\), the last coordinate is \(+3\). At \(x=3/2\), it is \(-5\). Since this entire rectangle path stays nondegenerate, continuity forces a square. The exact zero is \(x=1\); its side length is **2**.

This is a **controlled sign-crossing along an interval**, not an example of a nontrivial tensor line over a closed rectangle loop. It illustrates the extra information required by the conditional proposition: **the full rectangle equations remain satisfied, the two chords stay nonzero, and the perpendicularity sign is compelled to cross zero**.

## 6. The decisive shrinkout test: macroscopic antipodal chords, microscopic inscribed square

Reuse the exact [v219 microscopic-notch family](v219-equivariant-quotient-collision-link.md). For integer \(n\ge2\), write \(h=1/n\) and use the eight-vertex simple Jordan polygon

\[
P_n:\quad
(0,0),\ (1,0),\ (1,-h),\ (1+h,-h),\
(1+h,0),\ (2,0),\ (2,1),\ (0,1).
\]

The four consecutive vertices

\[
(1,0),\ (1,-h),\ (1+h,-h),\ (1+h,0)
\]

form an **exact inscribed square of side \(h\)**. These curves converge to a rectangle in the \(C^0\) sense after compatible reparameterization; the selected square shrinks:

\[
\boxed{s_n^2=1/n^2\to0.}
\]

But parameterize every polygon with **equal time \(1/8\) per edge** and form the *antipodal* secant \(v_n(t)=P_n(t+1/2)-P_n(t)\). Between its consecutive breakpoints \(t=i/8\), the vector is affine. Over the first four segments, its exact endpoint vectors are

\[
(1+h,0),\quad(1,0),\quad(1,1+h),\quad
(-1-h,1+h),\quad(-1-h,0).
\]

On each one of these four segments, **at least one Cartesian coordinate has magnitude \(\ge1\) throughout**. The second half of the path is their negative, by antiperiodicity. Consequently,

\[
\boxed{\min_{t\in S^1}|v_n(t)|^2=1
\quad\text{for EVERY }n\ge2.}
\]

The exact checker independently establishes this minimum by minimizing squared norm on **every rational affine secant segment**, and computes winding degree \(+1\) for \(n=2,3,4,8,64,128\). Thus we have simultaneously:

\[
\boxed{
\deg(v_n/|v_n|)=1,\quad
\min |v_n|=1,\quad
s_n=1/n\longrightarrow0.
}
\]

This is sharper than merely saying the winding is odd: **even a uniform positive lower bound for the chosen antipodal secants does not bound the size of *every* square**.

**Critical limitation:** The limiting rectangle has ordinary positive-size inscribed squares, and the \(P_n\) can have other large squares as well. This is **not a counterexample** to Toeplitz or to a potential *specially selected* large-square continuation. It only invalidates an insufficient implication from antipodal winding/length alone to all-square noncollapse.

## 7. Reproducible exact checks

The [v267 code](../code/secant_spin_rectangle_v267.py) uses only rational arithmetic and the existing exact Jordan polygon validator. The tests cover:

- exact antipodal reversal, odd secant winding, and curve-dependent positive secant length floors on CCW, CW, triangular, and **nonconvex simple polygons**;
- all three square-diagonal pairings, collapsed-square exclusion and wrong pairing controls;
- all **90** universally sufficient polynomial identity evaluations of the two dihedral symmetry generators, plus a separate **explicit nowhere-zero equivariant section** for this diagonal target representation;
- tensor line-sign monodromy, including the cancellation under reversing **both** diagonals;
- a strictly nondegenerate, complete one-parameter inscribed-rectangle family with a sign-forced square at \(x=1\);
- the six exact notch polygons with secant floor exactly **1** and chosen square side **\(1/n\)**;
- fail-closed responses to crossed/self-intersecting polygons and insufficient rectangle sign changes;
- retention of [v218](v218-global-configuration-symmetry-and-collision-audit.md) and [v219](v219-equivariant-quotient-collision-link.md) regression suites.

```sh
python code/secant_spin_rectangle_v267.py
python -m unittest -v tests.test_secant_spin_rectangle_v267
python -m unittest -v tests.test_configuration_symmetry_v218 tests.test_collision_quotient_v219
```

The v267 tests are exact **finite polygon** checks. The proof of odd secant degree for an arbitrary continuous Jordan embedding and the conditional twisted-line theorem are elementary **pen-and-paper topological arguments**; passing tests does not itself machine-certify these infinite-dimensional claims.

## 8. Relationship to the literature and prior repo attempts

- **Vaughan's rectangle theorem and Möbius pair geometry:** documented in B. Matschke, [*A survey on the Square Peg Problem*](https://doi.org/10.1090/noti1100), *Notices AMS* **61** (2014), 346–352. The unordered-chord twist is classical.
- **Configuration-space odd-square counts under regularity:** Cantarella, Denne & McCleary, [*Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*](https://arxiv.org/abs/2103.07506), 2021. Their smooth transversality approach requires boundary-disjoint configuration strata; it does not supply unrestricted \(C^0\) shrinkout control.
- **Symplectic/Floer noncollapse and rectangle spectral invariants:** Greene & Lobb, [*Floer homology and square pegs*](https://arxiv.org/abs/2404.05179) (2024), and [*Square pegs between two graphs*](https://doi.org/10.4171/CMH/619) (published online 2026). These provide genuinely additional geometric control in restricted settings.
- **v218–v219:** exact cyclic source/target orientation cancellation, the explicit globally nonvanishing absolute associated-bundle section, and the microscopic-notch selected-square shrinkout.
- **v262–v265:** nonzero *relative* quotient mod-2 degree for a certified complete moving convex pentagon atlas, and exact gluing of a square pair at an irrational polygon-vertex event.
- **v266:** independent convex two-corner and automatic size-cutoff results, unavailable for arbitrary nonconvex Jordan curves without further arguments.

**No-go:** The raw universal secant winding and Möbius chord twist **alone** cannot be substituted for a nonzero relative square residual on the whole configuration space. The cyclic sign cancellation remains, the absolute rank-four Euler obstruction studied in v219 still vanishes, and the chord twist does not select a size-protected family of rectangles.

**Potential v268:** Form a **rectangle-locus relative twisting ledger** (rather than a global four-point Euler class): identify the tensor line \(L_1\otimes L_2\) on genuine nondegenerate rectangle components, compute \(w_1\) on closed components or relative endpoints, and see whether any boundary class can be kept away from the total-collision stratum under \(C^0\) limits. Include all three chord pairings and test the v219 notch, v210 coils, and v266 hidden-corner constraints before claiming a universal invariant.

**Unrestricted Toeplitz square-peg conjecture: OPEN.**
