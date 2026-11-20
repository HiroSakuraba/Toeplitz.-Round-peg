# v256: Centrally symmetric Jordan curves inscribe a large centered square

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Result:** A short topological argument establishing an *exact centered square* on **every centrally symmetric Jordan curve**, with no convexity, smoothness, star-shapedness, rectifiability, or polygonality assumption. A disk **centered at the symmetry center** of radius \(r\) guarantees the square side is at least \(\sqrt2\,r\).

**Priority:** This is a natural classical-symmetry special case of the square-peg problem, not claimed to be a novel historical theorem. The contribution to this repository is a quantitative, precisely scoped proof and a reproducible exact polygon construction that is neither convex nor star-shaped about its center.

## The quarter-turn boundary-intersection lemma

Let \(\gamma\subset\mathbf R^2\) be **any** Jordan curve, and let \(c\) belong to its bounded interior \(\Omega\). Write \(R\) for the positive \(90^\circ\) rotation about \(c\).

**Lemma.** The curves \(\gamma\) and \(R\gamma\) intersect.

**Proof:** Both bounded Jordan domains \(\Omega\) and \(R\Omega\) contain \(c\). Suppose their boundary curves were disjoint. By the Jordan curve theorem, two disjoint Jordan boundaries with a point common to their bounded interiors have strictly nested bounded domains: either

\[
\overline{\Omega}\subset\operatorname{int}(R\Omega)
\quad\text{or}\quad
\overline{R\Omega}\subset\operatorname{int}\Omega.
\]

In the first case, apply the same rotation successively. Strict inclusion is preserved by rotations, and since \(R^4=\mathrm{id}\) we would get a strict cyclic chain

\[
\Omega\subsetneq R\Omega\subsetneq R^2\Omega
\subsetneq R^3\Omega\subsetneq\Omega,
\]

which is impossible. The other nesting direction yields the same contradiction in reverse. Therefore \(\gamma\cap R\gamma\ne\varnothing\). QED.

Because \(v\in\gamma\cap R\gamma\), the two points \(v\) and \(R^{-1}v\) both lie on \(\gamma\); they are equidistant from \(c\) and the two position vectors form a right angle. This **two-corner partial square** conclusion applies to *every* Jordan curve whose interior contains the center.

### The central-symmetry square theorem

Suppose in addition that \(\gamma\) is invariant under the half-turn about \(c\):

\[
R^2\gamma=\gamma.
\]

Choose \(v\in\gamma\cap R\gamma\). From the lemma, \(v\) and \(R^{-1}v=-Rv\) belong to \(\gamma\). By central symmetry, their opposites \(-v\) and \(Rv\) also belong to \(\gamma\), with all vectors measured from \(c\).

The four boundary points in cyclic order

\[
\boxed{v,\quad Rv,\quad -v,\quad -Rv}
\]

are exactly the vertices of a square centered at \(c\). Since \(c\) lies in the open Jordan domain, \(v\ne c\), so the square is nondegenerate.

**Quantitative corollary.** If the bounded Jordan domain contains the closed disk \(D(c,r)\), each boundary point has Euclidean distance at least \(r\) from \(c\). The square side satisfies

\[
\boxed{\operatorname{side}(S)=\sqrt2\,|v-c|
\ge\sqrt2\,r.}
\]

For a circle centered at \(c\) equality holds, so the constant \(\sqrt2\) is sharp for this disk-centered class.

**Crucial limit:** The hypothesis is a disk centered **at the symmetry point**. If the largest inscribed disk of a nonconvex centrally symmetric domain has a *different* center, this argument does **not** establish \(\operatorname{side}\ge\sqrt2\,r_{\mathrm{in}}\) for that arbitrary largest disk. Do not silently replace the centered disk radius by the true inradius.

## An exact NONCONVEX and NON-STAR-SHAPED example

Consider the following 12-vertex polygon, listed in counterclockwise order:

\[
\begin{aligned}
(-4,-4),&(-3,-4),(-3,-2),(-1,-2),(-1,-4),(4,-4),\\
(4,4),&(3,4),(3,2),(1,2),(1,4),(-4,4).
\end{aligned}
\]

It is a square with an *inward rectangular notch* in its upper-right boundary and a centrally symmetric opposite notch on the lower-left. Consecutive antipodal vertices satisfy \(P_{i+6}=-P_i\).

The polygon is **not convex** because its two notches create reflex corners. It is **not star-shaped about the symmetry center**: the boundary point \((4,4)\) is visible on the boundary, but the segment from \((0,0)\) to that point passes through \((5/2,5/2)\), which is **outside** the polygon in the upper notch.

Nevertheless, an independent exact test establishes that \((0,0)\) lies strictly in the bounded domain and every boundary edge is at squared distance at least \(4\) from the origin. Thus the entire closed radius-2 disk lies inside.

The certified quarter-turn intersection produces an exact centered inscribed square with

\[
\operatorname{side}^2\ge 2(2)^2=8.
\]

The [exact rational verifier](../code/central_symmetry_square_v256.py) computes every segment intersection between this polygon and its \(90^\circ\) rotation. It handles ordinary crossings, vertices and collinear overlap endpoints with rational arithmetic, checks that \(v\) and \(R^{-1}v\) really lie on the original polygon, and verifies all four square corners.

This computation is a reproducible example, **not** a substitute for the general proof: the nesting argument applies to arbitrary continuous Jordan curves, even with infinitely many corners.

## Relation to v255 and the broader Toeplitz program

v255 proved that **every square** inscribed on a convex curve containing a disk of radius \(r\) has side at least \(\sqrt2\,r\). This new symmetry class is different: it permits arbitrarily deep concavities and proves the existence of **at least one centered large square**, but makes no claim about the sizes of other squares.

The half-turn symmetry bypasses the hidden-corner obstruction by pairing the two orthogonal boundary points automatically. In the unrestricted Toeplitz problem, that antipodal pairing is unavailable.

The convex-polygon large-square conclusion in v255 also has known prior public discussion: David Bernier asked [whether every inscribed square on a convex polygon is large](https://math.stackexchange.com/questions/5140338/does-every-convex-polygon-have-only-large-inscribed-squares) on June 12, 2026, and Christophe Boilley supplied a positive answer that day. This manuscript **makes no historical priority claim** for v255.

## Reproduction

~~~sh
python code/central_symmetry_square_v256.py
python -m unittest -v tests.test_central_symmetry_square_v256
~~~

Every **algorithmic** input/output decision uses exact \`fractions.Fraction\`. The **topological proof** relies on the Jordan curve theorem; no discrete finite census is claimed to prove that theorem. See [v257](v257-near-central-symmetry.md) for a quantitative approximate-square result when exact half-turn symmetry is broken.

**Open problem:** Obtain an exact square with a polygon-complexity-independent positive side bound for *arbitrary* Jordan curves without relying on central symmetry or convexity. This is not established by either result.
