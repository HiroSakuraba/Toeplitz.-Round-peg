# v278 — Exact falsifiers for the atlas and continuation gates

The machine ledger records eight failed candidate statements, their gate, and exact evidence. The controls distinguish equality of geometric configurations from scalar coincidences, relative parity from raw counts, and a size estimate for one witness from a size estimate for another. [Exact replay](../code/adversarial_rectangles_v278.py) · [Machine ledger](../results/v278/adversarial-ledger.json)

## A nonconvex double-vertex birth

Consider the rational-affine quadrilateral
\[
P(t)=((0,0),(4,t),(1,1),(0,3)),\qquad |t|\le 1/100.
\]
The turns in cyclic order are \(12,4-t,-5-t,3\). The two nonadjacent edge pairs have coordinate gaps at least \(99/100\) in \(y\) and \(1\) in \(x\), respectively. These facts give a simple, consistently oriented concave polygon across the whole interval.

There are two selected square branches that meet at the unit square at \(t=0\). The common square touches two actual vertices of \(P(0)\), namely \((0,0)\) and \((1,1)\). For the first branch, set \(s=12/(12+t)\) and \(b=ts/4\). Its cyclic corners are
\[
(0,b),\ (s,b),\ (s,b+s),\ (0,b+s).
\]
For the second, set \(s=48/(48+4t+t^2)\). Start at \((ts/4,t^2s/16)\) and take successive side vectors \((s,ts/4)\) and \((-ts/4,s)\).

Both formulas solve their supporting-line equations identically. In each branch two contacts have inequalities of the form \(t h(t)\ge0\) with \(h(t)>0\) over the entire closed interval. The verifier recomputes these identities and uses rational Bernstein signs on the reparameterized interval. The other contacts remain strict. Both squared sides are greater than \(9/10\) as supporting-line squares throughout. Consequently neither branch lies on the required segments for \(t<0\), they give one common geometric square at \(t=0\), and both lie on the boundary for \(t>0\). In the latter case the first square has a horizontal side and the second has none, so they are different squares. This is a positive-size local birth, with two simultaneous vertex contacts at the shared square.

The independent v223 square solver enumerates **all positive geometric square orbits** at the three rational slices. This enumeration solves the square equations directly on edge products and is separate from v274's rectangle-angle charts.

| Time \(t\) | Whole-polygon geometric square count | Scope of the count |
|---:|---:|---|
| \(-1/100\) | 1 | All squares at this slice |
| \(0\) | 2 | All squares at the double-vertex event |
| \(1/100\) | 3 | All squares at this slice |

At \(t=0\), the second square has cyclic vertices \((0,8/11),(9/11,0),(17/11,9/11),(8/11,17/11)\) and squared side \(145/121\). Both square orbits have source order \((0,1,2,3)\) in their recorded counterclockwise labels and pairing \(02\mid13\). The unit square contributes 32 edge-label presentations and the other contributes eight. Neither the 40 presentations nor the two raw singular-time geometric points is a perturbation-invariant relative count. This is a counterexample to unrestricted **odd raw geometric parity** on polygons with singular vertex contacts. The continuously tracked result here concerns the two listed branches. The independent whole-polygon counts are claimed at the three displayed slices.

## Uniform embeddings with a small square and a large square

The v272 bump produces a square of side \(\epsilon=1/n\) on the top detour of a unit-square polygon. The four original unit-square vertices also remain boundary points of every approximant. They form a second square of side \(1\).

For convergence in a common parameterization, assign one time unit to each of the four sides. Use the unchanged affine parameterization on the bottom, right and left sides. Parametrize the modified top side by its normalized arclength \(s\in[0,1]\), including the detour, and use \((1-s,1)\) on the top of the limit square. The detour length is \(1+2\epsilon\). If \(V(s)\) is vertical length traversed so far, then
\[
x(s)-(1-s)=V(s)-2\epsilon s,\qquad 0\le V(s)\le2\epsilon,\quad
0\le y(s)-1\le\epsilon .
\]
Thus the uniform coordinate error is at most \(2\epsilon\) for actual Jordan embeddings on the same circle of parameters. The replay independently checks the piecewise-linear bound at every time breakpoint of a rational specimen. The chosen small square collapses at \((1/3,1)\), while the side-1 square stays exactly fixed. A size floor for the latter alone says nothing about continuation of a class represented by the former.

## Falsifiers and repairs

| Failed candidate | Exact control | Required repair |
|---|---|---|
| Every positive-size component has a finite square-zero count | A rational quarter-turn on the radius-\(\sqrt2\) circle is entirely squares of side \(2\). The finite-loop checker refuses a count. | Nonisolated zeros need a relative intersection class or a justified perturbation. |
| Every Jordan curve has an odd-monodromy rectangle loop | The complete noncircular ellipse locus is an open interval. The flat-end signs on semiaxes \(3,2\) are opposite and have relative parity one. | Admit intervals with compatible nonzero boundary sections. [v268 ellipse classification](v268-rectangle-monodromy-and-protected-square-criterion.md) |
| Matching center, diagonal radius and diagonal dot gives closure | A rational rotation of a \(6\times8\) rectangle preserves all three values and has different endpoint corner sets. The loop checker rejects it. | Compare corner coordinates and the actual cyclic lift. |
| Complete enumeration in one edge cell gives whole-curve completeness | The v271 bump chart has an exact interval of squares, and the original unit square lies outside its declared parameter box. The independent square census refuses a finite count for the whole polygon. | Audit all cells and geometric transitions independently of the restricted chart. |
| Every finite raw geometric square count is odd | \(P(0)\) has two squares, one at a double-vertex birth. | Prove a relative index on a stated admissible category and define contributions of singular events. |
| Positive squares can die in pairs only at a size cutoff or a matrix rank loss | The v254 convex pentagon has counts \(3,2,1\) before, at and after its isolated quadratic event throughout the v263 window, with all squares of squared side greater than \(1\). The joining Cramer charts are nonsingular. | Handle glued vertex events at positive size. [v254 geometry](v254-genuine-square-pair-death.md) |
| Uniform Jordan-embedding convergence protects every selected square | The parameterized moving bump above contains both a shrinking square and a persistent square of side \(1\). | Prove a size floor for a representative of the same class. |
| The bare collision-gap direction determines the normalized residual | On the same fixed v220 rectifiable curve, the equal-gap square and flat sequences give distinct exact normalized residual limits. | Track additional boundary data or exclude total collision for the protected class. [v220 limit proof](v220-fixed-rectifiable-collision-link-obstruction.md) |

These controls leave the unrestricted class definition, its nonvanishing with admissible collision boundaries, and a common size estimate under Jordan approximation open. The fixed rectifiable example addresses the chosen collision chart and belongs to an existing affirmative regularity class.

Run `python code/adversarial_rectangles_v278.py --check` from the repository root to recompute the ledger and `python -m unittest discover -s tests -p test_adversarial_rectangles_v278.py` for the failure controls.
