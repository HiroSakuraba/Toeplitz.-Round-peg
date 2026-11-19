# v247: A strictly smaller TRUE Big Peg ratio for a rational decagon

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Rigorous rational certificate with exhaustive maximum-square enumeration and **global** inradius coverings for the previous v245 and new v247 ten-gons. The numerical search proposes coordinates only. The mathematical result uses exact rational comparisons throughout.

## Main result

For a simple Jordan polygon \(P\), write

\[
J(P)=\frac{s_{\max}(P)}{r_{\mathrm{in}}(P)},
\]

where \(s_{\max}\) is the largest side of any boundary-inscribed nondegenerate square and \(r_{\mathrm{in}}\) is the *true globally largest* disk radius contained in the bounded polygonal region.

**Theorem.** The new rational ten-gon \(P_{247}\), recorded exactly in [v247-refined-ten-gon.json](v247-refined-ten-gon.json), has strictly smaller actual Big Peg ratio than the previously certified rational ten-gon \(P_{245}\):

\[
\boxed{J(P_{247})<J(P_{245}).}
\]

Numerical intervals below are outward-rounded displays, **not** inputs to a decision.

| Polygon | Rigorous interval for its true ratio |
|---|---|
| Existing v245 | \(1.4086952790854 < J(P_{245}) < 1.4086952826310\) |
| New v247 | \(1.4083545868226 < J(P_{247}) < 1.4083545868287\) |

Therefore the decrease in the *true* ratios is greater than **0.00034069225**. Since an actual example has \(J<1.40835458683\), any possible universal Big Peg lower-bound constant \(c\) must obey \(c<1.40835458683\). This is an *upper bound on a conjectural universal constant*, not evidence that no positive constant exists.

The new decagon descends from the exact ten-sided example introduced by **Grisha Pochuev** (see [his primary source](https://github.com/Grisha-Pochuev/sharp-big-peg-counterexample)). The new grid changes and certification are separate work and are **not** attributed to Pochuev.

## Exact rational candidate, derived from v245

The grid denominator is **1,200,000**. The following are all ten integer coordinate numerators in cyclic order:

| Vertex | X | Y |
|---:|---:|---:|
| 0 | 4114000 | -558000 |
| 1 | 2021390 | 1568001 |
| 2 | -1234100 | 1006000 |
| 3 | -868000 | 2516000 |
| 4 | -1089000 | 2813990 |
| 5 | -3767999 | 1167900 |
| 6 | -3459001 | 81100 |
| 7 | 140000 | -1348000 |
| 8 | 571334 | -1449110 |
| 9 | 1106401 | -485200 |

All coordinates are exact rational numbers \(X_i/1200000,Y_i/1200000\). The [canonical JSON data file](v247-refined-ten-gon.json) also stores each integer grid displacement from the *previous* v245 coordinates, such that \(P_{247,i}=P_{245,i}+\Delta_i/1200000\). This guards against any silent change to the candidate or its attribution.

### Search and reconstruction

The numerical search began from the already certified v245 polygon and evaluated neighboring rational-grid moves, ranking each by a full **floating-point** edge-assignment search and a candidate three-contact incenter approximation. There were successive grid resolutions \(1/12000\), \(1/120000\) and \(1/1200000\), with 24, 6 and 8 successful greedy coordinate steps respectively. The final exact grid displacement vector is preserved in the JSON.

The source [big_peg_search_v247.py](../code/big_peg_search_v247.py) is a reusable numeric-only neighborhood-ranker; it explicitly rejects its own output as a proof. Near-degenerate linear systems and changes in square incidence can fool floating-point ranking, so the final polygon is accepted **only** by a separately run Fraction-based certificate.

## First half of certificate: the largest square

The existing [v243](../code/big_peg_exact_v243.py) exact algorithm enumerates every \(2n^4\) ordered square-edge assignment/orientation chart. Every feasible chart is a rational compact polytope in two adjacent square-vertex edge parameters. The squared side is a convex quadratic, so its maximum lies at a polytope vertex, including in singular positive-dimensional families.

For \(n=10\), all **20,000** cases are exhausted. The exact attained maximal square side squared for the new ten-gon is

\[
\boxed{s_{247}^{2}=
\frac{407200306394224630512330297939647665150234737979752892405301}
{204751870350371073185004616046232753607402236547022400000000}}
\]

with decimal display \(s_{247}^{2}\approx1.988750118362407\). A maximizing edge assignment is \((1,4,5,6)\). The optional older [v223](../code/polygon_orbit_census_v223.py) four-parameter exact census independently verifies that the polygon has **exactly one geometric square orbit**, with precisely this maximal side length.

For comparison, the previous v245 exact squared side is

\[
s_{245}^{2}=
\frac{2263154862866037715637393248853633089}
{1137974336077532830058267809922880000}.
\]

Both maxima are recomputed, not loaded unverified from JSON.

## Second half: the true inradii

Let the new explicit **rational** disk center be

\[
C_{247}=\left(\frac{114334619}{100000000000},
-\frac{7507522}{100000000000}\right).
\]

The existing exact ray-crossing test proves that this center lies strictly inside the polygon. The minimum rational squared distance to every closed polygon edge equals

\[
\boxed{r_{\mathrm{lower},247}^{2}
=\frac{984883927824703700502571706090888809}
{982265423769090000000000000000000000}}.
\]

The entire closed disk of radius \(\sqrt{r_{\mathrm{lower},247}^{2}}\) is therefore contained in the polygon. This proves a **global inradius lower bound** without relying on a numerically estimated optimal center.

To prove a global **upper** radius bound, the existing [v246](../code/big_peg_inradius_v246.py) rational rectangle-cover algorithm is rerun on the *entire* new polygon's bounding rectangle with

\[
\boxed{r_{\mathrm{in}}(P_{247})\le
\frac{100133200316}{100000000000}=1.00133200316}.
\]

This cover classifies **all 1,257 rational subrectangles** (maximum subdivision depth at most 48). For every terminal box \(Q\), it proves one of:

1. **Boundary-close:** all points of \(Q\) have distance at most the claimed radius from the polygon boundary, using the 1-Lipschitz property of distance to a closed set.
2. **Exterior:** the midpoint is outside the polygon and the exact nearest-boundary distance exceeds a rational upper bound on the half-diagonal. The full box then contains no interior point.

If neither condition holds, the box is subdivided; the program *fails closed* if any box remains unclassified at its depth/node limit. No float decisions, assumed active edges, or unproved incenter formula enter this global certificate.

The new lower and upper inradius bounds differ by less than \(4.3\times10^{-12}\) in radius.

The old v245 polygon is checked afresh against its previously proved [v246](v246-true-inradius-global-cover.md) global radius upper bound

\[
r_{\mathrm{in}}(P_{245})\le\frac{100109166}{100000000}.
\]

The v245 rational cover processes **917 boxes**, also checked from scratch.

## Third half: an exact true-ratio comparison

From the complete square maxima and global radius bounds we obtain

\[
\frac{s_{247}^{2}}{
 (100133200316/100000000000)^{2}}
\le J(P_{247})^2\le
\frac{s_{247}^{2}}{r_{\mathrm{lower},247}^{2}},
\]

and

\[
J(P_{245})^2\ge
\frac{s_{245}^{2}}{(100109166/100000000)^{2}}.
\]

The decisive result is the fully **exact, strict rational inequality**

\[
\boxed{
\frac{s_{247}^{2}}{r_{\mathrm{lower},247}^{2}}
<
\frac{s_{245}^{2}}{(100109166/100000000)^{2}}
}.
\]

Because both ratios are nonnegative, it follows immediately that the actual \(J(P_{247})<J(P_{245})\). The inequalities do not merely compare convenient disks and are valid even if the true incenter is somewhere else.

## Reproduction

Pure Python standard library, from repository root:

~~~sh
python code/big_peg_true_ratio_v247.py
python code/big_peg_true_ratio_v247.py --crosscheck-v223
python -m unittest -v tests.test_big_peg_true_ratio_v247
~~~

Optional **numerical-only** diagnostics, requiring numpy:

~~~sh
python code/big_peg_search_v247.py
python code/big_peg_search_v247.py --start v245 --steps 4 --grid 12000
~~~

The exact CI job reruns both maxima, both global inradius coverings, all exact margin tests, and the v223 independent square census. A separate numerical job verifies the explorer can replay the candidate scores but does **not** contribute to the proof.

## Scope and next research questions

The result improves an **explicit upper bound** on any universal positive Big Peg constant and now also establishes a strict order of the true ratios of these two individual decagons. It does **not** establish the best possible decagon, any bound for all 10-, 12- or 16-gons, or the truth/falsity of the universal \(c>0\) assertion. Most importantly it does not prove or refute the unrestricted Toeplitz square-peg conjecture.

Next: build active-constraint and singular-chart-aware polygon optimization, expand the candidate space to 12 and 16 vertices, and require independent exact maximum-square and whole-polygon inradius audits for every new record.
