# v263: A complete all-real-time square-event ledger for the moving pentagon

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact-arithmetic, regression-certified **complete positive-square census at every real time** in the v254 pentagon's algebraic collision window, including the irrational collision itself. No time sampling is used for the exhaustion theorem. The unrestricted Toeplitz conjecture remains open.

## 1. Main theorem: all the squares, for all real times

Let \(P(t)\) be the explicitly rational-affine, strictly convex Jordan pentagon from [v254](v254-genuine-square-pair-death.md), and restrict to

\[
I=[a,b],\qquad
a=\frac{64722}{100000},\qquad b=\frac{64723}{100000}.
\]

Let \(\alpha\) denote the unique zero in \(I\) of

\[
Q(t)=72468080000000t^2
-496738335230000t
+291143052416293.
\]

The [new exact verifier](../code/pentagon_complete_ledger_v263.py) proves that **every positive-size inscribed geometric square** of \(P(t)\), at **every real time** \(t\in I\), belongs to exactly one of the three already known supporting-edge branch families:

\[
\begin{array}{lcl}
C=(0,1,2,4)&:&\text{persistent square, valid for all }t\in I,\\
A=(0,1,3,4)&:&\text{valid for }t\le\alpha,\\
B=(0,2,3,4)&:&\text{valid for }t\le\alpha.
\end{array}
\]

At \(t=\alpha\), branches \(A,B\) are **one and the same geometric square**, not two independent orbits. For \(t>\alpha\), their extended line-incidence solutions are off the actual polygon.

Consequently, the **complete whole-pentagon** count is

\[
\boxed{
N(t)=
\begin{cases}
3,& a\le t<\alpha,\\[2pt]
2,& t=\alpha,\\[2pt]
1,& \alpha<t\le b.
\end{cases}}
\]

This goes beyond v254, v260 and v262, which had all-real-time *selected branch* certificates but complete whole-polygon enumerations only at rational time slices. The new result is **exhaustive at every real time**, including possible untested algebraic instants, because every other edge assignment is excluded with a full-interval exact certificate.

Every actual inscribed square in this family has squared side **strictly greater than 1**. That stronger size statement follows from v254's exact proofs for each of the three surviving branch families, combined with the new exhaustive classification.

## 2. Why precisely 70 cyclic edge charts are complete

The polygon is strictly convex and counterclockwise for the entire original interval \([3/5,13/20]\), independently certified in v254. The four corners of a nondegenerate square on the boundary of a convex polygon occur in the same counterclockwise cyclic order along the polygon boundary as around the square. They are four distinct geometric points.

Choose the first corner to be the one with the smallest boundary parameter in \([0,5)\). Its ordered supporting-edge indices are then weakly increasing:

\[
0\le e_0\le e_1\le e_2\le e_3\le4.
\]

An actual polygon vertex belongs to both incident **closed** edges; choosing either valid edge label gives at least one such weakly increasing assignment. Therefore *including repeated edge labels and closed edge parameters* makes this enumeration complete.

There are precisely

\[
\binom{5+4-1}{4}=\boxed{70}
\]

multisets of four supporting edges. The enumeration includes all endpoint hits and repeated edges. It does not assume that all four corners touch four distinct open edges.

**Scope:** This simplification uses strict **convexity**. On a nonconvex Jordan polygon, inscribed square corners may appear in nongraceful boundary order; v218's six-sector problem then reappears. The 70-chart reduction is not a general simplification for all Jordan curves.

## 3. Complete classification of the 70 charts

| Class | Number | Why impossible / retained |
|---|---:|---|
| Three or more square corners assigned to one edge | **25** | Three distinct vertices of a nondegenerate Euclidean square cannot be collinear |
| Two *adjacent* supporting edges, with two corners each | **5** | Their nonparallel supporting lines cannot contain all four vertices of a positive square on the two incident closed segments |
| One strict infeasible edge parameter, for every real time | **37** | Full-interval exact rational Bernstein certificate |
| Three tracked square assignments | **3** | Full-interval v254 valid-branch certificates, including exact algebraic vertex gluing |
| **Total** | **70** | **No unresolved chart** |

The five two-adjacent-edge exclusions require a small geometric proof. If four distinct square vertices lie on just two lines, with exactly two on each, those pairs must be either (i) two opposite sides of the square, requiring **parallel** lines, or (ii) the two square diagonals, requiring perpendicular lines that intersect at the square center. Adjacent edges of the strictly convex pentagon have nonparallel supporting lines intersecting at their shared polygon vertex \(V\). In the diagonal case, the two points assigned to each *closed incident edge segment* lie on one ray from \(V\), so their midpoint cannot be \(V\). Thus they cannot form opposite vertices of a square centered at \(V\). A shared-vertex point assigned to both edges does not evade the argument: with four **distinct** corners one of the two edges would then have to carry three corners, already excluded.

The 37 off-edge exclusions are different: no numerical geometry threshold or sampled approximate root is used. Every corresponding square-defect incidence matrix has a fixed, strictly nonzero determinant sign across the full normalized time interval. Its unique extended solution has at least one parameter \(u_j\) for which one of

\[
\operatorname{sgn}(D)N_j(t),\qquad
\operatorname{sgn}(D)(D(t)-N_j(t))
\]

is **strictly negative at every real time**, as verified by rational Bernstein coefficients and rigorous subdivision. This excludes that candidate even at arbitrary irrational times and endpoint contacts.

Any chart not covered by one of the permitted certificates causes a hard failure. The code includes a separate **singular augmented-matrix inconsistency** verifier and explicitly classifies unknown determinant zeros or persistent singular charts as *unresolved* rather than pretending they contain no positive-size square. Neither singular fallback nor the positive-size cutoff is needed to discard any of the 37 off-edge charts in this particular example. Hence the result does not depend on an unproved generic-transversality assumption or on silently dropping singular charts.

The classification's exact machine-checked profile is

\[
\boxed{25+5+37+3=70.}
\]

## 4. Independent positive-size lower bound

The existing v254 whole-family Jordan certificate gives a strict global separation between every pair of nonadjacent polygon edges, throughout the larger time interval:

\[
\delta=\frac{3901}{10000}>0.
\]

A positive square on a simple polygon must have two corners on nonadjacent edges, a geometric fact audited in v250. Those corners are separated by at least \(\delta\), but no two vertices of a square of side \(s\) are farther apart than its diagonal \(\sqrt2\,s\). Therefore

\[
s^2\ge\delta^2/2
=\frac{15217801}{200000000}
>0.
\]

The verifier constructs a strictly smaller, independently safe threshold

\[
\eta=\delta^2/4
=\frac{15217801}{400000000}
\]

and verifies its exact equality to the v254 separation witness. This global geometric inequality is an **additional** positive-size guarantee and a safeguard for future cutoff-based chart exclusion. In v263's actual 70-chart classification, **zero charts are discarded as merely below \(\eta\)**: the other 67 charts are excluded by pure supporting-edge geometry or a strictly bad edge parameter.

Once the three tracked branches are proved exhaustive, v254 supplies the sharper family-specific statement \(s^2>1\) for **every** inscribed square of this pentagon over \(I\).

## 5. Reconciliation with the v260–v262 degree ledger

There are no hidden positive-square births, rank-drop events, repeated-edge square configurations, or untracked vertex events inside \(I\). The only actual geometric orbit event is the unique algebraic meeting and disappearance of \(A,B\) at \(t=\alpha\). This assertion follows from complete chart exclusion, not from a handful of sample counts.

| Time | Whole-pentagon geometric square orbits | Local index data on the oriented positive cyclic-order lift | Unlabeled quotient relative mod-2 degree |
|---|---:|---|---:|
| \(t<\alpha\) | 3 | \(-1,+1,+1\) | 1 |
| \(t=\alpha\) | 2 | Pair contributes **degree 0** at its nonsmooth interface, persistent root contributes \(+1\) | 1 |
| \(t>\alpha\) | 1 | Only the persistent \(+1\) root | 1 |

This now validates a **whole-zero-set interpretation of the selected v262 quotient atlas for this particular pentagon family and exact short interval**. It encompasses every positive square orbit, rather than only the three previously *known* branches. However, v262's quotient degree is still **relative to two explicitly selected zero-free-boundary neighborhoods**. This does not construct a canonical nonzero *absolute* Euler class on the full configuration-space quotient. v219 already established why the naive absolute obstruction vanishes.

At the singular event there are two distinct geometric square roots, an even raw count, but the degree remains odd. The shared nonsmooth pair root has mod-2 local degree zero. Counting singular roots without their local multiplicities would misstate the invariant.

## 6. Reproducibility and adversarial checks

The new standard-library-only program reads the exact v254 moving-polygon family and inherited v224/v249 polynomial Bernstein machinery. Its CI requires a **strict** all-real-time classification; an unresolved chart fails the test run, and the diagnostic never silently substitutes an endpoint census for proof.

```sh
python code/pentagon_complete_ledger_v263.py --diagnostic
python code/pentagon_complete_ledger_v263.py
python -m unittest -v tests.test_pentagon_complete_ledger_v263
```

The new [regression suite](../tests/test_pentagon_complete_ledger_v263.py) independently rebuilds and reconfirms each of the **37 strictly negative Bernstein exclusion polynomials**, verifies all 25 collinear and five adjacent-edge geometry exclusions, replays the exact irrational collision theorem, checks the global side-size gap and the independently complete v223 rational endpoint censuses, and forces failures for noncanonical geometry, altered time range, uncertified cutoffs and deliberately disabled sign proofs.

## 7. What remains open

The result is **global only within this explicit convex pentagon family and exact collision window**. It is not new general square existence for all convex polygons, a general Jordan-curve parity theorem, or a proof of Toeplitz. It does establish that the local-to-quotient degree construction and a complete positive-square geometric census can agree through a *genuine nonsmooth vertex-induced square-pair death*, without assuming missing squares away from the chosen local branches.

The next meaningful directions are:

1. **v264: reusable whole-family event ledger.** Generalize beyond this hand-selected pentagon and time window, removing reliance on hardcoded expected chart counts. Produce complete rational/algebraic event partitions for additional moving polygons, with strict rank-singular and repeated-edge handling.
2. **Multisector/equivariant topology.** Audit what survives when positive-square orbits are not gracefully ordered on the boundary, as in arbitrary nonconvex Jordan curves.
3. **Collision-size control under arbitrary \(C^0\) limits.** The central Toeplitz bottleneck remains uniform noncollapse: for approximating Jordan curves with arbitrarily fine features, the polygon-specific edge-separation \(\delta\) can tend to zero. v263's \(s^2>1\) is powerful but family-specific and does not remain uniform for arbitrary Jordan limits.

**The unrestricted Toeplitz square-peg conjecture remains open.**

## References

- [v223: exhaustive rational polygon geometric-orbit census](v223-complete-polygon-orbit-census-and-homotopy-audit.md).
- [v224: exact moving-polygon polynomial event candidates](v224-algebraic-moving-polygon-event-covers.md).
- [v250: positive square size separation on polygons](v250-noncollinear-corner-noncollapse.md).
- [v254: genuine positive-size pair annihilation](v254-genuine-square-pair-death.md).
- [v260: signed pair cancellation](v260-relative-degree-pair-death.md).
- [v261: oriented atlas with the surviving third square](v261-oriented-multi-region-square-atlas.md).
- [v262: cyclic quotient relative mod-2 degree](v262-cyclic-quotient-relative-parity.md).
- B. Matschke, *A Survey on the Square Peg Problem*, *Notices of the AMS* 61 (2014).
