# v264: A reusable exact moving-convex-polygon square-event atlas, stable under certified edge refinement

**Research program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Scientific status:** Exact rational algorithm and checked finite geometric theorems for affine motions of strictly convex rational polygons, including representations with permanently collinear edge subdivisions. One full moving pentagon, one independent moving quadrilateral, and positive affine similarity families are tested. A general full-time result is returned only when the candidate event set is empty; irrational event times are separately marked unresolved. **The unrestricted Toeplitz conjecture remains open.**

## 1. The new reusable algorithm

Versions [v260](v260-relative-degree-pair-death.md) through [v263](v263-complete-all-real-time-moving-pentagon-census.md) completed a rigorous 3-to-2-to-1 square-orbit ledger for one chosen moving pentagon. Crucially, v263's exhaustive classifier was tied to its five exact vertices and its hardcoded expected classes.

The new [implementation](../code/convex_square_event_atlas_v264.py) accepts **arbitrary rational endpoint vertex sequences** \(P_0,P_1\), with a common vertex count \(n\ge4\), and studies the entire linear family

\[
P_i(t)=(1-t)(P_0)_i+t(P_1)_i,\qquad 0\le t\le1.
\]

The algorithm has five independent gates, each failing closed when the stated hypotheses cannot be certified:

1. **Geometric path:** Validate the endpoints and prove strict convexity and consistent counterclockwise order for **every real** \(t\). Remove only provably harmless straight-edge subdivisions.
2. **Uniform noncollapse:** Construct explicit rational projection-separation witnesses for all nonadjacent edge pairs, uniformly over the time interval.
3. **Complete chart coverage:** Enumerate every weakly increasing assignment of four corners to the closed polygon edges, including repeated edges and vertex contacts.
4. **Exact algebraic classification:** Exclude charts with fully proved geometric or strict polynomial obstructions. Otherwise certify a stable square chart, or isolate every possible determinant/contact/order/size event.
5. **Geometric orbit census:** For every event-free time chamber, use v223's complete rational-polygon solver, deduplicating geometric square orbits. For rational event times, census the event slice itself. For *irrational* times, leave the event-slice count unresolved unless some other exact theorem handles that singularity.

It never asserts a global index merely from an odd count in one sampled polygon, or an all-real-time count merely from counts at rational endpoints.

## 2. Strict convexity and certified removal of artificial edge knots

An artificial new polygon vertex inserted into a straight edge creates consecutive collinear edge trajectories. Ordinary "no three collinear vertices" isotopy tests reject such a refinement even though it changes no geometric curve. A raw edge-index separation bound can also become arbitrarily small as subdivisions accumulate.

v264 resolves this for **permanently collinear, forward-pointing** subdivisions. Let \(A(t),B(t),C(t)\) be three consecutive moving vertices and let

\[
e_0(t)=B(t)-A(t),\qquad e_1(t)=C(t)-B(t).
\]

The middle vertex may be removed **only if**

\[
\operatorname{cross}(e_0(t),e_1(t))\equiv0
\quad\text{as a rational polynomial,}
\]
\[
\operatorname{dot}(e_0(t),e_1(t))>0
\quad\text{at every real }t\in[0,1].
\]

The strict dot-product certificate, checked by exact rational Bernstein coefficients or an equivalent Sturm root-free sign test, proves that \(B(t)\) is strictly between \(A(t)\) and \(C(t)\) for every real time. Removing that point therefore changes only the piecewise-linear **parameterization**, not the polygonal image or its geometric square set.

After all certified removals, the implementation checks

\[
\operatorname{cross}\bigl(P_j(t)-P_i(t),P_k(t)-P_i(t)\bigr)>0
\quad (i<j<k)
\]

for **all** \(\binom n3\) ordered vertex triples and all real \(t\). This is a sufficient exact certificate for a strictly convex, counterclockwise, non-self-intersecting polygon at every time. A nonconvex motion or a merely possibly convex family is rejected.

The 10-edge version of the moving v254 pentagon, formed by inserting a moving midpoint on each of the original five edges, is proved to compress to the **same exact five moving vertices**. The same holds with insertion at one-third of each edge. A tiny genuine bend is *not* silently compressed.

**Limit:** This handles exact persistent straight-edge subdivisions. A refinement whose new vertices follow only approximate straightness, whose collinearity is transient, or whose motion introduces a nonconvexity is not covered by this compression theorem.

## 3. A quantitative positive-square size floor, constructed from the polygon

Let \(E_i(t)\) and \(E_j(t)\) be nonadjacent closed polygon edges. A uniform projection certificate on a rational time slab \(J=[u,v]\) consists of a fixed rational nonzero direction \(w\in\mathbb Q^2\) such that

\[
\min_{x\in E_i(t)}w\cdot x-
\max_{y\in E_j(t)}w\cdot y\ \ge g>0
\]

for every \(t\in J\), or the same strict inequality with \(i,j\) exchanged. Because all endpoint trajectories are affine, the exact minimum and maximum of the linear projection on \((t,\text{edge parameter})\in J\times[0,1]\) occur at four rational vertices. The inequalities are checked **exactly**, not numerically.

Cauchy–Schwarz gives

\[
\|x-y\|_2^2\ge\frac{g^2}{\|w\|_2^2}
\]

for every corresponding pair of points. If one rational direction fails, the code tries additional rational directions (edge tangents, normals, endpoint differences, center difference) and, if necessary, subdivides the time interval into dyadic slabs. A failure to find sufficient witnesses within the declared finite budget yields **INCONCLUSIVE**, not a fabricated distance.

Taking the **minimum** of the finitely many rational lower bounds over **all nonadjacent edge pairs and time slabs** gives \(d^2>0\).

For a strictly convex polygon with at least four edges, every positive square has corners on some **two nonadjacent polygon edges**. Otherwise all four corners would have to lie on two adjacent edges, which is impossible for a positive Euclidean square. Their distance is at most the square's diagonal \(\sqrt2\,s\), hence

\[
\boxed{s^2\ge d^2/2>0.}
\]

The algorithm chooses a conservative cutoff \(\eta=d^2/4\), so no actual positive square is ever removed by a squared-side cutoff of \(\eta\).

**On the v254 pentagon window**, the independent v264 rational projection proof gives

\[
\boxed{d^2\ge51780781/100000000},
\quad
\boxed{s^2\ge51780781/200000000}.
\]

All five nonadjacent edge pairs are certified across the *entire* time interval without needing any dyadic subdivision. This quantitative witness is **stronger** than the earlier v263 axis-aligned rectangle certificate, which was conservative in a different way. The latter gave \(s^2\ge15217801/200000000\).

These bounds are statements about a chosen **fixed family**. The positive \(d\) may depend on its precise geometry and is not a universal bound for all Jordan curves. Certified straight-edge subdivision preserves the underlying geometric family and, after compression, gives exactly the **same** \(d^2\) bound.

## 4. Full positive-square chart coverage, including repeated edges

Every positive square inscribed in a convex counterclockwise polygon is gracefully ordered along its boundary. Number the corners in that order, beginning at the smallest positive boundary parameter. Its support edge tuple is therefore weakly increasing:

\[
0\le e_0\le e_1\le e_2\le e_3<n.
\]

There are \(\binom{n+3}{4}\) such quadruples, including repeated edges. A polygon vertex can belong to either incident closed edge and is therefore included without assuming genericity.

For each potential square the exact v224 four-dimensional affine incidence system gives

\[
F_t(u)=M(t)u+c(t)=0,\qquad
D(t)=\det M(t),\qquad N_j(t)=D(t)u_j(t).
\]

When \(D\ne0\), the contact conditions are \(0\le N_j/D\le1\), so the signed polynomial inequalities

\[
\operatorname{sgn}(D)N_j\ge0,\qquad
\operatorname{sgn}(D)(D-N_j)\ge0
\]

are necessary and sufficient for the root to lie on all four **closed edge segments**.

A repeated edge receives an additional constraint. If \(e_j=e_{j+1}\), the two corners must appear in increasing order along that edge:

\[
\operatorname{sgn}(D)(N_{j+1}-N_j)>0
\]

for a positive-size square. Without this check, weakly increasing *edge labels alone* would not ensure correct curve order.

The exact chart classifications are:

- **Geometrically excluded:** three assigned square corners on one line, or all four on two adjacent incident supporting edges.
- **Uniformly off-edge:** fixed nonzero determinant sign and an edge-contact/order polynomial strictly negative throughout the entire time interval.
- **Uniformly too small:** the algebraic candidate's squared side is below the independently certified positive-square cutoff; this cannot remove a genuine square.
- **Identically singular but inconsistent:** polynomial augmented-matrix minors prove \(\operatorname{rank}[M\mid -c]>\operatorname{rank}M\) for **every real time**.
- **Stable interior or vertex-locked chart:** fixed nonzero determinant, positive size, and all nonzero contact/order factors strictly signed; permanently zero vertex-contact factors are treated as identities, not ignored.
- **Potential event chart:** isolate real zeros of determinant, contact, repeated-order and size polynomials using exact rational Sturm intervals.
- **Unresolved:** if the preceding alternatives do not prove emptiness, stable occupancy, or an adequately isolated finite event set, the full-time request fails closed.

## 5. What an algebraic event bracket does and does not prove

Every change in a nonsingular square chart's support, side cutoff, or feasibility must happen at a zero of one of the enumerated rational polynomials. Determinant rank changes are also included. The verifier isolates all roots in \([0,1]\) as exact rationals or rational brackets enclosing the unique irrational root of a polynomial.

Brackets are conservatively merged when they overlap. This produces an exhaustive **cover of all possible exceptional event times**, not a claim that every candidate is a real geometric square birth or death.

On an event-free rational open chamber, the complete v223 polygon square census at one rational sample has the same *geometric orbit count* at every real time in that chamber, because no admissible chart can acquire or lose a positive square except through the covered event factors. The quantitative side floor excludes escape into zero-size configurations. The separate v223 census removes the duplication of square labels at vertices and supports.

For an **exact rational** event time, the polygon itself is rational and the exact v223 census determines the event-slice count. For an **irrational** time, that count is **not determined by the nearby rational samples**. This is not an incidental limitation: the v254 pair annihilation has **two** distinct geometric roots at the singular instant, even though adjacent chambers have three and one.

When there are no possible event roots anywhere in the full closed time interval, and no unresolved charts, the algorithm proves one **complete all-real-time** count and compares it against exact independent endpoint censuses.

For regular square configurations not on any polygon vertex, each rational chamber additionally reports the signed local Brouwer index sum in the fixed increasing source-parameter convention. Its **mod-2 orbit parity** is the number of geometric square orbits modulo two. These are **not** asserted to define a global absolute Euler class on the full quotient, and singular event slices do not inherit the raw root-count parity.

## 6. New independent examples

### A. Genuine moving pentagon, independently rediscovered

Without any hardcoded list of the expected branches or collision polynomial, the generic engine classifies all **70** charts as

| Category | Charts |
|---|---:|
| Three collinear corner assignments | 25 |
| Two adjacent support-edge assignments | 5 |
| Uniformly off-edge | 37 |
| Stable interior square | 1 |
| Potential event charts | 2 |
| **Total** | **70** |

The only candidate root factors are the second corner leaving edge 1 and entering edge 2:

\[
(0,1,3,4):\quad \text{corner-1-end},
\qquad
(0,2,3,4):\quad \text{corner-1-start}.
\]

They coalesce into **one** irrational candidate-event bracket in normalized time:

\[
\boxed{7991/65536<\tau<31965/262144}.
\]

The exact v254 quadratic root lies strictly within this bracket. Independent complete rational chamber censuses give

\[
N_\text{before}=3,\qquad
N_\text{after}=1,
\]

and both regular chambers have signed index sum \(+1\) and unlabeled orbit parity \(1\).

**Important:** v264 alone reports the irrational event-slice count as **unknown**. The earlier specialized v254/v263 proofs independently supply the missing fact \(N(\alpha)=2\) and prove that no other square exists. Combining those proven results gives the full 3-to-2-to-1 count, but the generic engine does not smuggle the v263 collision theorem into a universal event resolver.

### B. Translation and positive similarity, no event

An arbitrary rational translation of the pre-collision pentagon produces an **event-free, full-real-time** complete census of exactly **three** square orbits, odd parity and signed index sum \(+1\). A positive rational homothety of the post-collision pentagon gives exactly **one** square orbit for every real time, also odd parity and signed sum \(+1\).

These are genuine moving affine rational polygon families, distinct from v254's collision motion. The code reaches their full-time result from generic chart and gap certificates rather than from a named square branch list.

### C. An independent moving convex quadrilateral

A translating rational quadrilateral, with initial vertices

\[
(0,0),\quad(4,0),\quad(3,2),\quad(-1,3),
\]

is checked by exactly \(\binom{7}{4}=35\) complete edge charts. Its all-real-time count is compared directly to the independent v223 geometric-square census. The test imposes no pentagon-specific geometry or expected branch indices.

### D. Exact refinement invariance

Adding moving mid-edge vertices to all five original pentagon edges gives a ten-edge representation of the *same* moving curve. The new code removes exactly five permanently collinear forward-pointing vertices and certifies that all of the following match the unreﬁned polygon exactly:

- the same five coarse moving vertices;
- the same rational nonadjacent-edge separation and square side bound;
- the same complete 70-chart classification;
- the same event bracket, samples and geometric orbit counts.

One-third subdivisions reproduce the same exact compression. A perturbed bent vertex is not removed. This proves **geometric reparameterization invariance under the supported class of straight-edge refinements**, not an unlimited subdivision theorem for arbitrary rough Jordan approximations.

## 7. Verification

The dedicated v264 workflow runs the generic moving-pentagon example, **21 new exact tests**, and the retained **v262/v263 regression suites**. No floating-point decisions or root guesses are used for any certificate.

```sh
python code/convex_square_event_atlas_v264.py
python -m unittest -v tests.test_convex_square_event_atlas_v264
python -m unittest -v tests.test_pentagon_complete_ledger_v263 tests.test_square_quotient_parity_v262
```

The tests include a real irrational event check against the independently known exact v254 quadratic, general rational quadrilateral and affine similarity examples, the 10-to-5 edge compression equivalence, positive-size floors, singular and nonconvex fail-closed controls, and oriented chamber indices.

## 8. The next mathematical bottleneck

v264 removes the *hardcoded polygon* obstruction for one substantial class of affine convex motions. It also separates a false size-loss under **relabeling/refinement** from a genuine geometric loss of square size. But the Toeplitz problem requires a statement for arbitrary continuous Jordan curves, not just a fixed piecewise-linear convex family with finitely many edges.

The next high-value milestones are:

1. **Event-slice algebraic resolution:** build a general exact gluing verifier for irrational vertex-contact events. It should reason from Cramer polynomial identities and local cell orientations, not from hardcoded v254 data.
2. **Genuine approximation stability:** investigate whether any protected parity/relative-degree construction survives arbitrary \(C^0\) limits with a nonzero size floor **independent of the polygon mesh**. Merely subdividing a fixed polygon does not solve this.
3. **Nonconvex and nongraceful sectors:** extend source-chart coverage beyond the positive cyclic-order component without assuming convexity, retaining all six labeled cyclic-order sectors and the correct cyclic quotient/bundle geometry.

**No unrestricted Toeplitz solution is claimed.**

## Research lineage

- [v218: cyclic order and orientation obstruction](v218-global-configuration-symmetry-and-collision-audit.md).
- [v219: equivariant quotient and absolute Euler-class limitation](v219-equivariant-quotient-collision-link.md).
- [v223: exact geometric square-orbit census](v223-complete-polygon-orbit-census-and-homotopy-audit.md).
- [v224: chart polynomials and Sturm root isolation](v224-algebraic-moving-polygon-event-covers.md).
- [v254: genuine irrational moving-pentagon square-pair death](v254-genuine-square-pair-death.md).
- [v260: degree at singular kink](v260-relative-degree-pair-death.md).
- [v262: selected relative quotient-parity atlas](v262-cyclic-quotient-relative-parity.md).
- [v263: complete all-time family-specific square census](v263-complete-all-real-time-moving-pentagon-census.md).
