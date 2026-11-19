# v252: A COMPLETE size-filtered geometric square-orbit event ledger

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Scope:** A fully enumerated, exact-rational square-orbit ledger for one **genuinely changing, strictly convex nine-vertex Jordan polygon family**, restricted to the interval \(\lambda\in[1/2,501/1000]\) and the *strict positive squared-side cutoff* \(\eta=94257/10000\). It proves that the number of geometric squares **above this cutoff changes 1→2→1**, with unique algebraic entry and exit events. This is a concrete negative result for naive *size-filtered parity invariance*, not a disproof of unfiltered parity, positive universal Big Peg, or Toeplitz.

## 1. Why this is a step beyond v251

v251 selected two favorable square branches \(A\) and \(B\), proving that together they furnish some square with squared side above \(9.4257\) at every real time. However, it did **not** rule out *other, unobserved* large squares and so could not give a complete event count.

v252 does that exclusion. Every potential square with squared side strictly greater than the cutoff is considered, including exact polygon vertices, with no floating-point sampling of time or polygon shapes. Every general theorem invoked has a checked hypothesis, and ambiguous cases **fail closed**.

The family is the previously audited noncollinear moving-kink nine-gon of [v250](v250-noncollinear-corner-noncollapse.md); it is strictly convex and Jordan at **every** parameter. The exact upper-bound disk radius 3 is retained. We use the actual polygon, not a similarity scaling or collinear subdivision.

## 2. Why only 126 edge choices are sufficient

A square of side \(s\) has four distinct corners; the distance between **any two** of them is either \(s\) or \(\sqrt2s\), hence is at least \(s\).

For the whole time interval, exact Bernstein inequalities show every moving polygon edge has squared length **strictly below** \(\eta\). Therefore a square with \(s^2>\eta\) **cannot have two of its corners on the same closed polygon edge**.

Moreover, because each polygon is strictly convex, its four boundary-inscribed square corners must occur in the same counterclockwise cyclic order as the boundary. Orient the square counterclockwise, and begin at the square corner lying on the lowest-index polygon edge. Its four **distinct** supporting edge indices satisfy

\[
0\le i_0<i_1<i_2<i_3<9.
\]

Consequently, the exhaustive search has exactly

\[
\binom94=\boxed{126}
\]

ordered distinct-edge possibilities. Edge-endpoint duplications may lead to two *labels* for one square at an exceptional time, so an independent **geometric-orbit** census is additionally performed at the left endpoint.

## 3. Four exhaustive chart classifications

For each candidate edge quadruple, construct the exact time-dependent square-incidence matrix \(M(\lambda)\), Cramer denominator \(D\) and four coordinate numerators \(N_j\), as in v249. The resulting square, if \(D\ne0\), has rational parameter functions \(z_j=N_j/D\).

The exact classifier uses the following mutually exclusive tests, in this order.

| Exclusion or accepted class | Number of edge assignments | Exact justification |
|---|---:|---|
| Geometry | **87** | A pair of selected edges has global time-uniform rectangle bound \(\max\|x-y\|^2<\eta\); two square corners cannot both lie there |
| Outside a segment | **36** | With \(D\) of strictly certified sign, one of eight signed closed-segment inequalities is strictly false for **every \(\lambda>1/2\)** |
| Below size cutoff | **1** | Exact rational polynomial \(N_u^2+N_v^2-\eta D^2<0\) for **every \(\lambda>1/2\)** |
| Surviving branches | **2** | Only \(A=(0,3,5,7)\) and \(B=(1,3,5,8)\) remain; both have exact Cramer square roots, valid on their four edge segments throughout the window |

The total is **87+36+1+2=126**.

### Exact geometric exclusions

Since every family vertex depends affinely on \(\lambda\), each moving edge is contained, for every \(\lambda\in[1/2,501/1000]\), in a rational axis-aligned rectangle determined solely by the relevant vertex coordinate values at the two rational time endpoints.

The maximum possible squared distance between any point on rectangle \(R_i\) and any point on rectangle \(R_j\) is bounded above by

\[
(\max |x_i-x_j|)^2+(\max|y_i-y_j|)^2,
\]

where the maxima use their rational endpoint ranges. If this upper bound is **strictly below** \(\eta\), a square of side squared \(>\eta\) cannot have two corners on that pair of edges. This avoids even computing a singular Cramer determinant for 87 combinations.

### Polynomial exclusions cover *every real time*

For the remaining edge choices, the exact Bernstein coefficient test verifies \(D\) has constant **strictly nonzero** sign on the entire closed interval.

For an assigned square corner \(C_k=N^{\rm corner}_k/D\) on the segment \(A_i+tE_i\), the signed projected coordinate must satisfy two inequalities

\[
0\le \operatorname{sgn}(D)\,T_k,
\qquad
0\le \operatorname{sgn}(D)\,(D\|E_i\|^2-T_k),
\]

with \(T_k=E_i\cdot(N^{\rm corner}_k-D A_i)\). They are rational *polynomials*.

If the Bernstein coefficients of the **negative** of one required inequality are all nonnegative and the last coefficient is strictly positive, then the inequality is **strictly violated for every \(\lambda\in(1/2,501/1000]\)**. This works even if the failed inequality becomes equality at the left endpoint.

Exactly the same direct Bernstein argument excludes the single below-cutoff chart. **No possible chart is silently discarded because it is singular or numerically inconvenient.** Any unclassifiable survivor outside the two named branches stops the audit with an exception.

## 4. Independent exact endpoint orbit censuses

At \(\lambda=1/2\), square \(A\) places one corner precisely at the moving polygon vertex. The same geometric square is represented by **two different supporting-edge charts**, so chart count is not geometric orbit count at this instant.

We rerun the pre-existing **v223 complete geometric square census**, with all edge assignments, repeated-edge cases, both orientations and singular-chart handling. It deduplicates genuine geometric squares by exact corner coordinates and finds:

- At \(\lambda=1/2\): **exactly one geometric square above \(\eta\)**, namely \(A\), including the special vertex-on-edge labeling.
- At \(\lambda=501/1000\): **exactly one geometric square above \(\eta\)**, namely \(B\).

Their exact geometric corner coordinates are independently compared with the corresponding Cramer branch roots. The v223 audit is a second algorithm, not reuse of the 126-chart filter.

For every **interior** time the 126-classification proves the only possible high squares are \(A,B\). An exact sign of the difference of their first-corner coordinates proves they are geometrically **distinct for the whole closed window**. Together these checks make the cutoff-filtered census *complete*.

## 5. Two uniquely isolated size-threshold events

For each of the two surviving branches, its squared side minus \(\eta\) has the sign of

\[
F_E(\lambda)=N_{E,u}(\lambda)^2+
N_{E,v}(\lambda)^2-\eta D_E(\lambda)^2,
\]

because \(D_E^2>0\).

The exact derivative polynomials satisfy

\[
F_B'(\lambda)>0,\qquad F_A'(\lambda)<0
\]

**at every real time** in the observation interval, verified by direct rational Bernstein positivity. Endpoints have strict opposing signs, so each has exactly one threshold root by the intermediate-value theorem and strict monotonicity.

The exact rational isolating brackets are:

| Event | Certified unique root interval | Meaning |
|---|---|---|
| \(B\) enters | \((50053/100000,\;50054/100000)\) | \(B\)'s side squared crosses **upward** through \(\eta\) |
| \(A\) exits | \((500551/1000000,\;50056/100000)\) | \(A\)'s side squared crosses **downward** through \(\eta\) |

These event times are algebraic and are **not** replaced by rounded rational guesses. The brackets are disjoint. v251's independently verified **equal-size crossover** between \(A,B\) occurs strictly **between** the size-entry and size-exit events and, by itself, does *not* alter the number of squares above the cutoff.

### The full event ledger

Let \(\lambda_B\) and \(\lambda_A\) denote the unique entry and exit times. We have

\[
\frac12<\lambda_B<\lambda_A<\frac{501}{1000}.
\]

For the strict filter \(s^2>\eta\):

| Real time | Geometric square orbits above cutoff | Parity |
|---|---:|---|
| \([1/2,\lambda_B]\) | 1: \(A\) | Odd |
| \((\lambda_B,\lambda_A)\) | 2: \(A,B\) | Even |
| \([\lambda_A,501/1000]\) | 1: \(B\) | Odd |

At a threshold equality, that square is **not counted** because the filter is strict. No actual square necessarily vanishes at these times: the *membership in the chosen size-filtered set* changes.

\[
\boxed{\text{Filtered orbit count: }1\longrightarrow2\longrightarrow1.}
\]

**Consequence:** The parity of the raw number of squares **above an arbitrary positive cutoff** is **not an invariant under polygon deformation**, even for this strictly convex Jordan family. A proposed size-aware topological degree must account for contributions from crossings of its size-boundary hypersurface; one cannot infer constancy of the filtered parity from the ordinary generic parity arguments.

This is not a counterexample to the usual unfiltered odd-parity picture: the complete census here covers only squares above the specified \(\eta\), not all positive squares.

## 6. Reproducibility and scientific limitations

~~~sh
python code/big_peg_cutoff_ledger_v252.py
python -m unittest -v tests.test_big_peg_cutoff_ledger_v252
~~~

The dedicated CI also replays v251 and v250 and executes the existing repository-wide research checks. The 20 new tests check the exact 126-way partition, independent geometric endpoint counts, a third interior census with two distinct high squares, event uniqueness/ordering, rational event brackets, parity and fail-closed input/geometry/threshold controls.

**Trust base:** Exact rational arithmetic, finite combinatorics, convexity of polygonal Jordan curves, Cramer's rule, polynomial Bernstein-basis positivity and strict monotonicity, plus v223's independently implemented complete geometric orbit census at rational endpoint times. No numerical sampling or assumed approximate root locations enter the formal conclusions.

**Not proved:** A complete ledger for square orbits **below** this cutoff; arbitrary moving polygons; generically singular birth/death of unfiltered square orbits; a protected mod-two degree for the unrestricted square configuration space; any universal positive Big Peg constant; or the Toeplitz square-peg conjecture.

### Related mathematics

- Benjamin Matschke, [*A Survey on the Square Peg Problem*](https://www.ams.org/notices/201404/rnoti-p346.pdf), *Notices AMS* (2014), DOI 10.1090/noti1100. Classical parity/continuation arguments pertain to an appropriate full configuration space, not to a naive arbitrary strict side-length cutoff.
- [v223 independent complete orbit census](../code/polygon_orbit_census_v223.py) and [v251 size-aware witness atlas](v251-size-aware-multibranch-atlas.md), both in this repository.
- The original sharp Big Peg counterexample and its independently verified descendants are orthogonal to the present size-filtered parity phenomenon.

## Next meaningful step

Develop a **relative, size-boundary-corrected degree or mod-two event balance** for complete square-orbit configuration spaces, including singular edge charts. The v252 ledger provides an exact check that any such proposed invariant must assign **one boundary contribution at each upward/downward cutoff crossing**, otherwise it is immediately false.

A mechanism that then bounds those boundary contributions away from collapsing squares, uniformly over polygon complexity, is still missing. Without that uniform control, no universal Big Peg theorem or Toeplitz proof follows.
