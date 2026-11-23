# v213 — Exact PL periodic square certificates, stability, and the lift-coherence obstruction

**Research program:** Benjamin John Schulz  
**Prepared:** 18 November 2026  
**Status:** Elementary finite-model theorems with exact-rational executable tests. Not a universal Toeplitz proof, not a Floer-homology computation, and not a formal proof-assistant verification.

## Executive finding

A balanced four-strand periodic inscription theorem cannot be transferred to one Jordan curve merely by (i) proving a positive square-size floor and (ii) identifying four original subarcs in the *quotient* cylinder. There are **two independent localization obligations**:

1. **Tag confinement:** all four vertices must lie on *the original subarcs*, not on completion connectors (v212 provided a robust counterexample to a generic inference).
2. **Deck-lift coherence:** their respective integer translations in the universal covering plane must be the four *physical lifts* of those subarcs on the actual original curve. A square using other translates is not necessarily a square on that curve, even if every vertex projects to an original subarc.

This batch implements a complete finite exact decision method for PL periodic graphs to distinguish both obligations. It also gives a simple stable square-exclusion certificate and a complementary transversality-based persistence certificate.

## Theorem 1: finite exact decision on periodic PL graph pieces

Let `f_i : R/Z → R` for `i=0,1,2,3` be continuous piecewise-linear functions with **rational vertices**. Let each `W_i ⊂ [0,1]` be a closed interval, possibly a singleton, interpreted as a chosen subarc of the graph. Consider ordered Euclidean squares with one point `p_i` on the lift of the `i`th chosen graph arc, with either cyclic orientation. Each periodic strand admits arbitrary integer horizontal deck translations unless a fixed four-tuple of relative translations is imposed.

**Claim (complete finite decision):** There is a finite algorithm using rational linear algebra which either produces a nondegenerate labeled square with exact rational coordinates or proves that none exists in the chosen arcs, for all possible deck translates. In particular, the algorithm handles squares on segment endpoints and positive-dimensional families of squares, not just isolated transverse solutions.

### Proof of finiteness of the deck search

Write a square with first edge `e=(a,b)` and cyclic sign `σ=±1` as

```
p0,
p1 = p0 + e,
p2 = p0 + e + σ R90(e),
p3 = p0 + σ R90(e).
```

Its horizontal and vertical bounding-box lengths are **both** `|a|+|b|`. Therefore, if `H = max_{i,x} f_i(x) − min_{i,x} f_i(x)`, any such square satisfies `|x_i−x_0|≤H` for all four vertex coordinates in the covering plane. Translate the entire square by one integer horizontally so `x_0∈[0,1]`, and express each other x-coordinate as `t_i+k_i`, with `t_i∈[0,1]` and `k_i∈Z`. Then `|k_i|≤H+1`. It is sufficient to examine `k_0=0` and each `k_i` in `[-K,K]`, where `K=ceil(H)+1`. If `H=0`, no nondegenerate square can exist.

This is a genuine finite bound, not an unverified numerical cutoff. Fixing `relative_lifts=(0,k1,k2,k3)` further confines the test to specified physical lifts.

### Proof of exactness on four fixed segments and lifts

Represent the four graph segments as `p_i(u_i)=a_i+u_i d_i`, with `0≤u_i≤1`, including the chosen horizontal deck translations in the `a_i`. The ordered-square equations are

```
p2 − p1 − σ R90(p1−p0) = 0,
p3 − p0 − σ R90(p1−p0) = 0.
```

Their four scalar components form an **affine rational linear system** `Au=b` in four unknown segment parameters. Consequently, the feasible set `P={u∈[0,1]^4: Au=b}` is a compact rational polytope. If nonempty, it has at least one vertex. At a vertex, in addition to `rank(A)=r` independent equations, at least `4−r` independent box constraints `u_j∈{0,1}` are active. Enumerating every such combination, solving each determined `4×4` rational system, and checking the original equations finds every possible vertex of `P`.

A nondegenerate square exists on these segments **iff some vertex of `P` is nondegenerate**: the first-edge vector `p1(u1)-p0(u0)` is affine in `u`, and if it vanishes at every polytope vertex, it vanishes on their convex hull, i.e. on all of `P`. So discarding zero-side witnesses is complete even for rank-deficient families.

A finite product over the original PL segments, two orientations, and the bounded deck shifts gives the decision procedure. Necessary interval identities (`p0+p2=p1+p3`, `x0−x2=σ(y1−y3)` and `y0−y2=σ(x3−x1)`) safely prune cases but never justify a positive result by themselves.

**Important limitations:** The theorem addresses exact finite PL graphs, not arbitrary Jordan curves and not general non-graph arcs. Rational arithmetic avoids floating-point uncertainty, but the Python implementation is not formally verified. The potentially large number of segment combinations and deck shifts makes the method a *research auditor*, not a scalable general solver.

## Theorem 2: a stable exclusion margin

Let `W_i⊂R/Z` be the original parameter windows. Define

```
D = dist_circle(W0,W2),
M = sup{|f1(x)−f3(y)| : x∈W1, y∈W3}.
```

**Claim:** If `D>M`, there is no nondegenerate labeled square (of either orientation), even when the four vertices are allowed *arbitrary integer deck translates*.

**Proof:** Every labeled square has `|x0−x2|=|y1−y3|` in its planar lifts. The left side is at least the circular distance `D` between `W0,W2`; the right side is at most `M`. This contradicts `D>M`.

**Stability:** If every point of each selected piece moves horizontally by at most `δ` in the cylinder and vertically by at most `η` (not necessarily preserving linearity or even a graph structure), the corresponding left-hand lower bound decreases by at most `2δ`, and the right-hand upper bound increases by at most `2η`. Thus absence persists whenever `2(δ+η)<D−M`.

**Exact v212 strengthening:** For the four v212 hats, `D=7/16`. On the *selected* windows `W1,W3`, the height ranges are `[-1/16,0]` and `{0}` respectively. Hence `M=1/16` and the **exact strict exclusion margin is `3/8`**. The absence persists under any piecewise position perturbations satisfying `δ+η<3/16`. This improves v212's coarser global `|y1−y3|≤1/8` argument. It still does not give existence of a square on a general Jordan curve.

## Proposition 3: the deck-lift obstruction

Take the four horizontal period-one graphs with heights `(1,0,−1,0)`. Let each marked window be the same **nontrivial interval** `W=[7/16,9/16]`.

The alternating Liouville-area sum vanishes. The four graphs have no common intersection, and any labeled square has side `√2`. At `x=7/16`, the lifted points

```
(7/16,1), (−9/16,0), (7/16,−1), (23/16,0)
```

form an exact square. **Every projected parameter lies in `W`**, so it satisfies a superficial *all-original-window in the cylinder* test. Its relative deck shifts are `(0,−1,0,+1)`.

However, when all four actual physical subarcs are the ones in the **same** fundamental strip `x∈[7/16,9/16]`, no square exists: their total horizontal span is at most `1/8`, while the vertical span of any labeled square on these four horizontal levels is exactly `2`, and the x- and y-span of a Euclidean square are equal. The fixed-lift exact solver likewise finds no square when the shift tuple is `(0,0,0,0)`.

**Consequence:** Even obtaining a Floer-theoretic witness with *all projected vertices in original windows* is insufficient for the planar Toeplitz problem. We need to control the integer deck data as well as the tags. This is a conditional transfer obstruction, not a contradiction to any periodic-inscription theorem. The four horizontal strands are not asserted to arise from one Jordan curve.

## Theorem 4: a sufficient transverse persistence certificate

Fix four chosen linear graph segments with a nondegenerate square at strictly interior parameters `u*∈(0,1)^4`, and let `A` be the exact invertible `4×4` matrix of their square residual equations. Define

```
m = min_i min(u_i*,1−u_i*),
N = ||A⁻¹||_∞  (maximum absolute row sum),
g = dist(height-range(segment0),height-range(segment2)) > 0.
```

Here `g>0` means that opposite segments occupy disjoint **vertical** ranges. Consider *arbitrary continuous* perturbations of each selected graph segment, preserving its x-coordinate parameterization and changing its y-coordinate by at most `ε` uniformly.

**Claim:** For every `0≤ε<min(m/(2N),g/2)`, the four perturbed arcs contain a nondegenerate labeled square of the same orientation with four interior parameters.

**Proof:** The perturbed residual system is `A(u−u*)+E(u)=0`. Since each vertical perturbation is bounded by `ε`, every component of the four-vector square residual perturbation `E(u)` is bounded by `2ε` (the components are signed differences of two perturbed heights). The continuous map

```
T(u)=u*−A⁻¹ E(u)
```

sends the closed parameter cube of infinity-radius `m` around `u*` into its interior when `2Nε<m`. Brouwer's fixed-point theorem gives a zero of the perturbed residual in that cube. Since `2ε<g`, the two opposite vertices cannot coincide vertically, hence the square is nondegenerate. QED.

**Explicit exact certificate:** For four segments centered at `(7/8,1), (−1/8,0), (7/8,−1), (15/8,0)` with x-halfwidth `1/12` and graph slopes `(-2,3,7/2,-1/3)`, the exact rational computation yields

```
parameter margin m = 1/2
||A⁻¹||∞ = 819/55
opposite-arc height gap g = 37/24
ε threshold = min(55/3276, 37/48) = 55/3276.
```

In particular, **every** continuous vertical perturbation with sup-norm below `55/3276` preserves a nondegenerate labeled square on these four specified segments. This is an explicit finite open-family theorem, akin to the transversality strategy already used in v209, not a universal existence result.

## Tests and reproduction

Run from repository root:

```sh
python -m unittest discover -s tests -v
python code/periodic_square_cert_v213.py
```

The suite covers exact affine systems (including rank 0 and rank-deficient positive-dimensional families), both orientations, singleton parameter windows, periodic graph subdivisions, the v212 connector obstruction, exact non-square proofs on marked windows, tall-height lift bounds, the new deck-lift obstruction, 14 planted periodic squares with prescribed deck data, randomized planted line-segment configurations, and transverse persistence under randomly chosen rational PL vertical perturbations. CI (now `.github/workflows/ci.yml`, formerly `v212-exact.yml`) runs the entire test discovery suite.

## What is **not** proved and next target

1. We have not constructed four periodic completions of the **same arbitrary planar Jordan curve** satisfying exact area balance with a forced square on original pieces.
2. Neither the quotient tag test nor positivity of a separation gap controls *physical deck lifts*.
3. A positive certificate on finitely many PL arcs cannot be upgraded automatically to arbitrary rough Jordan curves; to pass to a limit, the vertex scale, parameters, source tags, and deck shifts must remain controlled.
4. A source-preserving Floer continuation or relative intersection invariant would have to preserve **both tag and deck data**. Establishing such an invariant is a separate, currently unproved theorem obligation.

**Next best constructive question:** For a specified family of Jordan-curve subarcs embedded into periodic completions, can a conserved local intersection index be *filtered simultaneously by source tags and prescribed deck shifts*? Use this exact PL auditor as the mandatory finite counterexample test, and do not promote a proposed filtered count until invariance through all relevant deformations (including boundary crossings) is proved.

## External basis

- Lorenzo Valerio, *A Proof of the Area Vanishing Square Peg Conjecture*, arXiv:2610.21076v2 (2026), https://arxiv.org/abs/2610.21076. Four periodic polygonal curves, signed-area balance, squares possibly degenerate.
- Ali Naseri Sadr, *The Four Periodic Inscription Problem*, arXiv:2610.24376v1 (2026), https://arxiv.org/abs/2610.24376. Related generalized four-periodic inscription theorem under distinct hypotheses.
- Internal primary sources: [v212 exact periodic bridge](v212-periodic-bridge-theorems.md), and original [v209–v211 bundles](../archive/).
