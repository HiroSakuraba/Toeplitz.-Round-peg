# v228 — Exact rank-deficient square feasibility at a rational or irrational polygon-event time

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** Exact, finite-dimensional **event-slice** theorem and a standard-library algebraic-feasibility solver. This is **not** a resolution of singular one-sided square branches, a complete polygonal cobordism, a global all-Jordan invariant, or a proof of the unrestricted Toeplitz square-peg conjecture.

## 1. The specific v227 obstruction

Versions [v223](v223-complete-polygon-orbit-census-and-homotopy-audit.md)–[v227](v227-exact-local-square-orbit-germs.md) enumerated fixed-polygon square orbits and identified algebraic event times. v227 uses Cramer's rule to match **nonsingular** square-branch germs but deliberately stops when `det A(alpha)=0`. That guard prevents an unobserved square supported at a rank-drop time from being lost. However, singularity by itself is not evidence of a genuine square, nor does it tell us whether the feasible square locus is discrete or continuous.

Here we solve the **instantaneous decision problem** exactly, including when the event time `alpha` is irrational. This is a substantial but carefully limited prerequisite to singular one-sided continuation.

## 2. Main theorem: bounded algebraic affine box slice

Fix a real algebraic time `alpha` specified by a nonzero rational polynomial `g(t)` and a *certified isolating interval containing exactly one specified root*. For a fixed ordered selection of four polygon edges and one square orientation `sigma`, the square residual equations at any time are

\[
 A(t) u=b(t), \qquad 0\le u_0,u_1,u_2,u_3\le1.
\]

The entries of `A,b` are affine rational polynomials when each polygon vertex travels affinely in time. At `alpha`, the feasible square-parameter set is the compact convex polytope

\[
 K_\alpha=\bigl\{u\in[0,1]^4 : A(\alpha)u=b(\alpha)\bigr\}.
\]

**Proposition 1 (exact complete vertex enumeration).** The rank `r` of `A(alpha)`, feasibility of `K_alpha`, and **every vertex** of `K_alpha` are decidable by finite polynomial determinant calculations and exact real-algebraic sign tests. Specifically:

1. Find a nonvanishing `r×r` minor of `A` at `alpha`, for maximal `r`, by testing the sign of every relevant minor at `alpha`. The selected `r` rows form an independent equation subsystem over `Q(alpha)`.
2. Every vertex of the bounded affine slice must have at least `4-r` linearly independent active unit-box constraints. Enumerate each choice of `4-r` coordinates set to 0 or 1.
3. Append these constraints to the `r` independent affine equations, giving a square `4×4` matrix `M(t)`. Where `det M(alpha)` is nonzero, Cramer's rule gives all four coordinates as **polynomial quotients** `N_j(alpha)/D(alpha)`.
4. Verify every original affine equation by zero tests of `A_i N-b_iD` at `alpha`. Verify each lower/upper bound by exact signs of `N_j D` and `(D-N_j)D` at `alpha`. Deduplicate candidate vertices using cross-multiplied polynomial equality.

The algorithm returns the empty set if and only if the affine box slice is infeasible, because a nonempty bounded polytope has a vertex. This is finite-dimensional rational linear algebra over a **real algebraic field**, not floating-point sampling or a general-purpose cylindrical algebraic decomposition. Its finite enumeration budget may cause a **reported inconclusive result** rather than a false negative.

**Proof.** Let `K` be a nonempty bounded polytope of affine dimension at most `4-r`. Any vertex must be cut out by an additional `4-r` independent active box hyperplanes, or else a nonzero null direction preserves all equalities and leaves the point on a feasible segment, contradicting extremality. Conversely, every candidate from an invertible augmented square system satisfying the original equalities and box bounds is a vertex. All entries lie in `Q(alpha)` and the defining polynomial root is isolated, so signs and equalities are decidable by polynomial gcd/Sturm isolation (the exact sign machinery already audited in v226–v227). QED.

**Proposition 2 (positive-size square family classification).** Let the squared side function of the ordered square be

\[
 q(u)=|P_{e_1}(\alpha)+u_1 d_{e_1}(\alpha)-P_{e_0}(\alpha)-u_0 d_{e_0}(\alpha)|^2.
\]

This is a **convex quadratic** in `u`. Every value is the image of a convex combination of polytope vertices, and Jensen's inequality gives `q(u) <= max_vertex q(v)`. Therefore the condition `q>kappa` anywhere in the affine feasible box can be decided solely from its enumerated vertices, even when the feasible system is singular. The exact test uses the polynomial numerator of `q-kappa` and the square of the nonzero Cramer denominator.

The checker classifies each single chart as:

- **`infeasible`**: no solution in the closed parameter box;
- **`no-positive-size-square`**: the feasible set exists, but `q<=kappa` throughout;
- **`isolated-square`**: the feasible set has one vertex (hence is a singleton) with `q>kappa`;
- **`positive-dimensional-square-family`**: the feasible set has two or more vertices and some point has `q>kappa`.

For an embedded Jordan polygon whose edges are nonzero, distinct chart parameters represent distinct ordered physical points on that chart, so the last case genuinely produces a continuum of physical square configurations. The solver **does not** assign a finite orbit parity to that continuum. These statuses refer to one ordered edge chart; global geometric-square orbit identification remains a separate step.

### Algebraic witness representation

All candidate parameters and square coordinates are returned as polynomial quotient pairs in `Q(alpha)`, with denominators proved nonzero. Equalities are certified by cross-multiplication and root-sign evaluation. In particular, no decimal approximation of an irrational event is used as a proof step.

## 3. An actual irrational rank-drop control

Set

\[
 \alpha=\sqrt2-1\quad(\alpha^2+2\alpha-1=0),
\qquad A(t)=\begin{pmatrix}
 t+1&1&0&0\\2&t+1&0&0\\0&0&1&0\\0&0&0&1
 \end{pmatrix}.
\]

The top-left determinant is `(t+1)^2-2`, so `A(t)` has full rank away from its two algebraic roots and **rank 3** at `alpha`.

Exact variants of `b(t)` exercise three different behaviors at the **same irrational event**:

1. For `b=(0,0,1/2,1/2)`, the feasible slice contains **one** boundary point: `u=(0,0,1/2,1/2)`. Although the equations have rank 3, the nonnegativity constraints exclude motion in either nullspace direction.
2. For `b=((t+1)/2,1,1/2,1/2)`, the feasible slice is a **line segment** with two box vertices.
3. For `b=(1,0,1/2,1/2)`, the affine equations are inconsistent at `alpha` and the feasible box is empty.

These are **abstract 4×4 affine systems**, not claimed to arise from a simple polygon or to describe square birth/death. They are exact controls establishing that a determinant zero does not by itself classify a geometric event.

## 4. Independent whole-polygon event-slice audits

### A. Convex quadrilateral v227 pair-birth at `t=8/9`

Use the v223–v227 polygon `P(s)=((0,0),(2,0),(3-s,3-3s/2),(0,2-3s/2))`, with normalized `s=9t/8`, cutoff `kappa=1/100`, and event time `t=8/9`. Of the 408 admissible labeled edge charts, **28** have `det A(alpha)=0`. The new exact checker classifies:

| Singular event-slice status | Chart count |
|---|---:|
| Affine box infeasible | 20 |
| Feasible but no square with `side² > 1/100` | 8 |
| Singular chart carrying a positive-size square | **0** |
| Inconclusive under resource cap | **0** |

Thus **none of those 28 singular charts hides a positive-size event-time square**. The independent v223 static census finds precisely two geometric squares at the event, and v227's nonsingular-germ certificates account for them: one continuing square and one two-branch birth. This result **strengthens the event-slice verification** of v227, but does not prove all possible singular *one-sided* continuations through arbitrary polynomial rank drops.

### B. Irrational scaling cutoff at `alpha=sqrt(2)-1`

For the quadrilateral `((0,0),(2,0),(3,3),(0,2))` uniformly scaled by `1+t`, with cutoff `kappa=144/25`, an exact algebraic-time audit finds **four** singular edge charts at `alpha`. Each has a feasible box but **no strictly above-cutoff square**; no chart is left unresolved. This complements v227's nonsingular-charts certificate for the one genuine irrational size-cutoff entry. A square **exactly on the cutoff** remains excluded by the strict filter.

### C. A genuine rectangle continuum

For the static `2×1` rectangle and edge selection `(0,0,2,2)`, orientation `-1`, the algebraic square system has rank 3 and **two distinct feasible box vertices**. Its square side can be positive, so it contains a **positive-dimensional family of square solutions**. Raising the strict squared-side cutoff to `1` removes every square from this chart's filtered locus, which the checker detects without deleting the continuum of unfiltered solutions. A separate adjacent-edge chart is proved entirely infeasible.

This is the indispensable negative control: a general singular-event pipeline cannot count an infinite family as one regular orbit or silently drop it.

## 5. Reproduction and implementation

The standard-library-only code is:

- [`code/singular_event_feasibility_v228.py`](../code/singular_event_feasibility_v228.py): exact algebraic determinant rank, complete affine-box vertex enumeration, feasible-square witnesses, positive cutoff tests, and audit of every singular chart at a specified time;
- [`tests/test_singular_event_feasibility_v228.py`](../tests/test_singular_event_feasibility_v228.py): **26 new tests**, including 36 seeded independent exact rational comparisons with the v223 static polytope routine, rational rectangle and trapezoid controls, irrational rank-drop cases, positive-size decisions, and both whole-polygon event audits.

Run:

```bash
python -m unittest discover -s tests -p 'test_singular_event_feasibility_v228.py' -v
python code/singular_event_feasibility_v228.py
```

All **26 v228 tests passed locally**. A broad standalone rerun of all legacy test modules hit the local execution timeout and is **not** reported as passing; GitHub Actions is responsible for running the full repository suite independently. All algebraic sign checks and original square residual calculations are exact; the theorems above provide the justification rather than test counts.

## 6. What is still missing (v229)

v228 settles the **static event slice** `A(alpha)u=b(alpha)` even when rank drops, including true irrational algebraic instants. But the following remain different, unresolved problems:

1. **Singular one-sided branch germs**: if `D(alpha)=0` and Cramer numerators also vanish, `u_i=N_i/D` may admit a removable finite limit. One must cancel the maximal common factor at the relevant algebraic root, check the limiting parameters and orientation/size signs, and prove the resulting germs cover the local time neighborhood. If numerator orders are lower than denominator orders, the branch has a pole and exits bounded source charts; this too needs certification.
2. **Exceptional instant-only roots**: a rank-deficient chart may be feasible at one instant even though it is inconsistent at neighboring times. Such an isolated zero need not carry a regular continuation branch. It must be represented separately, not invented as a birth or death.
3. **Persistent positive-dimensional square strata**: the rectangle shows that a finite number of square-orbit generators is not always defined. A correct future cobordism needs a perturbation, stratified-chain model, or explicit separate continuum handling.
4. **Global non-collapse under `C⁰` Jordan limits**: all polygonal event machinery leaves untouched the central unrestricted Toeplitz obstruction. v220 proves that even a fixed rectifiable curve can have arbitrarily small square witnesses; an all-Jordan existence proof must protect a positive-scale class or supply another noncollapse mechanism.

A useful **v229** milestone is to certify removable Cramer poles via exact algebraic valuations and combine that with v228's event-slice vertex data, explicitly preserving all singular-in-time caveats.

### Research provenance

Internal prerequisites: [v223 fixed-polygon exact census](v223-complete-polygon-orbit-census-and-homotopy-audit.md), [v225 rank certificates](v225-persistent-edge-strata-rank-certificates.md), [v226 algebraic event clustering](v226-exact-algebraic-event-clusters-and-orbit-chambers.md), and [v227 exact local germs](v227-exact-local-square-orbit-germs.md).

Prior mathematical framework, **not claimed as proving our proposed universal invariant**: Cantarella, Denne and McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem* (arXiv:2103.07506); Matschke, *On the Square Peg Problem and Its Relatives* (arXiv:1001.0186). Our finite-dimensional exact rank and polytope arguments do not imply new general Jordan-curve existence results.

**Scientific conclusion:** This is a valid, testable refinement of the singular-event portion of the polygonal program. **The unrestricted square-peg conjecture remains open.**
