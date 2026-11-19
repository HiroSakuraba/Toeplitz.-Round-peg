# v225 — Canonical polygon-vertex ownership and exact moving singular-chart certificates

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** Conditional elementary algebraic/topological lemmas, a reusable exact-rational checker, and adversarial tests. **Not a complete general polygonal square-orbit cobordism, a global mod-2 invariant, or a proof of the unrestricted Toeplitz square-peg conjecture.**

## Executive result

The v224 event-cover program reduced the square equations for four moving polygon edges to polynomial constraints in the deformation parameter, but flagged a substantial number of charts as unresolved when an endpoint polynomial vanished *identically*, or when the four-by-four residual Jacobian was identically singular. The flagship strictly convex four-vertex deformation produced **408** product charts (including repeated-edge assignments), of which **128** were flagged. This prevented straightforward global event bookkeeping.

v225 resolves **all 128 flags in this specific example** by combining a canonical polygon-vertex ownership convention with a purely geometric exclusion for configurations involving two incident edges twice apiece. It also adds a general exact **augmented-polynomial-rank obstruction**, which certifies that some identically singular charts contain no square except possibly at finitely many isolated real algebraic parameter values. Crucially, it retains genuinely singular positive-dimensional square families as **inconclusive**, rather than assigning them spurious finite parity.

The work is an exact classification of the *chart types* in the control deformation. It is **not** yet a complete classification of all *geometric event orbits*: candidate algebraic roots still require clustering, singular-time interpretation, and geometric orbit identification.

## 1. Geometry and inherited square equations

Let `P(t)` be a simple polygon with `n` vertices whose coordinates move affinely in `t∈[a,b]⊂[0,1]`, with rational starting and ending positions. Assume the same polygon edge labels persist and the full motion remains Jordan. The v224 construction chooses four labeled directed edges `(i0,i1,i2,i3)` and square orientation `σ∈{±1}`. A point on edge `ij` is

\[
p_j(t,u_j)=P_{i_j}(t)+u_j\bigl(P_{i_j+1}(t)-P_{i_j}(t)\bigr),\quad 0\le u_j\le1.
\]

The square residual `Fσ=0` yields

\[
A(t)u=b(t),\qquad A(t)\in\mathbb Q[t]^{4\times4},\ b(t)\in\mathbb Q[t]^4.
\]

For `D(t)=det A(t)` not the zero polynomial, Cramer's rule yields `u_j(t)=N_j(t)/D(t)`. The v224 size cutoff `κ≥0` is `Hκ(t)=|E(t)|²−κD(t)²`, with `Hκ>0` equivalent to side² `>κ` wherever `D≠0`. The existing degree bounds (at most 4 for `D,N` and 10 for `Hκ`) remain unchanged.

Two *different* meanings of singularity must be kept separate: `D(t*)=0` for one special time and `D≡0` as a polynomial for the entire homotopy. The latter demands a rank analysis over `Q(t)` rather than Cramer's rule.

## 2. Lemma A: a canonical half-open polygon atlas

**Lemma.** For a fixed, oriented simple polygon with nonzero edges, every point on its boundary belongs to exactly one outgoing-edge chart of the form

\[
E_i^\circ=\{(1-u)P_i+uP_{i+1}:0\le u<1\}.
\]

**Proof.** An interior point of an edge belongs to that unique edge, because a simple polygon has no nonadjacent edge intersections or overlapping adjacent edges. A polygon vertex `Pi` belongs to its outgoing edge `Ei` with `u=0`, and not to its incoming edge `Ei−1` with `u=1`. The final edge excludes its terminal vertex `P0`, which is included as the start of edge zero. Thus the half-open edges partition the parametrized polygon boundary. The partition persists during any Jordan isotopy preserving vertex labels. QED.

Consequently, on **any nonsingular Cramer chart**, an identically vanishing endpoint polynomial is not necessarily an unresolved square event:

* `N_j≡0` means `u_j(t)=0` for all `D(t)≠0`. The parameter is **owned by this outgoing edge**; the chart may contain genuine squares and must be retained.
* `D−N_j≡0` means `u_j(t)=1` for all `D(t)≠0`. The entire branch is assigned instead to the **next outgoing edge**, so this incoming-edge chart is a redundant duplicate for geometric orbit counting.
* `Hκ≡0` on a nonsingular chart means the square side is **exactly** the cutoff for all `D(t)≠0`; the strict filter `side²>κ` excludes that chart without a parity event.

We isolate every *remaining nonzero* event polynomial using exact Sturm arithmetic. This is a correct status/event cover on each chart **away from** `D=0`; at isolated determinant roots a separate singular-time treatment is still required.

**Exact positive witness that the endpoint convention is not discarding a square.** For the rational trapezoid

\[
(0,0),(2,0),(2,3/2),(0,1/2),
\]

edge tuple `(0,0,1,2)` has a true square with parameter tuple `(1/2,1,2/3,1/2)`, including the second vertex at the end of edge zero. The incoming chart is excluded, but outgoing tuple `(0,1,1,2)` has the **same physical square** with parameters `(1/2,0,2/3,1/2)`. Exact fixed-time v223 census confirms two total geometric squares rather than an accidental loss.

## 3. Lemma B: two pairs on incident edges are impossible

**Lemma.** A nondegenerate Euclidean square cannot have precisely two distinct vertices on each of two *adjacent edges* of a simple polygon, regardless of labeling or orientation.

**Proof.** The four vertices of a square split into two pairs in one of two geometric ways:

1. Both pairs are endpoints of two **opposite sides**, in which case their supporting lines must be **distinct and parallel**. Two incident polygon edges share a vertex, so cannot lie on distinct parallel lines. If their supporting lines coincide, all four candidate vertices are collinear and cannot form a nondegenerate square.
2. Both pairs are **opposite vertices**, in which case their supporting lines are the two **perpendicular diagonals** through the common center of the square. The polygon's two incident edges lie on rays from their common polygon vertex. The midpoint of two distinct points on each such ray is a positive distance along that ray. For distinct noncollinear incident rays, the only common point of the supporting rays is their origin, so these midpoints cannot coincide. (If the rays were collinear, the square would again have four collinear vertices.)

These exhaust pairings of four square corners. QED.

In particular, 2+2 incident-edge charts can be discarded **at all times** as long as the polygon remains simple, including charts whose algebraic determinant is identically zero. The static v223 affine polytope solver independently verifies this exclusion across all such labelings and both square orientations on three rational simple polygons.

## 4. Theorem C: singular-chart rank obstruction over `Q(t)`

When `D≡0`, compute

\[
r=\operatorname{rank}_{\mathbb Q(t)} A(t),\qquad
s=\operatorname{rank}_{\mathbb Q(t)}[A(t)\mid b(t)].
\]

These ranks are determined exactly by searching for nonzero polynomial minors. If `s>r`, choose a nonzero `(r+1)×(r+1)` minor `M(t)` of the augmented matrix. Every `(r+1)`-minor of `A` vanishes **identically**, so for every real time `t` we have `rank A(t)≤r`. If `M(t)≠0`, the augmented rank is at least `r+1`, so the equations are inconsistent. Therefore

\[
\boxed{\{t:\exists u\in\mathbb R^4,\ A(t)u=b(t)\}\subseteq\{t:M(t)=0\}.}
\]

This is a **necessary algebraic support set**, even before imposing `[0,1]^4` or a positive side. Exact Sturm root isolation gives rational intervals for every possible exceptional time. If there are no roots in the interval, the chart is **provably empty throughout**. If there are roots, the chart is empty off them but their actual feasibility remains undecided without additional algebraic-number or semialgebraic checks. For example, the system `u0=u1=u2=0, 0=t−1/3` can be consistent only at `t=1/3`; the solver reports that necessary event rather than pretending to have a transversal branch. An analogous `t²−1/2` example gets a certified irrational root bracket without floating point.

If `s=r`, the rank obstruction proves **nothing** about real box feasibility: a genuine one-parameter or higher-dimensional family may exist. These charts are explicitly classified **genuinely unresolved** unless another geometric lemma excludes them. At fixed rational times the code calls the v223 rational-polytope solver to classify such a singular chart as empty above the size cutoff, isolated at its parameter boundary, or positive-dimensional. This pointwise tool is **not** incorrectly extrapolated to the full moving interval.

## 5. Exact 408-chart audit and independent tests

The v224 convex polygon control uses

\[
P(t)=\big((0,0),(2,0),(3-9t/8,3-27t/16),(0,2-27t/16)\big),\quad 0\le t\le1.
\]

(The deformation parameter here is the *normalized* interpolation from the v223 polygon at original parameter 0 to original parameter 9/8. The exact vertex square from v224 occurs at **normalized `t=8/9`**.) The prior independent v224 three-vertex orientation certificate verifies that this entire deformation remains simple. Its 408 candidate ordered edge charts include repeated edges but skip configurations having three points on one straight edge.

The statuses after v225 classification, using `κ=1/100` and time interval `[3/4,1]`, are:

| Status | Charts | Meaning |
|---|---:|---|
| Regular | 280 | No identically vanishing event polynomial; exact root covers |
| Vertex-owned | 40 | One or more parameters permanently at the **outgoing** endpoint `u=0`, still valid |
| Incoming-edge duplicate | 40 | A parameter permanently at `u=1`; excluded by canonical ownership |
| Structurally empty | 48 | Two pairs on incident polygon edges; square impossible |
| **Genuinely unresolved** | **0** | All prior 128 flags classified for **this example** |

Notably 44 of the structurally empty charts were not previously flagged by v224, while **four** contained identically singular Cramer matrices; the other previously unresolved charts carried persistent endpoint identities. There are still **224 raw algebraic candidate entries** across the retained chart event covers. These are **not 224 distinct geometric square events**. They may be duplicate polynomial roots, roots at which the solution is off the edge, or singularities without a square.

The new exact implementation uses a bounded cache of polynomial chart expressions across snapshot tests; no numerical approximation is introduced. The local v224+v225 suite runs **50 regression tests**. The v225 tests check:

- Exact reduction of v224's **128** unresolved flags to **zero** in the control deformation.
- Agreement with **independent v223 full orbit census** at rational times `3/4`, `17/20`, `8/9`, `11/12`, and `1`, including the true two-square vertex-alignment snapshot and the later three-square snapshot.
- Preservation of the trapezoid's genuine endpoint square under half-open ownership, rotation/reversal checks, and an exact random rational perturbation.
- Full two-pair adjacent-edge exclusion on three unrelated Jordan polygons via a separate affine feasibility enumeration.
- Synthetic rank-deficient polynomial systems with always-empty, rational-exceptional, irrational-exceptional and generically consistent cases.
- A rectangle's persistently consistent singular chart with a **positive-dimensional square family** remains unresolved; the program does not pretend its orbit parity is finite.

## 6. Why this is not a completed global polygonal cobordism

A canonical half-open atlas is a **bookkeeping device for geometric square orbits**. It does not make an actual square disappear when a parameter passes through a vertex, and it does not prove that its physical branch is transverse there. In addition:

1. **Event clustering:** The 224 remaining raw candidates have not been reduced to distinct real algebraic times with geometrically identified orbit changes.
2. **Special determinant roots:** Cramer formulas are unavailable at `D=0`. A global solver still needs exact local semialgebraic or perturbation analysis at these times, including positive-dimensional and simultaneous endpoint strata.
3. **Generic orbit parity:** A polygon may have a continuum of squares (rectangle), so raw orbit parity need not be defined at each time. Even strictly convex polygon vertex events can yield exceptional counts (v223), and a global mod-2 cobordism requires an appropriate chain/cycle treatment.
4. **Arbitrary Jordan curves:** Neither a nonzero global orbit class nor a uniform positive square size under C0 approximation has been established. The v220 fixed rectifiable curve with collapsing squares remains a mandatory negative control.
5. **External novelty:** These are elementary facts about polygon parameterizations, square geometry, and polynomial rank. No claim is made that they establish a new major mathematical theorem or solve the Toeplitz conjecture.

**Recommended v226:** Cluster algebraic candidates from the retained half-open charts using polynomial gcd and disjoint rational isolating intervals, then verify actual geometric square-orbit incidence on either side of each clustered event. For the controls, compare exact orbit sets with v223 at rational chamber witnesses. Singular instantaneous roots and positive-dimensional strata must continue to be reported rather than silently discarded. Only then attempt a genuine relative parity class.

## 7. Reproduction and sources

```bash
python -m unittest discover -s tests -q
python code/persistent_strata_v225.py
```

All new code is standard-library-only, with rational arithmetic. Regression tests check finite instances; the mathematical arguments above supply the general statements within their hypotheses. No Lean/Coq/Sage/Arb interval proof is claimed.

- Internal exact base: [v223 polygon orbit census](v223-complete-polygon-orbit-census-and-homotopy-audit.md), [v224 moving polynomial event covers](v224-algebraic-moving-polygon-event-covers.md), and [v220 collision-link obstruction](v220-fixed-rectifiable-collision-link-obstruction.md).
- Siniša Vrećica and Rade Živaljević, *Fulton–MacPherson compactification, cyclohedra, and the polygonal pegs problem*, [arXiv:0810.1439](https://arxiv.org/abs/0810.1439).
- Jason Cantarella, Elizabeth Denne, and John McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, [arXiv:2103.07506](https://arxiv.org/abs/2103.07506).
- Benjamin Matschke, *On the Square Peg Problem and Its Relatives*, [arXiv:1001.0186](https://arxiv.org/abs/1001.0186).
