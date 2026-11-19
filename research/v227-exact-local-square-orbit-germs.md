# v227: Exact algebraic event germs and geometric-square continuation through polygon vertices and cutoffs

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** Conditional exact-rational **local event-germ** theorem and executable geometric-orbit incidence certificates for one affine-moving simple polygon. This is **not** a complete all-polygon cobordism, an all-Jordan invariant, or a proof of the unrestricted Toeplitz square-peg conjecture.

## 1. What v226 left unresolved

[v226](v226-exact-algebraic-event-clusters-and-orbit-chambers.md) identified exact algebraic event times and counted the *geometric*, rather than labeled, inscribed squares in each open event-free time chamber. It could establish that a convex polygon went from one to three squares at `t=8/9`, or that a uniform scaling went from zero to one square past the irrational threshold `sqrt(2)-1`. But an endpoint count does **not** specify which square branches meet at the event. A `1 → 3` transition might be one surviving square plus a pair birth, a more complicated reconnection, or a mixture of events.

The present batch constructs the **incidence graph local to one event time**. It matches polynomial root branches from all canonical four-edge charts to individual geometric squares *at the event*, while independently verifying each side against v223's complete fixed-polygon square census. It uses no floating-point time guesses or approximate point matching, and refuses unsupported singular cases.

## 2. Exact algebraic structure

Let `P(t)=(P_0(t),...,P_{n-1}(t))`, `t∈[a,b]`, be one Jordan polygon throughout an affine-in-`t` vertex motion with **rational** initial and terminal coordinates. Fix a squared-side cutoff `κ>0`. Use both orientations `σ∈{±1}` and all allowed ordered edge selections `e=(e_0,...,e_3)`, *including repeated edges*. For each edge chart, let

\[
X_i(t,u_i)=P_{e_i}(t)+u_i\bigl(P_{e_i+1}(t)-P_{e_i}(t)\bigr),\quad 0\le u_i\le1.
\]

The labeled square equations give `A_e,σ(t)u=b_e,σ(t)` with rational polynomial entries of degree at most one. Define

\[
D(t)=\det A(t),\qquad N_i(t)=\det A_i(t),\qquad U_i(t)=D(t)-N_i(t).
\]

For `D(t)≠0`, `u_i(t)=N_i(t)/D(t)` is a uniquely determined analytic branch. The four vertex coordinates can be written as rational functions

\[
X_i(t)=\frac{C_i(t)}{D(t)},\qquad
C_i(t)=D(t)P_{e_i}(t)+N_i(t)\bigl(P_{e_i+1}(t)-P_{e_i}(t)\bigr).
\]

As in v224, a rational polynomial `H_κ(t)` satisfies

\[
\lVert X_1(t)-X_0(t)\rVert^2-\kappa=\frac{H_\kappa(t)}{D(t)^2}.
\]

Hence **all geometry and admissibility information is contained in finitely many rational polynomials**. The prior v225 outgoing-edge ownership convention `0≤u_i<1` is used on *open time chambers*, while the *closed* `0≤u_i≤1` rule is used to find all incident event roots; incoming-edge endpoint descriptions are then deduplicated. Vertex ownership applies to the *curve parameter*, not merely its image in the plane.

## 3. Lemma: exact signs and one-sided germs at an algebraic event

Let `α` be the unique real root of a polynomial `g∈Q[t]` in a certified rational isolating interval `I`. Every polynomial `p∈Q[t]` has an **exactly decidable sign** at `α`:

- `p(α)=0` if and only if the polynomial GCD `gcd(p,g)` has a real root in the same isolating interval, with endpoint handling as in v226.
- Otherwise, refine `I` using exact Sturm isolation until it contains no root of `p`; evaluate `p` at any rational interior point. This terminates mathematically because `p(α)≠0` and the number of nearby roots is finite; the implementation has a finite resource cap and fails closed.

If `p(α)=0` but `p` is not identically zero, let `m` be the least positive integer with `p^(m)(α)≠0`. Then Taylor expansion gives

\[
\operatorname{sgn}p(\alpha+0)=\operatorname{sgn}p^{(m)}(\alpha),\qquad
\operatorname{sgn}p(\alpha-0)=(-1)^m\operatorname{sgn}p^{(m)}(\alpha).
\]

This determines whether each branch lies in the physical half-open edge windows on the **left** and **right** of `α`, and whether its squared side is strictly greater than `κ`, without any floating-point approximation. The implementation checks the resulting side eligibility again at a rational sample of the adjacent certified event-free chamber.

**Rational times** are exact `Fraction` values. For **irrational times**, the output retains the algebraic polynomial and isolating interval and does *not* substitute a rational midpoint for the instantaneous event.

## 4. Lemma: geometric orbit grouping over an algebraic number field

Two surviving edge charts at `α` may describe the same physical square under cyclic or reversed labeling. Let their point-coordinate numerators be `C_i(t)` and `E_j(t)` and denominators `D(t)` and `K(t)` respectively. At the event time their corresponding coordinates agree exactly if and only if

\[
\boxed{\bigl(C_iK-E_jD\bigr)(\alpha)=0.}
\]

The preceding algebraic-sign lemma decides this equality. Test the finite permutations of the four vertices and group **unlabeled geometric squares**, not chart labels. This remains exact at irrational `α`. On each side, compute each active chart's square coordinates at a rational sample *within the adjacent event-free chamber*, and deduplicate by sorted rational vertex sets. Thus two distinct right-hand square orbits may meet a **single** event-time geometric square. They must not be merged into one branch merely because they coincide at `α`.

This explicitly distinguishes a **pair birth** `(left,right)=(0,2)` from a **continuing square** `(1,1)`, a **pair death** `(2,0)`, and a **size-cutoff entry** `(0,1)` whose event-time squared side equals `κ`. Other incidence patterns are reported without inventing a generic topological interpretation.

## 5. Conditional local incidence theorem

**Theorem (exact local geometric-square orbit germ ledger).** Assume:

1. One rational-vertex affine polygon family is certified Jordan throughout the selected interval, using v224's sufficient non-collinearity test.
2. v225 produces a complete eligible-chart classification, without any persistently feasible singular square family or cutoff-locked chart.
3. v226 produces a complete, exactly ordered candidate event cover, and the two rational side samples lie in adjacent event-free chambers.
4. Every **physical positive-size square at the rational event time** has at least one nonsingular `D(α)≠0` chart presentation, verified against the full v223 event-slice census; every physical side-chamber orbit has such an event-limit anchor. For irrational `α`, in addition require that no potentially feasible four-edge chart has `D(α)=0`, and no generically inconsistent singular chart has an obstruction minor zero there. If this condition fails, report **inconclusive** instead of silently dropping possible instant-only square roots.
5. The algebraic isolations and all exact root-sign determinations complete within their stated resource bounds.

Then the algorithm gives a complete **finite local incidence ledger within these assumptions**: for each unique physical square at `α` with `side²≥κ` supported by a nonsingular chart, it records the distinct geometric square orbits incident from the left and from the right, including their cyclic/reversed chart descriptions. The union of all left and right incident orbits equals the independently certified v226/v223 geometric chamber census. At rational event times, every `side²>κ` event-slice orbit is independently matched. This distinguishes genuine local branch births, deaths, continuations, and simple size-filter crossings in the covered model.

**Proof sketch.** On each nonsingular chart, Cramer's rule supplies a unique continuous analytic root branch through the event. The rational algebraic-sign lemma determines its exact one-sided half-open edge and size eligibility. Exact polynomial-coordinate equalities identify all chart descriptions whose limits are one physical square. No side branch can enter or leave inside a v226 event-free chamber; consequently a rational sample in that chamber identifies the orbit represented by the incident analytic root. Taking unions over all chart descriptions is valid because the v225 chart classification covers the polygon edges and v223 supplies an independent complete geometric-square census. The final equality with that census is **checked**, not assumed. The fail-closed singular guards prevent hidden root branches from invalidating the matching. QED, within the explicit hypotheses.

**Not claimed:** a general semialgebraic-cobordism theorem for all moving polygons, a canonical orientation of the full quotient, or the existence of a positive-size square on arbitrary continuous Jordan curves. This is a conditional incidence classification, not the global invariant required for the unrestricted Toeplitz conjecture.

## 6. Concrete exact examples

### A. Genuine two-square birth at a polygon vertex

Use v223's convex quadrilateral

\[
P(s)=\bigl((0,0),(2,0),(3-s,3-3s/2),(0,2-3s/2)\bigr),
\]

with normalized motion `s=(9/8)t`, restricted to `3/4≤t≤1`, and cutoff `κ=1/100`. The v225/v226 pipeline has **408** eligible chart classifications, **224** raw polynomial candidate records, and exactly **one** geometric event at `t=8/9` (`s=1`). The exact event polygon has two geometric squares; one is the unit square

\[
(1,0),\ (2,0),\ (2,1),\ (1,1).
\]

The new event-germ ledger identifies:

| Event square | Left orbits | Right orbits | Exact chart presentations | Classification |
|---|---:|---:|---:|---|
| Unit square at vertex | 0 | 2 | 16 | **Pair birth** |
| Other event square | 1 | 1 | 8 | **Continuation** |
| **Total** | **1** | **3** | **24** | Net `+2`, no mod-2 parity change |

The event snapshot contains **two** geometric squares, not two event *germ branches*. Reversing the polygon motion places the same event at `t=1/9` and changes **pair birth** to **pair death**, preserving the continuing square.

### B. One genuinely irrational size-filter crossing

Scale the single Jordan quadrilateral `((0,0),(2,0),(3,3),(0,2))` by `1+t`, `t∈[0,1]`, with `κ=144/25`. The complete v226 cover has **88** raw candidates and **six** distinct algebraic times. The v227 classification yields:

- First three candidate times: **no geometric square orbit at or above cutoff**; these are spurious candidates.
- Fourth time: `t=√2−1`, certified as the unique real root of `t²+2t−1`; **one geometric square crosses the size cutoff**, with **eight chart/label presentations** of one physical square. Left count zero, right count one. The event-time physical square has `side²=κ` and is correctly **excluded** by the strict `>` cutoff.
- Last two candidate times: the same geometric square **continues** from left to right, with no change in count. These events reflect chart algebra, not a physical square birth or death.

The corresponding algebraic event-time test uses signs of polynomial germs **at the irrational root itself**. No numerical approximation is used as a proof step.

### C. Negative controls

- The `2×1` rectangle has positive-dimensional inscribed square families; the regular event-germ procedure refuses it as **inconclusive**.
- Artificially deleting chart evidence at the rational event causes the full independent v223 event-slice census comparison to fail, rather than asserting a complete classification.
- Inserting a possible exceptional singular chart at the irrational cutoff event forces a **fail-closed** result, because a real algebraic event-time square cannot safely be ruled out by counting neighboring chambers.
- One-sided exact signs are checked at simple, double, and triple roots, at rational and irrational algebraic times. Time reversal transforms birth into death as required.

## 7. Reproduction and next obligations

Run:

```bash
python -m unittest discover -s tests -v
python code/local_orbit_germs_v227.py
```

`code/local_orbit_germs_v227.py` uses only the Python standard library and the existing v223–v226 modules. `tests/test_local_orbit_germs_v227.py` adds **29** exact-rational regression tests, all locally passing. GitHub CI runs those alongside all prior project regressions, audits and archived-asset checksum checks. The tests supplement, rather than constitute, the conditional mathematical arguments above.

The next research step (v228) should handle **singular-time semialgebraic germs** when `D(α)=0` but the static event configuration may still carry isolated physical squares. Doing this completely requires exact real-algebraic feasibility over a neighborhood of an algebraic event (e.g. subresultants/CAD or a carefully bounded local Puiseux method). It should start with the genuine v225 rectangle/parallel-edge countercontrols, not simply weaken the fail-closed guards. Separately, a global positive-size class under arbitrary `C⁰` approximation remains **unproved**. Even a finished polygonal cobordism would not by itself solve Toeplitz.

### Research background

- v223 complete fixed rational polygon square-orbit census; v224 polynomial event cover; v225 persistent edge ownership/rank certificates; v226 root clustering and chamber orbits. These are all in this repository.
- Cantarella–Denne–McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, [arXiv:2103.07506](https://arxiv.org/abs/2103.07506): smooth/generic configuration-space framework, not an all-Jordan noncollapse theorem.
- Benjamin Matschke, *On the Square Peg Problem and Its Relatives*, [arXiv:1001.0186](https://arxiv.org/abs/1001.0186): survey and parity precedents.

**Final scientific caveat:** The unrestricted Toeplitz square-peg conjecture remains open. The v227 result is finite-dimensional, exact, and conditional on explicit audit hypotheses; none of the local counts is asserted to define a universal odd square-orbit parity on arbitrary Jordan curves.
