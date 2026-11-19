# v229 — Exact singular Cramer germs: removable limits, poles, and instant-only support

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** An exact algebraic **one-sided chart-germ** lemma, a bounded implementation, and a whole-polygon rank-drop audit. These are local results. **This is not a complete geometric-square orbit cobordism, a positive-scale invariant, or a proof of the unrestricted Toeplitz square-peg conjecture.**

## 1. The gap between v227 and v228

[v227](v227-exact-local-square-orbit-germs.md) matched nonsingular square branches through algebraic event times. [v228](v228-singular-algebraic-event-slices.md) decided the **complete feasible square set at the event itself**, including rank-deficient affine systems and positive-dimensional square families. Neither proved whether a Cramer branch with determinant zero at the event could approach a finite configuration from the left or right.

The key distinction is between (i) **rank loss of the chosen chart matrix**, (ii) a **genuine pole** in its unique nearby square-parameter branch, (iii) a **removable `0/0`** representation of a bounded branch, and (iv) a chart whose matrix is **singular at every nearby time**. These situations require different treatments.

Version v229 classifies the first three exactly when the determinant polynomial is not identically zero; the fourth remains a separate semialgebraic problem. The algorithm only analyzes *one ordered four-edge chart at a time*, with a chosen square orientation and a fixed squared-side cutoff.

## 2. The main singular-continuation lemma

Let `alpha` be the isolated real root of a rational polynomial `g(t)` in a certified rational isolating interval. Consider the moving-polygon square equations

\[
A(t)u=b(t),\qquad u\in[0,1]^4,
\]

where entries of `A,b` are rational polynomials (affine for a polygon whose vertices move linearly). Write the Cramer determinant and four replacement determinants as

\[
D(t)=\det A(t),\qquad N_j(t)=\det A_j(t).
\]

**Proposition (meromorphic Cramer branch).** Suppose `D` is not identically zero and `D(alpha)=0`. Put `m=ord_alpha D`, the least nonnegative integer `m` such that `D^(m)(alpha) ≠ 0`. Every punctured neighborhood of `alpha` small enough to exclude the other zeros of `D` has one uniquely specified algebraic root branch

\[
u_j(t)=\frac{N_j(t)}{D(t)}.
\]

Each coordinate satisfies exactly one of these cases:

1. If `N_j` is nonzero and `ord_alpha N_j < m`, then `u_j` has a **genuine pole** and the unique nearby branch cannot remain within the bounded source-parameter box.
2. If `N_j ≡ 0`, or `ord_alpha N_j ≥ m`, the singularity is **removable** and

\[
\boxed{u_j^*=\lim_{t\to\alpha}u_j(t)=\frac{N_j^{(m)}(\alpha)}{D^{(m)}(\alpha)}.}
\]

The ratio on the right has a nonzero denominator. This formula works for rational and irrational algebraic event times, and it naturally gives `0` when the numerator vanishes to a strictly higher order. Each limit is an exact element of `Q(alpha)`.

**Proof.** Taylor-expand `D` and each nonzero `N_j` about `alpha`. For `D`, the first nonzero term is `D^(m)(alpha)(t-alpha)^m/m!`. Comparing the first nonzero powers establishes a pole in the first case and cancellation of `(t-alpha)^m` in the second. The factorials cancel in the ratio of the two derivatives. Since a nonzero polynomial has isolated real roots, `D(t) ≠ 0` in a sufficiently small punctured neighborhood, so the Cramer branch there is unique. QED.

This proof does **not** say every finite limit is a physical inscribed square. Its limiting parameters must be in the closed source box; the actual one-sided branches must satisfy the canonical edge ownership and the **positive-size filter**.

## 3. Exact one-sided physical admissibility

The v225 unique polygon-vertex ownership convention uses `0 ≤ u_j < 1` on each chart. On each side of `alpha`, this is equivalent to

\[
\operatorname{sgn}\!\left(\frac{N_j(t)}{D(t)}\right)\ge0,
\qquad
\operatorname{sgn}\!\left(\frac{D(t)-N_j(t)}{D(t)}\right)>0.
\]

For a rational polynomial `p`, the one-sided signs at `alpha` follow from the first nonvanishing derivative of `p` at `alpha`: if `p` has a root of order `r` there, its right sign is the sign of `p^(r)(alpha)` and its left sign is `(-1)^r` times that sign. For irrational `alpha`, the sign of `p^(r)(alpha)` is found from polynomial GCDs and exact Sturm root isolation, never from a decimal approximation.

Let `H_kappa(t)` be the v224 squared-side polynomial, satisfying

\[
\lVert e(t)\rVert^2-\kappa=\frac{H_\kappa(t)}{D(t)^2},\qquad \kappa>0.
\]

Its one-sided strict size eligibility is exactly `sgn(H_kappa)>0`, since `D(t)^2>0` off the event. A **finite physical** parameter limit implies the square edge extends continuously; in that case `H_kappa/D^2` cannot have a pole. The sign of the cutoff difference at the event is then obtained from the `2m`-th derivatives if its vanishing order is exactly `2m`, or is zero when the numerator vanishes to higher order.

Together these tests yield **left and right** physical branch eligibility, without mistaking a source-vertex crossing or a size-cutoff event for a pole. If a coordinate genuinely diverges, that chart has no bounded nearby square branch. Other charts are not thereby excluded.

### Independent event-slice membership check

For a removable limit `u*`, v229 independently verifies **every original row** of `A(alpha)u*=b(alpha)` by cross-multiplying the four rational-function coordinate representations and using exact algebraic equality. It does not require `u*` to be a vertex of the v228 feasible polytope: after rank loss, a genuine branch may select an **interior point** of a positive-dimensional event-slice family. It also checks membership in `[0,1]^4` and compares with v228's complete event-slice feasibility decision.

A single removable branch only selects one event-slice solution. It is **not** the complete set of instantaneous solutions when the matrix loses rank. This qualification is essential.

## 4. What an identically singular chart can and cannot prove

If `D ≡ 0`, there is no unique Cramer branch to cancel; it is invalid to divide all numerators by a factor and declare an inscribed square. This case remains `persistent-singular` until a separate rank/feasibility argument is supplied.

There is nevertheless a useful exact exclusion law. If `rank_{Q(t)}[A|b] > rank_{Q(t)} A`, the v225 polynomial rank routine produces a **nonzero augmented-rank minor** `M(t)` such that any solution of `A(t)u=b(t)` must satisfy `M(t)=0`.

- When `M(alpha) ≠ 0`, the chart is infeasible throughout a neighborhood of the event, including the instant itself.
- When `M(alpha)=0`, the nonzero polynomial has an isolated root there. Therefore the chart has **no solutions at sufficiently nearby times distinct from alpha**. v228 must separately determine whether the event slice contains no square, an isolated square, or a positive-dimensional family. An event-positive witness in this situation is **instant-only**, not a conventional pair birth.
- When the generic coefficient and augmented ranks agree, the rank obstruction gives no exclusion. A persistently feasible singular square family remains a valid **inconclusive** outcome, not a finite parity count.

This supports an *instant-only incidence marker* without silently inventing continuation through a genuinely rank-deficient family.

## 5. Concrete exact controls and independent comparisons

### A. A simple rational removable rank drop

Let `alpha=1/2`, `D(t)=t-1/2` and

\[
(N_0,N_1,N_2,N_3)=D(t)(1/2,1/3,2/3,1/4).
\]

Every coordinate has a finite limit in the open source box. With `H_kappa=D(t)^2/8`, both one-sided germs are above the cutoff and the event-size excess is `1/8`. The determinant still vanishes exactly at `alpha`, so a checker that rejected every `D(alpha)=0` would miss this continuing branch. This is a **synthetic Cramer polynomial chart**, not an additional inscribed-square existence theorem for any Jordan curve.

Variations prove exactly that a first-order numerator with a second-order denominator gives a pole, while a numerator with higher vanishing order gives limit zero. The ownership tests also distinguish a source-parameter limit `1` (incoming edge duplicate) and an event parameter crossing through `1` (a one-sided source exit).

### B. An irrational removable event and genuine poles

Set `alpha=√2-1`, uniquely isolated in `(0,1)` by `g(t)=t²+2t-1`. The quotient `g(t)^2(2/7)/g(t)^2` has a **removable order-two** singularity with exact limit `2/7`. The quotient `g(t)/g(t)^2` has a **true simple pole**. All signs are determined by exact rational polynomial calculations at this irrational `alpha`.

### C. The actual convex moving-polygon control

For v223–v228's single strictly convex quadrilateral motion at normalized event `t=8/9`, v228 found **28** charts with singular square matrices. An independent determinant/numerator valuation audit of those 28 charts gives:

| Local matrix/chart behavior | Charts | What it implies |
|---|---:|---|
| Nonidentically singular, genuine Cramer pole | **20** | The unique nearby chart solution exits the bounded source box |
| Nonidentically singular, finite removable parameters | **4** | The limits correspond to degenerate/non-positive-size source configurations |
| Identically singular matrix | **4** | Requires v228 slice feasibility and v225 rank support, not Cramer cancellation |
| **Total** | **28** | Zero previously hidden positive-size event squares, consistent with v228 |

This classification is **chart-level**, not a count of distinct geometric square orbits. In particular, these 28 charts do not represent 28 physical events. The geometric pair birth and continuing square found by v227 occur in the supported nonsingular charts.

The full-repository v229 integration test verifies this chart-count decomposition and checks v228's independent finding that **no singular chart has a square above cutoff `1/100`** at this event.

### D. The irrational scaling cutoff and the rectangle control

At the v227 irrational scaling cutoff `alpha=√2-1`, all **four** singular charts have identically zero determinants as polynomial functions of time, so the finite-Cramer cancellation rule is not invoked. v228 confirms no singular chart contains a square **strictly above** cutoff `144/25` at that instant. The valid cutoff crossing remains in the nonsingular v227 orbit-germ calculation.

For the fixed `2×1` rectangle, the `((0,0),(2,0),(2,1),(0,1))` chart selecting edges `(0,0,2,2)` with orientation `-1` has a **positive-dimensional square family**. Version v229 leaves this chart persistent-singular and defers to v228's family certificate. It does **not** turn an infinite continuum into one alleged geometric orbit.

## 6. Reproduction and limitations

Code: [`code/singular_continuation_v229.py`](../code/singular_continuation_v229.py)  
Tests: [`tests/test_singular_continuation_v229.py`](../tests/test_singular_continuation_v229.py)

```bash
python -m unittest discover -s tests -p 'test_singular_continuation_v229.py' -v
python code/singular_continuation_v229.py
```

The module is **standard-library only** (`fractions.Fraction`) and works independently for synthetic polynomial rank drops. Its full-polygon bridge optionally imports the existing v223–v228 modules for an exact event-slice consistency check. An optional test uses SymPy as an **independent exact development oracle** for seeded irrational polynomial signs; SymPy is **not required** to run the code or CI.

**Test status before PR:** 40 new tests, with 35 passing locally and 5 integration checks skipped in the standalone working directory. GitHub Actions runs those integrations in the full repository checkout, along with the complete archived research verification.

### Remaining work

1. **Orbit-level graph composition:** Deduplicate removable singular chart descriptions against nonsingular event germs *by the actual geometric square*; record instant-only and positive-dimensional strata in a global event incidence complex. The present result is deliberately **per-chart**.
2. **Persistently consistent singular germs:** For `D ≡ 0` and equal generic augmented/coefficient ranks, one needs a stratified feasibility argument over a time neighborhood, not just a single Cramer quotient.
3. **Exceptional irrational event fibers:** An instant-only solution at irrational `alpha` may require an exact algebraic orbit key, rather than identifying it with approximate rational snapshots.
4. **Most important:** Even a complete polygonal cobordism does not force a positive-sized square to survive arbitrary `C⁰` Jordan limits. The v220 infinite-notch countermodel prevents any universal noncollapse assertion based solely on rectifiability or parameter-gap compactification.

**The unrestricted Toeplitz square-peg conjecture remains open.** The proofs here concern local polynomial charts and certified finite models, not a universal Jordan-curve invariant. Existing configuration-space precedents include Cantarella–Denne–McCleary (arXiv:2103.07506) and Matschke (arXiv:1001.0186); neither is claimed to establish this remaining general noncollapse step.
