# v231: Blind spiral-square discovery followed by an exact positive-scale degree certificate

**Research program:** Benjamin John Schulz  
**Status:** An a posteriori, mathematically independent **four-arc Brouwer-degree existence certificate** for three specifically identified non-Liouville-regular Jordan double spirals. The numerical search supplies candidate chart locations only. **Every final certificate is recomputed using rational inequalities on the actual smooth spiral, not just its polygonal approximation.** No proof of unrestricted Toeplitz is claimed.

## 1. What changes relative to v230

v230 reconstructed the previous v99–v101 Poincaré–Miranda certificates and proved a uniform **conditional** square side floor `>1/8` for all power exponents `p>=2`. That proof inherited a specific old two-phase Miranda box. Version v231 deliberately **does not use the old box, any supplied square radius or orientation, or a stored root during numerical discovery**. Instead it:

1. Takes the smooth power spiral `gamma_p(t)=t exp(i t^(-p))` and samples its first arm at increasing angular resolution, with rational radial source parameters.
2. Puts the sampled **polygonal segments** in a spatial grid. For every admissible pair of diagonal segments, it uses the square-diagonal formula to query the two other corner locations, then solves the complete four-parameter affine square equations on a small set of candidate segment tuples. The search has an explicit pair budget; it can return **inconclusive**.
3. Refines a proposed polygonal square to a numerical square of the **smooth** spiral using its analytic derivative. **The Newton computation is discovery only, not a proof.**
4. Chooses four small *rational radial parameter intervals* around the numerical candidate, independently reconstructs certified enclosures for their eight true transcendental endpoints, and proves that the exact rational four-chord affine square residual has **nonzero Brouwer degree that cannot change** under interpolation from the polygon to the smooth spiral.
5. Replays the final proof from eight immutable rational source endpoints, without depending on floating-point search or Newton at all. This is the authoritative certificate. The blind discovery is rerun separately as an algorithmic regression.

The first arm is a subset of the closed Jordan double spiral `D_p`; therefore a square on this arm is also inscribed in `D_p`. Unlike an arbitrary global polygon-parity argument, the local four-arc degree gives an independent quantitative non-collapse certificate for these instances.

## 2. General four-arc perturbation theorem (not power-specific)

Let `Gamma_i:[0,1]->R^2` be four continuous parametrized planar source arcs. Let `L_i(u)=(1-u)a_i+u b_i` be their straight-chord approximations, with exact rational endpoints. Fix a prescribed orientation `sigma in {+1,-1}` and the quarter turn `J(x,y)=(-y,x)`. Define the four-component square residual

\[
F_\sigma(P_0,P_1,P_2,P_3)=
\left(P_2-P_1-\sigma J(P_1-P_0),\;P_3-P_0-\sigma J(P_1-P_0)\right).
\]

Since the residual is affine in the four chord parameters, write `F_L(u)=A(u-u_*)` whenever `A` is invertible and its uniquely determined rational zero `u_*` lies strictly inside the cube `(0,1)^4`. Define

\[
m=\min_{0\le i<4}\min(u_{*,i},1-u_{*,i}),\qquad
M=\|A^{-1}\|_\infty.
\]

Suppose each arc and its chord satisfy the uniform **coordinatewise** bound

\[
\|\Gamma_i(u)-L_i(u)\|_\infty\le E
\quad\text{for every }u\in[0,1],\quad i=0,1,2,3.
\]

Each scalar component of `F_sigma` contains at most four signed coordinates of the four points, so

\[
\|F_\Gamma(u)-F_L(u)\|_\infty\le 4E.
\]

**Theorem 1 (portable a posteriori degree certificate).** If

\[
\boxed{4E<\frac{m}{M}},
\]

then `F_Gamma` has at least one zero in `(0,1)^4`. For each chord, set \(D_i=\max(u_{*,i},1-u_{*,i})\|L_i(1)-L_i(0)\|_1+2E\). If the affine root's squared side and point radii satisfy

\[
s_L^2>(s_0+D_0+D_1)^2,\qquad
|L_i(u_{*,i})|^2>(r_0+D_i)^2,
\]

then the true smooth-curve square has **side length strictly greater than `s_0`** and **all four point radii strictly greater than `r_0`**.

**Proof.** On the boundary of the source cube, `||u-u_*||_infinity >= m`. Thus

\[
\|F_L(u)\|_\infty=\|A(u-u_*)\|_\infty\ge m/M.
\]

The straight-line homotopy `F_s=(1-s)F_L+sF_Gamma` has no zero on the source boundary because its perturbation is bounded by `4E`. Its Brouwer degree therefore agrees with that of the invertible affine map `F_L`, namely `sgn det A` in `{+1,-1}`. Nonzero degree supplies a zero on the true curve. The true and affine roots may have different source parameters. At the true zero each corner is within \(2E\) of its chord point and that chord point is at most \(D_i-2E\) from the affine root's corner. Triangle inequalities give the strict floors. For the power spirals one can instead verify the radius floor directly from the minimum source radius. A positive side floor makes the square nondegenerate. QED.

**Important:** The existence proof does not rely on a polygonal square *converging* under arbitrarily many refinements. It uses the topological degree of one **specified four-dimensional source box** and a validated perturbation margin. This method is applicable to other curves whenever four local source arcs admit suitable certified chord errors.

## 3. Rigorous smooth power-spiral error bounds

For `gamma_p(t)=t e^(i t^(-p))`, direct differentiation yields

\[
\gamma_p''(t)=e^{it^{-p}}
\left[-p^2 t^{-2p-1}+i p(p-1)t^{-p-1}\right].
\]

On any rational radial interval `[a,b] subset (0,1]` the Euclidean deviation from linear chord interpolation is bounded by

\[
\frac{(b-a)^2}{8}\sup_{a\le t\le b}|\gamma_p''(t)|
\le\frac{(b-a)^2}{8}\left[p^2 a^{-\lceil 2p+1\rceil}
+|p(p-1)|a^{-\lceil p+1\rceil}\right].
\]

The ceiling powers are conservative because `0<a<=1`, ensuring the **entire interpolation bound is rational** for rational `p` and rational radial endpoints.

For `p in {2,5/2,3}`, set each endpoint `t=(u/2^30)^2`, where `u` is a positive integer. Then the exact spiral phase is `t^(-p)=(u/2^30)^(-2p)`, a **rational number**, despite the transcendental `sin/cos` values. To enclose those endpoint coordinates using rational arithmetic only:

* Apply the exact Machin formula `pi=16 atan(1/5)-4 atan(1/239)`. Alternating arctangent partial sums provide rigorous rational lower/upper bounds.
* Reduce the rational phase by an arbitrary integer multiple of twice the **rational midpoint** of this pi interval. The error from using the midpoint is at most `|k|` times the pi interval width.
* Round the reduced rational argument to a fixed dyadic. Its rounding discrepancy is computed exactly, not inferred from floating point.
* Evaluate rational sine/cosine Taylor polynomials through degree 41/40. Taylor's theorem bounds the remainder by `|x|^42/42!` and `|x|^41/41!`, respectively; `|sin'|,|cos'|<=1` transfers the argument uncertainty. Quantizing the resulting endpoint coordinates to dyadic rationals adds a **known exact rounding error**.

The maximum of the eight endpoint errors plus the four chord interpolation errors gives `E`. **All proof comparisons are `fractions.Fraction` comparisons.** Floating point is used **only** to select the angular-reduction integer (any integer is mathematically valid), display diagnostics, and discover the candidate. It does not enter the bound itself.

## 4. Exact certified outputs

The eight rational source endpoint values are stored as integers `u_i` with denominator `2^30` in `STORED_WITNESSES` in the code. Separate test replays work even when the numerical search is disabled and no candidate square is supplied. All three certificates satisfy the **same, exact rational separator**

\[
\boxed{4E<\frac{1}{20000}<\frac{m}{\|A^{-1}\|_\infty}}.
\]

This is a strict rational inequality checked by the proof program, *not* a floating-point observation. The resulting bounds are stronger than v230's common `>1/8` bound:

| `p` | Certified side length `>` | Every vertex radius `>` | Diagnostic polygon side (approx.) | Diagnostic `4E` (approx.) | Diagnostic boundary margin (approx.) |
|---:|---:|---:|---:|---:|---:|
| `2` | `1/4` | `6/25` | `0.29143` | `0.00003414` | `0.00005996` |
| `5/2` | `3/5` | `3/10` | `0.66113` | `0.00003953` | `0.00005404` |
| `3` | `7/10` | `3/8` | `0.77037` | `0.00004493` | `0.00006473` |

The decimal columns are **informational only**. The asserted side/radius bounds and degree separator are proved by exact rational inequalities in `replay_stored_witness` and checked again independently by the corresponding test cases. The polygonal root itself is the exact rational solution of an invertible `4x4` system, verified to satisfy the square equations identically.

A successful certificate guarantees an actual square on the **smooth** `gamma_p`, and hence on the corresponding **non-Liouville-regular Jordan curve** `D_p`. No prior v100/101 face-sign certificate is invoked for these three instances.

## 5. Reproducibility, failure controls, and performance

Run:

```bash
python -m unittest discover -s tests -p 'test_spiral_pruned_v231.py' -v
python code/spiral_pruned_v231.py
```

The CLI independently replays all three **fixed rational proof witnesses**, then runs a blind numerical search from a common outer-to-inner radius ladder `(.55,.45,.36,.30,.22)`, using a uniform angular sample spacing around `0.22` and a spatial segment grid. The search never reads the stored witness table and never imports the v100/v101 phase-box code. Its numerical Newton phase merely proposes rational radial source boxes that are then checked from scratch.

In the development run it found the three certifiable roots after examining **736**, **807**, and **1121** diagonal pairs, respectively, with **90**, **88**, and **93** sampled first-arm polygon segments at the successful radius. Search counts are implementation-dependent diagnostics, **not** complexity guarantees or part of the mathematical certificate. The current discovery search still has worst-case quadratic work in the number of sampled diagonal segments, but spatial lookups and a hard `max_pairs` cap prevent unbounded runs. No claim of scalable all-square enumeration on very large spirals is made.

The 22 new regression tests cover:

* Exact rational pi, trigonometric endpoint enclosures and analytic interpolation error;
* Exact rational affine-root and inverse-matrix calculations, and no-boundary-zero proof inequalities;
* Stronger side/radius floors and replay from immutable exact source values;
* Sign reversal, impossible side floors, invalid exponents and uncertified source charts (must **not** return a positive proof);
* Reproducible blind numerical discovery for all three exponents, independently followed by exact certification.

Optional float/independent numerical sanity tests are **not** relied on for a theorem. The SHA256/provenance and original Arb workflow from v230 remain unaffected.

## 6. What v231 proves and what it does not

The new **general conditional degree criterion** is portable to other four-arc configurations. The three **special-case certificates** are independent of the old power-spiral Miranda box and, unlike a mere finite-polygon square, actually supply squares on the original `D_2`, `D_(5/2)`, and `D_3`. All three lie outside the rectifiable, finite-variation-under-two and Liouville-regular existence classes, by v230 and Asano–Ike [arXiv:2609.11598v3](https://arxiv.org/html/2609.11598v3).

**Not established:** a universal algorithm that always finds a protected square on any arbitrary nonrectifiable Jordan curve; a general positive-scale bound for arbitrary irregular curves; full geometric-orbit continuation through infinite accumulation; or the unrestricted Toeplitz conjecture. The v220 notch and v210 physical-lift/localization controls remain mandatory negative cases for future broad claims.

**Next research target:** develop adaptive spatial/pruned source-box selection with a verified local degree margin for a wider class of rough-curve segments. Study stability under perturbations of the phase profile (v102), and test whether the method can construct square-protecting boxes on irregular families *not already known* to admit an explicit power-spiral square.
