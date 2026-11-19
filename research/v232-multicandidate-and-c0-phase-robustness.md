# v232 — Multi-candidate spiral discovery and an exact C⁰ phase-robust square theorem

**Research program:** Benjamin John Schulz  
**Status:** A rigorously checked **conditional local four-arc degree stability theorem** and bounded **numerical-only** candidate-search improvement. The quantitative results below require executable rational verification against the v231 source endpoints; no unrestricted Jordan-curve square-peg proof is claimed.

## 1. What v231 left open

v231 deliberately separated a numerical polygonal search, Newton refinement, and a reproducible **rational four-arc Brouwer-degree certificate on the original smooth spiral**. It successfully certified three non-Liouville-regular closed Jordan double spirals: `D_2`, `D_(5/2)` and `D_3`. The discovered four-arc source boxes were stored as eight dyadic-square rational radial endpoints. Their rational proof transcripts establish a strict separation between an affine zero's source-boundary margin and the uniform chord-to-smooth-curve perturbation.

One numerical weakness remained: the search returned a single first suitable polygonal chart. A more finely sampled polygon may have many near-squares; the first numerically acceptable root can be poorly placed for the **local degree certificate** even if another positive-margin root exists nearby. Such a failed search is **inconclusive**, never evidence that the target curve lacks an inscribed square.

The deeper opportunity is analytic: once **one** true four-arc degree is known, it should persist under **arbitrary continuous perturbations** of those arcs, including oscillatory perturbations for which no derivative or finite-length bound is available.

## 2. Theorem: exact phase perturbation neighborhood without derivative assumptions

For `p` in `{2,5/2,3}` let

\[
\gamma_p(t)=t e^{it^{-p}}, \quad t\in(0,1],
\]

and let `I_0,...,I_3` be the four rational radial source intervals **actually certified** in the v231 `STORED_WITNESSES` manifest. Fix a *real continuous* phase function `phi:(0,1]->R`. Define

\[
\widetilde\gamma_{p,\phi}(t)=t\exp\!\left(i[t^{-p}+\phi(t)]\right).
\]

Suppose `|phi(t)| <= eta` for every `t` in the four source intervals. (The same bound on the entire `(0,1]` is sufficient but not necessary.) The elementary chord-length inequality on the unit circle gives

\[
\left|\widetilde\gamma_{p,\phi}(t)-\gamma_p(t)\right|
 = t|e^{i\phi(t)}-1|\le t|\phi(t)|\le\eta.
\]

In particular, each **coordinate** changes by at most `eta`. Let `L_i` be the v231 exact rational chords; their verified coordinatewise deviation from the unperturbed smooth arcs is at most `E_p`. It follows **uniformly in all source parameters** that

\[
\|\widetilde\Gamma_i(u)-L_i(u)\|_\infty\le E_p+\eta.
\]

Let `A`, `u*`, `m`, `M=||A^-1||_infty` denote the same rational matrix, rational root, box margin and inverse norm as v231. Then

\[
\boxed{4(E_p+\eta)<m/M}
\]

implies that the perturbed four-arc square residual has nonzero Brouwer degree **with the same sign** as v231. Therefore it contains a genuine square, without any smallness assumption on `phi'` (which need not exist). The original strict squared side and radius lower bounds remain valid whenever the **stronger** exact inequalities

\[
 |L_1(u^*_1)-L_0(u^*_0)|^2>(s_p+4(E_p+\eta))^2,
\]
\[
 |L_i(u^*_i)|^2>(r_p+2(E_p+\eta))^2 \quad(i=0,...,3)
\]

hold. All quantities to be compared are **rational** because v231's endpoint phases have explicit rational enclosures.

**Executable claim for v232:** these comparisons hold for the common rational allowance

\[
\boxed{\eta=10^{-6}\ \text{radians}}
\]

for all three witness sets, retaining the original v231 floors:

| `p` | Inscribed square side `>` | All four radii `>` | Phase condition on the four source intervals |
|---:|---:|---:|---|
| `2` | `1/4` | `6/25` | `|phi|<=10^-6` |
| `5/2` | `3/5` | `3/10` | `|phi|<=10^-6` |
| `3` | `7/10` | `3/8` | `|phi|<=10^-6` |

The machine **recomputes** all original v231 rational chord endpoints, four-by-four matrix data, exact inverse norm, chord-error bounds and strict side/radius floors using the original eight integers. It does **not** rely on rounded diagnostic margins, a new numerical root finder, or the historical v100 Miranda box.

### Jordan closure and roughness

For any continuous `phi` on `(0,1]` with bounded values near zero, set `\widetilde\gamma(0)=0`. Strict radial modulus `|\widetilde\gamma(t)|=t` proves that each arm is embedded. Its negative copy meets it only at the center; an outer semicircle joining the phase-shifted unit-radius endpoints closes an embedded Jordan curve. The certified square remains on the first arm.

The degree argument does not require the phase perturbation to be differentiable, rectifiable, or of finite `q`-variation. **However, the perturbation theorem alone does not classify the modified curve as non-Liouville regular:** special oscillatory cancellations can change area-primitive behavior, and that question must be checked separately for each family. Nor does a fixed local neighborhood settle the existence of a square on every Jordan curve.

## 3. Bounded search: from a single first guess to multiple proposals

`scan_candidates` uses v231's sampled power spiral, spatial segment grid, diagonal midpoint geometry and affine polygonal square solver. Instead of returning the first chart above a threshold, it now:

1. examines at most `max_pairs` diagonal-segment pairs;
2. checks the small set of locally nearby edges offered by the spatial grid;
3. filters trivially unsuitable numerical square size and minimum radius;
4. deduplicates chart presentations using a **coarse radial signature for search only**;
5. retains multiple numerically distinct candidate roots ranked by the available interior-parameter margin times square size, with a strict storage cap;
6. tries multiple angular half-widths for each candidate with `certify_on_smooth_spiral`, recording why attempts fail;
7. returns **inconclusive** if all proposals or the explicit evaluation budget are exhausted.

Search de-duplication is **heuristic**. Two close roots can collide in the coarse signature. This is not used for geometric-orbit counting, square existence, or any proof. The only affirmative mathematical outcome is the unchanged v231 exact rational **Brouwer-degree certificate**.

This is more robust than first-hit return and has explicit memory, pair, and proof-attempt budgets, but retains a **worst-case quadratic** diagonal-pair loop. It is not a certified exhaustive square search, and increasing the resolution does not guarantee finding a certificate within fixed budgets.

## 4. Tests and falsification controls

`tests/test_spiral_multicandidate_v232.py` checks exact phase margins for all three exponents, monotonic increase of the perturbation bound, invalid/oversized phase allowances, proof independence from the numerical scanner, capped spatial search, approximate signature invariance under relabeling, explicit failure ledgers, false first candidate followed by a successful second, exception handling and exhaustion behavior. The full GitHub CI repeats previous v212–v231 tests and the v232 rational report.

**Crucial negative control:** setting `eta=1/100` must **not** return a certificate from these particular v231 boxes. That is a limitation of this *certificate*, not proof that such a perturbed curve lacks squares.

The v220 accumulating-notch and v210 false localization examples remain crucial for future attempted *universal* scale-protection claims; the current theorem is intentionally specific to four verified arcs and makes no inference about their global non-collapse behavior under arbitrary `C0` limits.

## 5. Next genuine research test

Apply this **same** numerical proposal and exact a posteriori degree checker to new parametric irregular families not merely small bounded-phase perturbations of power spirals. That will require certified chord-error envelopes for each new class. A good first stress case is a scale-dependent phase modulation with no uniform global phase bound but a candidate source window on which the local deviation is quantified. Then compare against independent Liouville-regularity and variation criteria.

**Scientific status:** A locally robust square-existence theorem and a more resilient candidate-search strategy. No general positive-scale invariant or unrestricted Toeplitz proof.
