# v230 — An irregular Jordan benchmark with a uniform, explicitly protected square

**Research author:** Benjamin John Schulz  
**Status (19 November 2026):** Elementary analytic conclusions proved below; square-existence and quantitative size statements **conditional on independent correctness of the recovered v100/v101 computer-assisted spiral theorem**. No unrestricted Toeplitz theorem is claimed.

> **Audit corrections (19 November 2026, after PR #19).** Changes are marked *[audit]* in the text.
> 1. **Attribution.** The curve `D_p` and the equivalence "arm Liouville regular ⇔ closed curve Liouville regular" are Asano–Ike's construction (arXiv:2609.11598v3, §6.1). What is new here is the **square** on that curve for `p≥2`; Asano–Ike do not discuss inscription for these spirals.
> 2. **Truncation caveat.** Asano–Ike cut their spiral at angle `θ≥e`; this note cuts at `t=1` (`θ≥1`). For large `p` the Box-A square uses the piece they cut off (numerically from `p≈23.75`), so say "`D_p`", not "Asano–Ike's curve".
> 3. **Stronger, fully stated bounds.** Theorem 2 now gives side `>1/4` (was `>1/8`) and states the outer-vertex bound `radius<1`, citing v100 §9.
> 4. **Arb rerun.** The eight v101 face certificates were rerun locally with python-flint 0.9.0 at 160 bits: all certified, zero unresolved cells. The three Theorem 2 inequalities are now also Arb-certified on all of Box A (`code/spiral_box_a_audit_v230.py`).
> 5. **iv101 is not outward-rounded.** `Decimal.sqrt()` ignores directed rounding, and the atanh tail bound rounds the wrong way for negative arguments. The recovered file is kept byte-identical; `related/power-spiral-v99-v102/iv230.py` is the corrected copy. The Arb proof path never used iv101.

## 1. Strategic correction

The recent v223–v229 exact polygonal-event computations are useful algorithms but do not enlarge the classes of curves for which square inscription is known. Pak's earlier polygon results and Asano–Ike's rectangular-peg theorem already cover all polygons and rectifiable Jordan curves. The challenging remaining class contains irregular **non-Liouville-regular** curves, including the power spirals studied in Asano–Ike, *Liouville regular Jordan curves*, arXiv:2609.11598v3, Section 6.1 and Theorem 6.3. This batch shifts the primary benchmark to that class.

## 2. Exact double-spiral closure

Fix any real `p>0`. Set

\[
\gamma_p(t)=t\exp(i t^{-p}),\qquad 0<t\le1,\quad\gamma_p(0)=0,
\]

and define the closed planar set

\[
\boxed{D_p=\gamma_p([0,1])\ \cup\ [-\gamma_p([0,1])]\ \cup\ \{e^{i\theta}:1\le\theta\le1+\pi\}.}
\]

*[audit] Attribution.* This is the closed curve `C = γ ∪ (−γ) ∪ η` of Asano–Ike, arXiv:2609.11598v3, §6.1, built from their spiral `γ_{α,β}` with radius `θ^{−α}(log θ)^{−β}` at `α=1/p`, `β=0`, reparametrized by `t=θ^{−1/p}`, with `η` a semicircle joining the outer endpoint to its negative. The only difference is the truncation: they keep `θ≥e` (that is, `t≤e^{−1/p}`), and `D_p` keeps `θ≥1` (`t≤1`). Lemma 1 is their observation, recorded here for completeness.

**Lemma 1 (Jordan embedding; Asano–Ike §6.1).** `D_p` is a Jordan curve.

**Proof.** On either arm, modulus equals the parameter `t`, so its parametrization is injective away from `t=0`. Opposite arms have the same radius only when they have the same `t`, at which their points are exact negatives; hence their sole common point is the origin. Their outer endpoints are `e^{i}` and `-e^{i}`. The unit semicircle meets the two arms only at those endpoints because all other points on the arms have modulus strictly below one. Traversing the first arm outward, the semicircle, and the second arm inward gives a continuous one-to-one image of a circle. QED.

The closing arc is **not** an approximation: the entire first arm `gamma_p((0,1])` remains as a subset of `D_p`.

## 3. Exact regularity thresholds

### 3.1 Length

\[
|\gamma_p'(t)|=\sqrt{1+p^2t^{-2p}}.
\]

The inequalities `p t^{-p} <= |gamma'| <= 1+p t^{-p}` prove directly that the arm has finite length exactly when `p<1`. The other arm and closing semicircle add only two finite pieces. Thus

\[
\boxed{D_p\text{ is rectifiable iff }0<p<1.}
\]

### 3.2 Finite q-variation

For any variation exponent `q>=1`, the compactified curve has finite `q`-variation **exactly when** `p<1` (any `q>=1` works) or `q>p`. Here is an independent elementary proof.

**Necessity when `p>=1`.** Parameterize by angle `theta=t^{-p}` and choose `theta_n=n pi`. Consecutive spiral samples have opposite angular direction. Their distance is at least a fixed positive multiple of `n^{-1/p}`. The q-variation sum of these successive increments therefore dominates `sum n^{-q/p}`, which diverges for `q<=p`.

**Sufficiency when `p>=1, q>p`.** In the dyadic angular annulus `2^k<=theta<=2^{k+1}`, the arm has amplitude at most `R_k=C 2^{-k/p}` and total arclength bounded by `L_k=C_p 2^{k(1-1/p)}`. For any partition increments `Delta_j` within this annulus,

\[
\sum_j|\Delta_j|^q\le(2R_k)^{q-1}\sum_j|\Delta_j|
\le(2R_k)^{q-1} L_k
\le C_{p,q}\,2^{k(1-q/p)}.
\]

The **q-variation seminorm** on this annulus is at most `C' 2^{k(1/q-1/p)}`. Concatenation and Minkowski's inequality bound the seminorm of the entire spiral by the sum over dyadic annuli. This is a convergent geometric series if `q>p`. The compactifying endpoint and the other two finite pieces do not spoil it. For `p<1`, finite length implies finite q-variation for every `q>=1`.

Hence:

\[
\boxed{\exists q<2\text{ with finite }q\text{-variation of }D_p\iff p<2.}
\]

### 3.3 Liouville regularity

Write `x=t cos(t^{-p})`, `y=t sin(t^{-p})` and `lambda=y dx`. Direct differentiation yields the **exact identity**

\[
y\,dx=d\!\left(\frac{t^2}{4}\sin(2t^{-p})\right)
+\frac p2t^{1-p}\,dt.
\]

The total-derivative term tends to zero at the origin. The improper primitive therefore converges precisely when `p<2`. At `p=2` it contains `log t`, and for `p>2` it contains `(p/[2(2-p)])t^{2-p}`, both divergent as `t->0`.

Asano–Ike, arXiv:2609.11598v3, Theorem 6.3, proves that the compactified corresponding power spiral is **Liouville regular iff `p<2`** (their angular exponent `alpha=1/p`, logarithmic exponent `beta=0`). Their Lemma 6.2 also proves that any Liouville primitive must agree with the smooth subarc integral. Since Liouville regularity restricts to subarcs, the divergence proves that the **closed** curve `D_p` is non-Liouville regular for `p>=2`. Conversely for `p<2`, both individual arms are Liouville-regular arcs by Theorem 6.3 and the semicircle is smooth; their finite concatenation is Liouville regular by their Proposition 2.7. *[audit]* Asano–Ike state this equivalence for the closed curve directly in §6.1 ("the arc γ is Liouville regular if and only if C is"), proving one direction from their Lemma 2.2 and the other from Theorem 5.2, Corollary 3.5 and Proposition 2.7; the argument above follows theirs. Therefore

\[
\boxed{D_p\text{ is Liouville regular iff }p<2.}
\]

The opposite arms can cancel *total* oriented-area integrals, but that does not restore the required **continuous local primitive on each arm**.

## 4. Previously certified power-spiral squares, and the NEW uniform margin

The saved project handoffs `v99`, `v100`, `v101` claim and document a computer-assisted Poincaré–Miranda theorem that **every** `gamma_p((0,1])`, `p>0`, contains a nondegenerate inscribed square (radial order `0132`). Their original final interval verifier is the **v101 Arb** code, *not* the defective v100 Decimal layer. The corresponding sources are restored under `related/power-spiral-v99-v102/`. This report does **not** independently rerun the Arb proof and does not elevate old source claims beyond their checked computational provenance. *[audit] The Arb proof has since been rerun: `verify_v230_arb_strict.py` with python-flint 0.9.0 certifies all eight faces at 160 bits with zero unresolved cells. That reruns the computation; it is not an audit of the v100 analytic reduction.*

**Conditional Theorem 2 (new explicit uniform non-collapse bound for `p>=2`).** Assume the v101-certified face inequalities are valid and supply the Box-A zero specified in the v100 handoff. Then **for every real `p>=2`, the curve `D_p` contains an inscribed square whose four vertices have radius strictly greater than `1/5` and strictly less than `1`, and whose side length is strictly greater than `1/4`.** *[audit: the original statement had side `>1/8` and left the outer radius bound implicit.]* The same square remains on each truncated-and-reclosed double spiral obtained by cutting the two arms at any `0<epsilon<=1/5` and connecting their inner endpoints by a radius-`epsilon` semicircle.

**Proof.** In v100 §5, Box A uses `q=1/(1+p)` in `(0,1/3]`, `c∈[1,19/10]` and `d∈[9/20,18/25]`. At a face-certified zero, two *opposite square vertices* are normalized to `z0=1`, `z2=-rho exp(i delta)`, where `rho=exp(cq)>1`, `delta=dq<6/25`, and the desired lifted phase drop is `Delta=3pi-dq`. The physical first radius `r0` satisfies

\[
r_0^p=\frac{1-\rho^{-p}}{\Delta}
=\frac{1-\exp(-c(1-q))}{3\pi-dq}.
\]

For `p>=2`, `c(1-q)>=2/3`. The elementary `exp(x)>=1+x` gives `1-exp(-x)>=x/(1+x)` for `x>=0`, hence the numerator is at least `2/5`. By `pi<22/7`, the denominator is less than `10`. Consequently `r0^p>1/25`, so

\[
\boxed{r_0>25^{-1/p}\ge1/5.}
\]

v100 §4 proves that the other two forced square radii fall strictly between `r0` and `rho r0`. Thus **all four vertices** have radius greater than `1/5`.

The square side length `s` is the diagonal length divided by `sqrt(2)`:

\[
s=\frac{r_0}{\sqrt2}\,|1+\rho e^{i\delta}|.
\]

*[audit: replaces the original `>1/8` step.]* Here `rho>=1` and `0<delta<=6/25<pi/2`, so `cos(delta)>=1-delta^2/2>0`. The function `rho↦1+rho^2+2rho cos(delta)` is increasing for `rho>=0`, so

\[
|1+\rho e^{i\delta}|^2=1+\rho^2+2\rho\cos\delta\ \ge\ 2+2\cos\delta\ \ge\ 4-\delta^2\ \ge\ 4-\left(\tfrac{6}{25}\right)^2=\tfrac{2464}{625}>\tfrac{25}{8}.
\]

With `r0^2>1/25` this gives `s^2 = r0^2|1+rho e^{i delta}|^2/2 > (1/25)(25/8)/2 = 1/16`, so

\[
\boxed{s>1/4.}
\]

(Even the crude step `|1+rho e^{i delta}|^2>1+rho^2>2` already gives `s>r0>1/5`.)

*[audit] Outer radius (v100 §9).* The outermost vertex has radius `rho r0` (v100 §9 calls its parameter `t3`), and

\[
(\rho r_0)^p=\rho^p r_0^p=\frac{e^{c(1-q)}-1}{3\pi-dq}.
\]

Since `c(1-q)<=19/10<2` and `e<272/100`, the numerator is below `(272/100)^2-1<37/5-1=32/5`; since `pi>3`, the denominator exceeds `9-6/25=219/25`. So `(rho r0)^p<(32/5)/(219/25)<1`: **all four vertices have radius strictly between `1/5` and `1`**, so they lie on the first arm `gamma_p((0,1))` of `D_p`.

*[audit] Machine check.* The three inequalities used above — `r0^p>1/25`, `(rho r0)^p<1` and `|1+rho e^{i delta}|^2>25/8` — are also certified in Arb ball arithmetic on **every point of Box A** (not only at the unknown zero) by `code/spiral_box_a_audit_v230.py`. A separate mpmath continuation (numerical, not a proof) locates the actual square for `p=2,…,100`: at `p=2` its smallest vertex radius is `0.249` and its side is `0.443`.

For any `epsilon<=1/5`, remove the two interiors `0<t<epsilon` and connect the truncated endpoints by a semicircle on the circle `|z|=epsilon`, which does not intersect the retained spiral arms. The four certified vertices all lie on the retained first arm and remain unchanged. QED, **conditional solely on the old Box-A face existence certificate**.

**Why this is useful:** unlike generic polygon approximations, the test family comes with an *a priori and uniform* quantitative size/position floor even across every `p>=2`. This is specific to the power-law phase and the Box-A proof, not a general irregular Jordan-curve theorem. The lower bound is conservative and not a sharp size estimate.

**[audit] Truncation caveat.** The theorem is about `D_p`, cut at `t=1`. Asano–Ike's §6.1 curve keeps only `θ=t^{-p}>=e`. The outermost Box-A vertex sits at angle `θ = Δ/(e^{c(1-q)}-1)`, which tends to about `2.4` as `p→∞`, below `e`. Numerically (mpmath continuation, not certified) the square stays on their truncated curve for `p` up to about `23.75` and leaves it after that. Do not describe Theorem 2 as a statement about Asano–Ike's exact curve for large `p`.

## 5. Recovery status and falsification requirements

Recovered file-library originals:

- `handoff_v99_small_exponent_log_spiral.md`, `verify_v99_small_exponent_log_spiral.py`;
- `handoff_v100_power_spiral_all_exponents.md`, **historical** `verify_v100_power_spiral_all_exponents.py`;
- `handoff_v101_interval_layer_audit.md`, corrected interval library `iv101.py`, independent `verify_v101_arb.py`;
- `handoff_v102_radial_phase_robustness.md`.

The saved v101 handoff also describes `verify_v101_faces.py`, `test_iv101.py`, `probe_bugs.py`, and a captured output file. **Do not assert that those files are recovered or rerun unless they are actually retrieved.** v101 explicitly found an unsound negative-argument enclosure in v100's `Efun_scalar`. The historical v100 script is not an acceptable proof certificate. We include it solely to document the issue.

The independent v101 Arb verifier requires `python-flint`. The working container could not install this dependency because external package resolution was unavailable. Its GitHub Actions reproduction must pass independently before labeling the recovered computational claim *freshly reverified*. *[audit] It has now been rerun outside that container (python-flint 0.9.0): all eight faces certified, zero unresolved cells, exit status 0. The CI job also runs the new Box-A inequality certificate and the audit tests with missing dependencies treated as failures.*

*[audit] Interval-layer correction.* The recovered `iv101.py` is **not** fully outward-rounded. Python's `Decimal.sqrt()` always rounds half-even and ignores `ROUND_FLOOR`/`ROUND_CEILING`, so its interval square root returned a lower endpoint above the true root, or an upper endpoint below it, for about half of all inputs; this reaches `atan_pt` and the Machin π. Separately, the atanh tail bound rounds `term*t*t` toward `+∞` with `t<0`, which can undershoot. `iv101.py` stays byte-identical for provenance; `iv230.py` fixes both (square roots verified by exact squaring; the tail bound uses `|t|`) and is tested in `tests/test_iv230_interval_layer.py`. None of this affects the Arb proof, which never used `iv101`.

A full independent publication audit should still inspect the theorem's analytic derivation, phase branch conventions, endpoint containment, every v99 face and the `v101` Arb transcript. Neither this manuscript nor the unit tests are substitutes for that analytic audit.

## 6. Research impact and negative controls

- **Irregular benchmark:** `p=2, 5/2, 3` are nonrectifiable, have no finite `q<2` variation, and are non-Liouville regular. They nonetheless inherit the power-spiral witness, conditional on v101. Use them as positive controls for any proposed non-collapse certificate.
- **Intermediate controls:** `p=1, 3/2` are nonrectifiable but **Liouville regular**. Their squares fall within Asano–Ike's theorem, so they are not frontier examples.
- **Finite-variation controls:** `p<1` is rectifiable, already in the established rectangular-peg class.
- **v220 notch curve:** rectifiable with arbitrarily tiny notch squares, useful for falsifying claims of a size floor for **all** squares; it is not a benchmark outside existing existence theorems.
- **Asano–Ike finite-coil arc:** their §6.2 gives a differentiable but non-Liouville-regular example. A closed version and a square certificate would require **separate** geometry and must not be inferred from its name.
- **Polygonal reconstructions:** approximate the two spiral arms and connecting semicircles, but avoid claiming certified square locations from plotted near-squares. Prove exact segment filters and independent intervals before asserting a rigorous limit.

## 7. Primary sources

1. Asano & Ike, *Liouville regular Jordan curves*, 17 November 2026, arXiv:2609.11598v3, §§2 and 6.1, Theorem 6.3. https://arxiv.org/html/2609.11598v3
2. Asano & Ike, *The Rectifiable Rectangular Peg Problem*, arXiv:2412.21057v4 (12 November 2026, per the arXiv version listing). https://arxiv.org/abs/2412.21057
3. Li & Pan, *Rectangular Pegs on Jordan Curves of Finite p-Variation*, arXiv:2609.04918. https://arxiv.org/abs/2609.04918
4. Greene & Lobb, *Square pegs between two graphs*, Comment. Math. Helv. (2026). https://ems.press/journals/cmh/articles/14299791
5. Schulz, recovered original v99–v102 power-spiral mathematical handoffs and independent Arb verifier; see `related/power-spiral-v99-v102/` and original file-library provenance.

**Final scientific limit:** The unrestricted Toeplitz conjecture remains open. A uniform square-size floor on this **special, explicitly solvable family** does not force a square on arbitrary continuous Jordan curves.
