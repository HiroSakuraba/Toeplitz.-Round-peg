# v212 — Periodic inscription bridge: a sharp size lemma and a localization obstruction

**Research program:** Benjamin John Schulz  
**Developed:** 18 November 2026  
**Status:** Two elementary proved statements, one exact-rational adversarial model, and an explicit list of unresolved global implications. The Toeplitz conjecture remains open.

## 1. Why this batch exists

Our v209 local theorem forces a nondegenerate square from four *actual crossing arcs* and a strict matrix/face-margin condition. The v210 rotating-scaling spiral packet disproves the attempted universal extraction of those arcs from action concentration. The v211 shape-sphere compactification keeps geometric shape but allows physical scale to tend to zero; its tangent estimate does not rule out moving-center escape.

Two **October 2026 arXiv preprints**, Lorenzo Valerio [V] and Ali Naseri Sadr [S], prove strong existence results for *four periodic strands*, conditional on signed-area balance. We ask what those theorems actually give when transported toward a **single arbitrary planar Jordan curve**.

The results here deliberately do not assume the desired transfer is available.

## 2. Exact external theorem hypotheses: do not conflate them

| Item | Valerio, arXiv:2610.21076v2 [V] | Naseri Sadr, arXiv:2610.24376v1 [S] |
|---|---|---|
| Domain | Four essential oriented polygonal simple closed curves on `(R/Z) × R` | Four smooth periodic embeddings of `R` into `C`, with `γ_i(t+1)=γ_i(t)+i`; quotients lie on `R × S¹` |
| Original four-strand balance | `∫σ₁ y dx − ∫σ₂ y dx + ∫σ₃ y dx − ∫σ₄ y dx = 0` | `t A(γ₁,γ₂)+(1−t) A(γ₃,γ₄)=0`, `0<t<1`, where `A(γ_i,γ_j)=∫γ_i x dy−∫γ_j x dy` |
| Additional common-intersection exclusion | **Not assumed** for the main polygonal existence theorem | **Assumed**: all four lifted images have empty common intersection (periodic quotient interpretation must be tracked) |
| Inscribed object | Ordered periodic square, possibly degenerate | Ordered type `(t,θ)` isosceles trapezoid for every `0<θ<π`; at `t=1/2`, rectangles including a square |
| Rough curves | Main polygonal theorem; its proof uses smooth approximation. This does *not* assert an arbitrary-rough-curve transfer on one Jordan curve. | Theorem 1.4 treats continuous periodic embeddings of **zero planar Lebesgue measure**, with common-intersection exclusion and a well-defined balance via cylinder region indicators. |
| Implication for one planar Jordan curve | **Not immediate** | **Not immediate**, and the zero-area-image condition cannot be imposed universally on arbitrary Jordan curves. |

**Coordinate caveat.** [V] uses horizontal period and `y dx`; [S] uses vertical period and `x dy`. Their alternating area conventions must be mapped by an explicit orientation-preserving coordinate rotation and sign check, not copied symbol-by-symbol. The exact-rational model below uses [V]'s horizontal-period `∫ y dx` convention. Its `t=1/2` equal-weight form illustrates [S]'s square case after the corresponding rotation.

**Square convention.** The ordered points `z₀,z₁,z₂,z₃` form a square when `z₁−z₀ = e`, `z₂−z₁ = i e`, etc., with `e≠0`. In a periodic problem, choose real planar lifts of the cylindrical points; their periods need not be the same.

## 3. Lemma A — four-strand non-collapse, with quantitative stability

Let `C_L = (R/LZ) × R` with its flat quotient metric and period `L>0`. Let `K₀,…,K₃ ⊂ C_L` be nonempty compact sets, and assume

`K₀ ∩ K₁ ∩ K₂ ∩ K₃ = ∅`.

Define the *collective intersection gap*

`d(K₀,…,K₃) = min_{p∈K₀} max_{j=0,1,2,3} dist_C(p,K_j)`.

**Lemma A.** We have `d>0`. Every Euclidean square `(z₀,z₁,z₂,z₃)` in the universal covering plane with `π(z_j)∈K_j` satisfies

`s >= d / √2 > 0`,

where `s` is its side length. Furthermore, if `K'_j` are nonempty compact sets with `d_H(K_j,K'_j) <= ε` for all `j` (same cylinder metric), then

`|d(K₀,…,K₃)−d(K'₀,…,K'₃)| <= 2ε`.

**Proof (positive gap).** The function `p ↦ max_j dist_C(p,K_j)` is continuous on compact `K₀`. Its minimum exists. A zero would be a point belonging to all four closed sets, contrary to the hypothesis.

**Proof (size).** For any square, each labeled vertex is at Euclidean distance at most `√2 s` from `z₀`, since the diagonal is the greatest vertex-to-vertex separation. The quotient projection `π:R²→C_L` is 1-Lipschitz. Thus for `p=π(z₀)`, `dist_C(p,K_j) <= dist_C(π(z₀),π(z_j)) <= √2 s` for each `j`. Taking the maximum and then the minimum gives `d <= √2 s`.

**Proof (stability).** Let `p∈K₀` minimize the first gap and choose `p'∈K'₀` with `dist_C(p,p')<=ε`. For every `j`, the Hausdorff bound and triangle inequality give `dist_C(p',K'_j)<=dist_C(p,K_j)+2ε`. Hence `d(K')<=d(K)+2ε`. Interchange primed and unprimed sets.

**Sharpness.** For four horizontal unit-period graphs at heights `(1,0,-1,0)`, the first-strand gap is `d=2`: every first-strand point is distance 2 from the third strand. The lifted square with vertices `(x,1), (x−1,0), (x,−1), (x+1,0)` has side `√2`, so equality holds.

**Consequence for approximations.** If a sequence of balanced periodic systems has `d_n >= d_* >0` (or converges in Hausdorff distance to a compact four-strand system having no common intersection), the corresponding square witnesses **cannot shrink to zero**. No fixed `d_*` follows merely from the word “nonintersecting” at each approximation level; `d_n` could tend to zero.

### An essential qualification: gap at each stage is not a uniform gap

There is a separate exact degenerating model. For each integer `n>=1`,
use four horizontal unit-period strands at heights `(1/n,0,-1/n,0)`.
They are embedded, satisfy the alternating area equality, and have
empty fourfold intersection for every individual `n`. Their collective
gap is `d_n=2/n`, while the square

`(x,1/n), (x−1/n,0), (x,−1/n), (x+1/n,0)`

has side `√2/n`. Thus both the geometric lower bound and the square
side collapse to zero. **Per-approximation nonintersection is not a
uniform non-collapse theorem.** A usable bridge must prove `inf_n d_n>0`
or an independent size invariant. The exact checker tests this family.

## 4. Proposition B — balanced periodicization does not localize witnesses

We construct four *piecewise-linear graph* embeddings in the horizontal-period cylinder (so the hypotheses of [V] hold) which:

1. obey exact alternating area balance;
2. have no fourfold common intersection, even with a uniform positive vertical gap;
3. jointly inscribe an explicit, nondegenerate square of side `√2`;
4. admit four specified small 'original-piece windows', **with no labeled square having all four vertices in those windows**.

This refutes any general inference of the form “area balance + common-intersection gap + square on completed periodic curves implies square on the four selected local pieces”. It does **not** refute [V] or [S], and the selected windows are not claimed to come from one Jordan curve.

### Construction, entirely over rational numbers

Let the periodic coordinate be `x mod 1` and set

`(h₀,h₁,h₂,h₃)=(1,0,−1,0)`.

For `i=0,1,2,3`, define a triangular hat `H_{c,w}(x)=max(0, 1−|x−c|/w)` on one period, periodically extended, and put

`c_i^+=(2i+1)/16,  c_i^-=(2i+2)/16,  w=1/128,  a=1/16`,

`f_i(x)=h_i+a [H_{c_i^+,w}(x)−H_{c_i^-,w}(x)]`.

All supports lie in `(0, 1/2+w)`; each `f_i` extends to a continuous, 1-periodic, piecewise-linear function. The entire lifted curve `γ_i(t)=(t,f_i(t mod 1))` is an embedding of `R` into the plane and has period vector `(1,0)`.

Since each positive and negative triangle has identical area `w`, exact integration gives

`∫₀¹ f_i(x) dx=h_i`, so `∑(-1)^i∫₀¹ f_i dx=1−0−1−0=0`.

Also `f₀>=15/16` and `f₂<=−15/16`, so the four graphs cannot have a common point. Their first and third strands are separated by at least `15/8` vertically. Lemma A therefore gives a positive size floor of at least `15/(8√2)` for every labeled periodic square.

At `x=7/8`, all four graphs have their unperturbed heights. The lifted ordered points

`z₀=(7/8,1),  z₁=(−1/8,0),  z₂=(7/8,−1),  z₃=(15/8,0)`

form an **exact square** of side `√2`. Their cylindrical x-coordinates all equal `7/8`, far outside the marked source windows.

Take the marked windows for the four strands to be, respectively,

`W₀=[0,1/16]`, `W₁=[1/4,5/16]`, `W₂=[1/2,9/16]`, `W₃=[3/4,13/16]`.

For any chosen planar lifts `x₀ in W₀ mod 1` and `x₂ in W₂ mod 1`, the absolute difference satisfies `|x₀−x₂|>=7/16`. However, for any labeled square we have the **exact identity**

`x₀−x₂ = y₁−y₃`.

The right-hand side has absolute value at most `2a=1/8`, since `h₁=h₃=0` and `|f₁|,|f₃|<=a`. As `1/8<7/16`, no square can have all four vertices in the marked windows. This is an inequality proof over continuous whole arcs, not a grid search.

The model remains valid for any rational amplitude `0<=a<7/32` with the graph separation maintained (`a<1`); the supplied checker uses `a=1/16`.

### Interpretation

This is a **localization obstruction**, distinct from v210's rotating-scaling obstruction and v211's scale collapse. Even when periodicization solves the latter, it can produce squares on purely artificial portions of the completions. One must prove *witness confinement* or an invariant that cancels connector witnesses before extracting a square on a single planar Jordan curve.

## 5. v212 theorem obligations / falsifiers

| Obligation | Current status | Falsifier / required evidence |
|---|---|---|
| Four periodic essential embeddings from local portions of arbitrary Jordan curve | **Not constructed** | Compact Jordan curve is not itself an essential periodic curve; added connectors may dominate. |
| Exact area balance in compatible orientation | **Constructible in artificial examples**, not from arbitrary curve | A balancing offset may move or create new connectors; proof must track it. |
| No-common-intersection gap uniform in approximation | **Lemma A proved conditionally** | For rough limits, `d_n` may tend to zero. |
| Guarantee a periodic witness uses actual Jordan pieces | **False without extra structure** by Proposition B | The exact four-graph model has no four-window square. |
| Continuous rough-curve extension | **Conditioned in [S] on zero planar measure** | A general Jordan curve can have positive planar measure; cannot assume otherwise. |
| Make v209 strict local crossings from periodic witnesses | **Not implied** | v210 orbit-packet obstruction still applies. |
| Relate full square-envelope winding to periodic Floer intersection | **Open** | Need a well-defined nontrivial global invariant insensitive to added connectors. |

## 6. Next research experiment

**Priority 1:** Design a *localization-preserving periodicization lemma*, not a generic periodicization. Specify original windows `W_i`, connector regions `C_i`, an explicit continuous deformation of four balanced strands, and a conserved/filtered algebraic count that cannot be discharged entirely by connector squares. The v212 graph model is the first mandatory failing test.

**Priority 2:** Study whether the Floer spectral invariant's physical-size information can be tied to this filtered count, while respecting the v211 moving-center obstruction.

**Priority 3:** Compare with the Nov. 17 updated Asano–Ike Liouville regularity criteria; keep a separate 'irregular' lane instead of implicitly assuming rough curves are zero-measure or action-regular.

There is not yet a demonstrated transfer from the Floer theorem to a single arbitrary Jordan curve.

## 7. Sources and provenance

- [V] Lorenzo Valerio, *A Proof of the Area Vanishing Square Peg Conjecture*, arXiv:2610.21076v2 (4 Nov 2026), https://arxiv.org/abs/2610.21076. Theorem 1.2, definition of joint periodic square, and the note that degenerate squares are permitted.
- [S] Ali Naseri Sadr, *The Four Periodic Inscription Problem*, arXiv:2610.24376v1 (17 Oct 2026), https://arxiv.org/abs/2610.24376. Theorems 1.2 and 1.4; common-intersection exclusion, area relation, and the zero-measure continuous case.
- Tomohiro Asano and Yuichi Ike, *Liouville regular Jordan curves*, arXiv:2609.11598v3 (17 Nov 2026), https://arxiv.org/abs/2609.11598. Assumption audit for rough curves.
- Internal v209, v210, and v211 primary research archives in `archive/`.

**Validation:** `python3 -m unittest discover -s tests -v` and `python3 code/periodic_bridge_v212.py`. All algebra in the counterexample uses `fractions.Fraction`; tests do not certify [V] or [S] themselves. For reproducibility, the archived v209–v211 checkers can be rerun independently.
