# v214: Source- and lift-filtered Brouwer degree for periodic square configurations

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** New finite-dimensional, conditional topology result and exact-rational implementation. **Not** a relative Floer homology theory. **Not** a proof of the Toeplitz square-peg conjecture.

## Executive result

The v212–v213 counterexamples require controlling **both** which parts of a periodic extension a square touches and which *integer planar lifts* it uses. The following construction does that at the finite-dimensional topological-degree level. For a fixed choice of four real source arcs, fixed relative deck shifts, and one square orientation, a signed count of transverse squares is constant under any homotopy that introduces **no zero on the original arc-parameter boundary**. It is possible to certify that condition with exact rational inequalities and thereby carry a square onto genuinely specified original arcs under a controlled deformation.

This is a valid conditional local theorem. It does **not** show that an arbitrary Jordan curve supplies a domain with nonzero degree, that periodic area balance implies nonzero *filtered* degree, or that this degree is compatible with a Floer chain differential. Those are open.

## 1. Four-arc square map and domain

Let each arc be a graph chart `p_i:[0,1] → R²`, oriented by increasing horizontal coordinate, with rational PL coordinates. The four lifted physical arcs are `p_i(u_i)+(k_i,0)`, with **fixed** relative integers `k=(0,k_1,k_2,k_3)`. Normalizing `k_0=0` removes a simultaneous horizontal deck translation. This data is fixed *before* counting squares; it cannot be changed merely because another periodic lift contains a square.

Let `R(a,b)=(-b,a)` denote 90° rotation. For orientation `σ∈{+1,-1}` and points `q_i=p_i(u_i)+(k_i,0)`, define

```
e=q1−q0,
F_σ(u)=(q2−q1−σR(e), q3−q0−σR(e)) ∈ R⁴.
```

Then `F_σ(u)=0` holds exactly when the four labeled points, in order, form a square of orientation `σ`, **possibly degenerate**. The algorithm additionally assumes the vertical ranges of source arcs 0 and 2 have strict positive separation `g>0`. Every zero is therefore nondegenerate, and this nondegeneracy survives vertical perturbations of size `<g/2`.

The parameter domain is the product `Q=[0,1]^4` with the endpoints representing the four **specified source-window endpoints**. The interior refers to the original arc windows, not the periodic completion, and all lift tags are held fixed.

## 2. The conditional relative-degree theorem

**Theorem.** Suppose `F_σ(u)≠0` for every `u∈∂Q` and `g>0`. Define `D_σ=deg(F_σ,Q,0)` by Brouwer degree. Then:

1. If the roots are isolated and transverse, `D_σ` is the integer sum of their local signs `sign det DF_σ(u)`.
2. If `D_σ≠0`, the four **physical** source arcs contain a nondegenerate inscribed labeled square of orientation `σ`.
3. For any continuous deformation of these four source arcs **holding the chosen windows and deck shifts fixed**, if the square map stays nonzero on the full boundary `∂Q` and the opposite-arc separation stays positive, the integer `D_σ` is unchanged, even if intermediate intersections become non-transverse.

**Proof.** The definition and homotopy invariance are the standard existence and homotopy properties of Brouwer degree for maps on a bounded open domain with the target point absent from its boundary. Every interior zero is a square by the exact residual identities. The fixed source charts/lifts ensure the square lies on the designated physical pieces rather than on connectors or other integer translates. The positive opposite-arc vertical gap rules out `q0=q2`, which a zero-size square would require. The transverse signed formula follows from the local orientation degree of an invertible affine map and additivity over isolated roots.

**Caveat:** This is a topological degree of a map between *finite-dimensional real parameter spaces*. It is not a claim about filtered Floer homology, Lagrangian non-displaceability, or a relative invariant in Valerio's or Naseri Sadr's symplectic setup.

## 3. Exact PL evaluator

For each combination of straight segments from the selected windows, use oriented local parameters `u_i∈[0,1]`. Since these parameters increase with the original arcs, chart transitions preserve orientation. The residual is affine `F_σ(u)=Au−b` with rational `4×4` matrix `A` and rational `b`.

- If `det A≠0`, solve `Au=b` exactly. A root strictly within `(0,1)^4` contributes `sign(det A)`. A root on **any** cell face, including an internal segment knot, is rejected rather than silently counted twice. The latter restriction is stronger than the theorem's global outer-boundary requirement; it keeps the evaluator straightforward and auditable.
- If `det A=0`, enumerate vertices of `{Au=b, 0≤u≤1}` exactly using the existing v213 rank-deficient feasibility routine. A feasible cell is rejected as non-isolated. An infeasible singular cell contributes zero but cannot be assigned the inverse-matrix quantitative robustness margin described next.
- If all cells pass, the total signed sum is the Brouwer degree of the chosen four-arc domain. It may be zero even if squares exist, since opposite signs can cancel.

The evaluator returns a positive *quantitative boundary margin* only when all cells are invertible. It never invents a margin for a singular cell. Invalid or non-transverse cells are explicit errors rather than silently successful certificates.

## 4. A purely rational sufficient boundary-avoidance bound

For an invertible PL cell `F(u)=A(u−u*)`, let `r_j` be the sum of the absolute values in row `j` of `A⁻¹` (so `r_j>0`). Then

`|u_j−u*_j|≤r_j‖F(u)‖∞`.

For a box face `u_k=c∈{0,1}`, we obtain

```
‖F(u)‖∞ ≥ max( |c−u*_k|/r_k,
               max_{j≠k} dist(u*_j,[0,1])/r_j ).
```

Take the minimum of this rational lower bound over **all eight faces of every cell**. Call the result `m>0`. If this succeeds, `‖F(u)‖∞≥m` on the full four-arc domain boundary, indeed on the union of every PL cell face.

Now perturb each selected arc vertically by a continuous function `h_i` with `‖h_i‖∞≤ε`. Each component of the square-residual perturbation has absolute value at most `2ε`: the perturbations enter as differences of *two* heights or as the 90-degree rotation of such a difference. Therefore for every time along the linear homotopy,

`‖F_t(u)‖∞ ≥ m−2ε` on the boundary.

The opposite-height gap decreases by at most `2ε`. Consequently a **single exact, sufficient robustness threshold** is

```
ε < min(m/2, g/2).
```

Under this inequality the full homotopy is boundary-free, degree is preserved, and any nonzero degree forces a nondegenerate square on the same four labeled physical lifts. This includes targets with arbitrary continuous vertical perturbations; they do not need to remain PL, differentiable, or transverse. The evaluator computes a conservative threshold, not a sharp threshold.

**Example.** The planted four-arc square has `D_+=−1`, one isolated square of side-squared `1/32`, opposite height gap `5/64`, cellwise rational boundary bound `m=55/4368`, and strict vertical perturbation threshold `ε<55/8736`. These values are reproducible by `python code/relative_degree_v214.py` and the exact tests.

## 5. Three adversarial controls

### A. Correct source tags with the wrong deck shifts

The v213 deck-lift configuration can have nonzero degree for the **correct** tuple `(0,−1,0,+1)` while having degree zero for the source windows interpreted in the **wrong** tuple `(0,0,0,0)`. The integer shift data are part of the domain, not an incidental detail recoverable from the projected cylinder.

### B. Zero degree with two actual squares

Hold arcs 0, 1, and 3 fixed in the planted positive example. Replace arc 2 by a two-segment tent whose endpoint heights sit `1/100` **below** the affine square-compatible line `y=1/4−2(x−1/2)` and whose central knot sits `1/100` **above** it. The exact solver finds **two** interior transverse squares, with determinant signs `−` and `+`; hence `D_+=0`. Lowering the tent by `1/50` removes both squares without changing the signed degree. Thus **degree zero does not prove absence**; a nonzero degree is sufficient, not necessary.

### C. A zero escaping through the source boundary

Translating a source arc vertically can carry its square parameter to the edge of a marked window. At the precise crossing the evaluator refuses to assign a transverse relative-degree certificate. Just beyond the crossing the degree can change from `±1` to `0`. This is why the fixed original-arc parameter boundary, not merely the avoidance of degenerate square sizes, is a genuine required hypothesis.

## 6. Reproduction and next step

```
python -m unittest discover -s tests -v
python code/relative_degree_v214.py
```

The complete suite covers v212, v213, and v214, with exact `fractions.Fraction` arithmetic and a seeded cross-check against v213's independent exhaustive PL square solver. The Python code is not proof-assistant verified, but the continuous-theorem statements and rational inequalities are set out explicitly above.

**Next research obligation:** We need to produce, from **one arbitrary Jordan curve**, four source windows and *their physical lift tuple* with nonzero relative degree, or find an alternate invariant that cannot cancel. The periodic Floer existence theorem supplies no such localization-filtered degree. To connect with symplectic methods, one would have to define a relative Floer chain complex or continuation map respecting source and lift tags and prove compactness, boundary avoidance, action control, and nonvanishing. None has been established.

## Sources

- Lorenzo Valerio, *A Proof of the Area Vanishing Square Peg Conjecture*, arXiv:2610.21076 (2026). https://arxiv.org/abs/2610.21076
- Ali Naseri Sadr, *The Four Periodic Inscription Problem*, arXiv:2610.24376 (2026). https://arxiv.org/abs/2610.24376
- Earlier v212–v213 notes and archived v209–v211 research bundles in this repository.
