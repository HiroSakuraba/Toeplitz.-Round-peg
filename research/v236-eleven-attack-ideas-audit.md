# v236 — Audit of *Eleven ways to attack the square peg*

**Input:** User-provided 19 November 2026 HTML artifact *Toeplitz attack ideas.html*. This is an independent brainstorm, not a peer-reviewed or certified theorem. The notes below respond to its exact six near-term and five moonshot proposals. **Do not silently treat the article's experimental figures or conjectures as theorems.** The user's v234 and v235 GitHub PRs were inspected on 19 November: **#24 and #25 were OPEN at review time**, with #25 targeting #24's branch. Their merged status must be checked again before reusing their calculations.

## Near-term proposals

| Source idea | Technical finding | Action / falsifier |
|---|---|---|
| **Put the roughness where the square is** | The right stress test: perturbation on the four actual certified source intervals, with amplitude comparable to the square side. PR #25's reported base-32 Takagi result is still **under review**, not a merged result. | Continue amplitude ramp *and recertify*; mere failure of an existing source-box certificate is **not absence of squares**. Try a genuinely rough four-arc instance not globally close to the v234 template. |
| **Template library with a wide C0 margin** | Portable v231/v232 Brouwer degree gives a quantitative margin, but a library of thin local neighborhoods does not automatically cover arbitrary Jordan curves. Conditioning and chord-error budgets compete. | Store the *maximum rigorously certified perturbation margin* in dimensionless units `epsilon/side`, not heuristic distance-to-template; distinguish nominal vs worst-case condition numbers. |
| **Open log spirals look square-free** | Corrected and sharpened in [theorem](v236-open-logarithmic-spiral-square-exclusion.md): `S_a` is **provably square-free for every a>=4**. However v99 already certifies a **square-bearing special pitch** `a≈0.1002747`. Self-similarity makes each square part of a family on a fixed pitch and defeats naive isolated-root IFT arguments. | Treat large-a arc as rigorous **negative** and v99 pitch as rigorous **positive** controls. Classify 0<a<4 without inferring a global no-square result for a Jordan closure. |
| **Follow the odd square down a canonical flow** | Analytic level curves of a Riemann map are useful, but a generic *odd number* of squares does not pick one canonical branch; births, tangencies, and reconnections must be addressed. | Monitor all branches / robust degree classes; do not designate an 'essential branch' without a proved selection rule. |
| **An action-pinning lab** | Potentially high-value, but the action used in rectifiable Floer arguments may **not be defined canonically** on a nonrectifiable curve. | Compute a precisely defined polygonal action proxy and its refinement dependence; do not call that a limit invariant until proved. |
| **Hunt for curves with only small squares** | Search errors are asymmetric: finding a large square gives a lower bound on the **maximum**, but finding no large square numerically does **not** give an upper bound. This is crucial for testing a universal inradius floor. | Use complete v223 census for manageable rational polygons or *certified exclusion* of all squares above a cutoff. Keep inradius disk itself certified. |

## Moonshots

| Source idea | Current audit | Primary obstacle |
|---|---|---|
| **Quantitative Toeplitz: square >= c × inradius** | This is the historically recognized **Big Peg conjecture**, discussed by Matschke (2014), and reported as Schwartz's Conjecture 3.5 (2022). It is not a new proposal. A 2026 blog reports a 10-vertex counterexample to the **sharp choice `c=sqrt(2)`**, but this does not disprove existence of some smaller universal `c>0`. | Certified upper bounds on the largest possible square and exact interior-disk containment. |
| **Few bad points theorem** | A valuable path if action/degree can be localized away from singularities, but 'finite bad points' does not automatically imply the extension to countably many of arbitrary Cantor–Bendixson rank. | Construct a global invariant with controlled singular supports and prove exhaustion, not merely pointwise regularity away from a finite set. |
| **Blow-up limits and collapse cost** | Useful as a negative/compactness laboratory. General continuous curve blow-ups need not be precompact as curves under scale normalization. | Specify topology and an energy compactness theorem before making claims about limit objects carrying square witnesses. |
| **C0 symplectic topology for gamma×gamma** | Aspirational: for an arbitrary continuous curve, calling the product a 'C0 Lagrangian' already involves nontrivial choice and rigidity hypotheses. | Verify precisely which Lagrangian spectral invariants remain defined and continuous under the proposed approximations. |
| **Dimension 2 benchmark** | Hausdorff dimension 2 excludes finite `p<2` variation, but **does not imply outside all square-existence criteria**: C0-near-circle / thin-annulus theorems can include extremely irregular curves. Some positive-area Jordan examples may be covered by such open criteria. | Build examples proven to fail **each actual theorem's hypotheses**, rather than using dimension alone as a proxy. |

## Revised sequence

1. **Now (v236):** establish large-openness square-free logarithmic-spiral theorem; use v99 special pitch as positive control and reject numerical 'no exact square' as proof. The local exclusion proof is new to this research program and may overlap prior literature; originality has not been established.
2. **Next:** test the claim that a globally continuous Jordan curve cannot hide behind a square-free arc: attach several different closures to the same large-`a` spiral and certify where the full-curve squares actually occur. A negative local result has no effect on global Toeplitz existence by itself.
3. **In parallel:** review/merge v234 and v235 only after their own checks pass, then quantify the margin-vs-amplitude/side ratio and test the v220 notch / v210 lift controls. Avoid treating unmerged PR results as mainline evidence.
4. **Only with certified maximum-square upper bounds:** run the inradius small-squares hunt, compare with the reported sharp-Big-Peg counterexample, and avoid optimizing based on an incomplete heuristic maximum.
5. **Longer-term:** action/conformal-flow/blow-up experiments should state their invariant and its compactness requirements before broad claims.

## Primary and supporting literature

- Benjamin Matschke, *A survey on the square peg problem*, Notices AMS (2014), https://www.ams.org/notices/201404/rnoti-p346.pdf — explains why a positive inradius-size bound would imply Toeplitz.
- Matschke, *On the square peg problem and its relatives*, Trans. AMS 375 (2022), https://doi.org/10.1090/tran/8649 — open sets of curves already known to contain squares.
- Asano–Ike, *Liouville regular Jordan curves*, v3, https://arxiv.org/abs/2609.11598 — precise regularity restrictions.
- Greene–Lobb, *Square pegs between two graphs*, Comment. Math. Helv. (2026), https://doi.org/10.4171/CMH/619.
- External 2026 sharp Big Peg counterexample report **(secondary; original coordinates not independently audited here)**: https://dbernier.ca/2026/06/16/the-sharp-big-peg-conjecture/
- The user's restored exact v99 proof, `related/power-spiral-v99-v102/handoff_v99_small_exponent_log_spiral.md` and `verify_v99_small_exponent_log_spiral.py` in this repository.

**No unrestricted Toeplitz proof is claimed.**
