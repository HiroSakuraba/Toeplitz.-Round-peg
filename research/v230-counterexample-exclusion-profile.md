# v230 — Counterexample exclusion profile: what a genuinely new Toeplitz test must escape

**Status:** Source-checked classification of sufficient hypotheses, not a counterexample or proof. See [v230 double-spiral note](v230-irregular-frontier-double-spiral.md).

| Known theorem / criterion | Applies to | `D_p` for p<1 | `D_p` for 1≤p<2 | `D_p` for p≥2 | Why / scope limit |
|---|---|---|---|---|---|
| Classical polygon square-peg results (Pak) | Simple finite polygons | Not literally finite polygon | Not polygon | Not polygon | PL approximations satisfy theorem, but this does not force non-collapse in limit |
| Asano–Ike rectangle theorem v4 | All Liouville-regular Jordan curves; class includes rectifiable | YES | YES | NO | v3 Thm 6.3 sharp for compactified power spiral; failure means theorem not applicable, *not* absence of square |
| Li–Pan finite `q`-variation theorem, some q<2 | Jordan curves with finite q-variation | YES | YES | NO | Sharp threshold for this spiral q>p (q≥1) |
| Greene–Lobb two Lipschitz graphs | Curves that are union of two uniformly Lipschitz graphs with specified slope bound | Not established | NO (would imply finite length) | NO (would imply finite length) | Do not assume generic spiral has any two-graph representation |
| Asano–Ike Hölder-domain/quasicircle criterion | Certain domain classes implying Liouville regularity | Not independently checked | Not independently checked | NO | Contrapositive: D_p non-Liouville for p≥2 cannot be in these sufficient subclasses |
| Chambers near-C² theorem | Specific geometric `C²`-closeness hypotheses | Unchecked | Unchecked | Unchecked | Ordinary Hausdorff closeness to a circle is not the stated certificate; audit the theorem before using it |
| v100/v101 power-spiral computer-assisted square theorem | Exact `gamma_p` for all p>0 | CONDITIONAL YES | CONDITIONAL YES | CONDITIONAL YES with uniform >1/8 side | Separate direct existence construction, not an all-Jordan sufficiency criterion |

**Important logic:** 'NO' in the third column means **not covered by the sufficient hypothesis**, not square-free. `D_p` is expected to contain squares for all p>0 by the recovered spiral theorem. Thus it is a *positive control outside some theorem classes*, not a plausible counterexample to Toeplitz.

## Backlog triage

- `v220` fixed infinite-notch curve: **solved existence class** by rectifiability; retain as a collision/noncollapse falsifier, not frontier existence research.
- `v221–v229` rational polygon cobordism pipeline: **solved existence class**, useful for exact oracle, chart deduplication and approximation diagnostics. Freeze as supporting infrastructure until a specific irregular-curve question requires it.
- `v99–v102` power spirals: **priority frontier lane** for p>=2. Restore the corrected certificates first, then use as protected-square benchmark.
- `v206` Liouville/Abel audit: revisit for true nonregularity assumptions; do not presume all irregular candidates violate the published conditions.
- `v207–v208` finite coils, `v210` selected-arc counterexample, `v211` shape-sphere: cannot classify their **whole-curve** regularity from short README summaries. Inspect original archive proofs and specific equations before assigning frontier status.
- Asano–Ike v3 nonregular finite-coil example: **potential new comparator**, but its authors establish nonregularity, not automatically square existence for a closed derived curve.

## Recommended v231 experiment

Use a fixed p>=2 (start p=2, 5/2 and 3). Take finite spiral truncations at `epsilon<1/5`; by the conditional v230 theorem they retain at least one original spiral square of side `>1/8`. Build geometric pruning for four-edge candidates before calling v223's rational algebraic solver. Produce a **separate a priori certificate** that finds a positive-square without using the original v100 point parameters; otherwise the experiment merely re-encodes the known witness.

Check its assumptions against v220's known tiny-square sequence and v210's source-window counterexample. An oracle that merely finds a square in every finite polygon has made no progress on the unrestricted limit.
