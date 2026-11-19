# Recovered power-spiral proof files (v99–v102)

These files are **exact original copies** recovered from the author's saved ChatGPT file library, not newly reconstructed or edited proof sources. Their SHA-256 checksums are recorded in `SHA256SUMS.txt`. The wrapper `verify_v230_arb_strict.py`, the corrected `iv230.py` and this README are NEW additions.

## Proof chain

1. `handoff_v99_small_exponent_log_spiral.md` and `verify_v99_small_exponent_log_spiral.py`: small `p` interval `0<p<=1/200`, with rational analytic inequalities and a **non-proof numerical locator**. Requires `mpmath` for the locator.
2. `handoff_v100_power_spiral_all_exponents.md`: analytic reduction, Poincaré–Miranda degree, compactified boxes for `1/200<=p<=2` and `p>=2`, and rigorous endpoint-containment claim. **The historical `verify_v100_power_spiral_all_exponents.py` must NOT be accepted as the final interval proof!**
3. `handoff_v101_interval_layer_audit.md`: identifies **an unsound interval enclosure** in historical v100 `Efun_scalar` on negative arguments and two other numerical issues, repairs the logic, and reports an independent second certified layer.
4. `iv101.py`: recovered Decimal transcendental layer, intended to be outward-rounded. **Audit: it is not fully outward-rounded** — `Decimal.sqrt()` ignores `ROUND_FLOOR`/`ROUND_CEILING`, and the atanh tail bound rounds `term*t*t` toward `+∞` with `t<0`. Kept byte-identical for provenance; use `iv230.py`. **The saved handoff mentions `verify_v101_faces.py`, `test_iv101.py`, and `probe_bugs.py`, which were NOT located in the available library.** Without its missing face driver, `iv101.py` alone is not a rerunnable full Decimal proof.
4a. `iv230.py` (**NEW**, audit): `iv101.py` with exactly those two fixes — square roots stepped outward and verified by exact squaring, and the tail bound computed with `|t|`. Tested in `tests/test_iv230_interval_layer.py`. Still no face driver, so it is a corrected library, not a proof by itself.
5. `verify_v101_arb.py`: recovered independent Arb ball arithmetic driver for **all eight face signs** using exact rational boxes. Requires `python-flint`.
6. `handoff_v102_radial_phase_robustness.md`: conditional robustness extension for a normalized radial phase-profile neighborhood; beyond the scope of v230's elementary inequalities.

## Run

```sh
python -m pip install mpmath python-flint
python verify_v99_small_exponent_log_spiral.py
python verify_v230_arb_strict.py
```

**Important:** `verify_v230_arb_strict.py` is the CI entry point. The recovered `verify_v101_arb.py` prints `FAILED` on unresolved cells **without setting a failing exit code**, so the new wrapper explicitly asserts every face has zero unresolved cells and positive margin. Neither driver should be confused with an independent inspection of the analytic derivation.

For a small independent standard-library regression with **no Arb dependency**, run `python code/irregular_spiral_frontier_v230.py` from the repository root and `python -m unittest discover -s tests -p 'test_irregular_spiral_frontier_v230.py'`.

## Provenance and verification boundary

- The original v101 audit handoff reported eight face signs independently certified in Arb ball arithmetic (160-bit precision), with zero unresolved cells.
- **Audit (19 November 2026):** the strict Arb wrapper was rerun with python-flint 0.9.0: all eight faces certified at 160 bits, zero unresolved cells, exit 0.
- **The original v230 working container could not install `python-flint`**, so v230's local tests certify the rational deductions and regularity classifications, **not the original eight transcendental face signs**. The dedicated GitHub workflow attempts to reproduce the eight-face proof independently and must finish successfully before saying the Arb path was rerun.
- The missing corrected Decimal face-driver and independent end-to-end analytical proof audit remain outstanding. A passing unit-test suite alone does not prove the geometric theorem.
- The v230 theorem `square side >1/4` (audit; originally stated as `>1/8`) for `p>=2` is conditional on the old **Box A** zero, regardless of the software test outcome.
