# Recovered power-spiral proof files (v99–v102)

These files are **exact original copies** recovered from the author's saved ChatGPT file library, not newly reconstructed or edited proof sources. Their SHA-256 checksums are recorded in `SHA256SUMS.txt`. The wrapper `verify_v230_arb_strict.py` and this README are NEW v230 additions.

## Proof chain

1. `handoff_v99_small_exponent_log_spiral.md` and `verify_v99_small_exponent_log_spiral.py`: small `p` interval `0<p<=1/200`, with rational analytic inequalities and a **non-proof numerical locator**. Requires `mpmath` for the locator.
2. `handoff_v100_power_spiral_all_exponents.md`: analytic reduction, Poincaré–Miranda degree, compactified boxes for `1/200<=p<=2` and `p>=2`, and rigorous endpoint-containment claim. **The historical `verify_v100_power_spiral_all_exponents.py` must NOT be accepted as the final interval proof!**
3. `handoff_v101_interval_layer_audit.md`: identifies **an unsound interval enclosure** in historical v100 `Efun_scalar` on negative arguments and two other numerical issues, repairs the logic, and reports an independent second certified layer.
4. `iv101.py`: recovered corrected outward-rounded Decimal transcendental layer. **The saved handoff mentions `verify_v101_faces.py`, `test_iv101.py`, and `probe_bugs.py`, which were NOT located in the available library.** Without its missing face driver, `iv101.py` alone is not a rerunnable full Decimal proof.
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
- **This working container could not install `python-flint`**, so v230's local tests certify the rational deductions and regularity classifications, **not the original eight transcendental face signs**. The dedicated GitHub workflow attempts to reproduce the eight-face proof independently and must finish successfully before saying the Arb path was rerun.
- The missing corrected Decimal face-driver and independent end-to-end analytical proof audit remain outstanding. A passing unit-test suite alone does not prove the geometric theorem.
- The v230 new theorem `square side >1/8` for `p>=2` is conditional on the old **Box A** zero, regardless of the software test outcome.
