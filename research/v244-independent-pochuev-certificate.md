# v244: Independently certify Grisha Pochuev's original ten-sided Big Peg counterexample

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact rational certificate replays the original ten-gon through our previously developed v243 two-parameter method, with an independent optional v223 geometric-orbit census. The original author's coordinates and proof sources are recovered and pinned.

## 1. Original primary source

- **Original author:** Grisha Pochuev
- **Original project:** https://github.com/Grisha-Pochuev/sharp-big-peg-counterexample
- **Pinned original commit:** b9e0c95cbab888405692c9efa7e592f2a0807ae4
- **Original rational data:** counterexample.json, Git blob 0f8b70e9aab72ca2c62e2e9713049e1558ba013f.
- **Published mathematical proof:** SOLUTION.md.
- **Three original exact verification programs:** sharp_big_peg_counterexample_verify.py, independent_verify.py, and independent_vertex_side_verify.py.
- **Prize organizer:** Daniel Bernier, https://dbernier.ca/2026/06/16/the-sharp-big-peg-conjecture/ (August 17, 2026 update).

The v243 note did not invent coordinates when these materials were unavailable. The exact input is now obtained directly from Pochuev's original public repository, not reconstructed from an image.

## 2. Ten exact rational polygon vertices

Each coordinate is an integer numerator divided by 1200, in cyclic order:

| i | X | Y |
|---:|---:|---:|
| 0 | 4114 | -558 |
| 1 | 2022 | 1568 |
| 2 | -1234 | 1006 |
| 3 | -868 | 2516 |
| 4 | -1089 | 2814 |
| 5 | -3769 | 1169 |
| 6 | -3459 | 81 |
| 7 | 141 | -1348 |
| 8 | 571 | -1449 |
| 9 | 1105 | -485 |

The original data are imported into [v244-pochuev-original-ten-gon.json](v244-pochuev-original-ten-gon.json) using exact rational-coordinate strings. The test suite compares all ten pairs to the published integral numerator list.

## 3. The sharp rational contradiction

The author's exact value for the smallest squared distance from the origin to any boundary segment is

$$
d_{\min}^2 = \frac{2403854890969}{2400326560000} > 1.
$$

The ten-gon is a simple polygon and the origin is strictly inside it. Therefore the closed unit disk lies inside its bounded region, giving a certified inradius lower bound of one.

The original certificate reports that the **largest** square whose four vertices lie on the polygonal boundary has squared side

$$
L_{\max}^2 =
\frac{56381275521625791352241234309}
     {28343589374868261752462760000}
\approx 1.98920732219.
$$

The exact strict gap to the proposed lower bound is

$$
2-L_{\max}^2=
\frac{305903228110732152684285691}
     {28343589374868261752462760000} > 0.
$$

Thus this polygon has a unit disk but **no boundary-inscribed square of side sqrt(2) or larger**. This refutes the sharp polygonal Big Peg conjecture, not the broader statement that some positive universal constant exists.

## 4. Independent local verification, using the existing v243 engine

The author's three original programs exhaust 10^4 four-corner edge assignments through exact four-dimensional linear systems. They use two square parametrizations and independent solvers.

The new local certificate deliberately relies on a different existing method, our v243 two-parameter convex-polytope verifier:

1. Choose two adjacent square corners P and Q on any pair of polygon edges, with rational parameters (t,u) in the unit square.
2. Construct the remaining corners by a +90 or -90 degree rotation of Q-P, and assign them to two other polygon edges.
3. Enforce both remaining edge line-membership equalities and their segment bounds as exact rational linear constraints in (t,u).
4. The feasible parameter region is a compact rational polytope of dimension zero, one or two. A convex quadratic function, the squared side length |Q-P|^2, maximizes at an extreme point, found exactly in all ranks.
5. Exhaust all 2*10^4 ordered chart/orientation combinations and retain the largest squared side with an attained witness.

The same v243 library checks polygon simplicity, center containment and the exact minimum squared distance to each polygon segment. Our runner asserts equality with all three original published rational invariants: boundary clearance, largest squared side and positive gap.

An optional second, pre-existing local engine, the v223 exact four-parameter geometric square-orbit census, checks that the polygon has a single positive square orbit and the same squared side. This adds a crosscheck independent both of the author's original implementations and of the v243 maximum-only algorithm.

This is a mathematical upper bound on **every** inscribed square, not a heuristic maximum obtained from sampled candidates. There are no floating point proof decisions.

## 5. Reproduce

From the root of this repository:

    python code/pochuev_big_peg_v244.py
    python code/pochuev_big_peg_v244.py --crosscheck-v223
    python -m unittest -v tests.test_pochuev_big_peg_v244

The GitHub Actions workflow v244-pochuev-exact.yml repeats all three commands in a fresh standard-library Python environment. The JSON output gives exact rational values, the maximizing square and complete audited-chart counts.

**Scope:** The sharp sqrt(2) polygonal Big Peg bound is false. The existence of a **smaller universal c>0** and the unrestricted Toeplitz square-peg problem remain open.

**Provenance and licensing:** The coordinate data and original numerical invariants are attributed to Grisha Pochuev and pinned by their Git commit/blob identifiers. The author's verifier code is not copied. His source code is released under the MIT License.
