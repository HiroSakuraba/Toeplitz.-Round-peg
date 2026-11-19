# v245: a new exactly certified rational-grid Big Peg decagon

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Scope:** An independently modified ten-vertex rational Jordan polygon with a *strictly smaller certified upper bound* on maximal-square-side / inradius than the original Pochuev origin-disk certificate. This is a concrete numerical-to-exact optimization result, **not** a proof that the true optimal ratio of this new polygon is lower than the true optimal ratio of Pochuev's polygon, and **not** a universal Big Peg theorem.

## What changed

v244 imported and independently confirmed Grisha Pochuev's exact rational ten-gon. The original has a contained disk centered at the origin with squared radius

\[
r_{0}^{2}=\frac{2403854890969}{2400326560000}.
\]

Its largest boundary-inscribed square has exact squared side

\[
s_{0}^{2}=\frac{56381275521625791352241234309}{28343589374868261752462760000}.
\]

Thus the explicit disk already gives the **rigorous upper bound** \(s_{\max}/r_{\mathrm{in}}\le s_0/\sqrt{r_0^2}=1.4093571594\ldots\). No claim is made that the origin-centered disk is maximal.

v245 began with numerical experiments on Pochuev's exact ten-gon. Unlike v242's Fourier-only search, the new candidate generator perturbs the actual counterexample on an integer coordinate grid and explicitly distinguishes numerical candidate **ranking** from exact acceptance. The search encountered narrow changes in square incidence, and several nearby integer-grid moves generated additional large squares, which is precisely why exact exhaustive checking is mandatory.

One improved candidate was rationally reconstructed at denominator \(12000\). The new vertices, in cyclic order, are

| i | 12000 x | 12000 y |
|---:|---:|---:|
| 0 | 41140 | -5580 |
| 1 | 20219 | 15680 |
| 2 | -12341 | 10060 |
| 3 | -8680 | 25160 |
| 4 | -10890 | 28140 |
| 5 | -37680 | 11679 |
| 6 | -34590 | 811 |
| 7 | 1400 | -13480 |
| 8 | 5711 | -14490 |
| 9 | 11050 | -4850 |

The complete exact [candidate JSON](v245-improved-ten-gon.json) records source attribution to Pochuev, his pinned original commit, and the exact integer-coordinate differences from his original vertices. These vertices are **new experimental modifications**, not coordinates attributed to him.

## The new exact certificate

Set the rational disk centre \(C=(3/10000,0)\). Exact segment-distance evaluation gives

\[
r_{\mathrm{cert}}^{2}
=\min_{e}\operatorname{dist}(C,e)^{2}
=\frac{1352506475516588161}{1349561502900000000}
=1.00218216999\ldots>1.
\]

The polygon is exactly checked as Jordan and the centre exactly checked as interior. Thus this entire closed disk is contained in its bounded component.

Our existing v243 exact *two-parameter convex feasible-chart* census checks all \(2\times10^4\) ordered edge and orientation cases, including singular nonisolated families, and returns an **attained exact maximum square** with

\[
s_{\max}^{2}
=\frac{2263154862866037715637393248853633089}
{1137974336077532830058267809922880000}
=1.98875738328\ldots
\]

The strict new scale-independent upper bound is

\[
J_{\mathrm{new}}^2=\frac{s_{\max}^{2}}{r_{\mathrm{cert}}^{2}}
=\frac{43286092375636809252886667578815515500510338125}
{21812891985778816953208198388629401097244304578},
\]

so

\[
\boxed{J_{\mathrm{new}}\le1.408696924883410}
\]

for the true maximum-side-to-*inradius* ratio (the true inradius can only be **larger** than the certified disk radius). This improves on the original published origin-disk certificate upper bound **1.409357159407647** by approximately 0.000660235 in certified upper-bound value.

Every comparison is made using exact Fraction arithmetic. The strict inequality between the two *squared rational upper bounds* is asserted directly; the decimal approximations are for display only.

### Important comparison caveat

We have **not computed the actual global inradius of Pochuev's original decagon**. Therefore we have established an improved *certified bound* over the published origin-disk baseline, **not** yet a theorem that this new polygon has a strictly smaller *true* \(J\) than the original. To establish that stronger statement, one needs an upper bound on the original decagon's true inradius (for example a complete exact medial-axis/disk-centre optimization).

This refinement is nonetheless a concrete new candidate with an exact maximum-square calculation and a strictly improved independently reproducible upper bound. It also disproves the sharp \(c=\sqrt2\) polygonal Big Peg inequality, as did the original.

## Numerical search versus rigorous proof

- [big_peg_search_v245.py](../code/big_peg_search_v245.py) provides numerical complete-*assignment* sampling of all nonsingular 4-by-4 edge systems, then ranks rational-grid neighbors. It is **not a proof**: floating-point solves near degeneracy, singular edge sets, and polygon validity can all invalidate a numerical candidate.
- The explicit candidate can be reconstructed from the original grid at denominator 12000 with the displacement list stored in the JSON. The numerical search was used only to *propose* it.
- [big_peg_improvement_v245.py](../code/big_peg_improvement_v245.py) uses the *previously merged exact v243 library*, not the numerical finder, to certify the maximum, disk, provenance, exact fractions, and strict comparison. The optional v223 crosscheck separately enumerates exact geometric-square orbits and compares its largest square.
- [tests/test_big_peg_improvement_v245.py](../tests/test_big_peg_improvement_v245.py) and dedicated CI fail if the candidate geometry, certificate values or comparisons change.

## Reproduction

Run from repository root:

    python code/big_peg_improvement_v245.py
    python code/big_peg_improvement_v245.py --crosscheck-v223
    python -m unittest -v tests.test_big_peg_improvement_v245

Numerical-only exploratory commands (requires numpy):

    python code/big_peg_search_v245.py
    python code/big_peg_search_v245.py --steps 40

**Limits:** This is an example of a *better certified upper bound*, not a determination of the smallest possible \(J\), nor an assertion that the general c>0 Big Peg conjecture fails. The unrestricted Toeplitz square-peg conjecture remains open. Future work should certify the actual inradius of both candidates and optimize polygon vertices while maintaining exhaustive square census.
