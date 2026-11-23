# v338. Robust terrain squares: the interior case is proved, and global parity is not enough

Benjamin John Schulz · November 2026

[v336](v336-robust-terrain-squares.md) conjectured that every grid terrain carries a robust inscribed square of side at least i/√2, where i is the side of its largest interior axis-aligned square. This note starts the proof for single terrains. It complements [v337](v337-terrain-component-degree.md), which settles two fixtures exactly: the robust size of (2, 7) and of the v334 glued witness A is exactly i/√2.

**Findings.**
1. **Theorem 338.1 (proved).** Suppose some largest interior square has a column of height less than i on each side. Then every registered perturbation of size ε < δ ≤ 1/8 carries a square of side at least i − 4δ. So the robust size is at least i, which beats 1/√2. The proof applies v316's bent-base theorem to the perturbed curve. In this case the v336 corner squares play no role.
2. **Cutoff parity.** N_σ, the number of squares of side greater than σ, has a parity that does not change under small perturbations, as long as σ is not the side of a square of the exact terrain. When that parity is odd, a robust square of side greater than σ exists. Across 282 terrains, 6 perturbations each, the parity agreed at every admissible σ, and the total count was always odd.
3. **Parity alone cannot prove the conjecture.** It certifies i/√2 on 279 of the 282 terrains. On the other three, every admissible cutoff at or above i/√2 has even parity, although their robust ratios are 1.12, 1.33 and 1.70. A proof therefore needs marked degrees on chosen arcs, as in Theorem 338.1, rather than the total square count.
4. **What remains is the end-wall case.** Here the excursion starts or ends on an end wall. It covers 252 of the 282 sampled terrains, including the sharp example (2, 7), which v337 settles exactly.

## Theorem 338.1: interior excursions

Fix a terrain with integer heights H(0), …, H(W−1), a largest interior square on the columns [a, a+i), and a height h just below i. v316's excursion leaves the roof at the wall x = l and returns at x = r.
- l − 1 is the last column left of the square with H < i.
- r is the first column right of the square with H < i.

If no such column exists on a side, the excursion runs to the end wall there: l = 0 or r = W. The theorem covers the case 0 < l and r < W.

**Statement.** In that case, let 0 < ε < δ ≤ 1/8. Every simple polygon obtained by moving each vertex of the terrain by at most ε in max-norm, keeping edges straight, has an inscribed square of side at least i − 4δ.

**Proof.** Put h = i − 3δ and apply the bent-base theorem of [v316](v316-roof-shear.md) to the perturbed curve. Take L to be the perturbed base edge and A the following arc.
- **The arc A.** The exact wall at x = l runs from height H(l−1) ≤ i − 1 up to H(l) ≥ i. Its perturbed copy is a segment whose endpoints are within ε of those, so it crosses height h exactly once, and the same holds at x = r. Let A run between these two crossings. Every perturbed roof vertex strictly between them comes from a vertex of height at least i, so it has height at least i − ε > h + 2δ. The edges joining such vertices therefore stay above h + 2δ too.
- **Hypothesis 1.** L lies in |y| ≤ ε < δ. It starts within ε of x = 0 and ends within ε of x = W. A lies in [l − ε, r + ε] with l ≥ 1 and r ≤ W − 1, so L extends beyond A by more than 2δ at both ends.
- **Hypothesis 2.** The part of A in the band h ≤ y ≤ h + 2δ lies on the two perturbed wall segments. Each projects to an interval of width at most 2ε < h − δ. Their separation is at least r − l − 2ε ≥ i − 2ε > h + δ, and A starts in one interval and ends in the other.

The bent-base theorem then gives a square with two vertices on A and two on L, of side at least h − δ = i − 4δ. ∎

Letting δ decrease to ε shows that perturbations of size ε keep a square of side at least i − 4ε. In the sample, the 30 terrains in this case have smallest robust ratio 1.0001. The replay also checks the bent-base hypotheses exactly on perturbed copies of three interior terrains, for example heights (1, 3, 3, 1, 4, 2).

The theorem does not reach the end-wall case. There one end of A lies on an end wall of the base, so the base cannot extend 2δ beyond A. This is the same corner geometry that makes the v316 square fragile.

## Cutoff parity

**Lemma 338.2 (sketch).** Let γ be a Jordan polygon and σ > 0 a number that is not the side of any inscribed square of γ. Then all generic polygons sufficiently close to γ have the same parity of N_σ.

*Sketch.* Squares of side at least σ/2 on curves converging to γ converge to squares of γ. So close enough to γ no curve has a square of side exactly σ. Join two nearby generic polygons by a generic path inside that neighbourhood. Along such a path squares are born and die in pairs, and a square vertex passing a polygon vertex does not change the count's parity. This is the mechanism behind the odd square count of generic curves, and no square crosses side σ. The transversality details for polygons are not written out here, so the lemma is used as a working tool.

**Measurement.** For each terrain the replay computes every exact square side, as isolated values and ranges along families (along a family the side vector is affine, so the range is exact). The gaps of width at least 0.02 between these ranges give the admissible cutoffs. The parity at each gap's midpoint is read from 6 perturbations of size 10⁻³. The certified robust side is the top of the highest gap with odd parity.

| Case | Terrains | Lowest robust M/i | Lowest parity-certified ratio | Below 1/√2 |
| --- | ---: | ---: | ---: | ---: |
| Interior (Theorem 338.1) | 30 | 1.0001 | 0.7071 | 0 |
| End wall | 252 | 0.70705 | 0.4714 | 3 |

On all 282 terrains the parity at every gap was the same for all 6 perturbations, and every perturbation had an odd total number of squares.

**Where parity fails.**
- **(8, 4, 3, 1, 1, 6), i = 3.** Every gap above side 2 has even parity, yet each perturbation keeps a square of side about 5.1.
- **(5, 1, 7, 4, 6, 4), i = 4.** The only gap above i/√2 = 2.83 has even parity. The robust ratio is 1.12.
- **(4, 3, 1, 6, 3, 8), i = 3.** Exact families cover every side from √2 to 4, so no admissible cutoff exists there, and the gaps above 4 have even parity. The robust ratio is 1.33.

In the first example the robust squares come in pairs that a larger, but still continuous, deformation could cancel. A small perturbation does not cancel them, but parity cannot see them. Even on interior terrains parity certifies less than Theorem 338.1 proves: for (1, 3, 3, 1, 4, 2) it gives only i/√2.

## What is left

The conjecture for single terrains is now reduced to the end-wall case: every largest interior square has, on at least one side, no lower column before the end wall.

[v337](v337-terrain-component-degree.md) supplies the method for this case on two fixtures. On (2, 7), the exact squares with one vertex on the base, one on each wall and one on the roof form a single family. It runs from the diamond, of side √2, to the corner square, of side 2. A perturbation can keep either end, and v337 shows that the whole family has odd degree in an isolating neighbourhood. v337 also gives a four-contact lemma: base, two walls and a roof step, all met in the interiors of their segments. That lemma certifies a robust square of side at least w/√2 directly.

Two steps remain for a general end-wall terrain.
- Select a large-square component whose contacts stay off arc endpoints, as v337 does by hand, or a four-contact configuration.
- Show that its degree is odd and that its squares have side at least i/√2.

The three terrains where parity fails show that the selection has to be local: the global count can be even while such a component survives.

## Replay

    python code/robust_parity_v338.py      # the sampled experiment, about 80 seconds on 4 cores
    python -m unittest discover -s tests -p 'test_robust_parity_v338.py' -v

[Recorded output](../results/v338/parity.json). The experiment is seeded and recorded, but not replayed in CI. The square searches are the v334 ones, which are seeded and heuristic.
