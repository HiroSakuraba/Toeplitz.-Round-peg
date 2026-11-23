# v336. Robust squares of grid terrains: the constant appears to be 1/√2

Benjamin John Schulz · November 2026

[v335](v335-even-loop-and-fragility.md) showed that the side-i square v316 guarantees on a terrain can rest on corner contacts and vanish under a perturbation of size 10⁻³. The noncollapse goal needs squares that survive approximation, so this note asks how large the **robust** squares of a terrain are.

**Findings.**
1. Over 282 random terrains and 194 random two-terrain unions, the lowest robust ratio sits just below 1/√2, short of it by an amount proportional to the perturbation size ε, and it approaches 1/√2 as ε → 0. So the evidence points to 1/√2 as a limiting constant, not a bound verified at fixed ε. This is Pettersson–Tverberg–Östergård's conjectured constant for grid polygons, here appearing as a robust constant for terrains.
2. The bound is attained by the terrain with heights (2, 7). There the exact largest square has side i = 2 but is fragile, and the surviving square is a 45° diamond of side exactly i/√2.
3. A natural degree argument has exactly the right constant, but in its naive form a fragile corner square can lie on the boundary of its parameter domain, namely when its roof vertex lies on the excursion. Cutting such squares out is the remaining step of a proof.

## Robust size, measured

For a polygon P, the estimate is the minimum, over 12 seeded perturbations moving every vertex by at most ε = 10⁻³ in max-norm, of the largest inscribed square. Each perturbed polygon is checked to be simple, and squares are computed exactly with the v334 machinery: every isolated square and both ends of every family. An adversarial perturbation could do worse than random ones, so this is an **upper** estimate of the robust size.

| Family | Shapes (i ≥ 2) | Lowest robust M/i | 10% / 50% / 90% quantiles | Robust below exact |
| --- | ---: | ---: | --- | ---: |
| Terrains (heights ≤ 8, widths 2–8) | 282 | 0.70700 | 0.90 / 1.27 / 2.06 | all |
| Two-terrain unions | 194 | 0.70712 | 1.00 / 1.25 / 2.00 | all |

**Stress test of the worst shapes.** For (2, 7), (6, 2) and (3, 3, 4), 200 perturbations each give robust ratios of 0.7069, 0.7068 and 0.7073 at ε = 10⁻³, and 0.7048, 0.7041 and 0.7059 at ε = 10⁻². The deficit below 1/√2 = 0.70711 is proportional to ε, so the robust ratio tends to exactly 1/√2 as ε → 0. The terrain (3, 6, 8) stays at 0.745.

**The sharp example.** The terrain with heights (2, 7) has boundary (0,0) (2,0) (2,7) (1,7) (1,2) (0,2) and i = 2.
- Its exact largest square is the axis-aligned (0,0), (2,0), (2,2), (0,2). It uses three convex corners and is lost under perturbation.
- The square that survives is the diamond (1,0), (2,1), (1,2), (0,1), of side √2 = i/√2. Its vertices lie on the base, the right wall, the **reflex** corner (1,2) and the left wall.

**Conjecture 336.1.** Every grid terrain, and every union of two terrains glued along a base segment, carries a robust inscribed square of side at least i/√2, meaning one present on every sufficiently small registered perturbation. The constant would be sharp.

This would repair the survival lemma of [v334](v334-two-terrain-gluing.md) in its robust form, with constant 1/√2 and no loss from gluing, at least on the families sampled.

## A degree argument with the right constant, and where it breaks

Consider squares with one vertex b on the base and the **opposite** vertex t on the v316 excursion above height h, and require the other two vertices to lie on the curve. Record their signed distances as G(b, t) ∈ ℝ². The diagonal satisfies |t − b| ≥ height(t) ≥ h, so every zero is a square of side at least h/√2. Letting h → i gives exactly the conjectured constant. If G had odd winding around the boundary of the parameter rectangle (base × excursion), the zero would be robust.

**Obstruction.** Let the roof first cross the diagonal y = x at (c, c), and consider the axis-aligned corner square (0,0), (c,0), (c,c), (0,c). It is a zero of G on the side b = (0,0) of the parameter rectangle exactly when (c,c) belongs to the excursion (c ≥ h) and (0,c) lies on the left wall (c is at most the first column's height). The mirror condition gives a zero on the side b = (W,0). The condition can fail: for heights (1, 8, 8, 8) the left corner square has c = 1 < i = 3, so its roof vertex is outside the domain and the left side carries no such zero, while the right corner square (c = 3) does give one. The replay checks the condition and the boundary contacts exactly on several terrains. Where the condition holds, which includes the sharp example (2, 7) on both sides, these boundary zeros are precisely the fragile corner squares of v335, the winding is undefined on the naive domain, and sampled windings (0, ±1, −2) carry no information.

**What a proof needs.** Remove the corner squares that meet the boundary from the parameter domain, for instance by restricting b to the interior of the interior square's base interval. Then show that the new boundary is zero-free and that the winding there is odd. Robustness of the resulting zero would then follow from the degree, giving Conjecture 336.1 for terrains. Gluing would need the same argument with the base replaced by the lower terrain's roof.

## Replay

    python code/robust_squares_v336.py      # the sampled experiment, about 5 minutes on 4 cores
    python -m unittest discover -s tests -p 'test_robust_squares_v336.py' -v

[Recorded output](../results/v336/robust.json). The experiment is seeded and recorded, but not replayed in CI.
