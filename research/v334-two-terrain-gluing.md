# v334. Gluing two grid terrains: the v316 bound fails, and the remaining lemma

Benjamin John Schulz · November 2026

**Question.** [v316](v316-roof-shear.md) proves M ≥ i for every grid terrain with a flat base. Here M is the largest side of a square with all four vertices on the boundary, and i is the largest side of an open axis-aligned square in the interior. Glue an upper terrain U and a lower terrain D along their common base segment. Does M(P) ≥ k·i(P) hold for the union P, with an absolute k > 0 independent of the mesh?

**Findings.**
1. **k = 1 fails.** Gluing a 2 × 1 block under the end of a terrain's base gives M/i = 5/7, which is confirmed by an independent exhaustive census.
2. **Straddle reduction (proved).** max(i(U), i(D)) ≥ i(P)/2, and this is sharp. So the gluing lemma reduces to a survival lemma about one terrain at a time.
3. **Search evidence.** Over seeded random, hill-climbing and structured searches, no ratio below 5/7 was found. That is consistent with the conjectured 1/√2 of Pettersson–Tverberg–Östergård for all grid polygons.

No positive constant is proved here.

## Setup and exact computation

The upper terrain has integer column heights over its base interval and the lower terrain integer depths under its own. The bases overlap in at least one column, the seam. The union is a vertically convex polyomino P = {g(x) < y < f(x)}, with g ≤ 0 ≤ f. The seam lies in the interior of P.

- **i(P)** equals the maximum over column intervals I of min(|I|, min_I f − max_I g). An open axis-aligned square spans columns whose heights admit it, and conversely each interval gives a lattice square. The largest such square may straddle the seam.
- **M(P)** is computed over *every* inscribed square, not only lattice ones. The repository has no lattice-vertex lemma for general grid polygons; v316 proves integrality only for base-supported squares of a single terrain. Isolated squares come from the v331 search with exact certification. Positive-dimensional families, which polyominoes have in abundance, contribute both endpoints of their parameter segment. That suffices, because the squared side is convex along a family.
- **Validation.** M ≥ i holds on random single terrains, as v316 requires. The witness below is confirmed by the exhaustive v223 census: with cutoff 25 no square exceeds side 5, and with cutoff 24 a positive-dimensional family reaches side² = 25.

## The constant 1 fails

**Witness A (ratio 5/7).** Take upper heights (8, 8, 8, 7, 7, 7, 7, 3) on columns 0–7, and lower depths (1, 1) on columns 6–7. The boundary is

    (0,0) (6,0) (6,−1) (8,−1) (8,3) (7,3) (7,7) (3,7) (3,8) (0,8).

- **Interior square.** i = 7, attained by (0,7) × (0,7), and i(U) = 7 as well.
- **What gluing destroys.** The v316 square of U has side 7 on the base [0,7]. Its corner (7,0) now lies in the open seam (6,8) × {0}.
- **The largest surviving square** is the tilted 3-4-5 square with vertices (4,0), (7,4), (3,7), (0,3), of side 5. Its vertices lie on the exposed base, the wall x = 7 above the short last column, the corner (3,7), and the left wall.

So M/i = 5/7. The lower terrain is only a 2 × 1 block, so the loss is caused entirely by the change of boundary at one base corner.

**Witness B (ratio 2/√5).** Take upper heights (6, 6, 6, 6, 5, 5) and lower depths (5, 1, 1) starting at column 3. Then i = 5, and the v316 square of U loses its corner (5,0). The largest surviving square is (2,0), (6,2), (4,6), (0,4), of side √20, which is again tilted.

In every witness found, the loss happens the same way. The base corner of a terrain's v316 square is swallowed by the seam, and the best replacement is a **tilted** square. An argument that only produces axis-aligned squares, like v316's shear, cannot reach constant 1 for glued terrains.

## Straddle reduction

**Lemma 334.1.** max(i(U), i(D)) ≥ i(P)/2.

*Proof.* Let Q = (a, a+m) × (c, c+m) be a largest interior square. If c ≥ 0, then Q ⊂ U and i(U) ≥ m, and similarly if c + m ≤ 0. Otherwise Q straddles y = 0, so each of its columns belongs to both terrains, and Q splits into parts of heights h₁ = c + m above and h₂ = −c below, with h₁ + h₂ = m. Lowering the upper part to the base by vertical filling gives U an interior rectangle of width m and height h₁. Likewise D gets one of width m and depth h₂. One of h₁, h₂ is at least m/2, and a rectangle of width m ≥ h contains a square of side h. ∎

The bound is sharp: 38 of 20,000 random glued pairs attain equality, with the largest square split evenly by the seam.

**Consequence.** The two-terrain question reduces to a survival lemma: is there an absolute k′ > 0 such that, for every terrain U and every lower terrain D glued to it, M(U ∪ D) ≥ k′·i(U)? Such a k′ gives M(P) ≥ (k′/2)·i(P). The witnesses show that any such k′ is at most 5/7.

## Search evidence

| Search | Instances | Lowest M/i found |
| --- | ---: | ---: |
| Random glued pairs (heights ≤ 6, widths ≤ 7) | 4,000 | 2/√5 ≈ 0.894 (2 instances below 1) |
| Hill-climbing, 24 seeded runs of 700 steps (heights and widths ≤ 12) | about 17,000 evaluations | 5/7 ≈ 0.714 |
| Structured family (three roof blocks plus an end block, scale up to 40) | 5,997 | 0.738 |

The structured family was a one-off exploration and is not part of the recorded searches. The searches are heuristic, so a lower ratio may exist. Nothing below 1/√2 ≈ 0.707 was found.

## What a proof needs

- **A tilted-square argument.** The witnesses force any proof of the survival lemma to produce tilted squares. A natural route is to deform the lower terrain from the flat base, where v316's marked shear problem has odd degree, to the actual D. Then show that a degree region loses odd parity only when a square crosses the seam, and that at that moment a replacement square with four contacts on exterior edges stays above the cutoff. Witness A gives a concrete target: its replacement uses the exposed base, a wall, a roof corner and the opposite wall.
- **Many seams.** A general simple polyomino decomposes into a tree of terrains. Repeating Lemma 334.1 loses a factor of 2 at every seam, so a proof for arbitrarily many seams needs an argument that does not compound the loss.

## Replay

    python code/terrain_gluing_v334.py              # both witnesses (CI replay)
    python code/terrain_gluing_v334.py --search     # seeded random and hill-climbing searches (about 30 minutes)
    python -m unittest discover -s tests -p 'test_terrain_gluing_v334.py' -v

[Witness report](../results/v334/witnesses.json) · [Search record](../results/v334/searches.json). The search record is not replayed in CI.
