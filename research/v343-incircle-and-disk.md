# v343. The incircle reading of 1/√2, and a disk version of the robust bound

Benjamin John Schulz · November 2026

[Conjecture 336.1](v336-robust-terrain-squares.md) asks for a robust inscribed square of side at least i/√2, where i is the side of the largest interior axis-aligned square. This note gives the constant a geometric meaning and tests a stronger, coordinate-free version.

**Findings.**
1. **i/√2 is the side of the square inscribed in the incircle of the interior square.** In the sharp example (2, 7), the surviving diamond is exactly that square: its four vertices lie on the incircle of [0, 2]².
2. **Conjecture 343.1 (disk version).** Let r be the curve's inradius, the radius of its largest inscribed disk. Then the curve carries a robust inscribed square of side at least √2·r, the side of the square inscribed in that disk. Since i/√2 ≤ √2·r ≤ i, this implies Conjecture 336.1. It also needs no choice of axes.
3. **The sampled evidence is consistent with it.** Over 576 sampled shapes (282 terrains, 194 two-terrain unions, 60 random star polygons and 40 circles with slots) and 12 named shapes, the robust estimate never falls below 0.9998·√2r. The deficit is of the size of the perturbation. The near-equality cases are regular polygons and terrains with i = 2r, where the robust square is the square inscribed in the incircle.
4. **What is already proved.** Theorem 338.1 gives the disk version on interior terrains, and v337 gives it, with equality, for (2, 7) and the glued witness A. The [v341](v341-ledge-and-component-certificates.md) certificates reach i/√2 on all 476 v336 rows, which is the disk version wherever i = 2r. [v342](v342-unmarked-components-and-square-frames.md) proves that a square boundary of side S has robust size exactly S/√2 = √2·r, so the disk version holds there with equality. Nothing is proved in general.

## The incircle reading

An open axis-aligned square of side i in the interior contains its incircle, of radius i/2. Every square inscribed in that circle has side i/√2. So Conjecture 336.1 says: the curve keeps, robustly, a square at least as large as the squares of the incircle of its largest interior square.

In the (2, 7) terrain the interior square is [0, 2]², with incircle centred at (1, 1) and radius 1. The robust square from [v336](v336-robust-terrain-squares.md) and [v337](v337-terrain-component-degree.md) is the diamond (1,0), (2,1), (1,2), (0,1). Its vertices are the four points where the incircle touches the square's sides. In the large-square family of v337, which runs from the corner square to the diamond, the two ends are therefore the circumscribed square and the inscribed square of one circle. A perturbation keeps one end or the other. The replay checks the incircle contacts exactly, and that the incircle is also the largest disk in the terrain.

## The disk version

**Conjecture 343.1.** Every Jordan curve with inradius r carries a robust inscribed square of side at least √2·r. Robust means here that for each s < √2·r, every sufficiently close registered Jordan curve has an inscribed square of side greater than s.

**Relation to Conjecture 336.1.**
- A disk of radius r contains an axis-aligned square of side √2·r, so i ≥ √2·r.
- An interior square of side i contains a disk of radius i/2, so r ≥ i/2.

Together these give i/√2 ≤ √2·r ≤ i. The disk version is therefore at least as strong as 336.1. They coincide exactly when i = 2r. In the sample, i/r ranges over [1.60, 2] on terrains and [1.58, 2] on unions.

**What is proved.**
- [Theorem 338.1](v338-robust-terrain-proof.md) gives a robust square of side at least i ≥ √2·r on interior terrains, so the disk version holds there.
- v337's two fixtures have i = 2r, and their robust size i/√2 equals √2·r. So the disk version holds with equality there.
- A square boundary of side S has r = S/2, and [v342](v342-unmarked-components-and-square-frames.md) proves its robust size is exactly S/√2 = √2·r, with the finite-error form (S − 2ε)/√2. It is a second equality case besides the circle.
- The v339–v340 certificates, and the [v341](v341-ledge-and-component-certificates.md) certificates on all 476 v336 rows, reach i/√2. That gives the disk version only where i = 2r. On the remaining rows, where i < 2r, the disk version asks for more, between i/√2 and √2·r.

For a circle the bound is attained, and it cannot be improved by any constant factor. An ellipse with semi-axes a ≥ b has r = b and an inscribed square of side 2ab/√(a² + b²) ≥ √2·b, with equality only for the circle.

## Test

The test compares the v336 robust estimate with √2·r. The robust estimate is the minimum over 12 seeded perturbations of size 10⁻³ of the exact largest square, and it is an upper estimate of the robust size. The inradius is bracketed by a branch-and-bound over square cells. The signed distance to the boundary is 1-Lipschitz, so a cell of side s whose centre has value d holds no value above d + s/√2. This gives a certified interval [r_lo, r_hi], up to floating point, and the ratio is computed with r_hi, so an inradius error can only lower it. On the recorded rows the interval has relative width at most about 10⁻⁵. A ratio clearly below 1 would point to a counterexample, which would then need an exhaustive census to confirm.

| Family | Shapes | Lowest robust / √2r | Below 0.99 |
| --- | ---: | ---: | ---: |
| Terrains (heights ≤ 8, widths 2–8) | 282 | 0.9998 | 0 |
| Two-terrain unions | 194 | 1.0000 | 0 |
| Random star polygons (20 vertices, radii in [2, 4]) | 60 | 1.0916 | 0 |
| 32-gon circles with 3–6 inward slots | 40 | 1.2043 | 0 |

**Named shapes.**
- **Near the bound.** (2, 7) at 0.9998, regular 8-, 12-, 16- and 20-gons at 1.0002–1.0015, a 16-gon with a thin stick attached at 1.0076, and the glued witness A at 1.0101.
- **Far above it.** Polygonal ellipses at 1.18–1.39, and squares with a slot at 2.68–3.82.

The lowest terrains, (2, 7), (6, 2), (2, 5) and (3, 3, 4), all have i = 2r. Their robust square is the square inscribed in the incircle.

**Shapes tried as counterexamples.** A disk with a thin stick attached interrupts the circle's family of squares at one point only, and the squares rotate past it (ratio 1.0076). Inward slots cut the circle in several places, but they lower r more than they lower the squares (ratios at least 1.20). No sampled shape came close to violating the bound except where equality is expected.

## Caveats

The robust estimate uses random perturbations and the v334 square search, which is seeded and heuristic. The inradius interval is computed in floating point, not exact arithmetic. The disk version for all Jordan curves would be a strengthening of the square-peg conjecture with an explicit size bound, so it is at least as hard as the open problem.

## Replay

    python code/incircle_disk_v343.py      # about 10 minutes on 4 cores
    python -m unittest discover -s tests -p 'test_incircle_disk_v343.py' -v

[Recorded output](../results/v343/disk.json). The experiment is seeded and recorded, but not replayed in CI.
