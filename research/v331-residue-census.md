# v331. A complete residue census of small squares on spiral polygons

Benjamin John Schulz · November 2026

**Question (from the [roadmap](positive-window-noncollapse-roadmap.md)).** On the retuned spiral, can the action residue alone decide which small square could carry the upper spectral value at the square angle? Or, failing that, does it reduce the question to a few explicit generators?

**Answer: no.** Once the spiral winds enough that the area primitive swings by more than twice the spectral gap, the number of small squares the residue cannot exclude grows explosively. It is 0 up to N = 16 turns, then 6, 2,617 and 158,721 at N = 20, 24 and 30 (all certified exactly), and about 3.4 million at the retuned N = 40 (a floating count). The marked v324 square is one of these at N ≥ 24 but is not distinguished. Deciding the carrier on such curves needs information beyond the action, such as gradings or the Floer differential. No Floer homology is computed here.

## Setup

**Fixture.** The boundary of [0,4]² with one v324 gadget of N turns at scale ε = 1/2, attached at the middle of the left side. Placed there, the gadget leaves the diameter at 4√2, so R² = 8 exactly, while the area is A ≈ 16.711. The spectral gap

    G = min(A/2, A − πR²/2) ≈ 4.1447

is therefore positive already at this coarse scale. At the corner placement of v330 it is positive only for k ≥ 3. The marked square's residue is (19N/24 − 1/(384N) − 1/2)/4, which exceeds G once N ≥ 22. N = 40 is the retuned value at which it approaches A/2.

**Reduction.** By spectrality, the upper spectral value at θ = π/2 is the action of some inscribed square (Greene–Lobb, [arXiv:2604.17116](https://arxiv.org/abs/2604.17116)). By subadditivity and the Lipschitz cap, its distance from AZ is at least G. By (323.1), the action residue of a square is −[(F(c₀)−F(c₁)) + (F(c₂)−F(c₃))] modulo A, where F is the boundary area primitive at its counterclockwise contacts. So a square whose residue lies within distance G of AZ cannot carry the spectral value. This is Lemma 330.1 at a single angle, applied square by square rather than through an oscillation bound. "Small" means half-diagonal below 1.

## Method and completeness

1. **Floating search.** For every ordered pair of edges carrying two adjacent corners, and both orientations, the third corner's locus is intersected exactly with every edge. Vectorized over edges, the fourth corner's segment is then intersected with every edge. Tolerances are widened (10⁻⁹), so the search over-reports. Each candidate's size and residue are evaluated, and only squares with q < 1 and residue distance at least G − 0.05 are kept.
2. **Exact certification.** For each kept edge quadruple the two line conditions are solved in rationals. All four contacts are checked on their closed edges, and the residue is recomputed exactly from the contact edges. Survivors are those with exact distance at least G. No certified candidate had a singular system.
3. **Agreement with the exhaustive census.** On an irregular 12-gon the search finds exactly the 7 squares of the repository's exhaustive v223 census (a regression test). On a 21-vertex polygon built from a one-turn gadget it finds exactly the same 128 squares, in 0.3 s instead of 71 s.
4. **Cases the search skips.** Collinear configurations are flagged rather than solved. `degenerate_audit` re-solves every flagged case exactly. At every N tested, exactly one positive-size family of small squares exists: the tilted squares inscribed in the unit box formed by the gadget's outer frame and the left side, with parameter λ ∈ [0, 1/2] and sides √½ to 1. Its residue is constant along the family, as it must be for a connected family of critical points, and equals 0.2110, far below G. Every flagged single-point configuration with residue at least G (0 at N = 20, 3 at N = 24 and at N = 30) is already among the certified survivors. All other flags are empty or have side zero.

## Results

| Turns N | Vertices | Marked-square residue | Max residue distance found | Small-square survivors (distance ≥ G ≈ 4.1447) |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 45 | 0.6665 | 1.128 | 0 |
| 8 | 77 | 1.4583 | 1.919 | 0 |
| 16 | 141 | 3.0416 | 3.503 | 0 |
| 20 | 173 | 3.8333 | 4.294 | 6 |
| 24 | 205 | 4.6250 | 5.086 | 2,617 |
| 30 | 253 | 5.8125 | 6.273 | 158,721 |
| 40 (retuned) | 333 | 7.7917 | 8.253 | ≈ 3.45 million (floating count) |

Survivor counts for N ≤ 30 are exact and certified. At N = 40, certifying every survivor exceeds 16 GB of memory, so only raw floating hits at or above G are counted: 27,585,728. A square with all contacts in edge interiors is hit exactly 8 times (4 starting corners × 2 orientations); the measured ratio is 8.27 at N = 24 and 8.01 at N = 30. The estimate is 27,585,728 / 8.

The six survivors at N = 20 lie between the gadget's outer frame and the spiral center, with sides 0.37 to 0.63. The marked square is not among them, since its residue 3.833 is below G.

**The survivors are robust.** At N = 24, 2,553 of the 2,617 survivors have all four contacts in the interiors of edges. Each solves a nonsingular exact system, so it persists under sufficiently small smoothing (implicit function theorem), and its residue changes continuously. The explosion is not an artifact of polygon corners.

**Genericity check.** The N = 24 census was rerun on three modified polygons, each verified to be a Jordan polygon, with every survivor certified exactly. The modifications are of size 10⁻⁴, well below the spacing of about 1.3·10⁻³ between spiral arms.

| Polygon | Vertices | Survivors | Base survivors with a survivor within 2·10⁻³ | New survivors |
| --- | ---: | ---: | ---: | ---: |
| Base fixture | 205 | 2,617 | — | — |
| Every vertex moved by a seeded rational offset in [−10⁻⁴, 10⁻⁴]² (seed 1) | 205 | 2,489 | 2,598 | 2 |
| The same, seed 2 | 205 | 2,610 | 2,612 | 0 |
| Every corner cut at max-norm distance 10⁻⁴ along both edges | 410 | 1,940 | 2,577 | 0 |

Random perturbations change the count by at most 5%. The losses come from survivors just above the gap whose residues shift slightly below it. Cutting every corner removes every contact at an original vertex, and about three-quarters of the survivors remain. The matching column is approximate: the radius exceeds the arm spacing, so a match may be a neighbouring square. The survivor counts do not depend on the matching.

## Interpretation

- **Threshold behaviour.** Below the threshold, the census agrees with Lemma 330.1: every small square is excluded, and the largest residue distance stays below G. Above it, residues fill the spectral window. By N = 40 the largest found distance, 8.253, is near A/2 ≈ 8.356, where the retuned family's spectral value lies.
- **The roadmap's v331 decision.** "If both large and small branches remain, isolate a specific collision or action ambiguity" assumed the ambiguity would be one explicit generator. It is not. It is a large family of generic, robust candidates, and the marked square is not distinguished among them.
- **Consequence for the program.** v324 showed that a single collapsing generator can keep an action gap. This census shows that above the threshold this is typical, not exceptional, among small squares. Arbitrary Jordan curves can carry such gadgets at every scale. So no argument that uses only actions and residues can select spectral carriers in general. A carrier-selection argument needs more structure: for instance, the Maslov grading of the grading-two class, or which generators the Floer differential lets represent it.
- **What this does not show.** Survivors are squares the residue cannot exclude. The spectral value may still be carried by a large square (for example, one of the squares inscribed in the base square). The census shows the residue test is inconclusive, not that the carrier is small.

## Mod-2 Floer grading of the survivors

The spectral value at π/2 belongs to the grading-two class. The integer grading of a generator is the Maslov index of its preferred capping, which Greene–Lobb do not give in computable form and which is not computed here. Its parity, however, is determined by the local intersection sign of L₀ = γ×γ and L₁ = R_{π/2}(L₀) at the generator, up to one universal constant. At θ = π/2 the rotation is rational, so the sign is an exact 4×4 determinant built from the edge directions at the four contacts. The product orientation of γ×γ does not depend on the direction of γ, so the sign is intrinsic.

- **Which generator.** A square (c₀, c₁, c₂, c₃), listed counterclockwise, gives two inequivalent generators: x = (c₁, c₃) = R(c₀, c₂), with action residue r, and x′ = (c₂, c₀) = R(c₁, c₃), with residue A − r. Each also has a swapped twin (w, z) with the same sign and action. The spectral value lies in [A/2, πR²/2], so exactly one generator of each surviving square can carry it. The sign convention of (323.1) is the one that gives ℓ(θ) = ρ²θ on the circle.
- **Calibration.** On an ellipse at θ < π/2, the major-axis rectangle continues as θ → 0 to the diameter, the maximum of H, so it has grading 2. The minor-axis rectangle continues to a saddle of H, so it has grading 1. Their signs are −1 and +1 at every angle tested and for aspect ratios 1.2, 2 and 5. So grading 2 has sign −1.
- **Check.** The homology is F₂ in gradings 2 and 1, so the signs over all generators of any generic polygon sum to zero. They do, exactly, on an irregular 12-gon (7 squares) and a convex pentagon (1 square), in both orientations. The pentagon's single square has one generator of each parity.

| Turns N | Survivors | Window generator even (may carry grading 2) | Odd (excluded) | Contact at a polygon vertex (parity undefined) |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 6 | 0 | 0 | 6 |
| 24 | 2,617 | 1,280 | 1,273 | 64 |
| 30 | 158,721 | 79,278 | 79,283 | 160 |

Parity removes about half of the candidates and leaves tens of thousands at N = 30. The near-exact balance between even and odd is consistent with the Euler-characteristic constraint. So the grading modulo 2 does not select the carrier either. The integer grading could do more, since only grading exactly 2 qualifies among the even ones, but it requires the Maslov index of each preferred capping. Beyond that lies the Floer differential, which decides which grading-two generators actually represent the class.

## Replay

    python code/residue_census_v331.py              # N = 4, 8, 16, 20 (CI replay, about 15 s on 4 cores)
    python code/residue_census_v331.py --large      # N = 24, 30 exact and N = 40 counted (minutes)
    python code/residue_census_v331.py --genericity # perturbed and corner-cut N = 24 (about 2 minutes)
    python code/floer_parity_v331.py                # calibration, Euler checks, parity at N = 20, 24 (CI replay)
    python code/floer_parity_v331.py --large        # parity at N = 30 (about 9 minutes)
    python -m unittest discover -s tests -p 'test_floer_parity_v331.py' -v
    python -m unittest discover -s tests -p 'test_residue_census_v331.py' -v

The search needs numpy. [CI replay output](../results/v331/census.json) · [parity replay](../results/v331/parity.json) · [Large runs](../results/v331/large-runs.json) and the [genericity check](../results/v331/genericity.json) and the [N = 30 parity run](../results/v331/parity-large.json), recorded and not replayed in CI.
