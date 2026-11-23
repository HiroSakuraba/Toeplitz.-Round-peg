# v335. A localized half-square proposal fails; v316's large squares are fragile

Benjamin John Schulz · November 2026

This note has two parts, each correcting a premise before it is built on.
- **Part A.** A localized "even-loop certificate" based on the half-square map, proposed in an outside idea sheet, is false as stated. A two-loop version is correct, and it is proved here with an explicit perturbation margin.
- **Part B.** The 5/7 loss in [v334](v334-two-terrain-gluing.md) is not caused by gluing. The terrain's square of side i rests on non-generic corner contacts, and a perturbation of size 10⁻³ already destroys it. The survival lemma should therefore be stated for robust squares.

## Part A. The even-loop proposal

**The proposal.** For x, y on the parameter circle, put

    S₁ = γ(x) + J(γ(y) − γ(x)),  S₂ = γ(y) + J(γ(y) − γ(x)),  h(x, y) = (f(S₁), f(S₂)),

where f is the signed distance to γ (negative inside) and J is the quarter turn. A zero of h is a square with adjacent vertices γ(x) and γ(y). Let M_a be the band y − x ∈ [a, 1 − a]. The proposal asserted:
1. if h has no zeros on M_a, every core loop in M_a has odd winding;
2. hence a core loop in M_a with even winding certifies a square of side at least m_γ(a), the least chord between parameters at distance at least a.

**Counterexample.** Take the unit circle parametrized in turns, and 1/4 < a < 1/2.
- Adjacent vertices of an inscribed square are 1/4 or 3/4 of a turn apart, so h has no zeros on M_a.
- On the antipodal core loop y = x + 1/2, S₁ = γ(x)(1 − 2i) and S₂ = −γ(x)(1 + 2i) both lie at radius √5. So h ≡ (√5 − 1, √5 − 1), and its winding is zero. This contradicts claim 1.
- At a = 3/10, claim 2 would certify a square of side 2 sin(3π/10) ≈ 1.618. But every square inscribed in the unit circle has side √2 ≈ 1.414.

**The same happens on a generic curve.** The circle has a continuum of squares, so the replay repeats the computation on the ellipse with semiaxes 1.05 and 1, which has a single transverse square. Every core loop y = x + c with c ∈ {0.3, 0.4, 0.5, 0.6, 0.7} has winding zero, and |h|_∞ stays at least 0.257 along it. That margin is 1.196 on the antipodal loop.

**Lemma 335.1 (two-loop version).** Let ℓ₀ and ℓ₁ be disjoint core loops in the parameter torus minus its diagonal, bounding an annulus A. Suppose the winding of h around ℓ₀ differs from its winding around ℓ₁. Then h has a zero in A, and γ inscribes a square of side at least m_A, the least chord |γ(x) − γ(y)| over (x, y) ∈ A.

Moreover, let μ be the minimum of |h|_∞ on ℓ₀ ∪ ℓ₁. Every Jordan curve γ′ whose registered parametrization is within δ < μ/5 of γ in sup norm inscribes a square of side at least m_A − 2δ.

*Proof.* If h had no zero on A, it would map A into ℝ² \ {0}. Since ∂A = ℓ₁ − ℓ₀, the two windings would then be equal.

For γ′, the corners move by at most 3δ: |S₁′ − S₁| ≤ |γ′(x) − γ(x)| + |(γ′(y) − γ′(x)) − (γ(y) − γ(x))| ≤ 3δ, and likewise for S₂. The signed distances satisfy |f′ − f| ≤ 2δ everywhere. Their absolute values differ by at most the Hausdorff distance δ, and their signs can differ only within δ of both curves, by straight-line interpolation of the parametrizations. Since f′ is 1-Lipschitz, |h′ − h|_∞ ≤ 3δ + 2δ = 5δ < μ. So straight-line interpolation from h to h′ has no zeros on ℓ₀ ∪ ℓ₁, and the windings are preserved. Then h′ has a zero in A, and its chord is at least m_A − 2δ. ∎

The localized single-band statement does not follow from the even-loop margins alone. Two loops whose windings are actually different are needed. As the attachment itself observed, m_A collapses on spirals, so the certificate is only as good as the large-scale geometry.

## Part B. Fragility of the v316 square

Take the v334 witness terrain, with upper heights (8, 8, 8, 7, 7, 7, 7, 3) and boundary (0,0) (8,0) (8,3) (7,3) (7,7) (3,7) (3,8) (0,8). Then i = 7, and its largest boundary square is (0,0), (7,0), (7,7), (0,7) of side 7, exactly as v316 guarantees. Two of its contacts, (0,0) and (7,7), are polygon corners, and (7,0) touches the base only as the end of the square's lower side.

All values below are exact, using the v334 computation of every inscribed square including families.

| Deformation | Largest square |
| --- | --- |
| Bump of depth t under the base on [6, 8], for t = 1/1000, 1/100, 1/10, 1/2, 1 | side 5 for every t > 0 (side 7 at t = 0) |
| Depth-1 bump on [c, 8] with c = 7, which leaves (7,0) on the boundary | side 7 |
| … with c = 6.9, 6.5 or 6 | side 5 |
| … with c = 4 | side √32 ≈ 5.66 |
| … with c = 1 | side √50 ≈ 7.07 |
| Terrain alone, each vertex moved by a seeded offset of size at most 10⁻³ (12 seeds) | side ≈ 5 for 6 seeds, ≈ 7 for 6 seeds |

Removing the single contact (7,0) by a change of registered size 10⁻³ drops M from 7 to 5. Generic perturbations of the terrain itself, with no gluing, do the same about half the time. So the v334 loss reflects how unstable the side-i square of a terrain is, not something special about seams.

**Consequences.**
- *The survival lemma should concern robust squares*: squares that persist under every sufficiently small registered perturbation, for instance squares of odd local degree. On this terrain the robust size is at most 5 = (5/7)·i, which equals the glued value. A natural target is therefore that every terrain carries a robust square of side at least k·i, and that robust squares survive gluing up to a constant.
- *v316's own noncollapse statement is consistent with this.* It promises odd local degree only above the cutoff i/2, which lies below 5.
- *Exact grid results overstate robust size.* Exact grid-polygon results such as the Pettersson–Tverberg–Östergård bound measure M on one exact polygon. For the noncollapse goal only robust squares count, and the gap can be a constant factor.

## Replay

    python code/even_loop_fragility_v335.py
    python -m unittest discover -s tests -p 'test_even_loop_fragility_v335.py' -v

[Recorded output](../results/v335/report.json).
