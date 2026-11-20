# v266: a two-corner proof of the hidden-corner lemma, stability, an automatic convex cutoff, and graceful-or-off-disk

**Research program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Five short lemmas with complete pen-and-paper proofs, plus exact-arithmetic support and one reusable exact certificate. Proposition 5 is classical for convex curves and is recorded only as the convex-family consequence for the v260–v264 degree program. **The unrestricted Toeplitz conjecture and the weak Big Peg problem remain open.**

## Summary

1. **Lemma 0 (two near corners).** The v255 hidden-corner lemma has a two-line proof that uses only two corners, and it says more: when a square is smaller than the bound, one of the two corners nearest the disk centre is hidden.
2. **Lemma 2 (stability).** If a convex curve around a unit disk has a square with half-diagonal at most 1 + ε, the square's centre is within √2·ε·(1 + 7ε) of the disk centre. The constant √2 is attained. So a convex curve that does not nearly touch its incircle at four orthogonal points beats the √2 ratio by an explicit margin.
3. **Lemma 1 (automatic cutoff).** A moving convex polygon that contains a disk of radius r0 at all times has every square of side² ≥ 2r0² at all times. `certify_moving_disk` proves the containment exactly for v264-style linear motions. On every v264 test family it gives a side floor 1.8 to 7.6 times larger than v264's projection floor.
4. **Lemma 3 (graceful or off the disk).** In a region where every point of a disk sees the whole boundary, a square with a diagonal through the disk is both graceful (corners in square order along the curve) and big. So ungraceful squares and small squares have both diagonals off the disk. On twin spikes, every small square is ungraceful, as predicted.
5. **Lemma 4 (collapse needs nonconvexity).** Inscribed squares cannot shrink to a point where the region is locally convex with room for a fixed disk. At the scale of each shrinking square, Lemma 0 hides a near corner behind every nearby disk.
6. **Proposition 5 (convex families).** On convex families, squares are graceful and never collapse, so v262's relative quotient parity is a global invariant (odd, classically). The genuinely open part of the program is nonconvex.

**Provenance.** The v255 convex theorem was discussed publicly before v255: David Bernier asked [whether every inscribed square on a convex polygon is large](https://math.stackexchange.com/questions/5140338/does-every-convex-polygon-have-only-large-inscribed-squares) on 12 June 2026, and Christophe Boilley answered positively the same day (recorded in v256). I have not checked whether that answer uses the two-corner argument below.

## Setting

D = D(c, r) is a closed disk. A square has centre p, half-diagonal h, and perpendicular unit diagonal directions u, v. Its corners are p ± h u and p ± h v.

A corner V is **hidden** if it lies in the interior of conv(D ∪ {the other three corners}). Equivalently, V is hidden if it lies in the interior of conv(D ∪ {its two neighbours}). The equivalence holds because any direction n supporting V against its two neighbours lies within 45° of the outward diagonal at V, and then automatically n·V ≥ n·(opposite corner).

The **near corners** are A = p − h·sgn((p−c)·u)·u and B = p − h·sgn((p−c)·v)·v: the corner on the disk-centre side of each diagonal.

## Lemma 0: two near corners

**Lemma 0.** If neither near corner is hidden, then h ≥ r, with equality only if p = c.

*Proof.* Scale and translate so that c = 0 and r = 1. Put X = |p·u| and Y = |p·v|, and use reflections to take u = e1, v = e2 and p = (X, Y) with X, Y ≥ 0. Then A = (X − h, Y), with neighbours (X, Y ∓ h). The directions supporting the square at A are n = (−cos f, sin f) with |f| ≤ π/4. A is not hidden iff some such n has n·A ≥ max over D of n·q, which is 1:

    (h − X) cos f + Y sin f ≥ 1.

Likewise, B = (X, Y − h) is not hidden iff for some |g| ≤ π/4:

    (h − Y) cos g + X sin g ≥ 1.

Add the two inequalities and regroup:

    2 ≤ h (cos f + cos g) − X (cos f − sin g) − Y (cos g − sin f).      (*)

For |f|, |g| ≤ π/4 we have sin ≤ 1/√2 ≤ cos, so both brackets multiplying X and Y are ≥ 0. Hence 2 ≤ h(cos f + cos g) ≤ 2h.

Equality forces cos f = cos g = 1, and then (*) reads 2 ≤ 2h − X − Y, so X = Y = 0. QED.

No hypothesis on where the corners sit is needed: a corner inside the open disk is automatically hidden. This replaces the three-case proof in v255 and sharpens its conclusion from "some corner" to "a near corner". The convex theorem of v255 follows as before. On a convex curve no boundary point is hidden, because conv(D ∪ corners) lies inside the region.

## Lemma 2: stability of the equality case

**Lemma 2.** Suppose neither near corner is hidden and h ≤ r(1 + ε) with 0 ≤ ε < 3 − 2√2 ≈ 0.1716. Then

    |(p−c)·u| + |(p−c)·v| ≤ T(ε)·r,   T(ε) = 2ε(1 + ε) / (1 − ε − 2√ε).

*Proof (r = 1).* From (*), cos f + cos g ≥ 2/h ≥ 2/(1 + ε). Since each cosine is at most 1, each is at least c_ε = (1 − ε)/(1 + ε), and so each |sin| ≤ s_ε = 2√ε/(1 + ε). Then cos f − sin g ≥ c_ε − s_ε = (1 − ε − 2√ε)/(1 + ε) > 0, and similarly for the other bracket. So (*) gives (X + Y)(c_ε − s_ε) ≤ 2h − 2 ≤ 2ε. QED.

**Convex form.** If no corner is hidden (for example, a square inscribed in a convex curve containing D) and ε ≤ 1/100, then

    |p − c| ≤ √2·ε·(1 + 7ε)·r,

and √2·ε·r is attained.

*Proof.* T(1/100) = 101/3950 < 2.557/100, so X, Y < 1. v255's chord step (which needs only "no corner hidden") gives h ≥ X + √(1 − Y²) ≥ X + 1 − Y², so X ≤ ε + T(ε)² ≤ ε(1 + 6.54ε), and the same for Y.

For sharpness, take p = (ε, ε) with u = e1 and h = 1 + ε. The corners (−1, ε), (ε, −1), (1 + 2ε, ε) and (ε, 1 + 2ε) are all unhidden, and |p| = √2·ε exactly. QED.

**Corollary (a convex Big Peg margin).** Let K be convex with D(0, 1) ⊆ K. Let G(K) be the least t such that some square inscribed in the unit circle has all four corners within distance t of ∂K. Then every square inscribed in ∂K has

    h ≥ 1 + min(1/100, G(K)/2.52).

*Proof.* If h = 1 + ε with ε ≤ 1/100, the unit diagonal directions u_k are the corners of a square inscribed in the unit circle. The corners of the actual square are q_k = p + h·u_k on ∂K, with |q_k − u_k| ≤ |p| + ε ≤ ε(1 + √2(1 + 7/100)) < 2.52ε. So G(K) ≤ 2.52ε. QED.

In words, a convex curve beats √2 by an amount proportional to how far it is from touching its incircle at four orthogonal points. Triangles, regular pentagons, and the v263 pentagon all have G > 0.

## Lemma 1: an automatic noncollapse cutoff for moving convex polygons

**Lemma 1.** Let P(t), for t in an interval, be a continuous family of convex polygons. Suppose each P(t) contains a closed disk of radius r0. Then every positive inscribed square of every P(t) has side² ≥ 2r0².

*Proof.* v255 at each time, or Lemma 0. QED.

**Exact certificate.** For a v264 motion P(t) = (1−t)P0 + tP1 whose strict convexity for all t is certified by v264, `certify_moving_disk(P0, P1, c0, c1, r²)` proves that D((1−t)c0 + tc1, r) ⊆ P(t) for all t in a slab. Convexity reduces this to each edge line. Its signed distance L(t) = cross(E(t), c(t) − A(t)) is a quadratic in t, and the edge length satisfies |E(t)|² quadratic. The certificate shows L > 0 and L² − r²|E|² > 0, polynomials of degree at most 4, by exact Bernstein coefficients with dyadic subdivision. `auto_disk_cutoff` chooses the centres and radius from floating-point Chebyshev centres, then certifies exactly. It halves the time slabs or shrinks the radius if needed, and fails closed otherwise.

| v264 test family | Certified floor 2r0² (Lemma 1) | v264 projection floor d²/2 | Ratio |
|---|---|---|---|
| v254 collision window | 0.909216 | 0.258904 | 3.51 |
| v254 full interval [3/5, 13/20] | 0.908845 | 0.258904 | 3.51 |
| translated pentagon | 0.909218 | 0.119304 | 7.62 |
| scaled pentagon (×3/2) | 0.909216 | 0.245410 | 3.70 |
| translated quadrilateral | 3.087892 | 1.724490 | 1.79 |

Each family needed one time slab. The true minimum side² over the v254 window is about 1.113, so the floor is safe. Unlike a projection gap between nonadjacent edges, the disk certificate is unchanged by subdividing edges, so it needs no collinear-refinement compression.

## Lemma 3: graceful or off the disk

**Lemma 3.** Let Ω be a Jordan domain and z an interior point of its kernel (z sees the whole boundary). If a diagonal line of an inscribed square passes through z, the corners occur along ∂Ω in the square's cyclic order.

Suppose further that the open disk int D(c, r) lies in the kernel. Then every inscribed square with a diagonal line meeting int D is graceful and has side ≥ √2·r. Hence every ungraceful square, and every square with side < √2·r, has both diagonal lines at distance ≥ r from c.

*Proof.* Every ray from an interior kernel point meets ∂Ω exactly once. So ∂Ω is a radial graph about z, and the order along ∂Ω is the angular order around z. The two corners on the line through z are distinct boundary points on that line, so they lie on opposite rays from z, at angles θ and θ + π. The other two corners lie on opposite sides of that line, one in each open half-turn. So the angular order is the square's order.

For the second part, int D lies in the interior of the kernel. Every point of the open chord of int D on that diagonal therefore lies between the two corners, which gives v255's chord inequality h ≥ x + √(r² − y²) for that diagonal. If both diagonals meet int D, v255's algebra gives h ≥ r. If only one does, the other is at distance x ≥ r, so h > r. QED.

**Twin spikes (numeric, circumscribed construction from v255).**
- Each twin spike (tip gaps 0.3, 0.6, 1.0 and 1.6) has five squares.
- The one big square is graceful, with both diagonals crossing the disk.
- All four small squares are ungraceful, with both diagonals off the open disk.
- Two of those small squares have a diagonal exactly tangent to the disk, lying along a spike edge, which is the boundary case of the lemma.

So the ungraceful sectors flagged in v262 (item 2) and v263 (direction 2) are real, and in disk-kernel regions they are confined to squares whose diagonals avoid the disk.

## Lemma 4: collapse needs nonconvexity at every scale

**Lemma 4.** Let Ω be a Jordan domain and B an open disk. Suppose U = Ω ∩ B is convex and contains a closed disk of radius ρ. Then every square inscribed in ∂Ω with all four corners in B has side ≥ √2·ρ.

*Proof.* A point of ∂Ω ∩ B is a boundary point of U. The closure of U is compact and convex and contains the disk. Apply v255, or Lemma 0. QED.

**Corollary (collapse locus).** Suppose Jordan curves γ_n converge uniformly and carry inscribed squares S_n that shrink to a point x. Then there are no ρ, ρ' > 0 such that, for infinitely many n, the region of γ_n inside B(x, ρ) is convex and contains a disk of radius ρ'. Moreover, for every disk D' inside the region of γ_n whose radius exceeds the half-diagonal of S_n, a near corner of S_n relative to D' is hidden behind D' and its two neighbours (Lemma 0).

So the uniform-noncollapse bottleneck named in v262 (item 3) and v263 (direction 3) only has to be fought at points that are nonconvex at every scale, and there a pocket must hide a near corner at the scale of each shrinking square.

## Proposition 5: convex families globalize the degree program

**Proposition 5.** Let K_t, for t in a compact interval, be a continuous family of convex Jordan curves with inradius ≥ r0 > 0. Then:
- every inscribed square is graceful, since four points in convex position appear in the same cyclic order along a convex curve;
- every inscribed square has side ≥ √2·r0 (Lemma 1);
- consecutive corners are therefore at least √2·r0 apart along the curve.

So the zeros of the square-defect map lie in one compact region of the positive cyclic-order sector, away from the collision boundary, for all t. The mod-2 degree on that region is defined and independent of t, and any v262-style relative quotient parity on an atlas containing all squares equals it. For generic members it is the parity of the number of squares, which is odd. That oddness is classical and not claimed here.

This is consistent with v263's complete count: 3 squares before the collision, 2 at it, and 1 after.

The proposition closes, for convex families only, the collapse and sector items listed in v262 §6 and v263 §7. **The open part of the program is nonconvex.**

## Computational support

Exact checks use the standard library (`python code/convex_lemmas_v266.py`):

- **Lemma 0.** All 17,328 rational grid squares with h < 1, centred anywhere, have a hidden near corner.
- **Lemma 2.** T(1/100) = 101/3950 and T(1/400) = 401/71800, exact because √ε is rational. On a rational grid of near-equality squares (h in [1, 1 + ε], nine exact unit directions), all 1,917 squares with unhidden near corners satisfy X + Y ≤ T(ε). All of them also satisfy the convex centre bound, and the largest ratio |p|²/(2ε²) is exactly 1, which is the sharpness case.
- **Sharpness example.** p = (ε, ε), h = 1 + ε has no hidden corner, and |p|²/ε² = 2.
- **Pochuev's square.** Its single hidden corner is a near corner.
- **Lemma 1.** Moving-disk certificates for all five v264 families, with the table above.

Numeric (`--numeric`, needs numpy): Lemma 3 on the twin spikes. All results are in [the JSON record](v266-convex-lemmas-checks.json).

## Limits

- Lemma 0 is a new proof of a known convex statement (Bernier/Boilley; v255). Its near-corner sharpening and Lemma 2's stability are, as far as I can tell, new, but the convex area is classical enough that a literature check is worthwhile.
- Lemma 1's certificate handles v264 linear motions using one affine moving centre and one radius per slab, and fails closed otherwise. Strict convexity for all t must be certified separately, which `auto_disk_cutoff` does by calling v264.
- Lemmas 3 and 4 constrain where small or ungraceful squares can be. Neither gives existence for nonconvex curves.
- Proposition 5 is the convex case only.

## Reproduce

    python code/convex_lemmas_v266.py                       # exact part (standard library)
    python code/convex_lemmas_v266.py --numeric --write     # add the twin-spike sectors, rewrite the JSON
    python -m unittest -v tests.test_convex_lemmas_v266
