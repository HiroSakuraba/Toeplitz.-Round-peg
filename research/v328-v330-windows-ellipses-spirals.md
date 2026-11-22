# v328–v330. Angular windows, ellipse calibration, and a carrier exclusion on the spirals

Benjamin John Schulz · November 2026

**Scope.** This is the first implementation tranche of the [positive-window roadmap](positive-window-noncollapse-roadmap.md). The three parts are:
- **v328:** a certified enclosure for the window average D_m from interval data.
- **v329:** an exact calibration on ellipses, which pins D_m between two explicit integrals that both converge to the squared radius of the inscribed square.
- **v330:** an exclusion of every small rectangle as a spectral carrier near the square angle. It holds on the original v324 spiral family, for every k ≥ 7, and on a fixed nonrectifiable target. The retuned family fails the same test, as the roadmap predicted.

Every computation is exact rational arithmetic. pi and arctan enter only through alternating-series enclosures. **No Floer homology is computed.** The spectral statements rely on Greene–Lobb's Theorem 5.1 ([arXiv:2604.17116](https://arxiv.org/abs/2604.17116)), the cap identity (323.1) and the localization statement v325 from [v323–v327](v323-v327-cap-winding.md). The curves treated here all inscribe large squares already, so nothing below is new information about square existence for them. The point is to test the window and carrier machinery on families where the answer is known.

Notation follows v318 and the roadmap: w_m = sin^(2m) θ, Z_m = π C(2m,m)/4^m, D_m(L) = Z_m^(-1) ∫ w_m dL, A is the enclosed area, R is half the diameter, and q is the half-diagonal (circumradius) of a rectangle.

## v328. Certified window enclosures

Use the half-angle variable t = tan(θ/2). Then sin θ = 2t/(1+t²) is rational in t, the square angle is t = 1, and θ = π is t = ∞. The kernel increases in t on [0,1] and decreases on [1,∞).

**Lemma 328.1.** Let L be nondecreasing and K-Lipschitz in θ. Take a mesh 0 = t_0 < … < t_N = ∞ containing 1, with L(θ(t_i)) ∈ [l_i, u_i]. On each cell let w⁻_i and w⁺_i be the smaller and larger endpoint kernel values. Then

    Σ w⁻_i max(0, l_(i+1) − u_i)  ≤  Z_m D_m(L)  ≤  Σ w⁺_i min(u_(i+1) − l_i, K (θ⁺_(i+1) − θ⁻_i)),

where θ± are rigorous bounds for θ(t) = 2 arctan t.

*Proof.* On each cell the kernel lies between its endpoint values, because the mesh contains the only interior maximum, t = 1. The Stieltjes increment of L is nonnegative and bounded by the interval data and by the Lipschitz bound. ∎

The checker `window_enclosure` implements this with a rigorous enclosure of Z_m. It recovers three things on the reference models:
- the circle model L = a²θ, whose enclosure contains a² and tightens under refinement;
- the v319 plateau, whose upper bound falls below 10⁻⁵ by m = 16;
- a deliberately coarse three-point mesh, whose answer [0, >2a²] is uninformative, as it should be.

The roadmap's order of limits applies unchanged: |D_m(ℓ_n) − D_m(L)| ≤ 2‖ℓ_n − L‖_∞ / Z_m.

**Limitation.** No certified values of ℓ_n are available for any actual curve in this repository. Lemma 328.1 is therefore machinery, checked here only on models. v329 and v330 avoid pointwise spectral values altogether by bounding L′ directly.

## v329. Ellipses: the window average converges to the square radius

Let E be the ellipse with semiaxes a > b > 0. Put α = a² + b², β = a² − b², and

    q₁²(θ) = 2a²b² / (α − β cos θ),   q₂²(θ) = 2a²b² / (α + β cos θ).

**Lemma 329.1 (all inscribed rectangles).** For 0 < θ < π, the inscribed rectangles in E whose diagonals meet at angle θ are exactly two. Both are centered at the origin and symmetric about the axes, with squared half-diagonals q₁²(θ) and q₂²(θ). At θ = π/2 both give the square, with q² = 2a²b²/(a² + b²).

*Proof.* Write the vertices as c ± u, c ± v with |u| = |v| = q, and H = diag(a⁻², b⁻²). Subtracting the conic equations at opposite vertices gives cᵀHu = cᵀHv = 0. Because u and v are independent, Hc = 0, so c = 0. Then uᵀHu = vᵀHv. For directions φ and φ + θ this reads (b⁻² − a⁻²)(sin²(φ+θ) − sin²φ) = 0. Since a ≠ b, sin²(φ+θ) = sin²φ, which forces 2φ + θ ≡ 0 (mod π). The two solutions φ = −θ/2 and φ = π/2 − θ/2 give the two rectangles. Substituting cos²(θ/2) = (1 + cos θ)/2 gives the two radii. ∎

**Theorem 329.2.** Let γ_n be registered, area-convergent generic smooth embeddings converging to E in C¹. Let L be any uniform subsequential limit of their upper spectral functions. Then for almost every θ ∈ (0, π),

    min(q₁², q₂²)(θ) ≤ L′(θ) ≤ max(q₁², q₂²)(θ).

Consequently, for every m,

    Z_m⁻¹ ∫ w_m min(q₁², q₂²)  ≤  D_m(L)  ≤  Z_m⁻¹ ∫ w_m max(q₁², q₂²),

and both bounds converge to 2a²b²/(a² + b²) as m → ∞. So lim D_m(L) exists and equals the squared radius of the inscribed square, whichever branch carries the spectral value.

*Proof.* The upper bound is v318 with M_E² = max(q₁², q₂²) from Lemma 329.1.

For the lower bound, fix a band [δ, π − δ]. C¹ convergence with nonvanishing derivative gives the uniform secant-cone condition of [v277](v277-localized-spectral-slope.md) with an arbitrarily small cone angle. So for large n every inscribed θ-rectangle on γ_n with θ in the band has radius at least a fixed c_δ > 0.

Next suppose that for some ε > 0 there were rectangles on γ_(n_j) at angles θ_j in the band with q² < min(q₁², q₂²)(θ_j) − ε. A subsequence converges to a nondegenerate θ*-rectangle on E. Its q² is q₁²(θ*) or q₂²(θ*), which is at least the minimum, a contradiction. So, uniformly on the band and for large n, every rectangle has q² ≥ min − ε.

By Greene–Lobb's Theorem 5.1, ℓ_n′(θ) is the q² of such a rectangle for almost every θ. Integrating gives L(y) − L(x) ≥ ∫ (min − ε) on every subinterval of the band. Let ε → 0, then δ → 0.

For the limit in m: both q_i² are continuous, they agree at π/2, and the normalized kernel concentrates there. ∎

The circle (a = b) is immediate: every inscribed rectangle has q = a, so D_m(L) = a² for all m. The theorem calibrates both sides of the roadmap's test without choosing a branch. The single-branch integrals ∫ w_m q_i² are equal for i = 1, 2 by the symmetry θ ↦ π − θ, and they reduce to exact rationals. The branch-free bounds involve arctan(b/a) and are enclosed rigorously.

| Semiaxes | Square q² | m = 0 | m = 4 | m = 64 | m = 256 |
| --- | --- | --- | --- | --- | --- |
| 2 : 1 | 8/5 = 1.6 | [1.181, 2.819] | [1.397, 1.930] | [1.537, 1.672] | [1.567, 1.635] |
| 8 : 1 | 128/65 ≈ 1.969 | [1.267, 14.733] | [1.604, 2.845] | [1.848, 2.120] | [1.905, 2.040] |
| 32 : 1 | 2048/1025 ≈ 1.998 | [1.273, 62.727] | [1.620, 2.944] | [1.872, 2.156] | [1.932, 2.072] |

Entries are the rigorous [lower, upper] bounds of Theorem 329.2 with b = 1. Multiplying both semiaxes by 3 multiplies everything by 9, and the replay checks that exactly. The single-branch value at m = 0 is ab = A/π, the mean slope.

Thin ellipses show the roadmap's point concretely. A/π grows like the aspect ratio, while q²/r² = 2a²/(a² + b²) stays in [1, 2). The window needs larger m to localize as the ellipse thins: at aspect ratio 32 the upper bound is still 1.47 times q² at m = 4.

## v330. Small rectangles cannot carry the spectral value on the original spirals

**Lemma 330.1 (oscillation exclusion).** Let γ be a generic smooth Jordan curve with enclosed area A and radius R. Put G = min(A/2, A − πR²/2). Suppose s > 0 and ω ≥ 0 satisfy the following: the four contacts of every inscribed rectangle with half-diagonal at most s lie in some arc σ of γ on which, for some origin o_σ, the primitive of ((x−o_x)dy − (y−o_y)dx)/2 has oscillation at most ω. If 2ω < G − R²h for some h > 0, then ℓ_γ′ > s² almost everywhere on |θ − π/2| ≤ h.

*Proof.*
1. By (323.1), the action residue of such a rectangle is ΔF₁ + ΔF₂ modulo A. Here ΔF_i are differences of the primitive along the two return arcs inside σ, and the sum does not depend on the origin. So its distance from AZ is at most 2ω.
2. By subadditivity, the endpoint value and the Lipschitz cap, A/2 ≤ ℓ(π/2) ≤ πR²/2. On the window, ℓ moves by at most R²h. Hence η := G − R²h ≤ ℓ(θ) ≤ A − η there.
3. Since 2ω < η, the hypothesis (325.1) of v325 holds, and v325 gives ℓ′ > s² almost everywhere on the window. ∎

This is the quantitative form of v326's regularity criterion. That criterion asked for small-scale oscillation tending to zero. Lemma 330.1 needs only oscillation below half the spectral gap.

**Applying it to the v324 square.** Let P_k be ∂[0,4]² with one v324 gadget of N turns at scale ε = 1/k (N = k² for the original family). Every rectangle with q < 1 has its contacts within distance 2 of each other. Contacts that close lie either on one edge, which allows no rectangle, or in the arc around a single corner. At the three plain corners the primitive about that corner is constant on both edges, so ω = 0 there. At the corner (0,0), the primitive about (0,0) is constant on the two edges, so ω is the oscillation Osc_k across the gadget.

The checker computes Osc_k, A_k and R_k² exactly. It also proves the following bounds for every k:
- Osc_k ≤ 19/24 + (19/4)ε². With the origin at the gadget center, the primitive rises monotonically along the inbound spiral to ε²(19N/24 − 1/(384N)) ≤ 19/24. It never falls below −5ε²/2 afterward, because every outbound term is dominated by the corresponding inbound term. Moving the origin to (0,0) costs at most √5 ε².
- A_k ≥ 16.
- R_k² ≤ 8 + 6ε + 5ε²/4.

The resulting margin, G − 2·Osc, is positive at k = 7 and decreasing in ε. So Lemma 330.1 applies to sufficiently fine generic smoothings of every P_k with k ≥ 7. The v324 marked smoothing changes the primitive, the area and the diameter by arbitrarily little. k = 6 also passes on its exact values.

| Family | k | Turns N | 2·Osc_k | G_k (lower bound) | Window half-width h < |
| --- | ---: | ---: | ---: | ---: | ---: |
| original | 4 | 16 | 1.970 | 1.133 | fails |
| original | 6 | 36 | 1.755 | 1.887 | 0.0146 |
| original | 8 | 64 | 1.680 | 2.269 | 0.0672 |
| original | 12 | 144 | 1.626 | 2.654 | 0.1208 |
| original | ∞ | k² | 19/12 ≈ 1.583 | 16 − 4π ≈ 3.434 | ≈ 0.2313 |
| retuned | 8 | 647 | 16.10 | 2.269 | fails |

**Conclusion for the original family.** For every s < 1 and every k ≥ 7, every small rectangle is excluded as a carrier on |θ − π/2| < h_k. So ℓ′ ≥ 1 almost everywhere on that window. Any subsequential limit L along k → ∞ therefore has L′ ≥ 1 almost everywhere on |θ − π/2| < 0.2313, and in particular D_m(L) ≥ 1 − o(1). This activates v325 on the family without any Floer computation, which is the outcome the roadmap's v330 test asked for. The marked shrinking square of v324, with action gap 19/24, is one of the excluded rectangles.

**The retuned family fails.** With N_k the nearest integer to (192/19)k², the tagged return integral tends to 8 = A/2, and 2·Osc_k tends to 16. The bound gives nothing. Carrier selection on this family is exactly the open question of v331.

**A fixed nonrectifiable target.** Put gadgets with ε_j = 2^(−j−3), N_j = ε_j⁻² at heights 2^(−j) on the left side of [0,4]². Every truncation satisfies Osc ≤ 0.9014, by bounding each gadget about its attachment point and adding the area drift of the earlier gadgets, each at most 4ε_j². It also has A ≥ 16 and R² = 8. So 2·Osc ≤ 1.803 < 16 − 4π, and every truncation, hence every approximation built from them, has ℓ′ ≥ 1 on |θ − π/2| < 0.2038. The exact oscillation for one to three gadgets is 0.8257, set by the largest gadget.

This target has infinite length, so the length-based shrinkout argument of Greene–Lobb (their Lemma 5.1) does not apply to it. Its area primitive also oscillates by about 19/24 on arbitrarily short arcs, so a vanishing-oscillation hypothesis does not hold. We have not checked this example against the precise definitions of Asano–Ike ([arXiv:2609.11598](https://arxiv.org/abs/2609.11598)).

## What this does and does not settle

- v329 shows that the window functional measures the right quantity on ellipses, with explicit rates and no branch selection.
- v330 gives the first families where the v325 hypothesis is verified, including a nonrectifiable curve. It also identifies the precise failure mode: the oscillation across a single gadget exceeds half the spectral gap, as in the retuned family.
- Neither part gives a bound for general curves. Lemma 330.1 needs both an interior spectral gap at π/2 (here supplied by A > πR²/2) and bounded small-scale oscillation of the area primitive. Arbitrary Jordan curves satisfy neither.

## Replay

    python code/angular_windows_v328_v330.py
    python -m unittest discover -s tests -p 'test_angular_windows_v328_v330.py' -v

[Recorded output](../results/v328-v330/windows-ellipses-spirals.json). The single-branch ellipse values were cross-checked against independent floating-point quadrature at m = 0, 1, 4, 16 and 64 for five ellipses.
