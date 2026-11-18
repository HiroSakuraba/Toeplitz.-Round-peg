# Literature review and v212 research plan (18 November 2026)

**Author/research program:** Benjamin John Schulz  
**Status:** Research note, not a proof. Primary papers are linked; several are arXiv preprints, not peer-reviewed final publications.

## The present obstacle

The v211 bundle separates the *shape* of a four-point configuration from its physical scale. Normalizing the center and size does not rule out a sequence of approximate squares with side length tending to zero. Its quantitative tangent estimate forces tangent-scale escaping, not a universal non-collapse theorem. v210 independently demonstrates that large signed area, nested spiral lobes, and small gaps do **not** suffice to extract a strict v209 four-arc crossing certificate.

Any successful extension must therefore either (A) force a square at a definite scale or (B) use an invariant that survives rough limits and explicitly exclude collapse.

## Directly relevant new work

### 1. Lorenzo Valerio, *A Proof of the Area Vanishing Square Peg Conjecture* (15 October 2026)

Primary source: https://arxiv.org/abs/2610.21076

For four essential polygonal simple closed curves on the cylinder, vanishing of the alternating Liouville-area sum
`∫σ₁ y dx − ∫σ₂ y dx + ∫σ₃ y dx − ∫σ₄ y dx = 0`
forces a joint periodic square inscription. Lagrangian Floer homology provides the existence mechanism. The periodic formulation explicitly allows degenerate squares, and it does **not** solve the arbitrary planar Jordan-curve conjecture.

**Possible use:** replace a proposed four-local-crossing extraction by a *four-strand area-balance* construction, but still prove admissible periodicization, localization of the returned witness, and a nonzero size bound.

### 2. Ali Naseri Sadr, *The Four Periodic Inscription Problem* (17 October 2026)

Primary source: https://arxiv.org/abs/2610.24376

A generalized balanced-area identity `t A(γ₁,γ₂)+(1−t)A(γ₃,γ₄)=0`, for `0<t<1`, gives prescribed isosceles-trapezoid inscriptions for four periodic embeddings under the stated nonintersection condition. In particular the square result follows in the balanced case. A continuous-curve extension is proved for zero-planar-measure periodic images under the prescribed hypotheses.

**Possible use:** Sadr's explicit exclusion of fourfold intersection is exactly the sort of noncollapse condition missing in an uncontrolled limit. Derive quantitative bounds when the four sets retain a positive common-point separation; determine whether our construction can ensure that separation uniformly as scales shrink.

### 3. Asano–Ike, *Liouville regular Jordan curves* (10 September 2026)

Primary source: https://arxiv.org/abs/2609.11598

Liouville regularity covers finite-p-variation, appropriate piecewise variation conditions, and conformal maps whose radial-tail image lengths uniformly vanish. Its consequences include rectangular inscriptions for quasicircles, suitable Hölder/John domains, and other controlled Jordan curves. Crucially, the authors construct compactified spirals and even arcs differentiable at a point whose action primitive fails to extend continuously.

**Possible use:** classify test curves into a regular lane already covered by known theorems and an irregular lane requiring new ideas. Include the finite-coil non-example as an adversarial regression against naive area cancellation.

### 4. Asano–Ike, *The Rectifiable Rectangular Peg Problem* (revised 11 November 2026)

Primary source: https://arxiv.org/abs/2412.21057

Uses microlocal sheaf theory to prove the prescribed-rectangle result for Liouville-regular continuous Jordan curves (including rectifiable and locally monotone cases). Re-read the November v4 revision as the current technical reference before citing or extending earlier arguments.

### 5. Li–Pan, *Rectangular Pegs on Jordan Curves of Finite p-Variation* (4 September 2026)

Primary source: https://arxiv.org/abs/2609.04918

For finite `p<2` variation, Young integration and controlled smooth Jordan approximation give the needed Liouville primitives and thus every rectangle similarity class. Independent reconstruction of a regularity criterion, not a solution for arbitrary Jordan curves.

### 6. Greene–Lobb, *Jordan curves inscribe a positive measure of rectangles* (18 April 2026; earlier but important for v211)

Primary source: https://arxiv.org/abs/2604.17116

For a Jordan curve of enclosed area `A` and half-diameter `R`, the set of diagonal angles with inscribed rectangles has measure at least `A/R²`. In the smooth generic setting the derivative of a Jordan-Floer spectral invariant at almost every angle equals the *squared circumradius* of a corresponding rectangle.

**Possible use:** this invariant's derivative directly measures physical size, unlike normalized v211 shape data. But the proven positive-measure result does **not** force the exact square angle `π/2`.

### 7. Binder–Ivrii–Prause–Smirnov, *Brennan's Conjecture and the Universal Integral Means Spectrum* (16 October 2026)

Primary source: https://arxiv.org/abs/2610.23853

Following the September OpenAI preprints, this paper provides a shorter and stronger approach to conformal derivative integral means. It could inform estimates for the v203–v205 conformal/screening lane. **Low direct relevance:** global integral-means estimates do not by themselves imply the uniform radial-tail control, action regularity, or square-size floor needed here.

## A concrete v212 attack plan

1. **Build a theorem-assumption matrix** for the October periodic-inscription papers. Record orientation, cylindrical period, disjointness/non-common intersection, allowed degeneracy, area-balance normalization, and continuous-limit hypotheses exactly.
2. **Prove a simple quantitative anti-collapse lemma.** If four compact sets `Γ₁,…,Γ₄` have no common point, define `d=min_{p∈Γ₁} max_i dist(p,Γ_i)>0`. For any square with its labeled vertices `p_i∈Γ_i` and side `s`, we have `d≤√2·s`, hence `s≥d/√2`. Check periodic quotient/lift variants carefully. This is elementary and conditional; it is **not** a new proof of Toeplitz.
3. **Attempt area-balanced four-strand periodicization** of a bounded 2+2 or finite-coil packet. Specify actual continuous embeddings, signed-area cancellation, prevention of spurious inscriptions on artificial connectors, and preservation of a common-point separation independent of the regularization parameter.
4. **Test hostile families:** v210 rotating-scaling logarithmic packet, v211 tangent-escape graph, Asano–Ike's non-Liouville finite-coil example, and the previously retained v187 critical spiral case. These are *attempted falsifiers*, not cases to assume away.
5. **Compare the Floer-spectral size route:** seek estimates on near-square angles for the positive-size witnesses supplied by Greene–Lobb; explicitly identify the step not provided by their theorem at `θ=π/2`.
6. **Stop conditions:** If balanced periodicization forces artificial-connector witnesses, loses separation, or merely repackages the discredited v210 strict crossing extraction, mark the proposed lemma false/inconclusive and pivot. Do not promote numerical patterns or conditional statements to universal theorems.

**Working assessment:** October's area-vanishing proofs are the strongest new external inputs. The missing bridge is from *one arbitrary Jordan curve* to *four compatible, area-balanced periodic strands with a nondegenerate witness*. The spectral-size route is the second serious option.
