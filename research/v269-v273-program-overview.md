# v269–v273: Unified relative-rectangle research ledger

**Program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Five completed restricted research investigations, exact verifier and 31 tests. The global Toeplitz conjecture remains open.

| Version | Verified result | Remaining universal obligation |
|---|---|---|
| [v269](v269-relative-index-ellipse-stability.md) | Relative opposite-endpoint sign on ellipses; exact nonzero 4-parameter square Jacobian and local C0 persistence near noncircular ellipses | Construct a natural global relative line/degree class for the compactified rectangle locus on arbitrary Jordan approximations |
| [v270](v270-spectral-axiom-no-go.md) | Analytic strict monotone Lipschitz subadditive action countermodel with derivative 0 at pi/2 | Identify additional *Floer/geometric* data that protects a class at exactly pi/2 |
| [v271](v271-six-sector-equivariant-ledger.md) | All 24 orderings, six C4 sectors, three D4 classes, all three opposite-pair matchings, universal exact square-residual equivariance, and a complete two-branch rectangle locus in one four-edge cell of a nonconvex Jordan polygon | Build a real full-curve rectangle branch atlas, with pairings, vertex gluing, singularities and boundary strata |
| [v272](v272-degeneration-trichotomy-and-shrinkout.md) | Complete geometric rectangle-limit trichotomy (genuine, two-pair flat, or total collapse); explicit genuine Jordan polygons with a selected square shrinking to zero | Prove noescape of a *protected global class*, not of a selected square |
| [v273](v273-ellipse-angle-continuation.md) | Unique rectangle for every diagonal angle on noncircular ellipse, sharp positive radius/side floors in any compact band around pi/2; conditional size-protected angle compactness theorem | Establish a uniform protected angle sequence approaching pi/2 on every Jordan curve |

## Integration and controls

All five versions share the [pure-standard-library exact verifier](../code/relative_rectangle_v269_v273.py), [31 independent tests](../tests/test_relative_rectangle_v269_v273.py) and a [dedicated workflow](../ci/replays.json).

The code uses rational arithmetic to verify the restricted ellipse and polygon examples, source/target group symmetries and geometric degenerations. The spectral countermodel's universal claims are analytic trigonometric identities in v270; any floating-point sampling is regression only, **not proof**.

**Recommended next actual theorem target:** a relative class defined on the true compactified 1D binormal/rectangle locus that survives generic birth-death events and has a geometric or Floer filtration preventing escape via total four-corner collapse. Exact sector bookkeeping and protected angle data are necessary checks but do not create such a class by themselves.

**Nonclaims:** No universal nonzero relative index, no full all-sector geometric branch census, no filtered Floer continuation theorem, no universal angle-neighborhood theorem, and no resolution of the unrestricted Toeplitz conjecture.

Primary literature: Greene–Lobb, [arXiv:2604.17116 (2026)](https://arxiv.org/abs/2604.17116), especially §2 and §3. The source's genericity result concerns a compactified locus of rectangles **and binormals** and does not automatically imply an odd loop of nondegenerate rectangles. The ellipse counterexample to universal odd-loop existence is in [v268](v268-rectangle-monodromy-and-protected-square-criterion.md).
