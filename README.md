# Toeplitz square-peg research

**Research program:** Benjamin John Schulz (GitHub: HiroSakuraba)  
**Status as of 2026-11-18:** Research in progress. The general square-peg conjecture for arbitrary Jordan curves is **not proved** here.

The [Toeplitz square-peg problem](https://en.wikipedia.org/wiki/Inscribed_square_problem) asks whether every Jordan curve contains four distinct vertices of a nondegenerate Euclidean square. This is a versioned research archive of constructive half-square approaches, local certificates, geometric obstructions, analytic energy methods, and supporting computation.

## Repository contents

- **[archive/](archive/)**: preserved original ZIP research bundles, versions **v202–v211**, with original internal source code, notes, experiments, and selected verifiers.
- **[archive/early/](archive/early/)**: independently developed local rough-Toeplitz recognizer and preconditioned certificate bundles.
- **[related/](related/)**: *separately scoped* exploratory Toeplitz Gram-matrix/SOS work in Extremal Geometry Lab. This is not a result about Jordan-curve square inscription.
- **[research/2026-11-literature-and-v212-plan.md](research/2026-11-literature-and-v212-plan.md)**: current literature review and research plan, explicitly distinguishing established results from proposed lemmas.
- **[archive/SHA256SUMS.txt](archive/SHA256SUMS.txt)**: integrity hashes of the archived ZIP files.

## Reading order

1. **v211** — shape sphere, quantitative tangent escape, and the remaining non-collapse obstruction.
2. **v210** — counterexample to a proposed extraction of strict v209 four-arc certificates.
3. **v209** — valid local four-arc square-forcing theorem under directional compatibility.
4. **v208–v207** — controlled finite-coil models.
5. **v206** — proof-route critique and Liouville/Abel-curve audit.
6. **v202–v205** — earlier energy, boundary response, stress localization, and fjord/mouth analysis.
7. **November 2026 literature plan** — how recent Floer-homology and Liouville-regularity results may alter the next attack.

Each research bundle preserves its own assertions, proof audit, and checker limitations. A passing numerical or exact-rational checker establishes only the statements it actually tests. Analytic and topological arguments are **not** automatically proof-assistant-certified.

## Scientific caveat

None of these manuscripts, finite checks, examples, or literature findings constitutes a proof of the unrestricted Toeplitz conjecture. Not every proposed intermediate implication survived scrutiny: in particular, v210 refutes a universal strict crossing-certificate extraction hypothesis while leaving the valid v209 sufficient criterion intact.

## Reproduction

Download the relevant ZIP, extract it, read its `README.md` and `notes/`, and run its verifier exactly as documented. Toolchains vary across versions. Where methods use Arb or external dependencies, verify those tools are installed before claiming certified output.

## External research

See the [November literature review](research/2026-11-literature-and-v212-plan.md) for primary sources including October 2026 proofs of periodic area-balanced inscriptions and the November update of Asano–Ike's rectifiable rectangular-peg work.

## v212 research batch (18 November 2026)

- [Periodic-bridge theorem and obstruction note](research/v212-periodic-bridge-theorems.md) proves a sharp collective-intersection-gap lower bound on square size (with a Hausdorff stability estimate).
- [Exact-rational periodic-graph model](code/periodic_bridge_v212.py) demonstrates that balanced periodic extensions can have a nondegenerate square while **none** has all four vertices in the specified original-piece windows. [Eight tests](tests/test_periodic_bridge_v212.py) and GitHub Actions cover the construction.
- **Open bridge:** To derive the Toeplitz conjecture, one still must construct four area-balanced periodic strands from an *arbitrary* Jordan curve, preserve non-collapse, and force a square wholly on its original pieces. The October 2026 Floer theorems do not supply this localization automatically.

## v213: exact PL square-localization and lift-coherence audit

- [Mathematical note and proof obligations](research/v213-exact-PL-localization-and-lift-audit.md): complete finite decision for labeled squares on rational piecewise-linear periodic graphs, including both square orientations, all bounded integer deck translates, and rank-deficient families.
- [Exact-rational solver](code/periodic_square_cert_v213.py) distinguishes **periodic quotient witnesses** from squares using the **prescribed physical lifts** of four original planar arcs. The latter distinction adds a second transfer obstruction beyond v212's connector-only witnesses.
- A strict local exclusion-margin certificate is stable under small coordinate perturbations; a separate exact transverse-square certificate yields a quantified persistence threshold via Brouwer's theorem.
- [New regression suite](tests/test_periodic_square_cert_v213.py) joins the existing tests under CI. **Neither these finite PL results nor the periodic Floer theorems prove the general Toeplitz square-peg conjecture.**

## v214: source/lift-relative signed degree and controlled continuation

- [Conditional Brouwer-degree theorem and proof](research/v214-relative-degree-source-lifts.md): exact signed count of oriented square configurations on four designated physical source arcs with a fixed deck-translation tuple. A **nonzero degree** forces a nondegenerate square *on those arcs* under the stated positive vertical-separation guard.
- [Exact rational degree evaluator](code/relative_degree_v214.py): sums local signs of isolated PL squares, rejects boundary or singular roots explicitly, and returns a sufficient perturbation threshold for **arbitrary continuous vertical deformations** preserving the marked windows and physical lifts.
- [Regression tests](tests/test_relative_degree_v214.py) cover lift mismatch, a root exiting the source boundary, two transverse squares with cancelling degrees, and 48 perturbations independently cross-checked against the v213 exhaustive solver.
- Example: signed degree **−1**, one square of side-squared **1/32**, certified perturbation threshold **ε < 55/8736**. **Degree zero does not imply absence of squares.**
- **Outstanding:** Neither the periodic Floer theorems nor this degree calculator provide a nonzero relative degree for four source arcs extracted from an arbitrary Jordan curve. This is *not* a relative Floer-homology construction or a proof of the unrestricted Toeplitz conjecture.

## v215: knot-robust boundary-certified source/lift degree

- [Mathematical proof and falsifiers](research/v215-knot-robust-relative-degree.md): exact Brouwer degree for fixed physical lifts of four rational piecewise-linear graph arcs, even if true squares lie on **internal PL knots** or isolated-chart methods fail on singular affine cells.
- [Exact degree evaluator](code/relative_degree_v215.py) checks for zeros on the **external source-arc boundary**, obtains an explicit rational boundary residual margin, and computes degree via a nearby rational **regular target**. Regular-target preimages are *not* actual squares; their signed sum computes degree.
- [16 new tests](tests/test_relative_degree_v215.py): knot/refinement invariance, singular empty cells, positive and negative degree cases, source/lift mismatch, degree cancellation, boundary escape at exactly `t=55/624`, and randomized continuation checks. The full test suite is also run by CI.
- For the original planted square: degree **−1**, certified external-boundary margin **1/384**, and guaranteed persistence under continuous vertical arc perturbations **ε < 1/768**. This is more generally applicable but more conservative than the earlier single-chart v214 threshold.
- **Outstanding:** There is still no theorem producing a nonzero relative degree from every Jordan curve or transporting the October periodic Floer existence results through connector and deck-lift localization. The universal Toeplitz square-peg conjecture remains open.

## v216: boundary-only degree and a priori polygonal arc selection

- [Proof and research audit](research/v216-boundary-degree-global-arc-selection.md): computes the relative Brouwer degree from **outward-oriented boundary-ray intersections**, avoiding interior square enumeration, for four arbitrary rational **planar polygonal** source arcs (not necessarily graph arcs).
- [Exact-rational implementation](code/boundary_degree_v216.py): validates simple input polygons, chooses vertex-half-edge and full-edge windows **without prelocating a square**, computes physical-lift Euclidean separation, and returns either a nonzero certificate or **inconclusive**, never a claim of universal square absence.
- [25 new tests](tests/test_boundary_degree_v216.py): 72 independently computed affine-root sign crosschecks, 15 two-dimensional convex-quadrilateral perturbations, nonmonotone arcs, source-boundary failures, orientation reversal, deck isolation, internal knots, degree cancellation, and a **known-square/no-certificate** control.
- **New conditional convex example:** for polygon `((0,0),(3,1),(2,3),(-1,2))`, preselected vertex arcs yield degree **−1**, external boundary margin **1/8**, opposite separation squared **13/4**, and an **arbitrary continuous planar perturbation** guarantee `ε<1/32` in coordinatewise sup norm. Another quadrilateral yields a full-edge certificate with `ε<1/64`.
- The original unit square itself defeats these specific finite windows because some squares hit their boundaries. **There is still no global source-arc selection theorem for arbitrary Jordan curves; Toeplitz remains open.**
