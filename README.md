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

## v217: exact multiscale signed-degree ledger

- [Mathematical proof and examples](research/v217-multiscale-relative-degree-ledger.md) establish a **telescoping boundary-chain identity** for subdivisions of the four-arc square-configuration domain. The oriented contributions of every shared internal wall cancel exactly.
- [Standalone exact-rational implementation](code/multiscale_degree_v217.py) computes parent degree, signed child face fluxes, physical opposite-arc separation, clean child degree certificates, and internal-wall pair cancellation. A child whose boundary contains a square **does not have a defined child degree**; only its formal ray flux is reported.
- [Regression suite](tests/test_multiscale_degree_v217.py) checks nonzero parent degree, 16-cell grids, rational perturbations, physical lift restrictions, a square precisely on an internal cut, and **two opposite-sign squares** which split into degrees `−1` and `+1` despite parent degree zero.
- **Limits:** This is a finite-dimensional Brouwer-degree calculation, not Floer homology or a global existence proof. Refinement conserves an existing total degree; it does not produce a nonzero total for an arbitrary Jordan curve.

## v218: cyclic orientation cancellation and collision-stratum audit

- [Mathematical note](research/v218-global-configuration-symmetry-and-collision-audit.md) analyzes the full labeled cyclic four-point configuration space. Cyclic relabeling has source Jacobian sign `−1`, whereas the induced linear map on the four square residuals has determinant `+1`. Consequently, a regular square's four cyclic labelings have alternating local signs and **cancel in the naive global integer count**.
- [Exact-rational symmetry and polygon code](code/configuration_symmetry_v218.py) checks both orientations, fourfold equivariance, all **six** curve-order sectors (only two graceful), and collision boundary degeneracies. A generic convex quadrilateral yields one geometric square with eight labeled descriptions in the finite distinct-open-edge model.
- [18 regression tests](tests/test_configuration_symmetry_v218.py), including 200 exact rational equivariance cases and an independent four-rotation comparison to the v216 degree calculator.
- **Limitation:** The C4 quotient is nonorientable, and the collapsed all-four-point diagonal remains on the naive compactification boundary. The code **does not construct** an all-Jordan mod-2 or equivariant invariant, or prove Toeplitz.

## v219: equivariant quotient and collision-link obstruction

- [Mathematical report](research/v219-equivariant-quotient-collision-link.md) proves that one cyclic-order component modulo the free cyclic `C4` action deformation-retracts onto a circle. For both square orientations the associated rank-four target bundle has an **explicit nowhere-zero equivariant section**, so the *absolute* top Euler/`w4` obstruction vanishes. This does **not** eliminate possible relative invariants or known generic smooth-curve mod-2 orbit parity.
- For fixed regular `C1` curves, the normalized total-collision square residual admits a **uniform nonzero first-order link** controlled by tangent magnitude and modulus of continuity. This statement fails to provide uniform control over arbitrary `C0` approximations.
- [Exact-rational code](code/collision_quotient_v219.py) and [21 tests](tests/test_collision_quotient_v219.py) verify the group action, retraction, target equivariance, collision tangent formula, and a simple outward-notch polygon family with inscribed squares of side `1/n` converging to zero while its boundaries converge to a fixed rectangle.
- **Next missing result:** construct meaningful relative boundary data (potentially in a cyclohedral compactification) that survive *arbitrary rough Jordan curves*. Neither an absolute parity count nor the constructed shrinking-square sequence solves the Toeplitz conjecture.

## v220: rectifiable fixed-curve collision-link obstruction

- [Proof and research audit](research/v220-fixed-rectifiable-collision-link-obstruction.md) construct **one fixed rectifiable Jordan curve** of perimeter `37/6`, formed from a rectangle with countably many disjoint square-shaped notches of size `s_n=4^(-n-1)`. Every notch is an exact gracefully inscribed square, so **even a single rectifiable curve need not have a uniform positive square-size floor**.
- Two four-point sequences on this same curve have **identical normalized successive parameter gaps** `(1/3,1/3,1/3)` and approach the same parameter collision chart point, but the residual divided by their common arclength parameter span has incompatible limits: **0** on the notch squares and **(37/18,−37/18,37/6,−37/18)** on nearby collinear controls. Thus that normalized residual **does not extend continuously on the parameter-only collision link**; retaining ambient shapes still leaves actual boundary zeros.
- [Exact rational constructor/checker](code/accumulating_notches_v220.py) and [13 regressions](tests/test_accumulating_notches_v220.py) check finite polygon simplicity, explicit square size, perimeter/area series, arclength gap ratios and the conflicting boundary limits. CI additionally runs the independent v216 polygon validator.
- **Scope:** This disproves a proposed boundary-regularity assumption, **not Toeplitz**. The constructed curve has infinitely many squares. The next approach must address collision-boundary zeros via a true relative cycle or macroscopic size selection.

## v221: positive-size-filtered incidence and exact wall crossing

- [Conditional mod-2 cobordism theorem and explicit examples](research/v221-positive-size-cutoff-incidence.md): for a compact transverse one-dimensional zero locus, the parity difference of **positive-size squares on four fixed source arcs and physical lifts** equals the parity of *source-boundary* plus *size-cutoff* crossings. This is standard topological bookkeeping, **not** a global all-Jordan parity invariant.
- [Exact rational single-branch solver](code/cutoff_incidence_v221.py): four translated PL segments yield one affine square branch and a quadratic squared-side function. All rational source-exit times and all *simple* cutoff crossings are certified using rational signs, including irrational root times. Simultaneous events, tangencies, and singular Jacobians are rejected as nongeneric.
- [21 tests](tests/test_cutoff_incidence_v221.py) cover one and two cutoff crossings, a source-wall exit at `55/19968`, a radical-time crossing, time reversal, rigid transformations, and v220-tagged corner-square cutoffs. The two standalone-skipped integration cases run on full GitHub CI.
- **Limits:** The v220 notch count includes only *explicitly tagged notch squares*, not all square orbits; its alternating cutoff parity does not refute a possible global mod-2 invariant. A uniform positive scale or controlled relative boundary class remains an open requirement for unrestricted Toeplitz.

## v222: exact multibranch PL knot and size-cutoff incidence

- [Conditional theorem and detailed proof](research/v222-multibranch-cutoff-and-orbit-audit.md): the mod-2 endpoint count on **four specified translating polygonal source arcs** changes only through exposed source-wall and positive-size-cutoff crossings. At an ordinary **internal PL knot**, two chart half-branches pair and cancel, even when two squares meet and disappear at a corner. The general Toeplitz conjecture remains **open**.
- [Rational multi-chart solver](code/multibranch_incidence_v222.py) enumerates regular square branches across all four-arc segment products, matches internal endpoints by **original arc parameters** rather than geometry alone, counts quadratic cutoff roots without irrational approximation, checks singular-empty charts by exact five-variable affine feasibility, and refuses unresolved nongeneric cases.
- [25 regression tests](tests/test_multibranch_incidence_v222.py): 23 passed locally, two optional crosschecks with v218/v220 skipped outside complete checkout. Tests include the two-to-zero tent example, a through-knot square, mixed cutoff/knot, source exit, wrong physical deck shifts, reflection, and 36 seeded independent comparisons against v221.
- A separate **D4 orbit-key utility** maps all eight ordered descriptions of one square to a single geometric vertex set. This is **not** a global orbit-level parity count. An abstract smooth fold normal form is treated separately, not claimed as an implementation of smooth square-residual bifurcations.
- **Next:** an actual geometric-square orbit cobordism for a homotopy of **one embedded Jordan polygon**, with a fixed positive side cutoff and all sectors/chart endpoints accounted for. No positive scale is guaranteed under arbitrary C0 approximation.

## v223: complete polygonal geometric-square census and a one-polygon cutoff model

- [Exact rational decision theorem and proof](research/v223-complete-polygon-orbit-census-and-homotopy-audit.md): exhaustively enumerate **all square orbits on one rational simple Jordan polygon**, including four points on repeated edges, square vertices on polygon corners, both orientations, and all six curve-order sectors. A singular affine chart is solved via rational polytope vertices; if positive-size square solutions form a **continuum**, the program gives a genuine rational witness and explicitly **refuses** a spurious finite parity count.
- [Implementation](code/polygon_orbit_census_v223.py) and [25 tests](tests/test_polygon_orbit_census_v223.py): exact geometric-orbit deduplication of all eight D4 labels and polygon-edge endpoint charts. Generic four-vertex polygon: exactly one square `q=72/25`. A trapezoid: two squares `q=45/49` and `q=1`. An ordinary square and a 2×1 rectangle both have **positive-dimensional square families**.
- A **genuine single-polygon homotopy** by positive uniform scaling and translation gives an exact **geometric-orbit size-cutoff event ledger**, including irrational crossing times represented by rational isolating intervals. It counts all geometric square orbits of the input polygon, not marked source-arc labels.
- A strictly convex affine-vertex deformation exhibits exact snapshot counts `1 → 2 → 3` as a square passes through a nongeneric polygon-vertex alignment. Snapshot counts are **not** claimed to be a complete homotopy incidence computation.
- **Outstanding:** an exact algebraic event ledger for arbitrary moving polygon vertices (determinant singularities, source-vertex crossings, fold events), then a relative global class with a nonzero size bound in arbitrary `C0` Jordan limits. Neither is established here; Toeplitz remains open.

## v224: exact algebraic event candidates for linearly moving Jordan polygons

- [Mathematical derivation and scope](research/v224-algebraic-moving-polygon-event-covers.md): on **one ordered four-edge chart** of a polygon with independently moving rational vertices, the square parameters are Cramer ratios `N_i(t)/D(t)` of rational polynomials of degree ≤4. Potential source-edge crossings occur at `N_i=0` or `D−N_i=0`; the positive squared-side cutoff has one polynomial of degree ≤10. All change times of a **nonpersistent nonsingular chart** lie in these candidate zero sets.
- [Exact standard-library checker](code/moving_polygon_events_v224.py) uses **rational Sturm sequences** to isolate every real candidate root, including irrational and repeated roots, with explicit resource limits; it does not mistake spurious candidates for actual events. A separate sufficient certificate rules out **every triple-vertex collinearity** over the deformation to guarantee the polygon stays Jordan throughout, not just at its endpoints.
- [23 exact regression tests](tests/test_moving_polygon_events_v224.py) passed locally, cross-checking moving Cramer solutions against v223's independent static solver. A complete convex-quadrilateral vertex-alignment candidate occurs at normalized time `t=8/9` with exact square `(1,0),(2,0),(2,1),(1,1)`; a scaling example isolates cutoff time `√2−1` without floating point. The new tests also confront repeated-edge charts, permanent endpoint identities, and paths that cease to be simple polygons.
- **Critical limit:** A 408-chart four-edge scan of the convex model has **128 persistently singular/endpoint-polynomial charts**. These are explicitly marked **unresolved**; v224 gives an exact *candidate cover*, **not** a global geometric-square orbit cobordism, cutoff-independent parity invariant, or unrestricted Toeplitz proof.

## v225: canonical polygon-edge ownership and exact singular-chart audit

- [Mathematical proof, quantitative audit, and limitations](research/v225-persistent-edge-strata-rank-certificates.md): in v224's embedded convex quadrilateral deformation, **all 128 previously flagged persistent/singular chart identities are classified**. The total **408** product charts now split into **280 regular, 40 vertex-owned, 40 incoming-edge duplicates, 48 geometrically impossible adjacent-edge double pairs, and zero genuinely unresolved**. This is chart classification, **not** yet a complete geometric-orbit cobordism or a Toeplitz proof.
- [Exact-rational classifier](code/persistent_strata_v225.py): one half-open outgoing-edge convention assigns a polygon vertex uniquely to its outgoing edge, retaining real corner squares while eliminating duplicate incoming-endpoint labels. A geometric lemma excludes squares with two vertices on each of two incident polygon edges. Identically singular charts are tested by the **polynomial rank of [A|b] versus A over Q(t)**; a nonzero augmented minor certifies that squares are impossible except possibly at its exactly isolated algebraic zeros. Generic consistently singular families remain explicitly **unresolved**.
- [27 new exact regressions](tests/test_persistent_strata_v225.py) (50 v224+v225 tests passed locally): independent v223 orbit-census crosschecks at five rational times (including the true square on a polygon vertex), two square orbits on a trapezoid, rational and irrational singular-support controls, moving Jordan checks, and a rectangle **positive-dimensional square family** that is correctly *not* assigned finite parity.
- **Still missing:** exact clustering of the **224 raw candidate entries** into real geometric orbit events, local analysis of exceptional determinant roots, and a protected nonzero-size global class under C0 Jordan limits. These limitations prevent any claim of a complete general-polygon cobordism or the unrestricted Toeplitz theorem.
