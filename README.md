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
- [Exact-rational periodic-graph model](code/periodic_bridge_v212.py) demonstrates that balanced periodic extensions can have a nondegenerate square while **none** has all four vertices in the specified original-piece windows. [Seven tests](tests/test_periodic_bridge_v212.py) and GitHub Actions cover the construction.
- **Open bridge:** To derive the Toeplitz conjecture, one still must construct four area-balanced periodic strands from an *arbitrary* Jordan curve, preserve non-collapse, and force a square wholly on its original pieces. The October 2026 Floer theorems do not supply this localization automatically.
