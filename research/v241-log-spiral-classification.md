# v241: Inscribed squares on open logarithmic spirals for a ≥ 3/50

**Research program:** Benjamin John Schulz  
**19 November 2026**  
**Status:** Computer-assisted theorem (Arb ball arithmetic) plus elementary lemmas. Extends the hand proofs of v239 (a ≥ 9/10) and v240 (a ≥ 17/20) downward and classifies the square pitches. Independent review welcome. Nothing here is a statement about arbitrary Jordan curves.

## Statement

Let S_a = {e^{(−a+i)t} : t ≥ 0} ∪ {0}, an open spiral arc compactified by its centre.

**Theorem.** For real a in [3/50, 9/10], S_a contains four points forming a nondegenerate square **if and only if** a = a\*, where

a\* ∈ [0.1002746560179786808441, 0.1002746562179786808441]  (a\* ≈ 0.10027465611797868).

At a = a\* the inscribed squares form exactly **one** orbit under the spiral's self-similarities z ↦ e^{(−a+i)c} z, c ≥ 0.

**Corollary (with v239).** v239 proves S_a is square-free for every a ≥ 9/10. So for all a ≥ 3/50, a\* is the only pitch at which S_a has an inscribed square; in particular **a\* is the largest such pitch**. The same holds for the full two-sided spiral, since a rotation-dilation moves any four of its points onto S_a.

v99 certified a square at a special pitch "near 0.1002747". By the uniqueness above, v99's square is this square.

## Proof

### N. Normalisation

Distance from the centre strictly decreases along the spiral, so a square on S_a has a unique farthest vertex, at some parameter t₀ ≥ 0. The similarity z ↦ e^{−(−a+i)t₀} z maps the spiral to itself, shifting parameters by −t₀. It sends the square to a square whose farthest vertex is v₀ = 1 = z(0) and whose other vertices have parameters > 0 (or are the centre), so it lies on S_a again. Conversely every such normalised square is a square on S_a, and it generates the orbit in the theorem. So it suffices to find all squares on S_a with v₀ = 1 as farthest vertex.

Label the square v₀v₁v₂v₃ in cyclic order, with v₁ a neighbour of v₀ and orientation σ ∈ {+1, −1}:

  v₃ = 1 + σi(v₁ − 1),  v₂ = v₁ + σi(v₁ − 1).

The same square also appears with the other neighbour as v₁ and the opposite σ. Writing v₁ = z(t₁), the unknowns are (σ, t₁, a). A square exists iff v₂ and v₃ lie on S_a with |v₂|, |v₃| < 1.

### R. Both neighbours of v₀ have radius at least √2 − 1

Since v₁ lies in the closed unit disk, σi(v₁ − 1) lies in the unit disk centred at −σi, so v₃ lies in the unit disk centred at 1 − σi. Hence |v₃| ≥ √2 − 1. The same argument with the roles of v₁ and v₃ swapped gives |v₁| ≥ √2 − 1. Therefore

  t₁ = −ln|v₁|/a ≤ ln(1+√2)/a < 0.8814/a,

where ln(1+√2) < 0.8814 is checked in Arb. This makes the search space compact.

### L. Tiny squares are impossible

Fix a band a ∈ [a₀, a₁]. Arb certifies inf over t ≥ 6/5 of |z(t) − 1| > M, where M = min(3/5, ⌊(4/5)(1 − e^{−2πa₀})⌋₁₀₀₀). The range 6/5 ≤ t ≤ T is covered by adaptive boxes, and the tail t ≥ T by |z(t)| ≤ e^{−a₀T} < 1 − M. The centre is at distance 1 > M.

Pick a rational ε with √2·√(1+a₁²)·ε < M; the code uses √2 < 1415/1000 and √(1+x) ≤ 1 + x/2. If t₁ ≤ ε, then the side is s = |v₁ − 1| ≤ arc length ≤ √(1+a²)·t₁. Every vertex is then within √2·s < M of v₀, so v₂ and v₃ have parameters below 6/5. All four vertices then lie on z([0, 6/5]).

Along a log spiral, the tangent direction turns by exactly Δt. Its tangents on z([0, 6/5]) therefore stay within 3/5 < π/4 of a fixed axis, which makes that arc a strict 1-Lipschitz graph in rotated coordinates. A square has two perpendicular sides, and one of them has slope at least 1 in absolute value. Contradiction. So t₁ ∈ [ε, 0.8814/a₀].

### B. Branch and bound

A point w with 0 < |w| ≤ 1 lies on the spiral iff u = w·exp(i·ln|w|/a) is a positive real. For a box of (t₁, a) values and a fixed σ, Arb encloses v₂, v₃ and their u-values. The box is excluded if, for v₂ or v₃, Arb shows any of:

- **outside:** |w|² > 1, contradicting maximality of v₀;
- **off the spiral:** Im u ≠ 0;
- **wrong arm:** Re u < 0.

The test on v₂ is skipped where |v₂| may be 0; v₃ never is, by Lemma R. Boxes that survive are bisected, and a box still alive at depth 96 makes the run fail. All inequality decisions are certain Arb comparisons; floating point only picks split directions.

The run covers [3/50, 9/10] in 168 bands of width 1/200, for both orientations, at 128-bit precision: 39,626 boxes, well under a second. The smallest ε is 0.1769. The only boxes not excluded are the ones inside the two Krawczyk boxes below.

### K. Exactly one square

Write F(t₁, a) = (Im u₃, Im u₂). For each orientation, the Krawczyk operator

  K(X) = c − Y·F(c) + (I − Y·J(X))(X − c)

is evaluated on the box X = [t₁ ± 10⁻⁹] × [a\* ± 10⁻¹⁰]. Here J(X) is an Arb enclosure of the Jacobian, from forward-mode differentiation over balls, and Y approximately inverts J(c). Arb shows K(X) lies inside the interior of X, so F has exactly one zero in X.

On all of X, Arb also shows 0 < |v₂|, |v₃| < 1 and Re u₂, Re u₃ > 0. So that zero is a genuine square on S_a, and there is no other zero in X.

Outside the two boxes, B excludes everything. So for each σ there is exactly one solution (t₁, a) in the whole range. The relabelling (σ, v₁) ↦ (−σ, v₃) maps solutions to solutions, so the two solutions are one square seen from its two neighbours, at a single pitch a\*. Numerically the σ = −1 parameter t₁ = 5.0789146… equals the σ = +1 parameter of v₃. ∎

## Controls

- Without the Krawczyk boxes, the band [0.1, 0.105] gets stuck exactly at (t₁, a) ≈ (1.0607064138, 0.1002746561), the real square.
- The Krawczyk test fails on a box shifted off the root, and on the wrong orientation.
- Running B on [0.9, 2.0] independently re-derives v239's hand-proved range; the main run already covers v240's range [17/20, 9/10].
- Below the range, B stalls near a ≈ 0.05507 and a ≈ 0.05134, and more densely further down. That is consistent with resonances crowding together as the spiral tightens. These stalls are **not analysed**: they may be true squares or places where the test loses resolution because v₂ sits near the centre. That is why the theorem starts at 3/50.

## What it means

The dimension count from the brainstorm predicted this picture. Once the spiral's one-parameter symmetry is removed, a square has 3 free parameters against 4 conditions. So squares should appear only at isolated pitches, which is what the theorem shows on [3/50, ∞): exactly one such pitch. A curve that is exactly a log spiral with a > a\* near a point has no squares of its own shrinking into that point. This does not reach the Toeplitz conjecture.

## Trust base

Arb (via python-flint) for ball arithmetic; Python `Fraction` for exact rationals. Floating point is used only to choose split directions and the Krawczyk matrix Y, and any Y gives a valid test.

## Reproduce

~~~sh
pip install 'python-flint>=0.6,<1' 'mpmath>=1.2,<2'
python code/log_spiral_classification_v241.py
REQUIRE_SPIRAL_AUDIT_DEPS=1 python -m unittest -v tests.test_log_spiral_classification_v241
~~~
