# v216: Boundary-computable degree and an a priori polygonal arc-selection experiment

**Research program:** Benjamin John Schulz  
**Prepared:** 18 November 2026  
**Status:** Conditional finite-dimensional theorem, exact-rational implementation, two certified convex quadrilateral examples, and falsifying controls. **This does not prove the unrestricted Toeplitz conjecture.**

## 1. Purpose and relationship to v215

[v215](v215-knot-robust-relative-degree.md) computes the relative degree by counting *preimages of a nearby regular target* of the four-point square residual. Its `Arc` type requires increasing horizontal coordinates, and its sample source arcs are chosen around planted square configurations. Although that degree is valid under the stated conditions, these implementation choices limit its usefulness for testing a global arc-selection argument.

This batch makes three changes:

1. It permits **arbitrary oriented planar polygonal arcs**, including vertical segments, negative horizontal slopes and nonmonotone geometry.
2. It computes the same Brouwer degree **from the four-dimensional configuration-space boundary alone**, with no interior square/root enumeration in the production algorithm.
3. It applies the boundary test to **a priori chosen** four-arc windows from a rational simple polygon, using all four distinct vertex neighborhoods and, if needed, full polygon edges.

The second and third changes are important for avoiding a circular rule of the form *find a square, then select four arcs near that square and call their degree nonzero*. The candidate windows here are fixed solely by the input polygon.

## 2. The square residual on arbitrary planar arcs

Let `A_i` be compact oriented planar polygonal arcs, each parametrized by a continuous piecewise-linear map `p_i:[0,1]→R²`, `i=0,1,2,3`. The arc parameter is the **orientation of the original polygon**, not necessarily horizontal position. Fix integer deck shifts `k=(0,k1,k2,k3)` and an orientation `σ∈{+1,−1}`. Put `z_i(u_i)=p_i(u_i)+(k_i,0)`, `J(x,y)=(-y,x)` and

```
Fσ(u0,u1,u2,u3) = (
    z2 - z1 - σ J(z1 - z0),
    z3 - z0 - σ J(z1 - z0)
) ∈ R⁴.
```

Then `Fσ=0` if and only if the four labeled points form an oriented square, *possibly of side zero*. The parameter domain is the oriented box `Q=[0,1]⁴`.

Require the two **physically lifted** opposite arcs `A0+k0` and `A2+k2` to have strictly positive Euclidean set distance `D>0`. This excludes a zero-side square and ensures that any resulting square has side at least `D/√2`. Unlike v215, this separation condition is **not** based on vertical projections; it survives arbitrary rotations.

The implementation computes the exact rational squared distance `D²` from all pairs of source segments. It checks the actual physical deck lifts, not merely their cylinder projections.

## 3. Main theorem: a boundary-only, rational degree certificate

**Theorem 1 (boundary-ray certificate).** Assume the arcs have rational vertices, `D>0`, and

```
0 ∉ Fσ(∂Q).
```

Then the Brouwer degree `deg(Fσ,Q°,0)` is an integer computable by a terminating exact-rational algorithm. There is a computable rational `m>0` such that `||Fσ(u)||∞ ≥ m` on `∂Q`. If the degree is nonzero, a nondegenerate oriented square exists with one vertex on **each designated arc in its specified physical lift**.

**Proof A, boundary zero detection.** Every product of four linear arc segments is a 4-dimensional affine chart for `Fσ`. A global boundary patch fixes exactly one arc parameter at one endpoint; each of the three remaining parameters lies in `[0,1]`. Consequently `Fσ=0` on that patch is a rational linear feasibility problem with three box-constrained variables. Existence is decided exactly by rank reduction followed by enumeration of enough active box constraints. This also detects nonisolated and endpoint solutions. There are finitely many patches. If a boundary zero exists, the algorithm refuses to assign a relative degree.

**Proof B, positive margin.** On every zero-free boundary patch, `Fσ` is affine. Dyadically subdivide its three-dimensional parameter cube. The exact coordinatewise range over each subbox can be calculated from signs of the rational coefficients. Once one coordinate interval excludes zero, the subbox has a positive rational `l∞` residual lower bound. Since the patch is compact and zero-free, sufficiently fine uniform boxes all satisfy such a sign condition; the adaptive subdivision terminates without a resource cap. Taking the minimum of all positive subbox bounds gives `m>0`. The implementation imposes a work cap and explicitly reports **inconclusive** if exceeded.

**Proof C, boundary ray intersection formula.** Choose a rational vector `r∈R⁴\{0}` such that the positive ray `λr` (`λ>0`) intersects the images of the oriented global boundary patches transversely, away from their internal 2-dimensional seams. For a boundary patch fixing coordinate `i` to `s∈{0,1}`, let the affine map in its three free parameters be

```
F(v)=b + C0 v0 + C1 v1 + C2 v2.
```

At a ray intersection, solve the exact rational 4×4 system

```
[C0  C1  C2  -r] · (v0,v1,v2,λ)ᵀ = -b.
```

Only solutions with `0<vj<1` for all `j` and `λ>0` contribute. A singular patch that meets the ray or any ray crossing at a seam invalidates the direction; try another. The oriented boundary of `Q` has incidence sign `(-1)^i` for the upper face and `(-1)^(i+1)` for the lower face. The intersection contribution is

```
upper: (-1)^i       sign det[C0 C1 C2 -r],
lower: (-1)^(i+1)   sign det[C0 C1 C2 -r].
```

The sum is the degree of the normalized boundary map `F/||F|| : ∂Q→S³`, equal to `deg(F,Q°,0)`. This is the standard boundary formulation of Brouwer degree, specialized to rational affine patches. **The algorithm never calls an interior square solver in order to obtain this count.**

An admissible ray exists generically because all bad directions lie in the finite union of proper linear spans of 2-dimensional patch seams and lower-rank patch images. The deterministic candidates `r_n=(1,n,n²,n³)` avoid every fixed proper linear subspace for all but finitely many positive integers `n`: a nonzero linear equation in `r_n` is a nonzero polynomial of degree at most three. The mathematical search therefore terminates without a cap; the implementation returns inconclusive if its candidate budget is exhausted. Degenerate affine patches are handled by exact *bounded* feasibility before discarding them, rather than by pretending they cannot intersect the ray.

## 4. Arbitrary two-dimensional C⁰ continuation

**Theorem 2 (stability).** Under the hypotheses of Theorem 1, suppose the four marked parameterized arcs undergo **arbitrary continuous planar perturbations** moving each point by at most `ε` in each coordinate, while preserving its marked identity and physical deck shift. Each scalar component of the residual changes by at most `4ε`, so the boundary remains zero-free if `4ε<m`.

The perturbed opposite-arc distance satisfies `D_new ≥ D−2√2 ε`. Choose any positive *rational* `L≤D` (computable from exact `D²`). The conservative sufficient inequality

```
ε < min(m/4, L/4)
```

ensures boundary exclusion and nondegeneracy throughout the straight-line homotopy. Thus the **same nonzero degree** forces a nondegenerate square on the perturbed marked arcs. This covers wild continuous graph and nongraph perturbations. It is a conditional **open-neighborhood theorem** around a certified polygon, not a uniform claim about arbitrary Jordan curves.

Unlike v215's vertical-only perturbation bound, this result permits shifts in **both** spatial coordinates and does not require the arcs to remain locally graphs. It is generally more conservative.

## 5. A priori global arc selection on rational Jordan polygons

For a simple closed polygon `P=(P0,…,P_{n−1})` (implicitly cyclic, no duplicated terminal vertex), choose for each vertex `Pi` the oriented arc from midpoint `(P_{i−1}+Pi)/2` through `Pi` to `(Pi+P_{i+1})/2`. These `n` windows are determined before any square computation. Enumerate all 4-element index subsets in curve order, and check both square orientations with physical shifts `(0,0,0,0)`. If none yields a nonzero degree, repeat with the original full polygon **edges** as candidate windows. Validate the simple-polygon condition exactly before scanning.

This is a *finite sufficient criterion*, not a decision procedure for polygonal square existence. It can return **inconclusive** when the degree cancels, a square touches a selected source boundary, or no candidate windows capture a useful configuration.

### Positive concrete example A: preselected vertex windows

For the convex quadrilateral

```
P = ((0,0),(3,1),(2,3),(-1,2))
```

all four vertex-centered half-edge windows are chosen *without locating a square*. The exact boundary-only calculation, with counterclockwise orientation, returns:

| Quantity | Exact value |
|---|---|
| Signed degree | `−1` |
| Number of transverse positive-ray intersections | `1` |
| Boundary residual lower bound `m` | `1/8` |
| Squared opposite-arc distance `D²` | `13/4` |
| Certified arbitrary planar `l∞` perturbation limit | **`ε<1/32`** |

Therefore **every continuously parametrized closed curve with four corresponding marked arcs each uniformly within `1/32` in both planar coordinates of these four selected polygonal arcs has a nondegenerate inscribed square**, provided the four arcs remain portions of the same curve in the correct cyclic order. This includes sufficiently small rough (not necessarily polygonal or differentiable) perturbations. The conclusion is conditional on the source-arc correspondence and does not cover an arbitrary Jordan curve.

### Positive concrete example B: full-edge windows needed

For the convex quadrilateral

```
P = ((0,0),(2,0),(3,3),(0,2))
```

the vertex-window signed count cancels to zero, but the next predetermined window family (the original four edges) produces degree `−1` and planar perturbation guarantee **`ε<1/64`**.

### Negative control: an exact square can defeat these windows

For `P=((0,0),(1,0),(1,1),(0,1))`, the polygon is already a literal square, yet the tested vertex/edge windows have zeros on their *external boundaries*. The local degree calculator therefore reports a boundary failure and the global selection routine returns **inconclusive**. This directly falsifies any assertion that these two fixed window families must always produce a certified nonzero degree on every square-bearing polygon.

A zero-degree result is also inconclusive, as v214–v215 demonstrated with two nondegenerate square roots whose signs cancel. No general positive lower bound for the degree on arbitrary Jordan curves is claimed.

## 6. Verification and limits

Run:

```bash
python -m unittest discover -s tests -v
python code/boundary_degree_v216.py
```

The v216 suite includes exact checks of linear-algebra rank and singular feasibility, independent interior-root-sign crosschecks on 72 random planted one-cell cases, reversal and 90° rotation covariance, internal PL knots, source-boundary escape, physical deck-shift isolation, nonmonotone arcs, two-arc Euclidean separation, 15 rational perturbations of the positive polygonal case, a positive edge-only selection, and a *known-square/no-certificate* negative control. All arithmetic in the certificate and tests is exact `fractions.Fraction`; the tests do not constitute proof-assistant certification. The checker does **not** implement Floer homology.

The program does not yet process arbitrary continuous Jordan curves directly, guarantee a nonzero degree under an exhaustive topological selection, or derive a relative Floer invariant from the October 2026 periodic-inscription theorems. Its main research value is a **noncircular source-arc selection audit** and an exact boundary observable that can be compared with proposed global geometry or topological invariants.

## 7. Honest next obligations

1. Extend a priori windows to multiscale partitions without destroying a fixed global topological or spectral invariant. A naive finite scan cannot be a general theorem; the square-polygon control already proves that.
2. Seek a *global* invariant whose restriction to marked physical source arcs is known to have nonzero total contribution, rather than selecting arcs after an unknown square is presumed.
3. Treat cancellation and the diagonal/collapse strata explicitly. Nonzero local degrees are sufficient but do not imply a universal nonzero global sum.
4. Compare the boundary observable with the October Floer theories only after specifying a precise continuation or pairing theorem. Neither theorem is presently known to imply this boundary count for arbitrary Jordan curves.

### Source context

- Internal mathematical base: [v215 relative degree](v215-knot-robust-relative-degree.md), [v214 source/lift-specific degree](v214-relative-degree-source-lifts.md), and [v213 lift audit](v213-exact-PL-localization-and-lift-audit.md).
- Lorenzo Valerio, *A Proof of the Area Vanishing Square Peg Conjecture*, arXiv:2610.21076 (2026), https://arxiv.org/abs/2610.21076 .
- Ali Naseri Sadr, *The Four Periodic Inscription Problem*, arXiv:2610.24376 (2026), https://arxiv.org/abs/2610.24376 .