# v217 — Multiscale source-arc degree: exact internal-wall cancellation and an obstruction to naive global counting

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Scientific status:** Conditional finite-dimensional theorem, exact-rational implementation, and reproducible adversarial checks. **Not a proof of the unrestricted Toeplitz square-peg conjecture and not a Floer-theory construction.**

## 1. The question inherited from v216

The v216 boundary-degree certificate works for *four already chosen physical arcs*, fixed relative integer lifts, and a square orientation. On some polygons it gives a nonzero signed Brouwer degree; on other polygons the degree cancels or a square lies on the chosen product's external boundary, making that particular relative certificate undefined. A proposed next step was to subdivide arcs repeatedly, look for locally nonzero counts, and search for a conserved quantity. The key issues are:

1. How do degree contributions combine when original arcs are split?
2. What if a genuine square lies on an **internal subdivision wall**?
3. Can refinement of a parent with **zero** signed degree reveal squares that the parent count hides?
4. Does conservation of local contributions imply a global nonzero square count on an arbitrary Jordan curve?

The first three questions receive an exact answer here. The answer to the fourth is **no implication has been established**.

## 2. Setting and exact definitions

Let \(A_i:[0,1]\to\mathbb R^2\), \(i=0,1,2,3\), be oriented rational polygonal arcs. Fix physical deck shifts \(k=(0,k_1,k_2,k_3)\in\mathbb Z^4\), a square orientation \(\sigma\in\{+1,-1\}\), and \(J(x,y)=(-y,x)\). With \(z_i(u_i)=A_i(u_i)+(k_i,0)\), put

\[
F_\sigma(u)=\left(z_2-z_1-\sigma J(z_1-z_0),\quad
z_3-z_0-\sigma J(z_1-z_0)\right)\in\mathbb R^4.
\]

Its zeros are labeled, oriented Euclidean squares, potentially collapsed. The **noncollapse guard** is the exact positive distance \(D=\operatorname{dist}(A_0+k_0,A_2+k_2)>0\). The algorithm determines \(D^2\) rationally from all pairs of segments. For a true square of side \(s\), \(D\le\sqrt2\,s\), hence \(s\ge D/\sqrt2>0\).

Let \(Q=[0,1]^4\). We require **zero to be absent from the external boundary image** \(F_\sigma(\partial Q)\); this is checked exactly via affine-box feasibility on all product boundary facets. A finite set of rational cut parameters in each arc subdivides \(Q\) into a grid of closed parameter boxes \(Q_\alpha\) whose interiors are disjoint, with inherited original arc orientations and **unchanged** shifts. Splitting is by the equal-segment PL parameterization and preserves every original knot.

For a generic nonzero target ray \(R_r=\{\lambda r:\lambda>0\}\subset\mathbb R^4\), write

\[
I_r(C)=\#_{\mathrm{oriented}}\bigl(F_\sigma(\partial C)\cap R_r\bigr).
\]

Every affine boundary-face intersection is determined by solving four rational linear equations in three face parameters and \(\lambda\), with positivity and strict face-interior tests. Singular faces are tested separately for **positive-\(\lambda\)** feasibility, not confused with a possible \(\lambda=0\) root. Thus \(I_r(C)\) remains an algebraic **face flux** when a *child* boundary contains a square; it is **not** called the child's Brouwer degree then.

## 3. Theorem A: telescoping multiscale boundary chain

**Theorem.** If \(0\notin F_\sigma(\partial Q)\) and a ray \(r\) is generic for all relevant piecewise-affine boundary faces of \(Q\) and its finitely many children, then

\[
\boxed{\deg(F_\sigma,Q^\circ,0)
= I_r(Q)=\sum_\alpha I_r(Q_\alpha).}
\]

Furthermore, for each *child* with \(0\notin F_\sigma(\partial Q_\alpha)\), its ray flux satisfies

\[
I_r(Q_\alpha)=\deg(F_\sigma,Q_\alpha^\circ,0).
\]

**Proof.** The first equality is the standard boundary-ray formulation of Brouwer degree (v216) for a zero-free external boundary. As oriented cubical chains, the child boundaries satisfy

\[
\sum_\alpha \partial[Q_\alpha]=\partial[Q].
\]

Each shared internal three-dimensional face has the **same** restriction of \(F_\sigma\) from both neighboring children but the **opposite** induced boundary orientation. Its generic positive-ray intersection contribution therefore occurs twice with opposite signs. Sum all face contributions: the internal terms cancel exactly and only the parent external boundary remains. A child whose boundary avoids zero has a well-defined Brouwer degree represented by the same ray count. If a child contains a zero on a cut face, only the formal positive-ray face count survives; it is incorrect to identify that child's count as a degree. QED.

**Existence of a generic exact ray.** The product PL subdivision has finitely many affine faces and seams. Nongeneric positive-ray directions lie in finite unions of proper linear spans associated with singular images and face seams. The rational moment-curve candidates \(r_n=(1,n,n^2,n^3)\) evade each fixed proper linear subspace for all but finitely many integers \(n\); thus a suitable rational direction exists. The implementation searches a configurable finite budget and reports `Inconclusive` if it is exhausted. Its signed algebraic contributions are computed using exact determinants, linear solves, and rational box-feasibility certificates, not floating point or gridded sampling.

**Corollary (clean cuts).** When every child has zero-free boundary, the parent's degree equals the sum of the *actual child degrees*. If the parent degree is nonzero, at least one child has nonzero degree and therefore contains a square on its assigned original physical subarcs (under \(D>0\)).

**What this does not say.** If a cut passes through a square, some child degrees are **undefined**. Their formal ray fluxes can depend on which generic ray is chosen. The parent sum is nevertheless invariant. The calculation does not magically supply a nonzero parent degree from an arbitrary Jordan curve.

## 4. Exact experiments and counterexamples

The code in `code/multiscale_degree_v217.py` includes two named rational constructions and an independent boundary-chain evaluator.

### A. One isolated transverse square

Four short linear arcs support a square with parent degree \(-1\). A clean preselected cut at \(u_0=1/3\) yields two child degrees \(0\) and \(-1\). Cutting instead at \(u_0=1/2\) passes through the square: **both** children have square zeros on their internal boundary and neither has a valid individual degree, yet their formal fluxes still sum to \(-1\).

There is an exact demonstration of ray dependence for these wall children:

| Generic ray \(r\) | Formal child fluxes | Parent degree |
|---|---|---|
| \((-2,-2,-2,-2)\) | \((-1,0)\) | \(-1\) |
| \((-2,-2,-2,1)\) | \((0,-1)\) | \(-1\) |

The change in allocation across the wall must **not** be interpreted as a change in the square-existence theorem.

### B. Signed cancellation conceals two genuine squares

Replace opposite arc \(A_2\) by the two-segment rational tent described in the v214–v215 cancellation tests. Both roots are transverse, with opposing local signs. The full four-arc parent has degree \(0\) but two actual squares. A clean cut of \(A_2\) at its center PL knot yields two children of **degrees \(-1\) and \(+1\)**. Their sum is \(0\), as required.

Thus a multiscale scan *can* reveal local sufficient certificates when its parent degree is zero, but **does not produce a new nonzero global index**. Any future argument that derives a universal nonzero invariant merely from subdividing an algebraically cancelling signed count is invalid.

### C. Higher dimensional grids and perturbation controls

The regression suite exercises a two-axis grid (four children, four paired interfaces), a four-axis grid (16 children, 32 paired interfaces), 10 seeded rational perturbations with varying clean cuts, source-arc reversal, incorrect physical deck shifts, an external-boundary square (refused), and an opposite-arc collision (refused). It also cross-checks parent counts against the *independent* v216 boundary-degree implementation when that module is installed. These are finite-model checks, not a proof-assistant formalization.

## 5. What is proved, what remains unresolved

| Assertion | Status |
|---|---|
| Exact source/lift-aware degree conservation under finite PL subdivision | **Proved conditionally**; implemented and tested |
| Interior-wall contributions cancel even when a square lies on the cut | **Proved** as a chain identity; do not assign child degrees |
| Refinement can reveal opposite-sign local squares invisible to a zero parent degree | **Constructed exact example** |
| Arbitrary polygon admits a preselected four-arc partition with globally nonzero total degree | **Not proved; naive assertion is unsupported** |
| Arbitrary Jordan curve has a uniformly noncollapsed, nonzero-degree four-arc selection | **Open** |
| This finite cubical chain count is a relative or filtered Floer invariant | **Not established** |

The remaining central problem is **finding the right global object** that produces nonzero total information, rather than simply redistributing an index that may already be zero. Possible next steps are: (i) construct a space of ordered four-point configurations on a single Jordan curve with collapsed/diagonal configurations blown up, (ii) understand equivariant cancellation and orientations under cyclic relabeling, and (iii) see whether the October periodic-inscription Floer complexes admit genuine source-arc and deck-lift filtrations with a nontrivial continuation map. No such bridge is claimed here.

## 6. Reproduce

```bash
python -m unittest discover -s tests -v
python code/multiscale_degree_v217.py
```

The v217 module is standalone, standard-library-only, and independent of the v216 implementation; the GitHub test suite also compares their parent degrees. For a clean child, its formal flux is a degree. For a child with a square on its boundary, the count is only a ray-dependent algebraic flux. The physical lift tuple, original arc labels, and orientation are never silently changed.

**Sources:** Internal v213–v216 research notes and source in this repository. External area-balanced periodic results: Lorenzo Valerio, arXiv:2610.21076 (2026), https://arxiv.org/abs/2610.21076 and Ali Naseri Sadr, arXiv:2610.24376 (2026), https://arxiv.org/abs/2610.24376. These preprints do **not** establish the global nonzero-degree arc-selection step discussed above.