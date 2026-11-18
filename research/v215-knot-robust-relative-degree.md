# v215 — A chart-independent, boundary-certified source/lift degree for piecewise-linear arcs

**Research program:** Benjamin John Schulz  
**Prepared:** 18 November 2026  
**Status:** A finite-dimensional exact-rational theorem, algorithm, and adversarial regression suite. **No proof of the unrestricted Toeplitz square-peg conjecture.** No claim to construct Floer homology.

## Executive result

The v214 relative-degree algorithm detected the signed count of *transverse, isolated* squares on four prescribed source arcs and their fixed physical integer lifts, but it deliberately rejected configurations in which a square lay on an **internal PL knot** or a feasible affine cell was singular. These exclusions were properties of the numerical chart, not necessarily of the underlying topology.

v215 removes both restrictions under one intrinsic condition: **the square-residual map must have no zero on the external product boundary of the four marked arcs**. It detects that condition exactly, builds a certified positive boundary margin, and evaluates the Brouwer degree at a nearby rational regular target. Internal knot crossings and singular affine cells cannot invalidate this regular-target count.

The work makes two conclusions precise:

- **Positive result:** Under a nonzero relative degree and an opposite-arc separation condition, at least one nondegenerate square lies on the *specific four source arcs in their specified physical lifts*. The conclusion survives a quantitatively controlled class of continuous perturbations.
- **Negative result:** This does not produce a nonzero degree for arbitrary Jordan curves. The periodic Floer-inscription theorems do not themselves supply the necessary marked arcs, lifts, or boundary-avoidance homotopy.

## 1. Exact setting

Let \(A_i\), for \(i=0,1,2,3\), be compact polygonal **graph arcs** with strictly increasing rational horizontal coordinates, each oriented left-to-right. Set \(Q=[0,1]^4\) and parametrize each arc continuously by \(p_i:[0,1]\to\mathbb R^2\), with each affine segment traversed monotonically. The computational parameter uses local \([0,1]\)-coordinates on each segment; the orientation agrees with increasing original arc parameter.

Fix integers \((k_0,k_1,k_2,k_3)=(0,k_1,k_2,k_3)\), an actual choice of *physical deck shifts*, and let \(\widetilde p_i=p_i+(k_i,0)\). Choose a cyclic orientation \(\sigma\in\{+1,-1\}\) and let \(J(x,y)=(-y,x)\). Define

\[
F_\sigma(u)=\begin{pmatrix}
\widetilde p_2-\widetilde p_1-\sigma J(\widetilde p_1-\widetilde p_0)\\
\widetilde p_3-\widetilde p_0-\sigma J(\widetilde p_1-\widetilde p_0)
\end{pmatrix}\in\mathbb R^4.
\]

The zero set \(F_\sigma^{-1}(0)\) consists of ordered oriented squares, **including possible zero-side degeneracies**. To exclude them, require the vertical ranges of \(A_0\) and \(A_2\) to be disjoint by a positive amount

\[
g=\operatorname{dist}(\operatorname{pr}_y A_0,\operatorname{pr}_y A_2)>0.
\]

For any zero, the two opposite vertices differ vertically by at least \(g\), so the square side is at least \(g/\sqrt2\). This is a sufficient geometric guard, not a necessary property of all square configurations.

## 2. Theorem: computable relative degree despite internal knots

**Theorem (finite PL relative degree).** Assume \(g>0\) and

\[
0\notin F_\sigma(\partial Q).
\]

Then \(\deg(F_\sigma,Q^\circ,0)\in\mathbb Z\) is well defined and can be calculated by a terminating algorithm using **rational arithmetic only**. The algorithm tolerates true square zeros on internal PL chart interfaces and even non-isolated families of zeroes. It also provides a rational \(m>0\) satisfying

\[
\|F_\sigma(u)\|_\infty\ge m\qquad (u\in\partial Q).
\]

**Proof, part A: exact boundary-zero decision.** The source arcs induce finitely many affine product cells. On a global boundary face \(u_i\in\{0,1\}\), the corresponding vertex is fixed. For each product of the three remaining segments, solving \(F=0\) amounts to solving a rational affine linear system with three free parameters in a box. A nonempty affine slice of a compact box has an extreme point. At an extreme point, sufficiently many independent coordinate bounds are active. Exhaustive enumeration of these bound choices, followed by exact Gaussian elimination, decides feasibility, including singular families and endpoint roots. This check is **only** on the actual global source boundary; internal chart faces are not treated as exit boundaries.

**Proof, part B: exact positive boundary residual margin.** On each zero-free boundary PL cell, \(F\) is affine in three variables. For a rational subbox, the exact minimum/maximum of each residual coordinate lies at a corner and can be computed from coefficient signs. If one coordinate has a fixed nonzero sign throughout the subbox, its distance from zero gives an explicit positive lower bound on \(\|F\|_\infty\). Otherwise subdivide dyadically. Because the cell is compact and contains no zero, every point has at least one nonzero coordinate and hence a neighborhood on which that coordinate retains its sign. Uniform continuity and compactness imply that sufficiently small boxes all receive a sign certificate. Thus the adaptive dyadic procedure terminates in exact arithmetic and yields a positive rational lower bound \(m\). The implementation has a configurable work cap and raises `CertificateInconclusive` rather than guessing if that cap is reached.

**Proof, part C: regular-target evaluation.** For rational \(\delta>0\), let

\[
q(\delta)=(\delta,\delta^2,\delta^3,\delta^4).
\]

Choose dyadic \(\delta=2^{-n}<m/2\). Every singular affine chart and every internal chart face has an affine image in a proper affine subspace of \(\mathbb R^4\). The moment curve \(q(\delta)\) is not contained in any proper affine hyperplane: a putative hyperplane condition would give a nonzero polynomial in \(\delta\) of degree at most four. Since there are finitely many bad images, only finitely many candidate values of \(\delta\) can be bad. Therefore a sufficiently small dyadic target is simultaneously outside every singular chart image and every chart-face image.

On each nonsingular chart, the equation \(F(u)=q\) has at most one solution. Find it by rational Gaussian elimination, retain it iff it is strictly interior to the chart, and add the sign of the exact \(4\times4\) Jacobian determinant. The result is

\[
\deg(F_\sigma,Q^\circ,0)
=\deg(F_\sigma,Q^\circ,q)
=\sum_{F_\sigma(u)=q}\operatorname{sgn}\det D F_\sigma(u).
\]

The first equality follows from the straight-line homotopy \(F-tq\): the boundary norm is at least \(m-\|q\|_\infty>0\). The second equality is the regular-value formula for Brouwer degree. Notice that the **regular-target preimages are not squares**; their signed sum certifies the zero-target degree and hence square existence if nonzero.

## 3. Persistence under continuous perturbations

Suppose the four selected arcs undergo arbitrary continuous **vertical** perturbations \(h_i\), retaining their horizontal parameterizations and physical deck shifts, with \(\|h_i\|_\infty\le\varepsilon\). Each coordinate of the square-residual perturbation has size at most \(2\varepsilon\), and the opposite-arc vertical gap decreases by at most \(2\varepsilon\). Consequently, if

\[
\boxed{\varepsilon<\min\{m/2,g/2\},}
\]

then the whole straight-line homotopy is nonzero on the external parameter boundary; the degree cannot change, and any degree-nonzero witness remains nondegenerate. This establishes a **source/lift-preserving continuation theorem** for the stated perturbation class, not for arbitrary connector motions or arbitrary Jordan-curve regularizations.

The algorithm's \(m\) is **conservative**, not the exact minimum boundary distance. The representative four-line example yields degree \(-1\), margin \(m=1/384\), positive opposite vertical gap \(g=5/64\), and the conservative strict persistence threshold \(\varepsilon<1/768\). The earlier single-chart v214 proof supplies a different, sharper estimate for that example; the v215 goal is substantially broader admissibility, not optimal margin size.

## 4. Adversarial examples and conclusions

| Experiment | Exact finding | Significance |
|---|---|---|
| Transverse planted square | Degree \(-1\); target preimages: 1 | Agreement with v214 on its domain |
| Collinear refinement exactly through the square | Degree remains \(-1\) | Result does not depend on arbitrary extra PL knots |
| Genuine kink through the square | Degree remains \(-1\) | Handles a root on an internal chart interface |
| Four horizontal levels with wrong physical deck shifts | Degree \(0\), positive boundary certificate | Singular empty charts need not disable robustness |
| Two real squares of opposite sign | Degree \(0\), two regular-target preimages | Zero degree **does not imply no squares** |
| Endpoint square at a marked boundary | `BoundarySquare` | Refuses to assign invariant across an exit event |
| Varying one source arc vertically | Degree \(-1\) before \(t=55/624\), undefined at the boundary contact, degree \(0\) after | Explicit wall-crossing test; source-boundary escape is real |
| Internal knot translated through square | Degree remains \(-1\) | Knot movement is not an external-boundary event |

The wall-crossing example follows from an exact affine root calculation: in the seed configuration the four local parameters are \((1/2,1/2,1/2,1/2)\); shifting arc 1 upward by \(t\) moves the first root to the boundary at \(t=55/624\). We test rational values on either side and at the contact itself. This is **not** a theorem about all possible crossing events; it exhibits one mechanism clearly.

## 5. What this does *not* resolve

1. **Global source-arc selection:** We have not proved that an arbitrary Jordan curve admits four arcs with nonzero source/lift-relative degree.
2. **Connector filtering:** Floer-theoretic existence for four balanced periodic strands does not automatically confine a square to the selected original planar arcs or deck lifts. The exact v212/v213 counterexamples remain valid.
3. **Zero-degree invisibility:** Pairs of squares can have opposite signs; a zero count cannot be read as nonexistence.
4. **Irregular Jordan curves:** Arbitrary rough Jordan curves need not permit the finite polygonal graph/vertical-homotopy model used here. No uncontrolled limiting argument is assumed.
5. **Floer functoriality:** No relative Floer chain complex, action filtration, continuation morphism, or proof of a bridge from the October 2026 periodic theorems has been constructed here.

**Recommended v216 target:** test source-compatible *relative homotopies* of the four-arc residual, and look for a principle supplying a nonzero degree from global shape information. Make a proof or counterexample for each proposed transfer rule before attempting a Floer formalism. In particular, investigate whether an arc-selection process can avoid the v210 spiral and v213 deck-lift obstruction **without assuming the desired square's existence**.

## 6. Reproduction and references

```
python -m unittest discover -s tests -v
python code/relative_degree_v215.py
```

All new arithmetic uses Python standard-library `fractions.Fraction`. Exact algebraic checks do not constitute a machine-checked proof in Lean, Coq, or another proof assistant. The mathematical arguments above are independent of the implementation.

- Internal starting point: [v214 source/lift degree](v214-relative-degree-source-lifts.md), [v213 finite exact PL localization](v213-exact-PL-localization-and-lift-audit.md), and archived v209–v211 research.
- Lorenzo Valerio, *A Proof of the Area Vanishing Square Peg Conjecture*, October 2026, https://arxiv.org/abs/2610.21076 .
- Ali Naseri Sadr, *The Four Periodic Inscription Problem*, October 2026, https://arxiv.org/abs/2610.24376 .
- Tomohiro Asano and Yuichi Ike, *Liouville regular Jordan curves*, arXiv:2609.11598v3, 17 November 2026, https://arxiv.org/abs/2609.11598 .
