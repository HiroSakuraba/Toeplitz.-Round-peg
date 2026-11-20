# v279 — A protected square in a uniform neighborhood of the atlas polygon

**Theorem (restricted).** Let \(P:\mathbb R/4\mathbb Z\to\mathbb R^2\) trace the v274 concave quadrilateral \(((0,0),(4,0),(1,1),(0,3))\) affinely on the four unit parameter intervals. If a Jordan embedding \(\gamma\) satisfies

\[
\sup_t\|\gamma(t)-P(t)\|_\infty\le 1/1000,
\]

then \(\gamma\) contains an inscribed square with side strictly greater than \(1\). It has one parameter in each of four fixed intervals of half-width \(1/100\), centered on the four boundary parameters of the interior square selected by the second v275 witness. The bound is in the original units of the polygon.

## Geometric class and its local continuation

**V274-path.** Identify the relative path in the second quotient interval. Its endpoints have opposite transported signs and its only square is the orbit of

\[
((9/11,0),(17/11,9/11),(8/11,17/11),(0,8/11)).
\]

This is node 337 of the exact artifact. In its counterclockwise lift the four local edge parameters are \(t_*=(9/44,9/11,3/11,25/33)\) on edges \(0,1,2,3\). The first relative component has a different witness and is not used in this continuation.

**V279-degree.** Let \(C=\prod_{j=0}^3[t_{*,j}-1/100,t_{*,j}+1/100]\). Use the square residual with vector components \(p_2-p_1-J(p_1-p_0)\) and \(p_3-p_0-J(p_1-p_0)\). Each boundary point in \(C\) stays on four interior edge intervals in the reference polygon. Its residual is \(A(t-t_*)\), where \(A\) is invertible and \(\|A^{-1}\|_\infty=1\), as computed in the exact rational replay. Thus the residual has norm at least \(1/100\) on \(\partial C\), and has degree \(\operatorname{sgn}\det A=\pm1\) on \(C\). With the same source parameters the residual for \(\gamma\) differs coordinatewise by at most \(4/1000<1/100\). The straight residual homotopy stays nonzero on \(\partial C\) and gives a zero in the interior. This continues the local degree at the specific square carried by the geometric witness in V274-path.

**V279-size.** The reference side at \(t_*\) is \(\sqrt{145}/11>12/11\). The first two reference edge directions each have \(1\)-norm \(4\); source movement in \(C\) changes their side vector by at most \(8/100\) in Euclidean norm. Two vertex perturbations contribute at most \(4/1000\). The actual side is therefore greater than \(12/11-8/100-4/1000=2769/2750>1\). A zero of the oriented square residual with positive first side has four distinct corners. For uniformly converging embeddings within the bound, compactness of the four fixed source intervals gives the same protected square on any Jordan limit. Equivalently one can apply the theorem directly to the limiting embedding.

No part of this proof uses an index from an unrelated configuration component. It uses the local degree of the square on the exact relative path at the base polygon. Extending the full graph class over arbitrary perturbations is a separate question.

## Broader continuation and remaining gates

The v275 relative obstruction yields a square on its own representative whenever that graph and its boundary framing are admissible. The v277 secant-cone estimate controls **all** near-square configurations on the specified approximants. It can therefore bound the size of a geometric representative without identifying it with a spectral generator. In the separate filtered route, Greene–Lobb's derivative theorem and the v277 localized lower slope imply a protected square for the selected spectral carrier. Identifying that carrier with the geometric relative class would require a comparison that preserves both the boundary data and the filtration. [Greene–Lobb, Theorem 5.1](https://arxiv.org/html/2604.17116v1#S5)

| Gate | Established here | Remaining for unrestricted Toeplitz |
| --- | --- | --- |
| Definition | Complete atlas and compatible relative flat boundary frames for \(P\). A relative class on verified positive-radius graphs. | A compatible compactification and framing for arbitrary approximants, especially at total collision and artificial size boundaries. |
| Nonvanishing | Both components of the concave test polygon carry a nonzero relative class. | A nonzero protected representative for every required approximation. |
| Comparison and nonescape | The selected atlas square has a local degree and uniform size bound. Conditional secant, spectral and local cap estimates give other sufficient controls. | The relevant geometric and filtered carriers have no proved identification in general, and the uniform size or action bound is open in that category. |
| Limit | Uniform convergence of parameterized embeddings, compactness of the source circle, and a common size floor imply distinct limiting square corners. | Establish the hypotheses for arbitrary Jordan limits. |

The v276 counterexamples defeat global area as a substitute for local action control. The v278 exact polygon tests defeat odd raw orbit counts and size control of arbitrarily selected branches. The unrestricted conjecture remains open.

To reproduce the computational inputs, replay v274 in SageMath and run `python code/relative_index_v275.py`. The other exact controls and their commands are listed in the [tranche ledger](v274-v279-program-ledger.md).
