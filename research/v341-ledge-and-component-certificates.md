# v341. Large squares on all 476 terrains in the v336 sample

Benjamin John Schulz · November 2026

For each of the 282 single terrains and 194 glued terrains in the [v336 experiment](v336-robust-terrain-squares.md), a finite rational certificate proves a robust lower bound of $i/\sqrt2$, where $i$ is the largest contained open axis-aligned square side. Every strict cutoff below the certified floor holds for every sufficiently close **registered Jordan curve**. Its parametrization may introduce arbitrarily fine detail. The neighborhood generally depends on the reference and cutoff.

| Sample | Rows | Rows with a nonzero labeled degree and floor at least $i/\sqrt2$ |
|---|---:|---:|
| v338 single terrains with an interior excursion | 30 | 30 |
| v338 end-wall terrains | 252 | 252 |
| v336 glued two-terrain unions | 194 | 194 |

These counts refer to the original seeds, including repeated shapes. The smallest certified squared-side ratio in each of the two v336 samples is $1/2$. The [certificate ledger](../results/v341/component-certificates.csv) gives the generating data, labeled seed square, chart count, exact squared floor, regular direction, signed chart contributions and degree for every row. A degree of $+1$ or $-1$ is obtained in every row. This proves a finite collection of instances of Conjecture 336.1. A proof that every terrain supplies such a component remains open.

## What changing the marked ledge buys

[v339–v340](v339-v340-end-wall-selection.md) left 24 end-wall rows open. The original strict/one-corner test used the base and three vertical contacts. We tried two changes. Mark any boundary segment as the base in a quarter-turned frame. Then drop the restriction that the opposite contact must project between the two side walls, and allow two vertices to use the same segment. Both extensions preserve the linear contact equations. With the between-walls restriction and three distinct wall abscissas, only two of the 24 pass; the full [ledge selector](../code/terrain_ledge_v341.py) proves the target on 11. The component calculation below certifies all 24 independently.

For a horizontal marked contact at height $y$ and ordered wall roles at abscissas $a,k,b$, the square is
$$
(a+b-k,y),\quad(b,y+b-k),\quad(k,y+b-a),\quad(a,y+k-a).
$$
The three abscissas may occur in any order. Exclude the collapsed case $a=k=b$. Require strict interior contacts on the marked segment and opposite wall. For the two remaining contacts require either two strict hits or one strict hit and a corner crossing. With the right wall fixed, the free left point traces the line $X+Y=k+y$. With the left wall fixed, the right point traces $X-Y=k-y$. Testing the signs on both incident rays gives the one-corner degree of v340. The two points on a repeated segment have independent parameters.

### Both smallest misses have sharp robust size $\sqrt5$

| Heights | Counterclockwise vertices in original coordinates | Reference maximum | Robust size |
|---|---|---:|---:|
| (1,1,5,2) | (1,1), (3,0), (4,2), (2,3) | $\sqrt5$ | $\sqrt5$ |
| (1,4,2,2) | (0,1), (2,0), (3,2), (1,3) | $\sqrt5$ | $\sqrt5$ |

For the first square there is a short proof one can check without enumerating any charts. Put the first vertex on the ledge $y=1$, the second on the base and the fourth on the tower wall $x=2$. If the first edge is $(d,-1)$, completing the square forces the first vertex to be $(1,1)$ and the third to be $(2+d,d)$. This point traces $X-Y=2$. At $(4,2)$, the right wall below has $X-Y>2$ and the ledge to the left has $X-Y<2$. The crossing gives degree of absolute value one. The independent full chart census bounds every reference square's squared side by 5 on both fixtures (19 and 28 positive charts up to cyclic relabeling, respectively). The reference itself gives the matching robust upper bound.

**An explicit tolerance for fine detail.** On either polygon let $\varepsilon$ be the registered sup error measured in the coordinate max norm. For $0\le\varepsilon<1/24$, every such Jordan boundary has a square of side at least
$$
\sqrt5-26\sqrt2\,\varepsilon .
$$
In particular $\varepsilon\le 1/100$ gives side greater than 1.86 grid units, above the target $\sqrt2$.

To verify it, give each of the four short marked arcs its oriented unit-speed coordinate centered at the reference contact. Use the parameter cube of radius $1/2$. Each contact is at least one unit from the next different polygon vertex in each applicable edge direction. Only one of the four parameters straddles a corner. On each half of the cube the square residual below is $A u$. In the displayed labeling both matrices have determinant $+1$ and inverse maximum row-sum norm 3, checked from their four signed unit tangents. The determinants agree, so this isolated piecewise-linear zero has local degree one (the two images occupy opposite sides of their common face). The boundary residual is at least $1/6$ in max norm. Moving each marked point by at most $\varepsilon$ changes the residual by at most $4\varepsilon$. Homotopy invariance therefore gives a perturbed zero inside. There $\|u\|_\infty\le12\varepsilon$, and each vertex is within $13\varepsilon$ in max norm of its reference point. The ordinary Euclidean triangle inequality gives the displayed floor.

## A component certificate from an exact seed

Let $P$ be an embedded rational polygon. Mark four boundary parameters independently on the parameter torus and write $p,q,r,s$ for their evaluated points. With $J$ the positive quarter turn, use
$$
F=(r-q-J(q-p),\ s-p-J(q-p)).
$$
A positive-size zero is a counterclockwise labeled square. On an ordered choice of four closed polygon edges this map is affine in their four unit-interval parameters, with equations $Au=b$ at zero. We orient the four parameters along the fixed polygon traversal. We keep each label fixed during tracing and degree counting.

1. **Trace the connected component.** Start at an exact seed square. Enumerate all incident edge choices, including both edges at a corner. In each chart use the [v223 exact box solver](../code/polygon_orbit_census_v223.py) to find every vertex of the feasible polytope $Au=b$, $0\le u\le1$. Continue from each such vertex, including across the seam in the boundary parametrization. This exhausts the component: a transition between charts lies on a common box face, whose feasible slice is a face of either feasible polytope. Each nonempty face contains a vertex, so the traversal discovers every neighboring chart. Each visited slice is convex and connected. Simplicity of the polygon identifies each geometric vertex uniquely with its boundary parameter. The number of charts is finite.
2. **Bound the size everywhere on it.** All accepted chart slices are points or segments. Along a segment, the first square-edge vector is $e_0+t(e_1-e_0)$, $0\le t\le1$. Project zero onto that segment and clamp the parameter to obtain the exact minimum squared side, including an interior minimum. Take the least value over the entire component. The implementation declines traces meeting collapse or requiring a higher-dimensional minimum.
3. **Compute the degree.** For $k=2,3,\dots$, try a small positive residual target $\tau v$, $v=(1,k,k^2,k^3)$. On each rank-deficient chart require $v$ outside the column space of $A$, so the shifted chart contributes no zeros. On a full-rank chart the reference parameter $u$ shifts exactly to $u+\tau A^{-1}v$. At any active face require a nonzero normal velocity. Count the chart precisely when all active velocities point inward, with sign $\operatorname{sgn}\det A$. There are only finitely many forbidden proper linear subspaces or hyperplanes for $v$. No nonzero linear functional vanishes identically on the moment curve, so a suitable integer $k$ exists. In the recorded rows $k=2$ always suffices. This displacement is applied to the residual for counting; it need not be realized by moving a curve.
4. **Transfer to nearby boundaries.** The full zero set on the compact torus is a finite polyhedral complex. A traced component disjoint from collapse has a small open neighborhood whose boundary is zero-free and whose closure meets no other component. All small shifted zeros in that neighborhood lie on traced charts: otherwise a sequence tending to zero on one fixed closed chart would limit to a reference zero there. Thus the signed count is its local degree. Registered perturbations change $F$ uniformly. The boundary margin preserves degree, and a nonzero degree forces a zero. As the neighborhood shrinks to the component, every strict cutoff below its minimum side follows. This argument does not require a bound on edge count, derivatives, or oscillation frequency of the perturbation.

Integer lattice squares supply convenient seeds. A successful certificate follows the entire real component in rational edge charts. The search makes no claim that every component has a lattice seed. The code and directed signed ledger can be replayed independently. [Component implementation](../code/polygon_component_v341.py) · [Corpus replay](../code/terrain_replay_v341.py)

### A boundary case of the method

On the perimeter of the ordinary $2\times2$ square the marked component through its four corners has squared-side floor 2, spans 60 labeled charts and has signed degree **zero**. The regular target gives two positive and two negative contributions. Thus the proposed universal claim that the signed degree of this full labeled component must always be odd already fails on a square. This example concerns the certificate, and poses no difficulty for inscribing a square in the reference boundary. On the two v337 fixtures the same computation recovers the earlier odd degrees and squared floors $2$ and $49/2$, independently of their special shears. A cyclic change of the chosen seed label can change the sign convention.

The remaining research question is to select some configuration region carrying a large square from the geometry of an arbitrary terrain, or supply another persistence argument when labeled degrees cancel. Larger samples alone cannot settle that question. For an arbitrary limiting Jordan loop one must additionally keep a common positive size floor along the chosen approximation sequence.

## Replay

Install NumPy for the archived v336 generator imports. The component calculation uses standard-library rational arithmetic. From the repository root run

```sh
python code/terrain_replay_v341.py
python -m unittest discover -s tests -p 'test_terrain_component_v341.py' -v
```
