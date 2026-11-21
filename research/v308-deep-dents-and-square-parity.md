# v308 — Large squares on loops with deep dents

Benjamin John Schulz · November 2026

A Jordan curve can make arbitrarily deep excursions into a rectangular hull while a single, explicitly computable gap in square sizes survives. The permitted excursions lie in narrow, separated vertical strips; near the two horizontal walls the boundary can move across the full width. Within these regions it may fold or have infinite length. An odd number of runs of enclosed cells forces a large inscribed square. The parity can also be computed when that number is even.

## The directional question from v307

Write $B(O,r)\subset K\subset B(O,R)$ for a convex body and concentric disks. For **any** square of half-diagonal $q<r$ in $K$, take a vertex $p$ closest to $O$. In side coordinates based at the square center $C$, this is the corner in the quadrant of $O-C$, with either choice allowed at a tie. Its depth below the hull satisfies
$$
\operatorname{dist}(p,\partial K)\ \ge\ \frac{q(r-q)}{R+q}.
$$
Indeed the support estimate in [v307 §1](v307-convex-shell-square-gap.md) only uses the slack at that corner. With $d=\operatorname{dist}(p,\partial K)$ it gives $r-d\le q+dR/q$.

This gives a precise one-scale test. Let $A\subset K$ be a closed permitted region for all point zeros during the deformation from the convex reference to the target. For example $A$ may include the hull boundary and all of $K$ outside the target's bounded interior. One only needs to test candidate $p\in A$ at least the displayed depth. For $|v|=q$, put $C=p-v$. Keep the facing orientations $(O-C)\cdot v\ge |(O-C)\cdot Jv|$, where $J$ is a quarter turn. If at least one of $p-2v,p-v+Jv,p-v-Jv$ lies outside $A$ for every retained pair, a square at the cutoff is impossible. The odd geometric-square degree above it then transfers by v307 §2. This reduces the test to potential facing corners; geometry is still needed to rule out their mates. The next result does this by strip widths and separations and also permits references with several components.

## Strip theorem

Let $K=[L,U]\times[0,H]$ in some orthonormal frame, allowing a translation. Specify closed vertical strips $I_j\times[0,H]$ ordered left to right, including strips abutting both vertical walls. They have positive widths, maximum $w$, and successive intervals have minimum intervening gap $g>0$. Allow two full-width horizontal strips $0\le y\le d$ and $H-d\le y\le H$, with $0<2d<H$. Write $A$ for the union, and put
$$
a=\max(w,\sqrt2\,d),\qquad b=\tfrac12\min(g,H-2d).
$$
Assume $a<b$. Between successive vertical strips there is one open rectangular cell with $d<y<H-d$. For a Jordan curve $\Gamma\subset A$, mark each cell 1 if its midpoint is inside the curve and 0 otherwise. The mark applies to the whole cell by connectedness. A **run** is a maximal consecutive block of 1s.

**Result.** No square with vertices in $A$ can have half-diagonal in $(a,b)$. For each cutoff $q$ in that interval the mod-two degree of the square equations above $q$ equals the number of runs modulo two. If this number is odd, $\Gamma$ inscribes a square of half-diagonal **at least $b$**, equivalently side **at least $\sqrt2 b$**. This is about boundary vertices; the filled square may cross a dent. For a nonisolated family, degree means the parity of a transverse perturbation on the bounded configuration domain.

### Why the gap survives a deep dent

Let $a<q<b$. A square of half-diagonal $q$ has diameter $2q$, so it reaches at most one vertical strip and at most one horizontal strip. The least strip width containing three square vertices is $q$: the three form a right triangle, whose shortest altitude is $q$. The minimum width of a triangle is its shortest altitude (minimize its projection width over directions; between changes of extreme vertex the width is a positive cosine arc, and its minimum is at an endpoint). Each of the two strips has width strictly below $q$. If four vertices are covered, exactly two must be in each strip with no overlap.

If a pair is adjacent, the other pair is the parallel opposite side. Thus the common side vector has horizontal component at most $w$ in magnitude and vertical component at most $d$. Its length squared is at most $w^2+d^2<2q^2$, a contradiction.

If the two horizontal-strip vertices are opposite, their midpoint is the square center, at depth at most $d$ from the corresponding horizontal wall. All corners lie inside that wall. By central symmetry all four have depth at most $2d$. But the vertical span of a square of side $\sqrt2q$ is at least $\sqrt2q>2d$. This proves the gap. The fact that the horizontal strip adjoins an exterior wall matters: two perpendicular *interior* lines can carry the diagonals of a whole square.

### How the cell marks determine parity

For each run form a closed reference rectangle of full height $H$, with its two vertical sides down the midlines of the two vertical strips bounding the run. These disjoint rectangles have width at least $g$. Two distinct rectangles are separated horizontally by at least $g$, since between their runs lies an unoccupied cell. Their boundaries lie in $A$. Let $f$ be signed distance to their union boundary, negative on their interiors. For no runs take $f=1$. Let $f_\Gamma$ be signed distance to the Jordan boundary, negative inside. On every cell and outside $K$, the two signs agree strictly. Linear interpolation from $f$ to $f_\Gamma$ can therefore have point zeros only in $A$. The gap proved above excludes four simultaneous corner zeros at $q$ throughout.

Compute the reference degree as follows. Translate successive rectangles horizontally by successively larger nonnegative distances until every gap exceeds $H$, monotonically increasing the gaps along the way. Their signed-distance functions vary continuously. A square touching two components has diameter at least $g>2q$. A square on one component has half-diagonal at least its inradius by v307 at depth zero; its inradius is at least $\min(g,H)/2>q$. Thus the cutoff face stays zero-free, also during the translation. At the end a mixed square is impossible: its horizontal span would exceed $H$ while its vertical span is at most $H$, and a square has equal spans on perpendicular axes. The compact individual zero sets are disjoint. Near a pure zero all four signed-distance evaluations agree with those for its single rectangle (shrink the neighborhood using the positive separation). By excision the endpoint degree is the sum of the individual degrees. Each is odd by the ellipse-to-convex-body computation of v307 §2, using its own inball of radius greater than $q$. The empty reference has degree zero.

These degree statements use the geometric square angle modulo $\pi/2$ and the corresponding permutation of the four residuals as in v307. Take the outer bounds on center and half-diagonal sufficiently large for every compact zero-level set in both deformations. Homotopy invariance gives the claimed parity. In the odd case any zero above $q$ must already have radius at least $b$ by the gap. No choice or matching of squares across refinements is required. An even count leaves existence above the cutoff undecided.

## Examples and strict comparisons

Take $[-6,6]\times[-2,2]$, horizontal cap thickness $1/32$, vertical edge strips of width $1/32$ and internal strips of width $1/4$ centered at $x=-2,2$. Here $a=1/4$ and $b=15/8$. Cut a rectangular dent down from the top in the left internal strip and another up from the bottom in the right internal strip. Each stops just $1/64$ short of the opposite wall. The curve is simple; the cells read $111$. Every Jordan curve in this same strip arrangement with these marks has a square of side at least $15\sqrt2/8$. There is no lower bound on the remaining tip-to-wall distance needed for the specified Jordan curves. One may also roughen every portion inside the strips while preserving simplicity and the marks.

This example defeats even a change of convex reference in the **uniform-depth** v307 test. The four box corners remain on the curve, so every enclosing convex body contains the box. The curve contains a dent-wall point at height zero whose distance from the box boundary is 2, forcing any uniform layer thickness $\delta\ge2$. But the v307 sufficient condition implies $\delta<r/4$ since $R\ge r$ and $q(r-q)\le r^2/4$. Enclosure of $O$ and avoidance of the open core put the open disk of radius $r-\delta$ about $O$ inside the target (as in v307). It must fit in a strip of height 4, so $r-\delta\le2$. These inequalities conflict. This comparison concerns the v307 certificate; it makes no priority claim about other square existence theorems.

The replay also checks marks $101$ and $10101$ obtained by wide bottom indentations ending in the top cap. The first has even degree above the gap; in the actual fixture a side-four square uses the vertical lines $x=-6,-2$. The second has three occupied runs and the same guaranteed lower bound as the one-run example. A geometric negative control is a diamond of half-diagonal one on crossing interior lines. Moving the horizontal strip away from an exterior wall invalidates the diagonal exclusion.

**Extra perturbation allowance.** Suppose a second Jordan curve is within registered Euclidean sup distance $\eta$ of a certified input. Enlarge the outer rectangle by $\eta$ on all sides, every vertical interval by $\eta$ on each end, and each cap thickness to $d+2\eta$. Coordinatewise proximity puts the new curve in this enlarged arrangement. Its bounds are $a'=\max(w+2\eta,\sqrt2(d+2\eta))$ and $b'=b-\eta$. If $a'<b'$, the old horizontal and vertical cell midpoints are unchanged and have original clearance at least $b>\eta$ from the old allowed region. Straight interpolation of the two parametrizations avoids each midpoint, preserving its winding and hence the cell marks. With odd runs the perturbed curve therefore has side at least $\sqrt2(b-\eta)$. For the opposing dents, $\eta=1/8$ still gives side at least $7\sqrt2/4$ whenever the target remains Jordan. This allows displacement eight times the tip-to-wall distance and requires no derivative control.

For a uniformly convergent sequence of registered curves with the same parameters and odd run parity, a subsequential limit of the geometric squares has the same size floor on the limiting curve, since centers, radii and angles lie in a common compact set. For a general Jordan curve, the existence of coordinates and separated strips with an odd number of occupied runs remains an additional hypothesis.

## Independent square counts

For each fixture we also perturb the vertices by less than $1/1024$, check the cell certificate on the enlarged strips, and run the repository's full v223 square census at half-diagonal cutoff 1. This second computation solves the edge-incidence equations, including repeated-edge and rank-deficient charts and both orientations. Each retained geometric square has regular interior contacts and eight labeled representations.

| Perturbed input | Runs | Predicted parity | Actual squares above cutoff | Oriented edge charts checked |
|---|---:|---:|---:|---:|
| Opposing deep dents | 1 | odd | 19 | 40,392 |
| Two occupied runs | 2 | even | 2 | 7,728 |
| Three occupied runs | 3 | odd | 5 | 40,392 |

The positive counts on the even example illustrate why zero degree leaves existence undecided. These are counts on the specified perturbed polygons. The proof above applies to all Jordan curves satisfying the strip and winding conditions.

## Replay

Run `python code/deep_dents_v308.py` and `python -m unittest discover -s tests -p 'test_deep_dents_v308.py' -v`. The [rational verifier](../code/deep_dents_v308.py) accepts `--input input.json` with keys `polygon` (cyclic coordinate pairs), `lanes` (ordered horizontal-coordinate intervals), `ymin`, `ymax`, and `cap`. Fractions may be quoted. It verifies simplicity, covers every entire polygon edge by clipped strip intervals, computes the cell windings, and returns the [parity and physical size bound](../results/v308/deep-dents.json). The analytic gap and the parity argument establish the theorem for continuous Jordan curves.

Run `python code/deep_dent_census_v308.py opposite_dents` for the one-run independent enumeration, substituting `two_runs` or `three_runs` for the other two fixtures.
