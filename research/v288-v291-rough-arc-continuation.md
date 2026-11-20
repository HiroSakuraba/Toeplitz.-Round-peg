# v288–v291 — A full rectangle witness through rough arc replacements

The concave finite-coil curve has a full flat-to-flat rectangle path, with circumradius at least $1/2$. Another patch on the opposite edge can be added to that path. Two near-square patches can also overlap in the same rectangle. A coupled degree certificate then connects both flats with circumradius at least $499/1000$. With two specific finite-coil replacements, the same connected set contains a square at both rough accumulation points.

## Reference path

Number the vertices of $P=((0,0),(4,0),(1,1),(0,3))$ and their forward edges $0,1,2,3$. The persistent v280 path has one contact on each edge, in geometric CCW order. With $x=w/(w+h)$ for successive side lengths $w,h$, set $H=1+3x+x^2$. Its source parameters are

$$
 u(x)=\frac{1}{H}
 \left(1-\frac34x-\frac14x^2,\ 4x+x^2,\ 1-x^2,\ \frac{10}3x+\frac53x^2\right).
$$

The fixed-angle affine system has determinant $12H$. The verifier substitutes the four numerator polynomials and checks their derivative signs $-,+,-,+$ on $[0,1]$. Therefore every source patch has an unambiguous pair of cut points on the reference path. The end angles $x=0,1$ are flats and their diagonal-dot products are $25,-2$. Edge $0$ lies at height zero, and edge $2$ has height at least one. This gives circumradius at least $1/2$ along the entire reference path, including the flats. The square in the middle has vertices $(9,0),(17,9),(8,17),(0,8)$ divided by $11$ and squared side $145/121$.

## v288 — Replacing the bottom contact by the finite-coil arc

Replace the bottom segment from $(2,0)$ to $(3,0)$ by a Jordan arc in $[2,3]\times[-1/8,1/8]$ with the same endpoints. Require the replacement to meet the remaining boundary only at those endpoints. The v287 finite-coil arc and all its truncations satisfy this condition.

Let the first corner travel from $(3,0)$ to $(2,0)$ along the replacement, and write it as $q=(a,b)$. Setting $z$ to be the horizontal coordinate of the third corner, the midpoint equations give

$$
 p_0=q,\quad p_1=(4-3u,u),\quad
 p_2=(z,3-2z),\quad p_3=(0,b+(5+a-5z)/3),\qquad u=(4-a-z)/3.
$$

Perpendicularity is equivalent to

$$
 F(z)=(4-a-z-3b)(5+a-5z)-9az=0.
$$

Throughout the box, $F(0)\ge35/8$ and $F(1/2)\le-9/16$. For the latter bound maximize first at $b=-1/8$, then use that the resulting quadratic decreases on $[2,3]$. Also $F'(z)\le-145/8$ for $0\le z\le1$. There is a unique admissible root in $0<z<1/2$, continuous in $q$. All three other contacts lie strictly inside their assigned segments. For example $u>1/6> b$, and the fourth corner has height strictly between zero and three. The first side has positive horizontal component $z$, and the fourth corner is above $q$ and to its left. This gives the required CCW orientation at the perpendicular solution.

At each end of the arc the original path has the same anchor and satisfies the same contact equation with $z\in[0,1]$. Uniqueness identifies the full corner tuple at each join. The horizontal separation of $p_0,p_3$ along the inserted lift is $a\ge2$, so its circumradius is at least $1$. Joining to the two retained reference pieces gives a continuous flat-to-flat path with circumradius at least $1/2$, the same endpoint frames, and the exact reference square. The construction gives uniform convergence of the lifted paths for the compatible finite-coil parametrizations, because the lift is continuous on the compact box.

## v289 — A second patch at a different part of the path

On edge $2$ replace source interval $[9/100,11/100]$ by a Jordan arc in $[22/25,23/25]\times[117/100,123/100]$, retaining the endpoints and the stated boundary disjointness. Anchor $p_2=q=(a,b)$ on this arc. Put $p_0=(p,0)$. The midpoint equations then give

$$
 u=(4-a-p)/3,\qquad p_1=(4-3u,u),\qquad p_3=(0,b-u).
$$

The scalar equation is $G(p)=-ap+u(b-u)=0$. On the given box $G(3/20)\ge189/5000$, $G(1/2)\le-2471/40000$, and $G'(p)\le-16/75$ whenever $0\le u\le1$. The root in $(3/20,1/2)$ is unique in the feasible supporting chart. It gives $0<u<1$, $0<b-u<3$, and a positive CCW cross product $a(b-u)+pu$. Its circumradius is at least $117/200$ by the height of the diagonal $p_0p_2$. The feasible-root uniqueness gives both exact joins to the reference path.

The first patch cuts the reference path before $x=1/4$ and the second cuts it after $x=3/5$. During the bottom lift the edge-$2$ parameter is $1-z>1/2$. During the opposite lift the edge-$0$ parameter is $p/4<1/8$. Thus each lift stays on the retained portion of the other edge. The two joins can be made independently, leaving the reference square, both flats, and the whole $1/2$ radius bound intact. The test for composing such substitutions is strictly monotone source coordinates on the reference path, a continuous unique contact root on each enclosure, and clearance from the other removed intervals. A second finite-coil arc fits in the opposite enclosure by mapping its unit segment to the patch chord of length $\sqrt5/50$. Its disks extend at most $\sqrt5/400<1/100$ from that chord.

## v290 — Overlapping modifications with fixed end rectangles

Use the fixed-angle equations $F_x$ of v284 and the reference contact vector $u(x)$. Let $u_*=u(1/2)$. Independently replace the intervals of halfwidth $\eta=1/100$ about $u_{*,0}$ and $u_{*,2}$ on edges $0,2$. Assume compatible source parametrizations with coordinate error at most $\varepsilon=1/1000$ on those intervals and equality to the reference on their complements in the near-square region. Let $I=[2/5,3/5]$ and use the moving box $D_x=\{u:\|u-u(x)\|_\infty\le\eta\}$. The entire box stays in the open source cube. At either end of $I$ the minimum gap from a patched coordinate center is $205/3476>2\eta$. Thus each end box uses exactly the reference affine maps and has a unique zero.

Writing $A(x)$ for the affine reference matrix, the verifier bounds $\|A(x)^{-1}\|_\infty\le2$ on $I$. It checks every signed row sum of the adjugate against $2\det A(x)$ by rational polynomial signs. In coordinates centered at $u(x)$ the preconditioned change from the identity is at most $4\cdot2\varepsilon=1/125$, leaving a lateral margin of $1/500$. The straight homotopy avoids zero on the lateral boundary and the slice degree is $+1$. Apply the connected-continuation argument of v284 on this fixed-width cylinder. A compact connected zero set meets both end plates and hence contains both exact join rectangles. Joining it to the unmodified outer pieces gives a compact connected witness between the original flats. Within the coupled part, opposite corners have height difference at least $1-2\varepsilon$, giving circumradius at least $499/1000$ on the whole witness.

The same tube can be combined with the two separated substitutions. Its edge-$0$ and edge-$2$ windows miss the far source patches. The bottom lift stays at edge-$2$ parameters above $1/2$, clear of the near-square patch. The opposite lift stays at edge-$0$ parameters below $1/8$, clear of its near-square patch. The reference sections in between are available by strict monotonicity.

In all these constructions the chosen ordered labels define a lift of the witness to the cyclic labeling cover. A nowhere-zero extension of the flat frame would have a constant sign on the connected lift. The flat signs $+,-$ rule out that extension. The normalized side coordinate is continuous at the positive-radius flats and ranges from zero to one along the witness. Consequently every prescribed aspect ratio occurs there. The coupled theorem uses a compact connected set for general continuous replacements. The separate scalar lifts give an explicit path.

### A common square at two finite-coil accumulation points

On each patched edge insert a scaled copy of the arc in [Asano–Ike, Proposition 6.6](https://arxiv.org/html/2609.11598v3#S6.SS2) starting at $u_{*,i}$, using the forward half of its source patch. Retain only the coils indexed by $n\ge4$. The arc and its finite truncations retain the reference point at the accumulation parameter $s=0$. In the unit source interval the error on the $n$-th coil is at most $2R_n$, and $s\ge 2^{-n}-R_n$ there, where $R_n=2^{-3n}$. Hence the error divided by $s$ is at most $2/255$. Scaling along either edge gives

$$
 \|\widetilde p_i(u_i)-p_i(u_i)\|_\infty
 \le K|u_i-u_{*,i}|,\qquad K=8/255.
$$

Here the edge lengths are at most four. This bound also implies the required coordinate error bound on the whole patch. Both copies are embedded in disjoint neighborhoods of their own edges. The preserved non-Liouville obstruction follows from the primitive excursions on indefinitely late coils. The same estimates apply to the finite truncations, which converge uniformly in their source parametrizations.

At the square angle the original contact vector $u_*$ remains a zero. Every zero in the tube satisfies $\|u-u_*\|_\infty\le 4\cdot2K\|u-u_*\|_\infty$, and the coefficient is $64/255<1$. Thus the reference square is the unique square in the tube and belongs to its selected connected crossing on the limit and every truncation. In the limit, two opposite corners of this very square lie at the two finite-coil accumulation points.

## v291 — Failure of a sequential lift and the fold control

Replace each near-square patch by a two-segment tent with the same endpoint parameters. At the central point shift edge $0$ by $(0,-1/2000)$ and edge $2$ by $(1/2000,0)$. This is a simple rational polygon. Anchor the bottom contact at the displaced kink and try to keep the other three contacts on their original lines. The bottom equation from v288 still has negative derivative on $[0,1]$, and its only possible edge-$2$ contact has $z\in(789/1100,811/1100)$. This entire interval is deleted from the old edge-$2$ segment. The shifted tent meets that original line only at the two endpoints. Thus this sequential lift fails with straight contacts on the other three tagged edges.

The coupled certificate supplies a crossing here. The error of each tent is Lipschitz in its source parameter with constant $(1/2000)/\eta=1/20$. In each moving box, solving the preconditioned equation is a contraction with constant at most $4\cdot2/20=2/5$, and it maps the closed box strictly inside itself. Consequently the coupled solution is a unique continuous path in the tube. At its square both patched parameters are strictly inside their two tent intervals, so both contacts are on changed pieces. An independent exact piecewise-linear solve verifies this square and the unchanged join rectangles.

The v281 time-fold remains a separate control for turning back at positive radius. On the extra strip time is $t=|s|$. The transverse scalar equation $|s|-t=0$ has two regular roots with opposite degrees for $t>0$, no root for $t<0$, and a root of degree zero at $t=0$. The nonzero relative class within each rectangle interval does not supply a nonzero crossing degree in the time direction. The coupled proof requires the verified lateral isolation and end joins along its chosen angle direction.

## Smoothing the original vertices

The v286 rounding construction applies at both concave flats. Measured in the parameter pointing into the reference path, the short-side derivative at $x=0$ is $(3,4)$ and the unnormalized inward longitudinal midpoint projections on $(-4,3)$ and $(4,-3)$ are $60,30$. At $x=1$ the short-side derivative is $(-1,1)$ and the projections on $(1,1)$ are $1,2/5$. All four original vertices can therefore be rounded in arbitrarily small collars while keeping the respective witness, the flat frames, and the original lower bound on radius. Patch charts stay a positive distance from those vertices, so the rounding can be confined to the flat portions. In particular the finite-coil variants can have a smooth boundary away from their rough accumulation points.

## Replay and scope

Run `python code/rough_arc_continuation_v288_v291.py --check` for the rational certificate. The SageMath tests independently solve the patch equations at the enclosures and both pairs of algebraic join angles. The standard-library tests check the coupled rational tent fixture.

The bounds concern the selected witness on these tagged replacements of the concave path. A theorem selecting such a witness on arbitrary Jordan curves remains open.
