# v280–v283 — Moving relative classes and refinement

The concave quadrilateral has a complete positive-radius moving atlas on $|t|\le 1/100$. One relative interval persists through the deformation. A second moving component projects to time by an absolute-value fold, giving zero, one, or two interval fibers. The same framing construction yields a uniform positive-radius path on both v276 diamond spirals and every truncation. The shrinking-notch control has a square in its local rectangle family, where the relative class vanishes.

## v280 — Complete moving atlas for the concave quadrilateral

Let $P_t=((0,0),(4,t),(1,1),(0,3))$, with vertices and forward edges numbered $0,1,2,3$. Use CCW corner labels $p_0,p_1,p_2,p_3$, and let $x$ be the length of $p_0p_1$ divided by the sum of its length and the length of $p_0p_3$. For $0<x<1$, the equations are

$$
 p_0+p_2=p_1+p_3,\qquad
 x(p_3-p_0)=(1-x)J(p_1-p_0),
$$

where $J$ is the counterclockwise quarter turn. They give a four-by-four affine linear system in the four edge parameters. Its determinant and Cramer column determinants have degree at most two in $x$. The square angle has $x=1/2$.

The opposite edges $0,2$ have vertical separation at least $99/100$, and edges $1,3$ have horizontal separation at least $1$. The adjacent directions are nonparallel and the opposite edges are disjoint, so every $P_t$ in the window is Jordan. In fact all six pairs of directions are nonparallel, as checked by their determinants. A genuine rectangle uses at least three supporting edges. If two edges contain two corners each, the two corresponding pairs are opposite sides or the two diagonals. Opposite sides require parallel directions. The diagonals would make two polygon edges intersect in their interiors, which simplicity forbids. Three corners on one edge are impossible. Any set of at least three edges includes an opposite pair, so every rectangle and its limits in this atlas have circumradius at least $99/200$.

The checker enumerates the 168 ordered edge charts in 42 cyclic orbits. For each representative it factors the determinant and the eight Cramer wall polynomials $N_i,D-N_i$, with $x,x-1/2,x-1$ as additional walls. It projects leading coefficients, coefficient denominators of the monic factors, discriminants, and pairwise resultants onto the $t$-axis. Within each open time chamber, continuous simple roots of the distinct factors keep their multiplicity in the original polynomials, their order, and their incidences with each wall. The signs within root bands and on their root sections are fixed. This makes one exact fiber audit sufficient on each chamber. The only interior projection time across all charts is $0$. The two endpoints and $0$ are audited directly. At each representative interior singular angle some Cramer column determinant is nonzero, which rules out consistency. The projected root incidences preserve that certificate on the two time chambers. The checker also rejects isolated feasible points outside the listed bands.

The table gives all feasible bands. Here $\alpha=(t^2+16)/(t^2-4t+32)$ and $\beta=(4-t)/(8-t)$. The right column gives the closures of the angle bands for $t\ge0$. For $t<0$ only the first row occurs.

| Ordered contacts | Strip | Closed angle band |
|---|---|---|
| $(0,1,2,3)$ | Persistent | $[0,1]$ |
| $(0,0,2,3)$ | First | $[0,\alpha]$ |
| $(0,0,1,3)$ | First | $[\alpha,1]$ |
| $(0,1,3,3)$ | Second | $[0,\beta]$ |
| $(0,2,3,3)$ | Second | $[\beta,1]$ |

In each of the four repeated-edge charts a common factor $x$ cancels from all five Cramer determinants. The reduced denominators of the active charts are zero-free on their closed time and angle domains. This follows from the same projection audit with the marked flat angles. Consequently the formulas extend jointly to the displayed closed bands. Endpoints at $x=0,1$ are included by closure of genuine rectangles. At $\alpha$, the first strip changes the third corner's edge at vertex $2$. At $\beta$, the second strip changes the second corner's edge at the same vertex. Equality of the full ordered source tuples at these seams is a rational-function identity in $t$.

For $t>0$, all contacts in the open bands are strict and the only nonflat contacts on the band boundaries are the indicated seam contacts. Hence the support orbits exclude further identifications. At $t=0$, the first strip has its first corner fixed at vertex $0$, and the second strip has its fourth corner fixed at that vertex. Rotating the first strip's labels once and sending $x$ to $1-x$ gives the second strip identically, including the double-vertex contact at the square. The persistent strip has strict contacts at all genuine angles and stays separate. Normalize source positions in $\mathbb{R}/4\mathbb{Z}$. The six unordered pairs of distinct flat positions are

$$
 \{1,3\},\quad \{0,2\},\qquad
 \{3t/(16+t^2),3\},\quad \{0,1\},\qquad
 \{1,4-t/3\},\quad \{0,3\}.
$$

For $t>0$ these are distinct. At $0$ the first and second strips identify their endpoints in reverse order, exactly as the labeled transition predicts. This accounts for the gluing at both flat angles. The finite union of closed bands is compact and its physical quotient is Hausdorff, so the resulting continuous bijection from the glued model is a homeomorphism.

It follows for **every** $t$ in the window that the geometric square counts are $1,2,3$ for negative, zero, and positive time. The independent v223 affine square census agrees on both endpoints, both open-chamber samples, and the event. The certificate describes the closure of the geometric nondegenerate rectangle locus in time, angle, and corner positions.

## v281 — Relative class through the fold

On every strip there is a continuous choice of CCW label, with diagonal dot signs $+,-$ at $x=0,1$. At a positive-radius flat the two distinct corner positions are paired and cyclic rotation acts freely on the four labels. Thus the labeled sign line and its flat frame from v275 extend over the moving atlas. On any strip interval, comparing the two framed endpoints in its lift gives $w_1(L,\tau)$ evaluation $1$. This agrees with the two v275 evaluations on the central fiber.

Write the persistent component as $[-1/100,1/100]\times[0,1]$, with time $t=s$. The other component is also a strip $[-1/100,1/100]\times[0,1]$, with time $t=|s|$. Use $s>0$ and interval coordinate $x$ for the first chart strip, and $s<0$ and interval coordinate $1-x$ for the second. Rotate the latter's labels back once. The lift and the frames then match on $s=0$. Both strips carry the pullback of the nonzero class of an interval relative to its endpoints.

For $t<0$ only the persistent interval is present. For $t>0$ the folded component meets the time slice in two intervals, each with evaluation $1$, for a total of $0$ in $\mathbb{F}_2$. The same moving component meets $t=0$ in a single relative interval of index $1$. This interval is an interior fold of the time map, so it cannot be used once as the outgoing horizontal boundary of the cobordism from negative times. Its two time-positive sheets form the valid paired boundary contribution. This explains the obstruction to continuing the unit-square interval toward negative time despite the uniform radius bound. The persistent interval continues in both directions. The usual cobordism comparison between regular time levels holds with total parity $1$ on both sides.

## v282 — The shrinking notch has zero local relative class

In the v278 control let $\varepsilon=1/n$ and $a=1/3$. Consider the closure of rectangles with all corners on the notch and the upper rails, at height $y\ge1$. If a corner lies on a rail strictly outside $[a,a+\varepsilon]$, it is at height $1$. If both its adjacent corners are higher, their horizontal displacements have the same strict sign and their vertical displacements are positive, so the two side vectors cannot be perpendicular. If one adjacent corner is also on the rails, a nondegenerate rectangle has a horizontal lower side, putting a higher corner outside the allowed horizontal range. Thus all four corners are supported on the three notch sides. One supporting line contains at least two of the corners. It cannot contain a diagonal since all four corners are in the same closed half-plane of that line. The rectangle is axis-aligned and has width $\varepsilon$, with its two distinct heights between $1$ and $1+\varepsilon$.

In the fixed CCW labeling starting at the lower left, the compactified family is the triangle $0\le h_0\le h_1\le\varepsilon$. Its flat edge is $h_0=h_1$. The diagonal dot product is $(h_1-h_0)^2-\varepsilon^2$, so the flat frame has a fixed negative sign and extends across the triangle. Its relative obstruction is zero. The unique square in the triangle is the apex $h_0=0,h_1=\varepsilon$, and the entire triangle has squared circumradius at most $\varepsilon^2/2$. The separate unit square of the notch polygon persists. This calculation applies to the closure of the specified local family. Membership of the square in a component of the full rectangle locus remains open here.

## v283 — A uniform relative witness on the diamond spirals

Use the outer blocks of the v276 arms and the lower diamond cap. Allow any Jordan replacement of the inner arc from $(1/2,0)$ to $(-1/2,0)$ inside the closed diamond $|a|+|b|\le1/2$. This includes the straight radial truncation, every later monotone or counterwinding truncation, and both infinite limit curves in v276.

Keep the first three contacts on the following fixed segments, in order.

$$
 [(-1,0),(0,-7/8)],\qquad [(-1,0),(0,-1)],\qquad [(0,-1),(1,0)].
$$

Let the fourth corner $q=(a,b)$ traverse the complementary arc from $(1,0)$ inward on the first arm, through the diamond, and outward on the opposite arm until $(0,-7/8)$. For $p_1=(-\ell,\ell-1)$ and $p_2=(r,r-1)$, the parallelogram equation gives $p_0=(a-\ell-r,b+\ell-r)$. The first segment's supporting line yields $\ell=15r-L$, where $L=7+7a+8b$. Put $B=1+b$. Perpendicularity of $p_2-p_1$ and $q-p_2$ becomes

$$
 g(r)=-2r^2+(16a-14B)r+L(B-a)=0.
$$

Take the larger root. The common exterior arc is a finite union of segments. Checking the affine inequalities at their endpoints and at the four inner diamond vertices gives $0\le L\le14$, $B-a\ge0$, $1+a-b>0$, $a+L\ge0$, $b\ge-7/8$, and $|a|+|b|\le1$. Only the start has $B-a=0$, and only the terminal point has $L=0$. Hence in the open arc the root is unique and positive, and depends continuously on $q$, with endpoint values $1,0$.

The supporting-line solutions lie on the actual segments. At $\ell=0,1$ the comparisons are

$$
 g(L/15)=L(1+a-b)/225 \ \ge 0,\qquad
 g((L+1)/15)=-(B-a)(L+226)/225 \ \le 0.
$$

Thus $0\le\ell\le1$ and $0\le r\le(L+1)/15\le1$. For the range of the first corner's horizontal coordinate, the comparisons are

$$
 g((a+L)/16)=((8(a-b)-7)^2+16b+14)/128 \ \ge0,\qquad
 g((a+L+1)/16)=((b-a)^2-1)/2 \ \le0.
$$

This gives $-1\le a-\ell-r\le0$. Finally $g(B)=-(9-7a+8b)(B-a)<0$ on the open arc, so the side $q-p_2$ has strictly positive vertical component. The adjacent side $p_2-p_1$ has strictly positive horizontal component. Together with perpendicularity this proves that the four distinct corners form a CCW rectangle. The endpoints are positive-size flats with diagonal dot products $4,-1/64$.

For every rectangle or flat on the path, $p_0$ belongs to the half-plane $x-y\le7/8$ and $p_2$ to the line $x-y=1$. They are opposite corners, at distance at least $1/(8\sqrt2)$. The squared circumradius is therefore at least $1/512$, uniformly over the whole path and every allowed inner replacement.

In the full geometric rectangle closure on such a curve, restrict to squared radius at least $1/512$ and frame the sign line on the flat subset. The cyclic action is free also at a flat, since its tuple has two distinct paired positions. Pulling the line and flat frame back to the displayed path gives a nonzero relative evaluation. Any global nowhere-zero extension of the flat frame would restrict to one on the path, which its opposite endpoint signs rule out. With the compatible parametrizations of v276, the paths converge uniformly to the path on the limit by continuity of the larger root on the fixed compact domain. The outer portions and both flat frames are identical throughout the sequence.

There is an exact common square on the fixed first arm segment. For a square the extra equation $q-p_2=J(p_2-p_1)$ gives $\ell=a$, $r=(B-a)/2$, and $31a+b=1$. Its four corners, in order, are

$$
 \frac{1}{241}((-225,-14),(-1,-240),(225,-16),(1,210)).
$$

Its squared side is $101252/58081>1$. This square lies on the very path with opposite framed signs. Since the monotone spiral is non-Liouville-regular as proved in v276, the example gives a geometric relative witness on a specific curve beyond that regularity hypothesis.

## Verification and scope

Run `sage -python code/moving_rectangle_atlas_v280.py --check` for the moving atlas and `python code/relative_witness_refinement_v282_v283.py --check` for the exact notch and spiral controls. The latter audit checks sample truncations and the affine domain bounds, endpoint frames and rational square in the argument above. The uniform theorem uses the displayed inequalities for an arbitrary admissible inner arc. The moving atlas and spiral witnesses use different fixed outer geometries. A positive relative carrier for arbitrary Jordan embeddings and the identification with a filtered Floer or sheaf carrier remain open.

- [Fixed concave atlas](v274-full-rectangle-binormal-atlas.md)
- [Relative sign line and cobordism](v275-relative-class-and-cobordism.md)
- [Diamond spiral and non-Liouville proof](v276-liouville-sheaf-rough-frontier.md)
- [Adversarial controls](v278-adversarial-program.md)
