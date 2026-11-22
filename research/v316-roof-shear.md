# v316. A mesh-independent square bound from a base and an excursion

Benjamin John Schulz · November 2026

For a square-grid terrain polygon, let $i(P)$ be the largest side of an open axis-aligned square in its interior. Its boundary carries an axis-aligned lattice square of side **at least $i(P)$**, with its lower side on the base. The factor one is sharp and independent of the grid size. A terrain polygon here is a consecutive row of positive-height columns filled from a common horizontal baseline. Its roof may have any number of steps, peaks and valleys. Only the four vertices of the resulting square must lie on the curve.

For arbitrary grid polygons, Pettersson, Tverberg and Östergård proposed the bound $i(P)/\sqrt2$ and verified their conjecture through grid size 13. Here the stronger bound holds for the stated subclass at every grid size. A second result below permits a rough base and an overhanging upper arc and gives an explicit region of odd degree. [A Note on Toeplitz' Conjecture (2014)](https://doi.org/10.1007/s00454-014-9578-5)

## Why a horizontal shear detects a square

Put the base on $y=0$. Suppose an arc $A$ starts at $(l,h)$ and ends at $(r,h)$, stays strictly above $h$ between them, and $r-l>h>0$. Assume the vertical projections of its points belong to the base segment. Define $T(x,y)=(x+y,y)$. The endpoints of $A$ and $TA$ on the line $y=h$ occur in the order
$$
l,\quad l+h,\quad r,\quad r+h.
$$
Consequently the arcs intersect above that line. At an intersection $q=Tp$ the original arc contains $(x,y)$ and $(x+y,y)$. Projecting them down to the base gives the other two vertices of a square of side $y>h$.

More precisely, $(t,u)\mapsto A(u)-TA(t)$ has odd planar degree on the parameter square. Its boundary is zero-free since the only arc points at height $h$ are the four distinct endpoints of the two traces. Homotope the traces independently, keeping their endpoints fixed and their interiors in $y>h$, to the two upper semicircles on their respective endpoint intervals. Straight interpolation of parametrizations has these properties. The semicircles intersect transversely once. This proof works for continuous arcs with nonisolated intersections. At equality $r-l=h$, the two endpoints themselves provide a square of side $h$.

## From an interior square to the excursion

Let the Jordan boundary consist of a base $[0,W]\times\{0\}$ and a return arc $B$ in the strip $0\le x\le W$, with its interior in $y>0$. Orient $B$ from left to right. Suppose $(a,a+m)\times(0,m)$ lies in the bounded component for some $0\le a<a+m\le W$. Then the boundary carries a square of side at least $m$ on the base, even if $B$ folds back or crosses a horizontal line infinitely often.

Fix $0<h<m$. The compact set of parameters on $B$ with height at most $h$ splits into two disjoint compact sets, one having $x\le a$ and the other $x\ge a+m$. The arc starts in the first and ends in the second. Take the first visit to the second set and the last preceding visit to the first. Between these times the arc stays strictly above $h$. Both excursion endpoints have height $h$ by continuity, and their horizontal separation is at least $m>h$. Apply the shear argument. As $h$ increases to $m$, boundedness of the vertices gives a convergent subsequence of squares with limiting side at least $m$. The same argument passes to the supremum over admissible $m$.

In a terrain every contained open axis-aligned square can be lowered to the baseline, by vertical filling. This proves the stated bound. On the integer grid a largest base-supported square can be chosen with integer coordinates. Intersect the full polygonal roof, including the lateral end walls, with its shear. The maximum intersection height is attained at an endpoint of an intersection of two closed segments. The original segments have integer horizontal or vertical supporting lines. The sheared segments have integer horizontal supporting lines or lines $x-y=k$ for integer $k$. Every isolated intersection is integral, as are endpoints of collinear overlaps. The upper contact and its original partner therefore both have integer coordinates. Rescale for a different mesh size. A rectangle wider than its height gives equality because a square's vertical projection is at least its side. Since an inball of radius $r$ contains an open axis-aligned square of side $\sqrt2 r$, this also proves the inradius bound $\sqrt2 r$ for terrains.

**An approximation rule.** For a continuous positive roof $y=f(x)$ over $[0,1]$ with vertical end walls, take columns of width $1/n$ and round the left-sampled height down to the $1/n$ grid, with minimum height $1/n$. If $\omega$ is a modulus of continuity, the height error within each column is at most $\omega(1/n)+1/n$. Match each graph interval to its horizontal step and following vertical connector, using the same parameter interval and omitting a zero connector. On the last interval use the last horizontal step, and match the end walls linearly. The registered sup error is at most $\omega(1/n)+1/n$ in each coordinate. Every grounded interior square of side $m$ gives a grounded square of side at least $m-\omega(1/n)-1/n$ on sufficiently fine stages. Thus the largest boundary side has liminf at least $i(\Gamma)$. On each stage select the lattice square and apply the [marked polygon lemma](general-noncollapse-goal.md) with an additional displacement tending to zero. The resulting approximations have odd local degree above the **fixed** cutoff $i(\Gamma)/2$ on every sufficiently fine stage. In particular this meets the project's noncollapse goal on that class by an explicit mesh rule.

## A bent base with a computable odd degree

Use orthogonal coordinates in the plane with the same unit scale. Let $A$ be an excursion as above and $L$ another marked arc of the same Jordan curve. Suppose:

1. $L$ lies in $|y|\le\delta<h$. Its starting abscissa is strictly below $\min x(A)-2\delta$ and its ending abscissa strictly above $\max x(A)+2\delta$.
2. The part of $A$ in $h\le y\le h+2\delta$ projects into two closed intervals. Each has width less than $h-\delta$, and their separation exceeds $h+\delta$. The starting and ending points of $A$ lie in the left and right intervals respectively.

Then a square has two vertices on $A$ and two on $L$, with side **at least $h-\delta$**. The marked four-parameter problem has odd mod-two degree. Neither arc needs to be a graph. The projection test involves only the low portion of the excursion. Its height and total variation farther up are unrestricted.

**Proof.** For ordered upper points $p=A(t),q=A(u)$, put $(d_x,d_y)=q-p$. Complete clockwise to the two lower points $p+(d_y,-d_x),q+(d_y,-d_x)$. Match them independently to $L(v)$ and $L(w)$, using the two planar differences as a four-coordinate map on the parameter cube. Deform $L$ to the flat segment that linearly interpolates its two endpoint abscissas at height zero. Straight interpolation stays in the strip and fixes those endpoint abscissas.

At any zero, $|d_y|\le2\delta$ by the lower heights, and $d_x\ge h-\delta>0$ by the upper and lower heights. If $t$ or $u$ is on a parameter boundary, that upper point is at height $h$, so $d_x\in[h-\delta,h+\delta]$ and both upper points are in the low band. Their horizontal separation can be neither this large in the same interval nor this small in different intervals. If $v$ or $w$ is on the boundary, its abscissa cannot be an upper abscissa plus $d_y$, by hypothesis 1. Thus the entire homotopy has zero-free cube boundary. Every zero gives a genuine positive-size square with the claimed side bound.

At the flat end the lower horizontal coordinates solve uniquely for $v,w$ near each zero. The remaining vertical equations are $(p_y-d_x,q_y-d_x)$, an invertible linear transform of $(d_x-p_y,d_y)=q-Tp$. Eliminating the two horizontal equations preserves mod-two degree by the product rule. Equivalently perturb the two remaining equations transversely near their compact zero set, then solve the horizontal equations there with $v,w$ strictly interior. The shear calculation gives degree one. This is a degree for the indicated choice of contact arcs. Other squares of the whole curve may have any parity.

## What the exact checks measure

For integral column heights, the optimum contained side is $\max_{I}\min(|I|,\min_{j\in I} H_j)$ over nonempty consecutive column intervals $I$. Indeed any open square spans columns with height at least its side and is bounded in width by their combined width. Conversely such an interval supplies a square of the indicated size. At an integer abscissa the roof occupies exactly the height interval between its two neighboring column heights, using outside height zero at the two end walls. These facts give separate calculations of the contained side and the largest boundary square on the base. Across all **4,096** six-column roofs with heights from 1 through 4, the boundary square is strictly larger in **2,489**. For heights $(1,1,1,2,1,2)$ the contained optimum is 1 and a boundary square of side 2 starts at $x=3$. These finite checks exercise the code. The all-mesh result follows from the proof.

For an independently checked bent-base example, take the upper arc through $(1,2),(2,5),(7,11/2),(9,2)$ and the lower arc through $(-2,0),(3,0),(4,1/8),(15/4,1/5),(8,1/12),(11,-1/8),(13,0)$. The lower arc turns backward. Joining corresponding ends by straight segments gives a Jordan curve: the two connectors stay outside the upper span except at its endpoints and above the adjacent lower horizontal/sloping ends except at attachment. The low-band check gives a floor of $7/4$. The independent four-edge census checks all $3^2 6^2=324$ clockwise marked charts including repeated edges and singular systems. It finds exactly one regular square, with squared side $176409046145/6998326336>25$.

## Why scanning the longest straight grid edge fails

Consider the physical trapezoid with vertices $(0,0),(12,12),(4,14),(-2,8)$. At mesh $1/n$, trace each side by repeating first a horizontal run and then a vertical run. In grid units the successive run pairs and repetition counts are $(1,1)$ repeated $12n$ times, $(-4,1)$ repeated $2n$ times, $(-1,-1)$ repeated $6n$ times, and $(1,-4)$ repeated $2n$ times. Each chain is monotone in both physical coordinates. Consecutive chains are separated by the horizontal line through their common extremum or by the vertical line through it, with just their common endpoint shared (the final horizontal/vertical runs give the same conclusion at the origin). The opposite side chains lie in disjoint horizontal-coordinate ranges or opposite height bands in the coordinates below. The polygon is simple. Alternatively, at the origin the last chain has $x\le0$ and the first $x\ge0$, with only the final vertical run on $x=0$ and only the first point of the first chain there. The longest maximal straight run is $4/n$. A selector confined to one literal base edge therefore loses its size bound under refinement.

In coordinates $u=(x+y)/2,\ v=(y-x)/2$ lengths are divided by $\sqrt2$. The first staircase is the base, with $0\le u\le12$ and $|v|\le1/(2n)$. On the left slant $0\le u-3v/5\le4/(5n)$. On the right slant $-4/(5n)\le u+3v/5-12\le0$. Both are monotone in $v$; the top has $v\ge5$ and $3\le u\le9$. For $n\ge4$ choose the upper excursion cut at $h=5-2/n$ and $\delta=1/(2n)$. Its abscissas are in $[2,10]$. Its low band stays on the slants and projects into
$$
[3-6/(5n),3+1/(5n)]\quad\text{and}\quad
[9-1/(5n),9+6/(5n)] .
$$
These satisfy the strict inequalities for every $n\ge4$. Consequently the **physical** square side is at least
$$
\sqrt2\,\big(5-5/(2n)\big).
$$
It tends to $5\sqrt2$ even though every individual grid edge tends to zero. The same gate proof applies to any folded replacement arcs obeying the bounds. On the unperturbed staircase there are also directly visible squares at matched grid points. The gate calculation supplies a marked odd degree without selecting any one of them or requiring those coincidences.

**Remaining selection problem.** A flat terrain supplies its own baseline and a grounded interior square. The bent-base test expands the search to a long arc lying in a thin strip, paired with a separated upper excursion. Arbitrary Jordan curves can have overhangs, obstructions beneath an interior square, and no straight boundary segment. We have no general estimate that selects a direction and a pair of such arcs at positive scale. For an extension beyond these classes the question is whether failure of the two projection intervals forces a different positive-size crossing that can be used with [v315](v315-cutoff-payment.md).

[Exact verifier](../code/roof_shear_v316.py) · [Tests](../tests/test_roof_shear_v316.py) · [Rational results](../results/v316/roof-shear.json)
