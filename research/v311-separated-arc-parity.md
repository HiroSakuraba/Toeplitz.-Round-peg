# v311. Detecting squares on separated arcs through infinitely many refinements

Benjamin John Schulz · November 2026

On every nonempty stage of the [v220 accumulating-notch construction](v220-fixed-rectifiable-collision-link-obstruction.md), and on its limiting Jordan curve, the same four marked arcs carry degree one. Every registered Jordan perturbation of any of these curves by at most $1/512$ in coordinate sup norm has a square meeting all four arcs with side at least $255/256$. For a finite transverse census away from the cuts, this degree equals the sum modulo two of four ordinary mixed-square counts. A two-dimensional crossing calculation gives the degree even when those counts are individually infinite.

## The four-cut identity

Take two disjoint closed parameter intervals $I=[p_0,p_1]$ and $J=[q_0,q_1]$ on the oriented parameter circle, in that cyclic order. In a finite collection of distinct geometric inscribed squares assume no vertex parameter is a cut. For each cut pair $(p_i,q_j)$ call a square mixed if the two complementary open arcs both contain vertices. Let $M_{ij}$ be its mixed count modulo two. Then
$$
\sum_{ij} M_{ij} =
\#\{(S,e):e\text{ joins consecutive contacts of }S
\text{ in source order, one in }I,\text{ one in }J\}\pmod 2 .
$$
Consecutive in source order can include a geometric diagonal on a nongracefully ordered square.

**Proof.** For a single square give each of the four gaps between successive parameter contacts its own basis vector over $\mathbb F_2$, and write $g(t)$ for the vector of the gap containing $t$. The mixed indicator is $1+g(p_i)\cdot g(q_j)$. Its four-term sum is $(g(p_0)+g(p_1))\cdot(g(q_0)+g(q_1))$. The change across either interval is the sum, over its contacts, of the two neighboring gap vectors. Two distinct contacts share a gap precisely when consecutive. Sum over squares. Thus a bit of one certifies a square with side at least $\operatorname{dist}(\gamma(I),\gamma(J))/\sqrt{2}$. If the bit is one on arbitrarily fine registered approximations with the same marked intervals, compactness and injectivity on the limit give a positive-size square.

**Extremal detection theorem.** Let $H$ be the greatest distance between consecutive source contacts of any square in a finite nonempty collection of distinct geometric squares. Some contact pair of length **exactly $H$** occurs an odd number of times. Thus some four-cut test detects a pair at least as far apart as the largest square's side. Conversely no test can recover a consecutive pair longer than $H$.

**Proof.** Form the multigraph of consecutive parameter pairs. If all its pairs of length $H$ are geometric diagonals, each occurs once: its two ends uniquely determine the square. Otherwise a geometric side has length $H$. Every square of that side length visits its vertices in perimeter order, since a non-perimeter four-cycle includes a diagonal of length $\sqrt2 H>H$. Fix the orientation of one such square, modulo quarter turns. Among side-$H$ squares of this orientation choose the center that maximizes a generic linear functional. Choose the functional so that the maximizing center and the maximizing one of the four side-midpoint offsets are both unique. Adding these two maxima gives a midpoint $m$ where exactly **one** side contribution of this length and orientation occurs.

There are two perpendicular segments of length $H$ through $m$ as midpoint in these directions. A side from another orientation cannot equal either. A square supplying either segment as a *consecutive diagonal* must be the unique square having these two segments as its diagonals. In any non-perimeter traversal both diagonals occur, so it contributes equally to the two multiplicities. Their sum is therefore odd. One of the two segments survives cancellation. Injectivity identifies its geometric endpoints with unique source parameters. Choose disjoint intervals small enough to isolate them from all other distinct contacts. The four-cut identity now gives bit one. Finally, each traversal uses a side or a longer diagonal, so $H$ is at least the largest side. This proves the claim.

For a finite graceful collection all consecutive pairs are geometric sides; the midpoint argument applies separately at **each** side length present. The solver computes the odd edges and returns isolating intervals for a longest one. This finite statement controls separation of their center contacts. Their widths may shrink and their positions may change between approximations.

## Computing the bit without the microscopic census

Use the $2$ by $1$ rectangle of v220 with outward bottom notches at $x=2^{-n}$ of depth and width $4^{-(n+1)}$, for $1\le n\le N$ or for every $n\ge1$. Use these four disjoint closed boundary arcs in counterclockwise order:

| Arc | Endpoints and support |
| --- | --- |
| $I$, first marked interval | From $(0,1/8)$ down the left side, through the origin and along the notched bottom to $(1/2,-1/32)$, halfway down the first notch wall |
| $J$, second marked interval | Bottom, $7/8\le x\le13/8$ |
| $K$ | Top, $3/4\le x\le7/4$ |
| $L$ | From $(5/8,1)$ along the top to $(0,1)$, then down the left side to $(0,3/4)$ |

Use the same marked parameters at every stage, e.g. register $I,J,K,L$ and the intervening gaps on fixed parameter intervals by rescaling arclength on each. The cut inside the first notch is defined from $N=1$ onward. On the exact curves, the unit squares with lower left $(x,0)$ for $3/5<x<4/5$ are mixed for **each** cut pair. Thus even finite-stage polygons have infinitely many mixed squares. The degree below gives the combined bit meaning at these singular examples.

Write $P=(x,y)\in I$ and $R$ for counterclockwise quarter turn. For a square ordered $P,Q,B,D$ positively, $B=Q+R(Q-P)$ and $D=P+R(Q-P)$. With $Q=(z,0)\in J$ and $B=(u,1)\in K$, elimination gives
$$
z=1+x,\qquad u=1+x+y,\qquad D=T(P)=(x+y,1+y).
$$
The projected point $T(P)$ must meet $L$. The first endpoint maps to $(1/8,9/8)$ above the rectangle, and the last to $(15/32,31/32)$ inside it. For $P$ on $I$ one has $0\le x\le1/2$ and $-1/16\le y\le1/8$. Consequently $z,u$ stay strictly inside their assigned segments. No $T(P)$ can meet either endpoint of $L$. Its distance in sup norm from $(0,3/4)$ is at least $3/16$. For $(5/8,1)$ use $T_x-T_y+1=x\le1/2$ to get distance at least $1/16$. The images of both endpoints of $I$ miss all of $L$ by at least $1/32$.

More formally take the four-dimensional square residual $(B-Q-R(Q-P),D-P-R(Q-P))$ on the product of the four parameter intervals. These estimates keep the residual nonzero on its entire boundary. Straighten the parameterized first arc to the segment joining its endpoints inside the same coordinate box. The other arcs stay fixed. This is a homotopy of residuals and does not require a new global Jordan curve. On the straight segment $T(P)$ crosses the top transversely once, at $(2/5,1)$, and misses the left side. Elimination of $z,u$ has invertible linear part, so this single transverse intersection gives four-dimensional degree one modulo two. The argument permits arbitrary oscillation of $I$ inside its box with the given endpoints. It does not count any notch squares.

## All four arcs may move

Let $\eta=1/512$. Register the entire new Jordan curve to a reference stage or to the limit, with coordinate sup distance at most $\eta$. Take the straight registered interpolation. At any four-contact zero on the product, comparison with reference coordinates gives
$$
|z-(1+x)|\le4\eta,\qquad |u-(1+x+y)|\le8\eta,\qquad
\|D_{\rm ref}-T(P_{\rm ref})\|_\infty\le8\eta.
$$
Each original residual component has at most four coordinate errors of size $\eta$. Substitute the first estimate in the other equations. The reference $z,u$ predictions are at least $1/8$ from their segment endpoints. The margins above for $I,L$ are at least $1/32,1/16$, respectively. Since $8\eta<1/32$, there are no boundary zeros, and the degree stays one. The vertical separation of the actual $Q,B$ gives side at least $1-2\eta=255/256$. This product-degree argument only needs the stated registration on the four designated arcs.

**Identification with the four-cut test.** Under registration of the whole curve, *every* square with contacts in $I,J$ has exactly one vertex on each of $I,J,K,L$, in that order. Here are the exclusions. Write reference $P=(x,y),Q=(z,0)$, so $v=z-x\ge3/8$. For actual opposite corners the lower of the other vertices has height at most $(y-v)/2+2\eta\le-1/8+2\eta$, below the entire curve. For an adjacent pair the downward completion has a corner at height at most $-v+3\eta$. Only the upward completion remains. Its upper right vertex is within $3\eta$ of $(z+y,v)$, with reference horizontal prediction in $[13/16,7/4]$. It can meet only the top, since the reference curve away from the top is on the sides or at height at most zero. Thus $|v-1|\le4\eta$ and the refined horizontal prediction puts the reference contact strictly inside $K$. The last vertex has reference position within $8\eta$ of $(x+y,1+y)$, with horizontal error alone at most $4\eta$. It can meet only the top or the left. On the top, $|y|\le8\eta$ gives reference abscissa at most $1/2+12\eta<5/8$. On the left its height is at least $15/16-8\eta>3/4$. Both put it inside $L$. The face margins exclude equality at the marked endpoints. The designated contacts are consecutive in source order because the four disjoint arcs are cyclically ordered. Each such square contributes exactly once. For finite transverse configurations with cuts missing all contacts, the four-cut identity now identifies the bit with the degree just computed. For arbitrary targets the degree alone certifies existence.

## Checks and remaining question

A full independent edge census on the v310 triangle returns the four-cut matrix with rows $(1,1),(0,1)$. For two rational perturbations of the **same one-notch polygon**, the complete census gives three squares with rows $(0,1),(1,1)$ and one square with rows $(1,1),(0,1)$. A two-notch perturbation gives three squares with rows $(0,1),(1,1)$. Every one of these matrices has bit one. The [rational replay](../code/separated_arc_parity_v311.py) verifies the margins and compares direct counts to consecutive-contact counts. [Recorded output](../results/v311/separated-arc-parity.json).

The finite detection theorem rules out complete cancellation for any single finite nonempty square census. The notch theorem selects the **same** two separated arcs across infinitely many refinements, including a registered neighborhood of each stage. To carry this strategy further we need a way to select intervals of uniformly positive image separation and retain a nonzero bit as the approximations change. The extremal argument controls the distance on each finite census; persistence at a fixed scale between censuses remains open. No general Jordan-curve noncollapse conclusion follows.

