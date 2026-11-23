# v361. Finding three traversals on an unmarked boundary

The input is a rational polygon, a relative strip width and an optional minimum square side. The search chooses a frame, three stretches of boundary and two probes. A hit certifies an inscribed square by [v360](v360-three-windows-free-fourth-contact.md). The computed lower bound can also be used with any strictly smaller cutoff on sufficiently close registered Jordan redrawings. Every polygon used below is given to the search without its contact labels.

## Which frames and traversals are covered?

Project onto a rational vector $u$ and its counterclockwise quarter turn $Ju$. Subtract the projected left and bottom levels and divide both coordinates by the same positive width $L$. This is an orientation-preserving similarity. In the normalized plane the three lines are $y=0,x=1,x=0$, in bottom, right, left order. The automatic candidate list uses every oriented polygon-edge direction and the levels of actual parallel and perpendicular edges on these lines. The frame edges only propose lines. Any part of the boundary can supply a traversal. The separate frame-scanning entry point accepts other rational frames.

Fix a relative half-width $e$ and tangential padding $p$ with $0\le e<p<1/2$. On each of the three lines restrict the tangential coordinate to $[p,1-p]$. The resulting closed strip rectangles are disjoint. Clip every edge to each rectangle, join adjacent parameter intervals including across the polygon's starting vertex, and find the two tangential extrema on each connected run. The subarc between the extrema realizes the full range, including when it folds. Thus using only this one subarc per run loses no candidate within these rectangles: every shorter traversal's endpoint span is contained in the extrema span. The certificate records fractional positions on the original edges. The padding and candidate-frame rules restrict the search.

Write $[B_-,B_+],[R_-,R_+],[W_-,W_+]$ for the three spans. Both probes must lie strictly between
$$
\max(W_-,1-B_++4e,1-R_++8e)
\quad\text{and}\quad
\min(W_+,1-B_--4e,1-R_--8e).
$$
On the line $y=1$, compute exactly the complement of the closed max-norm $8e$ neighborhood of the boundary. For each edge, clip to $[1-8e,1+8e]$ in height, project the remaining segment to the horizontal axis and expand that interval by $8e$ on each side. The union is precisely the forbidden set. Each open complementary interval has one inside/outside sign. After intersecting with the admissible probe interval, a rational midpoint realizes that sign. Opposite signs yield Euclidean clearance greater than $8e$, as required in v360. The physical squared side floor is $L^2(1-2e)^2/(2\|u\|^2)$. This procedure exhausts these clipped traversals and max-norm-clear probes for each tested frame.

## Exact results

On the three contact-exchange stages of v360, both straight and folded, the automatic search at $e=1/2816$ selects the same width-12 frame. Its bottom reference is the original vertical line $x=15$. Its other two reference lines are $y=0$ and $y=12$. In normalized coordinates the outside probe is at $1791/2048$ and the inside probe is at $17/24$. All six squared side floors are $17816841/247808$, giving **side greater than 8.47**. The earlier hand-marked width-11 certificate gives side greater than 7.77 and follows a different family of contacts. These are lower bounds, not maximum-side calculations.

In the old width-11 frame, blind clipping also rediscovers successful arcs on the folded example after reversing traversal, subdividing every edge and moving the parameter seam. The tests reinsert the chosen parameter endpoints and independently run v360 on those arcs. A rational rotation, scaling and translation also passes the automatic search. The triangle with vertices $(0,0),(4,0),(2,3)$ supplies a control: the finite edge-frame search misses, although the horizontal-base square of side $12/7$ is explicit.

## A limit even unrestricted frame search cannot fix

**A round circle admits no v360 three-traversal certificate.** In the original three-parameter proof the residual has no zero on the boundary of the trimmed box, and a successful certificate gives a zero in its interior. For the circle centered at the origin, the noncollapsed solutions in the full three-parameter torus have the form $(A,JA,-JA)$ in geometric bottom/right/left roles, as $A$ runs once around the circle. Indeed opposite vertices of any inscribed square are antipodal. This solution set is connected. Disjoint marked arcs exclude all collapsed solutions in the box. A connected family meeting its interior and exterior must meet its boundary, yet the family reaches the exterior as soon as $A$ leaves the proper bottom arc. This contradiction holds for any orientation, widths or choice of traversals allowed by v360.

This example identifies a genuine reason to pursue another certificate. The [relabeling-path criterion](v362-relabeling-path-noncollapse.md) uses the motion of a whole family to force separated parameters.

[Search implementation](../code/blind_windows_v361.py) · [Exact replays](../tests/test_blind_windows_v361.py)
