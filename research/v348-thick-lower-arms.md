# v348. A square survives every upper-arm collision when the lower arm is thick

Use the six-sided L of [v346](v346-tilted-inner-wall.md). Its total width is $W$, the ledge height is $b$, the ledge runs from $(W,b)$ to $(a,b)$, and the upper arm ends at height $H>b$. The missing cases in the [v347 collision proof](v347-square-through-the-tilted-wall-collision.md) are $b\ge W$.

**Result.** For $b>W$ the robust largest-square side is exactly $W$. At $b=W$ a robust lower bound is
$$
\sqrt{d^2+(W-d)^2},\qquad d=\min(W-a,W/2).
$$
Both reach the disk target $\sqrt2 r$: the region contains the bottom $W$-by-$W$ square and lies in a strip of width $W$, so its inradius is exactly $W/2$. Here robust means that every strict cutoff below the stated bound survives sufficiently close registered Jordan redraws. In conjunction with v347 this covers the entire v346 positive-slope collision locus, including a neighborhood of each collision in the height parameter with $W,a,b,k$ fixed.

The upper wall may even be an arbitrary continuous simple arc from $(a,b)$ to $(c,H)$, $0<a,c<W$, whose interior stays in $(0,W)\times(b,H)$. It may turn back in either coordinate. Join its top endpoint horizontally to $(0,H)$ as usual. The following argument explains why this freedom does not cost any square size.

## Strictly taller end walls

The reference axis-aligned squares spanning the strip slide from bottom level $s=0$ to $s=b-W$. Label their corners starting at bottom left. Every contact is on an unchanged boundary piece, including the corner $(W,b)$ and its incident ledge. The varying upper wall stays a positive distance from these contact traces: its compact image is disjoint from both end walls, the base and the upper-right endpoint of the traces. The same holds throughout linear interpolation to the vertical upper wall at $x=a$, with the top endpoint moving along height $H$. The top edge remains above the traces. Thus this labeled sliding component has exactly the local geometry of the [v340 tall-wall proof](v339-v340-end-wall-selection.md), including at both ends. That proof gives odd degree and minimum side $W$. Every square in the reference strip has side at most $W$, giving equality.

## At exact equality the elbow belongs to the contact traces

When $b=W$, label the following interval of squares:
$$
(u,0),\quad(W,u),\quad(W-u,W),\quad(0,W-u),
\qquad 0\le u\le W-a.
$$
At $u=0$ this is the bottom axis-aligned square; at the other end the third vertex is precisely the elbow $(a,W)$. The shape of the outgoing wall can now affect the local edge chart, so we verify isolation through a deformation.

Parametrize each of the six boundary pieces on a fixed interval. Interpolate the upper arc pointwise to the straight segment from $(a,W)$ to $(a,H)$, while moving its top endpoint and the adjoining top edge continuously. Every interior point of the interpolated arc remains strictly above $W$ and strictly between the end walls. The degree homotopy uses only continuous parametrizations and is valid even if an intermediate arc crosses itself.

The displayed family uses the same four parameter traces throughout (choose linear parameters on the fixed straight pieces). Away from the elbow endpoint, all four boundary germs are unchanged from the v340 plateau model, reflected horizontally. This includes $u=0$, where the v340 corner-chart check isolates the family. Near the elbow endpoint the bottom, right and left contacts are all interior to straight edges. If the first and fourth contacts are $(x,0)$ and $(0,y)$, completing the square makes the second contact $(x+y,x)$. The right wall forces $x+y=W$, so the third contact must be $(y,W)$. On the marked elbow arc the points at height $W$ belong exactly to the ledge, including its endpoint. This gives just the displayed family locally, uniformly over the homotopy. In particular no tangent or nonzero slope at the elbow is required.

For precision, take the square residual on the product of four labeled parameter intervals as in v345, using circle charts across polygon vertices. The preceding local descriptions give an open neighborhood of the compact family times the homotopy interval in which every zero belongs to that product. Compactness gives one isolating neighborhood in parameter space with zero-free boundary for the entire homotopy. At the straight-wall end v340 gives odd local degree. Homotopy invariance transfers it to the given arc. Shrinking the neighborhood and using uniform continuity of the residual gives the claimed survival under arbitrary registered redraws above every strict cutoff below the minimum side. Minimizing $u^2+(W-u)^2$ gives the formula.

For $W=b=6$ and $a=4$, the surviving family has minimum side $\sqrt{20}$, exceeding the disk target $3\sqrt2$. For $a=2$ its minimum is exactly $3\sqrt2$. This comparison concerns the selected family; other squares may coexist.

## Exact checks and scope

The [tests](../tests/test_thick_lower_arm_v348.py) trace the labeled component on straight and folded upper walls, including a wall that reverses height. They verify the minimum squared side and degree on both sides of the minimum's switch at $a=W/2$, at the switch itself, and at lower-arm heights equal to and greater than the width. The theorem proves persistence for continuous arcs. The remaining selection problem for tilted L shapes lies in the thin lower-arm regime away from the collision neighborhoods already covered.
