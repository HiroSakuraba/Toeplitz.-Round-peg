# v356. A roof contact continues a square; a far notch wall can destroy two

Use the box $[0,W]\times[0,H]$ with two disjoint roof cuts $[p,q]$ and $[A,B]$ at floor heights $g,h$ respectively, with $0<p<q<A<B<W$ and $0<g,h<H$. Assume $W-B\le A$. Write $r$ for the inradius. A robust floor means that every strictly smaller side cutoff survives sufficiently close registered Jordan redraws.

## Following a square through the outer roof

Suppose $\max(A,g)<h+p<B$ and $H>A$. The following formula gives one robust square for all such heights, including $H=h+p$:
$$
 t=\min(H,h+p),\quad v=t-h,\qquad
 (h,0),\quad(t,h),\quad(v,t),\quad(0,v).
$$
The contact $(t,h)$ is strictly on the floor $(A,B)$. Below the critical height, $(v,t)$ is on the roof left of the cut; above it, it is on the wall $x=p$. At equality it is exactly the top corner $(p,H)$.

Indeed, with interior bottom and left anchors $(u,0),(0,v)$ the completions are $(t,u),(v,t)$ with $t=u+v$. The first floor forces $u=h$. Near the top-left corner of the cut the third contact requires $v\le p$, $t\le H$, with at least one equality. This gives the displayed solution uniquely and isolates it at the tie. Off the tie, the second active equation is either $v=p$ or $u+v=H$, independent of $u=h$. The remaining two contact parameters are recovered at nonzero speed. The regular solution therefore has degree of absolute value one. Changing $H$ slightly computes the same local degree at the tie: just one nearby regular solution lies in its isolator. Uniform residual continuity then transfers every strict cutoff below the reference side to arbitrary registered redraws.

For $t=H$ the square spans the box height and its side is at least $H/\sqrt2\ge\sqrt2r$. For $t=h+p$ the anchors obey $t>A$ and $u=h$; the size comparison of [v355](v355-tall-strip-and-grounded-disks.md) gives the disk bound at every strip height. This removes the strict outer-roof clearance condition on this contact pattern.

## A far-wall contact with a different outcome

Now suppose $A<h+q<\min(B,H)$ and vary $g$ near $g_*=h+q$. Keep the bottom and left anchors and the first-floor contact. Near $v=q$, the second floor or its right wall requires $v\le q$, $t\ge g$, with either $t=g$ or $v=q$. Since $u=h$, the candidates are:
- $v=g-h$, on the second floor;
- $v=q$, on its right wall.

For $g<g_*$ sufficiently close, both candidates lie inside their respective edges and are regular. At $g=g_*$ they coincide. Above the threshold neither is feasible in the marked neighborhood: $t=h+v\le h+q<g$. The collision is isolated with degree zero, as can be seen by raising $g$ slightly. Thus continuing the old marked square across this event is impossible.

**Replacement criterion.** If $h<g<H$, $g\ge A$ and $p<B-g<q$, use the square
$$
 (B-g,0),\quad(B,0),\quad(B,g),\quad(B-g,g).
$$
All contacts are interior, including the two separate bottom contacts. They force the first side horizontal, the upper-left floor fixes its length $g$, and the upper-right wall fixes its right abscissa $B$. These are nonsingular affine equations on the four edge parameters. Thus this square is robust with floor $g$. The largest disk fits in the region obtained by filling the second cut. By v351 that disk size is attained in the comparison L with tower width $A$ and strip height $h$. Any axis-aligned square inside the L has side at most $\max(A,h)$; a larger square would extend above the strip into the narrow section. Inscribing an axis-aligned square in the disk gives $\sqrt2r\le\max(A,h)\le g$. If $p<B-g_*<q$, these strict contact conditions persist throughout an open interval across the collision.

**Finite-error example.** Take $W=15,H=14$, cuts $[2,5]$ and $[8,13]$, and $h=4$. At $g=9$ the two tilted squares merge, while the replacement has side 9. For a Euclidean registered error at most $1/32$, the redraw still inscribes a square of side at least $9-6/32$.

To see the finite margin, parametrize the four reference contact edges by $(x,0),(y,0),(13,z),(w,9)$ near $(x,y,z,w)=(4,13,9,4)$. Use the residual comparing the last two corners with the positive quarter-turn completions from the first two. Its reference entries are $13-y,z-y+x,w-x,9-y+x$. A redraw of error $e$ changes each entry by at most $4e$, also in the straight residual homotopy. A zero therefore has $|y-13|\le4e$, $|x-4|\le8e$, $|w-4|\le12e$, $|z-9|\le8e$. A coordinate box of radius $1/2$ stays inside all four edges, and $12/32<1/2$ keeps the homotopy zero-free on its boundary. At a redrawn zero the first side has length at least $y-x-2e\ge9-6e$. The disk target of the reference is at most 8.

## When the replacement also hits the far corner

The replacement itself continues through $B-g=q$. Set $s=\max(g,B-q)$ and suppose $s\ge A$, $h<s<H$, $B-s>p$. Its lower corners are $(B-s,0),(B,0)$ and its side is $s$. The two bottom contacts and the right wall force the remaining upper point to be $(B-s,s)$. At the second cut's right foot the constraints are $B-s\le q$, $s\ge g$, with at least one equality. They give exactly this maximum formula. If $s>g$ both upper points lie on different vertical walls, whose abscissas fix both bottom parameters. At equality the isolated corner transfers degree from either regular side as $g$ varies. The same disk comparison holds since $s\ge\max(A,h)$.

For example changing the first cut to $[8,14]$ in a $16$ by $15$ box makes $g=h+q=B-q=9$ a simultaneous collision. The tilted pair disappears while the axis square switches contacts and retains odd degree. At the other end $B-g=p$, the adjacent wall and floor both have $x+y\ge p+g$ near their shared foot. For $g>B-p$ neither meets $x+y=B$ there: that replacement can genuinely disappear.

[Exact selector](../code/notch_contact_switch_v356.py) · [Independent degree and size checks](../tests/test_notch_contact_switch_v356.py). Thirteen geometries cover the roof transition, a far-wall collision, and a simultaneous corner switch. Finding a replacement when the axis square hits the left foot at the same time as the tilted pair disappears is still open here.
