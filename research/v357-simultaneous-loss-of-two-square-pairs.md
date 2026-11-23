# v357. When the tilted squares and their replacements both disappear

Use a $W$ by $H$ box with roof cuts $[p,q]$ and $[A,B]$ of floor heights $g$ and $h$, respectively. Assume $0<p<q<A<B<W$, $0<h,g<H$ and $W-B\le A$, as in [v356](v356-roof-switch-and-far-wall.md).

**Result.** The robust disk bound holds throughout
$$
 A<h+q<\min(B,H),\qquad g<B .
$$
This includes the simultaneous event $g=h+q=B-p$ left open by v356, as well as an open neighborhood of that event within the stated geometry. Robustness means that for every strict cutoff below $\sqrt2 r$, where $r$ is the original inradius, all sufficiently close registered Jordan redraws have a square above that cutoff.

## What actually disappears

Fix $g_*=h+q=B-p$, with $A<g_*<\min(B,H)$. Near the old tilted square, bottom and left anchors $(h,0),(0,v)$ give $t=h+v$ and the third point $(v,t)$. At the right foot of the second notch it satisfies $v\le q$, $t\ge g$, with at least one equality. For $g<g_*$ sufficiently close, this gives the two regular squares with $t=g$ and $t=g_*$. For $g=g_*$ the solutions coincide; above it there are none locally.

Near the old axis square the two bottom points and the upper right contact on $x=B$ force side $s$ and upper left contact $(B-s,s)$. At the *left* foot of the second notch both inequalities are reversed in their effect on feasibility: $B-s\ge p$, $s\ge g$, with at least one equality. Again $s=g$ and $s=g_*$ give distinct regular solutions below the threshold, one at it, and none above. All other edges are strictly clear for small changes. Regularity off the collision follows by eliminating the bottom coordinates and imposing the active horizontal or vertical edge equation, each with nonzero derivative. The two collision neighborhoods are disjoint because one square is tilted and the other axis aligned. Raising $g$ demonstrates zero local degree at both collisions.

## Where to find the survivor

Use bottom and left anchors $(u,0),(0,v)$ and their completions $(t,u),(v,t)$ with $t=u+v$. Let $T=\min(H,h+A)$.

1. **If $T<B$, use $t=T$, $u=h$, $v=T-h$.** The second point lies in the first notch floor since $A<h+q<T<B$. The third lies on the intervening roof ($q<v\le A$) or the first cut's left wall ($v=A$ and $t>h$). The floor fixes $u=h$, while the remaining local conditions are $t\le H$, $v\le A$ with at least one equality. They select the formula uniquely. For $H\ne h+A$ the two active affine equations are independent. At the roof corner $H=h+A$, changing $H$ slightly computes the same absolute degree one, since exactly one regular solution remains nearby. This square uses boundary pieces that the motion of $g$ leaves untouched.
2. **If $T\ge B$ and $H>B$, use $t=B$, $v=q$, $u=B-q$.** Both completed points are strictly on the right walls of the cuts: $h<B-q<H$ and $g<B<H$. The independent equations $u+v=B$ and $v=q$ give absolute degree one. In particular varying the second floor up or down while $g<B$ leaves the square fixed.
3. **If $T\ge B$ and $H=B$, follow the interval $t=B$, $q\le v\le B-h$, $u=B-v$.** Here $q<B-h\le A$. All squares retain height $B$. At its ends the wall and roof contacts can reach notch corners. The following local count includes $B-h=A$.

In the third case the marked arc for $(t,u)$ follows the first cut's right wall and floor; it requires $t\le B$, $u\ge h$, with $t<B$ forcing $u=h$. The marked arc for $(v,t)$ follows the intervening roof between $q$ and $A$, including the adjacent walls: $q\le v\le A$, $t\le H$, with $t<H$ forcing $v=q$ or $A$. These conditions are valid near the specified traces, since their endpoints stay clear of the nonincident edges. At $H=B$, any $t<B$ would equal $h+q$ or $h+A$. The former is separated below $B$ and the latter is at least $B$. Thus this interval is isolated. Increase $H$ slightly, raising $h$ slightly too if $h+A=B$. In the same neighborhood $t\le B<H'$ forces $v=q$ or $A$. The latter would give $t\ge h'+A>B$. At $v=q$, the floor solution $t=h'+q$ remains separated below the isolator, so only the regular solution $t=B$ remains. Homotopy invariance gives absolute degree one to the interval. The minimum squared side is attained at the value of $v$ closest to $B/2$.

In each case shrinking the isolator protects every strict cutoff below the minimum reference side. The marked residual changes uniformly under a registered redraw, so its nonzero degree forces a square there. For the disk comparison every selected square has $u\ge h$ and $t\ge A$. The [v355 comparison](v355-tall-strip-and-grounded-disks.md) used only these inequalities, $v\ge0$ and $W-B\le A$: filling the second cut bounds the disk by the wider comparison L, and the anchor inequalities give $u^2+v^2\ge2r^2$. Hence it applies to the entire selected family.

## A fifth square with finite clearance

Take $W=15,H=12,p=2,q=5,A=8,B=11,h=4$. At $g=9$ both pairs collide. The square $(6,0),(11,6),(5,11),(0,5)$ has side $\sqrt{61}$ and remains stationary for every $0<g<11$. At the collision the reference disk target is at most $4\sqrt2$ by v351 ($A=2h$).

At $g=9$, every registered Euclidean redraw of error $e\le1/32$ retains a square of side greater than $7.4$. Parametrize the four marked points by $(x,0),(11,y),(5,z),(0,w)$ and compare the second and third with the completed anchors. The four residual entries are $11-x-w,y-x,5-w,z-x-w$. Each changes by at most $4e$ in the straight residual homotopy. Thus at a zero, $|w-5|\le4e$, $|x-6|\le8e$, $|y-6|\le12e$, $|z-11|\le8e$. The parameter box of radius $1/2$ stays inside the four edges and has zero-free boundary. Degree persists there. The distance between the two actual anchors is at least $\sqrt{61}-(\sqrt{80}+2)e>7.4$.

[Exact selector](../code/double_loss_v357.py) · [Component tests](../tests/test_double_loss_v357.py). The independent tracer checks minimum side and degree on 45 parameter rows, including the two height ties. It also checks both merging pairs just below and at their common threshold. The result leaves other two-notch regimes open; in particular it uses strict separation $h+q<B$ and $g<B$ to keep the right-wall contacts available.
