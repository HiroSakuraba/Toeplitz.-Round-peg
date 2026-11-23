# v355. A taller strip forces a longer square

[v354](v354-when-the-second-cut-covers-the-endpoint.md) assumed $2h\le A$ to bound the disk diameter by the first tower width $A$. We remove this condition throughout that theorem, allow the outer height to vary, and include a third contact pattern. The disk really can exceed the old bound. The additional square size comes from its orientation.

## The two-notch theorem

In a box of width $W$ and height $H$, make disjoint interior roof cuts $[p,q]$ and $[A,B]$, in that order, with positive floor heights $g,h<H$. Assume $W-B\le A$. Put $t=\max(A,g,h+p)$ and start with the closed interval $I=[h,t-p]$. Intersect it with $\{h\}$ if $t>A$, and with $\{t-p\}$ if $t>g$. By the maximum formula this leaves a nonempty interval, possibly a point. Suppose
$$
 t<\min(H,B),\qquad t-\min I<q .
$$
Then the robust side reaches $\sqrt2$ times the inradius $r$. An explicit squared floor is the minimum of $u^2+(t-u)^2$ on $I$, attained at the point closest to $t/2$. The hypotheses hold for all the old endpoint parameters even after deleting $2h\le A$: there $W=B+d$, $H=B$, $p+h<A$, $g<B$, and $q\ge B-h$. If $t=A>g$ the third contact has abscissa $p<q$; otherwise $t-h<B-h\le q$.

**Contacts and survival.** Mark the bottom, first cut's left wall/floor at $A$, second cut's left wall/floor at $p$, and the outer left edge in geometric order. Two adjacent anchors $(u,0),(0,v)$ force the other contacts to be $(s,u),(v,s)$ with $s=u+v$. Locally their constraints are
$$
 s\ge A,\quad u\ge h,\quad (s=A\text{ or }u=h);
 \qquad s\ge g,\quad v\ge p,\quad(s=g\text{ or }v=p).
$$
Thus $s\ge\max(A,g,h+p)$; a strictly larger $s$ would force $u=h,v=p$ and hence $s=h+p$, impossible. Conversely the solutions at $s=t$ are exactly $u\in I$. The strict hypotheses give clearance from the far endpoints of the marked edges, so this is an isolated compact labeled family. At a unique maximum of $A,g,h+p$ the two active equations are respectively $(s=A,v=p)$, $(s=g,u=h)$ or $(u=h,v=p)$. Each pair is independent in $u,v$, and recovers the remaining two edge parameters at nonzero speed. The third case uses the first floor and the second wall.

At ties change the dimensions generically by arbitrarily small amounts, registering the edges. The interval formula puts every limit of the unique regular solutions inside the original family, and all solutions in a sufficiently small isolator satisfy the same contact equations. Therefore the family has odd local degree, also at a three-way tie. Uniform registered redraws preserve that degree. Shrinking the isolator establishes every strict cutoff below the stated floor.

**The size comparison.** Fill the cut at $[p,q]$. By [v351](v351-one-rectangular-roof-notch.md) its wider tower has width $A$, so the original radius is at most the radius of the comparison L with parameters $W,H,A,h$. If $h\le A/2$, use $r\le A/2$ and $t\ge A$. If $h\ge A$, an axis-aligned square inside the comparison L has side at most $h$: a larger square would extend above the floor into a section of width only $A$. A disk contains an axis-aligned square of side $\sqrt2 r$, and our boundary square has side at least $u\ge h$.

The missing range is $A/2<h<A$. By the [v345 inradius formula](v344-v345-convex-and-L-disk-bounds.md),
$$
r\le R=A+h-\sqrt{2Ah}.
$$
Every member of our family has $u\ge h$, $v\ge0$ and $u+v\ge A$. Minimizing over those three inequalities gives $u^2+v^2\ge h^2+(A-h)^2$: for fixed $u$ take $v=(A-u)_+$, then minimize over $u\ge h>A/2$. Write $z=\sqrt{2h/A}$. The exact difference factors as
$$
 h^2+(A-h)^2-2R^2=A^2(z-1)^2(2z-1)\ge0 .
$$
This proves the full comparison. Increasing the strip height moves the anchor away from the middle of its span and increases the minimum possible side.

**Example.** In a $10.5$ by $10$ box use cuts $[1,5]$ at height 9 and $[8,10]$ at height 6. Here $A=8$ but the region fits a disk of radius $4.2$. In fact its radius is $14-\sqrt{96}$: the comparison disk grounded at the left corner lies below the second floor since its diameter is less than 9. The selected square has $(u,v)=(6,3)$ and side $\sqrt{45}$, above the disk target of about $5.94$.

## Exact disk feasibility for any finite list of roof cuts

For a proposed rational radius $s$, translate a candidate disk vertically down until its center is at height $s$. The region contains the vertical segment below each of its points, so this preserves containment. The box requires $2s\le\min(W,H)$ and the center abscissa $x\in[s,W-s]$. A cut $[a,b]$ of floor height $k\ge2s$ imposes no further restriction, including at tangency. For $k<2s$, the supremum of the half-width of the disk above $k$ is
$$
\ell=\sqrt{s^2-(k-s)_+^2}.
$$
Avoiding the interior of the bite is equivalent to $x\le a-\ell$ or $x\ge b+\ell$. Sort active cuts by position. The center must be right of some initial group and left of the remaining cuts. For each split, compare every lower endpoint against every upper endpoint, including the box bounds. Each comparison has the form $\sqrt{U}+\sqrt{V}\le D$ and is decided by $D\ge0$, $D^2-U-V\ge0$ and $(D^2-U-V)^2\ge4UV$. Thus this test handles exact tangencies without approximating a radical. It determines disk feasibility, not square survival outside the contact theorem.

[Exact implementation](../code/tall_strip_v355.py) · [Tests](../tests/test_tall_strip_v355.py). Independent polygon traces verify the degree and squared floor on 72 shapes. Disk tests include two-cut tangencies, the larger-than-$A/2$ example and a separate single-notch inradius oracle. The remaining geometric restriction is contact clearance: when the proposed square reaches a far wall or the outer roof, additional boundary edges enter the degree calculation.
