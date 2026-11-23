# v349. The disk bound for every outward-tilted L

Use vertices $(0,0),(W,0),(W,b),(a,b),(c,H),(0,H)$ with $0<a<c<W$ and $0<b<H$. The upper arm widens from $a$ to $c$. Its right boundary is $f(y)=a+k(y-b)$ with $k=(c-a)/(H-b)>0$.

**Theorem.** If the largest interior disk has radius $r$, every strict side cutoff below $\sqrt2 r$ survives all sufficiently small registered Jordan redraws. The conclusion covers every positive slope and all real dimensions above. [v348](v348-thick-lower-arms.md) proves the result when $b\ge W$. Assume henceforth $b<W$.

## Short upper arm

If $H\le a$, the top and bottom extend from the left wall to lengths $c$ and $W$, both strictly greater than $H$. The side-$H$ squares slide horizontally from left position zero to $c-H$. Apply the [v340 tall-wall argument](v339-v340-end-wall-selection.md) with axes exchanged. Here are the local hypotheses explicitly. The whole polygon is in the height-$H$ strip. For positive slide positions the two left contacts are interior to its parallel supports. At position zero the two right contacts are interior to those supports and the left closing edge is straight. A square near the interval with its left pair on those supports must also have its right pair on them, by the height constraint, and thus belongs to the interval. Near zero the right pair gives the same isolation. In v340's taper the two supports converge in the slide direction, excluding a square with its left pair on the supports. Only the zero end contributes, with the same three transverse solutions at the unchanged straight closing edge. Thus the family has odd degree irrespective of the slope at its far end. Its floor $H$ clears the disk target since $2r\le H$.

## Anchor squares and span cutoffs

Suppose $H>a$. Mark adjacent contacts $(u,0),(0,v)$ on the bottom and left. The other two must be $(t,u),(v,t)$ where $t=u+v$. Use the right wall or ledge at the first of these contacts and the top or slope at the second. Near the selected family we have $t>\max(a,b)$ and
$$
t\le m=\min(W,H),\quad u\le b,\quad v\le f(t),\qquad
(t=W\text{ or }u=b),\quad(t=H\text{ or }v=f(t)).
$$
On the slope $v=f(t)$; on the top the inequality uses $f(H)=c$. Thus $g(t)=t-b-f(t)\le0$, with equality forced for $t<m$. Since $g(b)=-a$, this selects
$$
t=\begin{cases}
\min(m,b+a/(1-k)) & k<1,\\
m & k\ge1 .
\end{cases}
$$
Conversely, at this $t$ take $u=b$ if $t<W$, or $v=f(t)$ if $t<H$. When both apply these equations agree. When $t=W=H$, take any $t-c\le u\le b$. In every case $u,v>0$: at a free root $v=a/(1-k)$, when $u=b$ use $t>b$, and on the width cap use $0<f(W)\le c<W$. Also $t>\max(a,b)$ since both outer spans exceed $a,b$ and a free root exceeds them. All contacts therefore belong to the indicated open arcs, including their internal corners. At a double cap $t-c>0$.

**Persistence including ties.** In labeled parameter space this compact point or interval is isolated: the anchor equations just listed give every nearby zero. There is positive clearance from the ends of the marked arcs, so small changes in dimensions register these neighborhoods edge by edge. Perturb the lengths slightly to break cap and root ties. The formula then gives one interior solution in the neighborhood (every limit lies in the original family). Its two anchor equations have independent normals: they are $u=b$ and $v=f(u+v)$ at a free root (determinant $1-k>0$), or one cap $u+v=m$ paired with $u=b$ or $v=f(u+v)$ (absolute determinant one). The other two contact parameters are recovered uniquely at nonzero edge speed. Thus the local square residual has odd degree, including at the original ties by homotopy invariance. Uniform registered redraws preserve it, giving every strict cutoff below the minimum side on the family, as in [v345](v344-v345-convex-and-L-disk-bounds.md).

If $t=m$, every such square has side $\sqrt{u^2+v^2}\ge m/\sqrt2\ge\sqrt2 r$. This settles the cap regime.

## When the elbow and upper square are both needed

It remains that $t=b+v<m$, $v=a/(1-k)$ and $0<k<1$. The elbow square has robust side $E=\sqrt{b^2+v^2}$. There is also a regular upper square from [v346](v346-tilted-inner-wall.md): writing $D=1+k+k^2$, its side is $U=c\sqrt{1+k^2}/D$. Indeed $H-b>v>a(1+k)$, so $c=a+k(H-b)>aD$. This is exactly the strict condition that both of its sloping contacts lie above the elbow.

Consider any interior disk of radius $R$ with center height $y$.

- If the disk lies entirely above the ledge, it fits in the triangle cut out by the top, left and extended sloping lines. The upper square lies on that triangle's boundary. The zero-depth [convex lemma](v307-convex-shell-square-gap.md) gives $U\ge\sqrt2 R$.
- If the disk goes below the ledge but its center is above it, then $b<y<b+R$. Its horizontal diameter satisfies $2R\le f(y)<a+kR$. Hence $R<a/(2-k)$ and $E\ge v=a/(1-k)>\sqrt2 a/(2-k)$.
- If the center is at or below the ledge, inscribe an axis-aligned square of side $\sqrt2 R$ in the disk. If its top is at or below $b$, its side is at most $b$. Otherwise it crosses height $b$, and horizontal slices immediately above $b$ bound its side by $f(b+)=a$. Thus $\sqrt2 R\le\max(a,b)\le E$.

Consequently $\max(E,U)\ge\sqrt2 r$, and each of the two certificates is robust. This proves the theorem. In particular the conclusion does not require that the largest disk or any one labeled square be tracked while the shape changes.

The [rational selector](../code/outward_tilt_v349.py) returns both witnesses when needed. [Tests](../tests/test_outward_tilt_v349.py) compare their true component floors and degrees to the independent v341 chart tracer and check disk estimates across a grid. The zero-slope case is v345. Inward tilts narrow the upper arm; their thin-arm cases require a separate argument because the extended triangle and the free-root ordering change.
