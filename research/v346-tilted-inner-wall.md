# v346. Tilting the inner wall of an L can move the large square into the upper arm

Benjamin John Schulz · November 2026

The square selected at the elbow in [v345](v344-v345-convex-and-L-disk-bounds.md) can survive a tilt and nevertheless become arbitrarily small compared with the largest interior disk. A different square, with two corners on the tilted wall, proves the robust disk bound for an explicit family exhibiting this failure.

Use vertices $(0,0),(W,0),(W,b),(a,b),(c,H),(0,H)$, where $0<a,c<W$ and $0<b<H$. Above height $b$ the right boundary is $f(y)=a+k(y-b)$ with $k=(c-a)/(H-b)$. All guarantees of side $s$ below have the v343 meaning: every strict cutoff below $s$ survives sufficiently small registered Jordan redraws.

## What happens at the elbow

For $k<1$ put $v=a/(1-k)$ and $t=b+v$. When $a<t<W$ and $t<H$, the points $(b,0),(t,b),(v,t),(0,v)$ form a regular inscribed square provided $v<H$ (automatic here). Its side is $S_E=\sqrt{b^2+v^2}$. Indeed the bottom and left anchors $(u,0),(0,v)$ uniquely complete to $(u+v,u),(v,u+v)$. The ledge sets $u=b$ and the slope sets $v=f(b+v)$. Both constraints are independent since $1-k\ne0$. All four contacts are interior to their edges. The four-parameter square residual of v345 therefore has invertible linearization and local degree of absolute value one, which persists under uniform registered redraws.

A sufficient disk test for this square is $S_E\ge\max(a,b,c)$. An axis-aligned square of greater side extends above $b$ and is wider than every cross section there. Every disk of radius $r$ contains an axis-aligned square of side $\sqrt2 r$. For example $W=10,H=9,a=3,b=4,c=1/2$ gives an inward tilt and $S_E^2=20$, clearing this test.

For an outward tilt set $a=b=1$, $k=1/2$, $c=(H+1)/2$ and $W=H+9$ for $H\ge81$. The elbow square always has side $\sqrt5$. Meanwhile the rectangle $[0,H/4]\times[H/2,3H/4]$ fits in the upper arm, hence $r\ge H/8$. The ratio $S_E/(\sqrt2 r)$ tends to zero. Its odd degree alone gives no control of that ratio.

## A square supported by the sloping wall

For $0<k<2$ set $D=1+k+k^2$, $d=c/D$, $y_0=b-a/k$, and $P=(d,y_0+d/k)$. Starting at $P$, take side vectors $(kd,d)$ and $(-d,kd)$ in that order. These give a square of side $S=d\sqrt{1+k^2}$. If $d>a$, both endpoints of the first side lie strictly inside the tilted edge. The third vertex has height $H$ because $H-y_0=c/k$, and the fourth has abscissa zero. The third abscissa is $kd<c$, the fourth height is strictly between $b$ and $H$, and the upper slope contact lies below $H$ by $kd$. All four contacts are interior.

This square has local degree of absolute value one. In its marked edge chart, write the first contact as $(x,y_0+x/k)$ and the edge to the second as $(kh,h)$. Completing the square forces the left contact to have abscissa $x-h=0$ and the top contact to have ordinate $y_0+x/k+(1+k)h=H$. The two equations have determinant $1+k+1/k>0$. The four contact parameters are then uniquely recovered. This proves regularity and the same persistence as above.

**Disk comparison.** A sufficient condition for the full disk bound, in addition to $d>a$, is
$$
S^2\ \ge\ 2\max(b,a/(2-k))^2 .
$$
To prove it, split interior disks according to whether they go below the ledge. A disk entirely above height $b$ fits in the triangle bounded by the left wall, the top line and the extended tilted line. Our square has all four corners on this triangle. The zero-depth convex lemma of [v307](v307-convex-shell-square-gap.md) gives $S\ge\sqrt2$ times the triangle's inradius. For a disk going below the ledge, write $y$ for the center height and $r$ for the radius. If $y\le b$, the bottom edge gives $r\le b$. Otherwise $y<b+r$ and the horizontal diameter gives $2r\le f(y)<a+kr$, hence $r<a/(2-k)$. This proves the test.

For the entire outward family above, $d=4c/7>1$ and $S^2=20c^2/49>2$. Thus the replacement square clears the disk target for every $H\ge81$, even as the elbow square falls arbitrarily far short. At $H=81$, their respective side lengths are approximately 26.19 and 2.24.

## The corner transition has degree zero

Hold $a,b,k$ fixed and vary $L=H-b$ near $L_0=a(1+k)$, adjusting $c=a+kL$ and keeping $W>c$. At equality $d=a$ and only the first vertex is at a corner. Thus locally there are just two incident charts. The two-slope-contact chart exists on the side $L>L_0$. In the other chart the first point is $(x,b)$ on the ledge and the first side vector is $(h,j)$. Left and top contact give $x=j$ and $h+j=L$. Contact of the second point with the sloping line gives $L=a+kj$, so $j=(L-a)/k$. The first point enters the ledge precisely for $j>a$, again $L>L_0$. All other constraints remain strict near equality. The two distinct squares on that side approach the same corner square. On the other side neither exists in this neighborhood. Consequently the isolated corner zero has degree zero, and these registered redraws actually remove it locally. Another component is needed to carry a guarantee through this transition.

## Checks and remaining question

The [exact tests](../tests/test_tilted_wall_v346.py) trace both marked components on rational examples using the independent v341 polygon engine and check their minimum squared side and signed degree. The formulas prove the family statements. The tests also recover degree zero at the corner transition. A selection rule covering every tilted L remains open.
