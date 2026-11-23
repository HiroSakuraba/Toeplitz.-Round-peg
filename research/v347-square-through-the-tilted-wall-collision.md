# v347. A third square survives the tilted-wall collision

At the [v346 collision](v346-tilted-inner-wall.md), two squares of side $\sqrt{20}$ disappear together when the upper arm is shortened. A third square has side $\sqrt{40}$ at the collision and persists on both sides. Its corners touch the bottom, ledge, top and left wall. The sloping wall can move without disturbing any of those contacts.

## A full-height square

Use the v346 polygon with bottom width $W$, ledge height $b$, elbow at $(a,b)$ and top right at $(c,H)$. Put $L=H-b$. Whenever $a<H<W$ and $L<c$, use
$$
(b,0),\quad(H,b),\quad(L,H),\quad(0,L).
$$
The consecutive side vectors $(L,b)$ and $(-b,L)$ define a square of side $S=\sqrt{L^2+b^2}$. Each contact is strictly inside its edge. Since $2(L^2+b^2)\ge (L+b)^2=H^2$, its side is at least $H/\sqrt2$. Every interior disk has diameter at most $H$, proving the v343 disk bound once persistence is checked.

Write the four marked points in local coordinates as $(x,0),(y,b),(z,H),(0,w)$. Completing the adjacent bottom and left corners gives $(x+w,x)$ and $(w,x+w)$ for the other two. The four scalar residuals are therefore
$$
y-x-w,\quad b-x,\quad z-w,\quad H-x-w .
$$
Setting them to zero determines $x=b$, $w=H-b$, $y=H$, $z=w$ uniquely. This affine system is invertible, so the square has local degree of absolute value one. Uniform registered redraws preserve the degree in any sufficiently small parameter box around the solution. Shrinking the box gives survival above every strict side cutoff below $S$.

At a v346 collision of positive slope $k$, we have $L=a(1+k)$ and $c=a+kL$. Consequently $c-L=ak^2>0$ and $H>a$. Thus **every collision with $W>H$** satisfies the strict contact hypotheses on an open interval through the collision. For the original example $W=12,a=4,b=2,k=1/2$, the entire interval is $4<H<10$, with the collision at $H=8$.

## A finite redraw tolerance at the collision

For $H=8$ the four coordinates are $(2,8,6,6)$. Take the closed coordinate box of radius $1/2$ around them, which lies inside the four contact edges. A Euclidean redraw of size $e$ perturbs each completed-corner residual by at most $4e$ per coordinate, by the triangle inequality applied to the four point terms (including the quarter turn). The same holds throughout the straight residual homotopy. At any zero on this homotopy, substitution in the four original residuals gives
$$
|x-2|\le4e,\quad |w-6|\le8e,\quad |y-8|\le8e,\quad |z-6|\le12e.
$$
For $e\le1/64$ no zero reaches the box boundary. Degree therefore forces a square on the redraw. Its side, measured between the redrawn bottom and left points, is at least $\sqrt{40}-(\sqrt{80}+2)e>6$. Even the *redraw's* disk target is at most $(8+2e)/\sqrt2<6$, by its enclosing horizontal strip.

## When the outer width is shorter

Suppose $b<W<H$, $k>0$ and $c=a+k(H-b)<W$. Set $v=a+k(W-b)$ and $u=W-v$. If $W-b<v$, then $0<u<b$ and
$$
(u,0),\quad(W,u),\quad(v,W),\quad(0,v)
$$
are interior contacts on the bottom, right wall, slope and left wall. They span width $W$. Their side is at least $W/\sqrt2$, enough for every interior disk. Completing the anchors as above reduces the regularity check to $u+v=W$ and $v=a+k(u+v-b)$. The determinant in the anchor coordinates is one. At the collision, $a+k(t-b)-(t-b)$ is positive at both $t=b$ and $t=H$ (values $a$ and $ak^2$), hence at $t=W$. This proves persistence on both sides here as well.

At $W=H>b$ and $H-b<c$, $a<H$, with positive slope and $c<W$, the two constructions meet along an interval. Set $t=H$ and take $(u,0),(t,u),(v,t),(0,v)$ with $u+v=t$ and $t-c\le u\le b$. These squares span both dimensions. To isolate their marked family, use the open bottom and left anchor edges, the right-wall/ledge arc through $(W,b)$ and the slope/top arc through $(c,H)$, restricting around the contact traces. All anchor contacts are interior since $0<t-c<b<t$. The other contacts stay away from the elbow since $t>a,b$. Near the interval, $t'=u+v\le W,H$. If $t'<W,H$, the two middle contacts must be on the ledge and slope, forcing $u=b$ and $v=a+k(t'-b)$. This is impossible near $t'=H$ since $H-b<c$. Thus the displayed interval is the whole local zero set. Slightly increasing $W$ with $H,a,b,c$ fixed leaves exactly the regular full-height solution locally (the same strict inequality rules out $t'<H$). Degree is therefore odd on the original interval. Its minimum side is at least $H/\sqrt2$. Together the three cases cover collisions with $b<W$. The argument leaves the thick lower-arm regime $b\ge W$ outside this result.

[Exact tests](../tests/test_tilted_collision_v347.py) replay all three width orderings and check the finite-error inequalities. The universal tilted-L selection question remains open.
