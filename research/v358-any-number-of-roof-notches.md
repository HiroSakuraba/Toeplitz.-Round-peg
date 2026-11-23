# v358. Rectangles with any number of roof notches

Cut $k\ge1$ disjoint open rectangular bites $(a_j,b_j)\times(k_j,H]$ from the top of the box $[0,W]\times[0,H]$, with $0<k_j<H$ and positive roof length at both outer corners and between consecutive cuts. Write $\Gamma$ for the boundary, $r$ for the inradius, *towers* for the roof pieces and *walls* for the vertical notch edges $x=a_j,b_j$, each spanning heights $(k_j,H)$. As before, a robust floor $s$ means that every strict cutoff below $s$ survives every sufficiently close registered Jordan redraw. The target is the disk bound $s\ge\sqrt2 r$.

[v352](v352-all-rectangular-roof-notches.md) proved the disk bound for one notch, and [v353](v353-two-cuts-and-degree-transfer.md)–[v357](v357-simultaneous-loss-of-two-square-pairs.md) handled several two-notch regimes. This note does not choose two notches and track collisions. It describes whole classes of squares for every $k$ at once.

**Results.**
1. **Squares standing on the bottom edge.** For every $k$, every square with two vertices in the open bottom edge is of one of four kinds. Three are transverse (degree $\pm1$): an axis square between two walls, and a square with one top vertex on a wall and the other on a floor, in either order. The fourth is a slide along one horizontal level. Its mod-two degree is given by a four-letter interleaving rule. Development runs checked the rule against the independent v341 tracer on 1,757 roof slides, and the tests check every closed-form component they generate.
2. **Squares anchored on the bottom and a wall.** For every $k$, squares with consecutive vertices in the open bottom edge and in an open wall have four closed-form transverse types. A floor–roof square always reaches the disk target.
3. **The inradius is computable exactly.** It is $H/2$, half a floor height, or a root of one of finitely many explicit linear or quadratic equations. So $\sqrt2 r\le s$ is decided exactly, including the many exact ties of integer shapes.
4. **Earlier theorems for any $k$.** The tall-box, square-box, wide-end-tower, end-plateau and interior-slot theorems hold verbatim for every $k$.
5. **Computer-assisted proofs.** The certifier below decided every shape in five exhaustive integer families with $k=1,\dots,5$ (154,574 shapes) and in 5,000 seeded random rational shapes, with no failures. Every certificate is an exact proof for its shape.
6. **Five notches at a finite tolerance.** Consider a $46\times10$ box with five notches whose largest contained axis square has side 8. Every registered redraw within $1/32$ keeps a square of side greater than $8.18$. The disk target is $\sqrt2(15-4\sqrt7)\approx6.25$.

A selection rule proved for every notched rectangle remains open. The last section shows examples of where the closed forms stop.

## 1. Squares standing on the bottom edge

Two adjacent vertices on the line $y=0$ force the square to be axis-aligned: $(x,0),(x+s,0),(x+s,s),(x,s)$. It is inscribed exactly when $(x,s)\in\Gamma$ and $(x+s,s)\in\Gamma$. Put $S(x,y)=(x-y,y)$. In the $(x,s)$ plane these squares are the intersections of the two curves
$$
A=\Gamma,\qquad B=S(\Gamma).
$$
This is the shear of [v316](v316-roof-shear.md). Eliminating the two bottom parameters, exactly as in v316's bent-base proof, shows that **the mod-two local degree of an isolated family equals the local intersection number of $A$ and $B$** there. This needs both bottom vertices in the open bottom edge.

$S$ maps a horizontal edge to a horizontal segment and a vertical edge $x=c$ to the antidiagonal piece $\{(c-y,y)\}$. Vertical pieces of $A$ cross antidiagonal and horizontal pieces of $B$, and horizontal pieces of $A$ cross antidiagonal pieces of $B$, all transversally. This gives three transverse kinds:

| Kind | Square | Conditions (all strict) |
|---|---|---|
| wall pair | $[w_1,w_2]\times[0,s]$, $s=w_2-w_1$ | both walls' floors $<s<H$ |
| wall–floor | $[w,w+k_m]\times[0,k_m]$ | wall floor $<k_m$, $a_m<w+k_m<b_m$ |
| floor–wall | $[w-k_m,w]\times[0,k_m]$ | wall floor $<k_m$, $a_m<w-k_m<b_m$ |

The walls may belong to any notches, and the bites between them do not matter: only the four vertices are on $\Gamma$. In the five-notch example below, the wall pair on $x=29$ and $x=38$ has side 9 and passes over two bites.

The remaining intersections are overlaps of horizontal pieces at one height $c$: the roof ($c=H$) or floors of height $c$. Let the top-left vertex run on piece $[s_1,e_1]$ and the top-right vertex on piece $[s_2,e_2]$ at the same height. The maximal family is $x\in[x_0,x_1]$ with $x_0=\max(s_1,s_2-c)$ and $x_1=\min(e_1,e_2-c)$.

**Theorem 358.1 (slide parity).** Say $A$ *leaves* at an end of the overlap when the top-left vertex reaches an end of its piece, and $B$ leaves when the top-right vertex does. The boundary leaves a roof piece downward and a floor piece upward. $A$ leaves along a vertical ray and $B$ along the sheared image of a vertical wall, which points down-right or up-left. Read the exits of $A$ and $B$ counterclockwise around a thin neighbourhood of the overlap:
- **roof:** west, $A$ left, $B$ left, $A$ right, $B$ right, east;
- **floor:** east, $A$ right, $B$ right, $A$ left, $B$ left, west.

A curve that does not leave at an end exits there to the west or east instead. For a single square ($x_0=x_1$) the two ends merge, and the vertical exit comes before the sheared one. **The family has odd degree exactly when the exits of $A$ and $B$ alternate.** Exclude a roof slide that reaches a box corner with a wall at distance $H$ from the opposite outer wall: $x_0=0$ with a wall at $x=H$, or $x_1=W-H$ with a wall at $x=W-H$. Every other slide is an isolated family whose degree is given by this rule. A roof slide has side $H\ge2r$, so an odd one gives the disk bound.

*Proof.* Away from the box corners, both lower vertices of every nearby square lie in the bottom edge. So the $(x,s)$ picture is complete near the family. Near the overlap, $A$ and $B$ consist only of the two pieces and the rays leaving them, and those rays point in different directions. So the overlap is a component of $A\cap B$. Push $B$ off the line by a small amount. Its intersections with $A$ are where the rays cross, and their number is odd exactly when the exits alternate.

At a box-corner end, say $x_0=0$, $H$ lies in an open tower. The computation at v352's join applies to any square near $[0,H]^2$ with its lower left vertex $(0,y)$ on the outer wall and lower right vertex on the bottom. Its upper left vertex would be $(y,D+y)$ on the roof, and its upper right vertex $(H,H-y)$. Near $(H,H)$ the boundary is only roof, so $y=0$. Thus all nearby squares are axis squares. Now push the corner $(0,0)$ outward by extending the bottom edge slightly to the left before the wall rises. This creates no new zero, since the upper left vertex would leave $\Gamma$. After the change the bottom edge is straight beyond the corner, so the elimination applies. The zero set is constant along the straight homotopy, so the degree does not change. ∎

Two consequences are worth stating separately.

- **Corollary 358.2.** Suppose $H$ lies in an open tower $(s,e)$ with $e<a_1+H$. Then the slide from the left box corner has the pattern "$A$ left, $B$ right", so it is odd and the robust floor is at least $H$. The same holds after reflection. With $(s,e)=(b_k,W)$ this is the condition that all notches lie in $(W-H,H)$. That is v351's sliding test, now for any number of notches.
- In the $28\times17$ box with notches $[11,13]$ at height 16 and $[18,24]$ at height 5, the slide $x\in[0,1]$ runs from the box corner to the notch corner $(18,17)$. It has the pattern "$A$ left, $B$ right" and odd degree, as the tracer confirms. The largest contained axis square there is smaller than $H$.

## 2. Squares anchored on the bottom and a wall

Let $x=c$ be a wall with the region on its right: the outer left wall ($c=0$, spanning $(0,H)$) or the right wall of a notch ($c=b_j$, spanning $(k_j,H)$). Write a square with consecutive vertices in the open bottom edge and in the open wall as
$$
(c+u,0),\quad (c+t,u),\quad (c+v,t),\quad (c,v),\qquad t=u+v .
$$
The two anchors lie on perpendicular straight edges. Eliminating them, as in v345, identifies the mod-two degree with the local intersection number of the planar curves $\varphi_1(\Gamma-c)$ and $\varphi_2(\Gamma-c)$ in the $(u,v)$ plane, where $\varphi_1(x,y)=(y,x-y)$ and $\varphi_2(x,y)=(y-x,x)$. The second vertex lies on $\Gamma$ exactly when $(u,v)\in\varphi_1(\Gamma-c)$, and the third exactly when $(u,v)\in\varphi_2(\Gamma-c)$.

Consider the vertical line $x=c+t$ at height $u$ and the horizontal line $y=t$ at abscissa $c+v$. The second vertex can be on a floor ($u=k_j$) or a wall ($c+t=w$). The third can be on a wall ($c+v=w'$), a floor ($t=k_m$) or the roof ($t=H$). A roof position for the second vertex would force $v\le0$. The four transverse combinations are:

| Kind | $(u,v)$ | Second vertex | Third vertex |
|---|---|---|---|
| floor–wall | $(k_j,\ w'-c)$ | interior of floor $j$ | interior of wall $w'$ |
| floor–floor | $(k_j,\ k_m-k_j)$ | interior of floor $j$ | interior of floor $m$ |
| floor–roof | $(k_j,\ H-k_j)$ | interior of floor $j$ | open tower |
| wall–wall | $(w-w',\ w'-c)$ | interior of wall $w$ (or the right outer wall) | interior of wall $w'$ |

In each case the anchors must also be interior: $k_c<v<H$ and $0<c+u<W$. Each row is a transverse intersection of a vertical or antidiagonal piece with a horizontal or antidiagonal piece, so it is a regular zero with degree $\pm1$. Reflection gives the walls with the region on their left. The excluded combination puts the second vertex on a wall $x=c+t$ and the third on a horizontal edge at height $t$. It needs a coincidence and gives the sliding families of [v345](v344-v345-convex-and-L-disk-bounds.md) and v355. The squares of [v345](v344-v345-convex-and-L-disk-bounds.md), [v355](v355-tall-strip-and-grounded-disks.md)–[v357](v357-simultaneous-loss-of-two-square-pairs.md) with $c=0$ are instances of the first, third and fourth rows.

**Corollary 358.3.** A floor–roof square has squared side $k_j^2+(H-k_j)^2\ge H^2/2\ge2r^2$, so it always reaches the disk target. For the outer left wall it exists when $a_j<H<b_j$ and $H-k_j$ lies in an open tower. For the right wall of notch $m$ it exists when $a_j<b_m+H<b_j$, $b_m+H-k_j$ lies in an open tower and $k_m<H-k_j$. Either case gives the robust disk bound, for any number of notches, and the same holds after reflection.

## 3. The inradius, exactly

Lower a disk until it touches the bottom. Its centre $(x,s)$ must then satisfy $x\ge s$ or $x\ge b+\ell$ for each obstacle on its left, and $x\le W-s$ or $x\le a-\ell$ for each obstacle on its right. Here $\ell=s$ at a notch whose floor is at most $s$, and $\ell=\sqrt{2ks-k^2}$ at a notch whose floor $k$ lies in $(s,2s)$. This is the [v355 test](v355-tall-strip-and-grounded-disks.md). At $s=r$ the feasible centres stop existing. This happens in one of three ways:
- at $s=H/2$;
- at $s=k_j/2$, when a notch above the only feasible centres becomes active;
- at a pair $\ell_1+\ell_2=D$ that becomes tight.

The pair equations are $2s=D$, $s+\sqrt{2ks-k^2}=D$ (a quadratic after squaring once) and $\sqrt A+\sqrt B=D$ with $A,B$ linear in $s$ (a quadratic after squaring twice). So $r$ is one of finitely many numbers $(p+\sigma\sqrt q)/m$.

The decision procedure brackets $r$ by exact bisection with the v355 feasibility test. It stops when every candidate in the bracket lies on the same side of $\sqrt{F/2}$, where $F$ is the squared side being compared. Each comparison is the sign of $\alpha+\beta\sqrt q$. This settles exact ties such as the $10\times9$ box with notches $[1,2]$ at height 4 and $[5,9]$ at height 8. There $r=4$ and the best odd component has squared side exactly $32=2r^2$. The tests compare the procedure with the exact comparison-L inradius of [v351](v351-one-rectangular-roof-notch.md) on single notches.

## 4. Earlier theorems that hold for any number of notches

- **$H>W$.** Both outer walls are taller than the width, and the roof may dip ([v340](v339-v340-end-wall-selection.md)). The robust side is $W\ge2r$.
- **$H=W$.** [v353](v353-two-cuts-and-degree-transfer.md) already allowed any finite number of notches.
- **An end tower wider than $H$, or of width exactly $H$.** v351's comparison L uses only the first notch. Its sliding family, or the plateau family when the width equals $H$, has every contact in $x\le a_1$, so it transfers unchanged.
- **An interior slot.** Suppose some largest contained axis square of side $i$ lies between a wall descending below $i$ on each side. Then Theorem 338.1 gives the robust floor $i\ge\sqrt2r$. Its proof uses the bent-base theorem with $h=i-3\delta$. Real widths require only that $\delta$ be small compared with the distances of those walls from the outer walls and with $i$ minus their floors.

## 5. Certificates and sweeps

The certifier tries three levels in order. Each returns an exact proof for the shape, or passes it on.
1. **Theorem.** One of the cases in Section 4 applies.
2. **Closed form.** The best odd component from Sections 1–2 reaches the disk target in the exact test.
3. **Search.** An exhaustive enumeration of every transverse square: each vertex lies in an open edge and fixes one coordinate, and the $4\times4$ contact system is nonsingular. Failing that, a v341 trace from every closed contact solution above $\sqrt2$ times the lower radius bound. A component counts only if its traced degree is odd and its exact minimum side reaches the target.

The optional numpy pass only prunes the enumeration. Every reported square is recomputed in rational arithmetic.

| Family (integer shapes up to reflection, $H<W$) | Shapes | Theorem | Closed form | Search | Failures |
|---|---:|---:|---:|---:|---:|
| $k=1$, $W\le14$, $H\le13$ | 10,073 | 2,709 | 6,215 | 1,149 | 0 |
| $k=2$, $W\le10$, $H\le9$ | 20,093 | 1,948 | 14,154 | 3,991 | 0 |
| $k=3$, $W\le10$, $H\le8$ | 46,193 | 8,272 | 33,199 | 4,722 | 0 |
| $k=4$, $W\le11$, $H\le6$ | 27,115 | 10,230 | 16,717 | 168 | 0 |
| $k=5$, $W\le13$, $H\le5$ | 51,100 | 29,129 | 21,962 | 9 | 0 |
| random rational, $k=1,\dots,5$, denominator 4, 1,000 each | 5,000 | 3,191 | 1,782 | 27 | 0 |

The $k=1$ row only cross-checks v352. Integer shapes produce many exact ties, which is why the search level is visited so often there. In the random rational samples, where such ties are rare, the search level was needed 27 times in 5,000. The [results](../results/v358/) record the certificate kind for every family.

## 6. Five notches at a finite tolerance

Take $W=46$, $H=10$ and notches $[7,12]$, $[13,29]$, $[30,31]$, $[33,38]$, $[40,43]$ at floor heights $8,5,3,2,9$. No theorem of Section 4 applies, and the largest contained axis square has side 8. The inradius is exactly $15-4\sqrt7\approx4.417$, so the disk target is about $6.25$.

The wall pair $(29,0),(38,0),(38,9),(29,9)$ stands on the right walls of the second and fourth notches. Its top vertices have clearance 1 on walls spanning $(5,10)$ and $(2,10)$.

**Finite margin.** Let $z$ collect the four edge coordinates of the vertices. The four-point residual is $M(z-z^*)$, and the exact computation gives $\lVert M^{-1}\rVert_\infty=3$. A registered redraw of Euclidean error $e$ changes each residual coordinate by at most $4e$, also along the straight homotopy. So every zero has $\lvert z-z^*\rvert_\infty\le12e$, and the box of radius 1 keeps a zero-free boundary when $12e<1$. At a zero, the first two redrawn vertices are at least $9-24e-2e$ apart. For $e=1/32$ this is $131/16>8.18$, which clears the target.

The tilted floor–roof square $(8,0),(10,8),(2,10),(0,2)$ gives $\sqrt{68}-13/16>7.4$ under the same tolerance. The tests check both squares with the independent tracer and both margins with the exact disk test.

## 7. Where the closed forms stop

The search-level shapes we inspected show three mechanisms. This is a description of examples, not a classification:
- **Exact ties at notch corners.** A transverse square of Section 1 or 2 has a vertex exactly at a corner. Examples are the wall pair $[2,5]$ under a floor of height 3 in the $8\times4$ box with notches $[1,2]$ at height 3 and $[5,6]$ at height 1, and the tie $F=2r^2$ above.
- **Turning at a box corner.** Some wall lies at distance $H$ from an outer wall, and a roof slide reaching that corner turns along the wall. This is v352's collision $H=B$. With more notches the turning branch can be cut by further roofs, as in v353.
- **Tilted squares with no anchor on an outer wall.** Rarely, transverse squares with every vertex on notch edges or the right outer wall.

The traced degree settled each case. But no proof shows that some component of Sections 1–2, or a turning branch, always reaches the target. That is the open selection problem for this class.

[Code](../code/multi_notch_v358.py) · [Sweeps](../code/multi_notch_sweep_v358.py) · [Tests](../tests/test_multi_notch_v358.py) · [Results](../results/v358/)

```sh
python -m unittest discover -s tests -p 'test_multi_notch_v358.py' -v
python code/multi_notch_sweep_v358.py exhaustive 3 10 8
python code/multi_notch_sweep_v358.py random 5 1000 4 5
```
