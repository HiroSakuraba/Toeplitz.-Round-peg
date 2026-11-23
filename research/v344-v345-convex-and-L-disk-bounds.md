# v344–v345. A square survives redraws of convex curves and orthogonal L shapes

Benjamin John Schulz · November 2026

Draw the largest disk that fits inside a simple closed curve. The square whose four corners lie on that disk has side $\sqrt2 r$, where $r$ is the disk radius. Can we guarantee a square at least that large with its corners on the *curve*? Can we still get arbitrarily close to that size after redrawing the curve? [v343](v343-incircle-and-disk.md) asks this for every Jordan curve. We prove the bound for two families: every convex boundary, and every six-sided orthogonal L with arbitrary real dimensions. The redraw may have infinitely many folds. For the convex case we also bound the loss at a specified drawing error.

Let $M(\Gamma)$ be the supremum of inscribed square sides (zero if there are none). A **robust lower bound** $s$ means that every $s'<s$ admits a tolerance such that $M(\Gamma')>s'$ for every registered Jordan redraw within that tolerance. Registration pairs points via boundary parametrizations, and the error is the largest Euclidean distance of any paired points. A single redraw can have a slightly smaller square than $s$. The four corners may move to new points of the boundary.

## v344. A finite-error bound for any convex reference

**Theorem.** Suppose $B(O,r)\subset K\subset B(O,R)$ for a compact convex body $K$. Let $\Gamma'$ be a registered Jordan redraw of $\partial K$ with error at most $e$. If $e<r$ and $D=(r-e)^2-8e(R+e)>0$, then $\Gamma'$ has a square of half-diagonal at least
$$
\frac{r-e+\sqrt D}{2}.
$$
Equivalently, a sufficient test for a square of side strictly greater than $\sqrt2 q$, for $0<q<r$, is
$$
q(r-q)>e(q+2R)+2e^2.
$$
Taking $r$ to be the inradius and letting $e$ tend to zero gives the robust bound $\sqrt2 r$.

**Proof.** Thicken $K$ outward by $e$, forming the convex body $K_e=K+B(0,e)$. It contains $B(O,r+e)$ and lies in $B(O,R+e)$. For a registered point $p'$ near $p\in\partial K$, take any supporting line at $p$. The parallel supporting line of $K_e$ is $e$ farther out. The inward distance of $p'$ from this line is at most $2e$, and $p'$ lies in $K_e$. Following the outward normal from $p'$ reaches the boundary of $K_e$ within $2e$. The straight homotopy between paired curve points stays at distance at least $r-e>0$ from $O$, so the redraw still winds once about $O$ up to orientation.

Apply the [v307 large-square theorem](v307-convex-shell-square-gap.md) in this layer with radii $r+e,R+e$ and depth $2e$. Its gap polynomial becomes $q^2-(r-e)q+2e(R+e)$. If it is negative, v307 forces a square above $q$, on the large side of the gap. When $D>0$ and $e<r$, the choice $q=(r-e)/2$ makes it negative. Its larger root gives the asserted bound. Conversely, for each fixed $0<q<r$, the strict cutoff test holds for all sufficiently small $e$.

**Controls.** With $r=1$, $R=2$ and $e=1/100$, the guaranteed half-diagonal is $(99+\sqrt{8193})/200>9/10$. For a radius-$r$ circle the reference itself admits no square larger than $\sqrt2 r$, giving equality in the robust bound. For the side-$S$ square the inradius is $S/2$ and the robust bound is also attained, by [v342](v342-unmarked-components-and-square-frames.md). The finite-error estimate above is sufficient and can be conservative.

## v345. Every orthogonal L

After a rigid motion or reflection, the region is
$$
L=([0,W]\times[0,b])\ \cup\ ([0,a]\times[0,H]),\qquad
0<a<W,\quad 0<b<H .
$$
Write $m=\min(W,H)$. These four real lengths allow unequal arms, arbitrarily thin arms, and ties. A rectangle obtained by filling the missing corner is covered by v344.

**Theorem.** Every such $\partial L$ has robust largest-square side at least $\sqrt2\,r(L)$.

### The two arms meet at the square

First suppose $a,b<m$. Mark four boundary arcs in counterclockwise order: the open bottom edge, the right-hand boundary where $x>a,y>0$, the upper boundary where $y>b,x>0$, and the open left edge. The middle arcs each go around a convex corner and omit their endpoints. The conditions $x>a,y>0$ and $y>b,x>0$ select, respectively, the right wall plus horizontal ledge and the upper inner wall plus top. Mark adjacent square corners $(u,0)$ and $(0,v)$ on the bottom and left. In counterclockwise order the four corners must be
$$
(u,0),\qquad(t,u),\qquad(v,t),\qquad(0,v),\qquad t=u+v .
$$
On the indicated arcs this is equivalent to
$$
0<u\le b,\quad 0<v\le a,\quad \max(a,b)<t\le m,
\qquad (t=W\text{ or }u=b),\quad(t=H\text{ or }v=a).
$$
Indeed the second corner must be on either the right wall or the ledge, and the third must be on either the top or the inner wall. In particular $t=\min(W,H,a+b)$: it cannot exceed any of these numbers, and if it is below both outer spans, both arm constraints are equalities. Conversely, this value exceeds $\max(a,b)$ and the following choices give **all** marked squares.

| Outer spans reached by $t$ | Bottom coordinate $u$ | Left coordinate $v$ |
|---|---|---|
| Neither | $b$ | $a$ |
| Only $W$ | $t-a$ | $a$ |
| Only $H$ | $b$ | $t-b$ |
| Both | any $u\in[t-a,b]$ | $t-u$ |

The squared side is $u^2+(t-u)^2$. On an interval its minimum occurs at the point closest to $t/2$. This gives a computable robust floor that can exceed the disk target.

**Why the square survives.** Use four independent boundary parameters, one on each marked arc, and the four-coordinate residual
$$
(p_2-p_1-J(p_1-p_0),\ p_3-p_0-J(p_1-p_0)),
$$
where $J$ is counterclockwise rotation by a right angle. The table gives its entire zero set on this product of open arcs. That set is a compact point or interval strictly inside the product: $t>\max(a,b)$, $u,v>0$ and the two anchor points lie in their edge interiors. Work in a small neighborhood whose closure is inside the product and whose boundary has no zero.

Change the lengths slightly so $W,H,a+b$ are pairwise distinct, keeping $a,b<m$, and register old and new boundaries edge by edge. There is exactly one marked solution on the perturbed polygon, at interior points of its four active edges. It lies in the chosen neighborhood for all sufficiently small generic changes: the table shows every limit of these solutions belongs to the original family. Conversely a square in the neighborhood satisfies the same marking conditions by compact containment and edge registration. For example, near the middle two contacts the inequalities $x>a$ and $y>b$ defining their respective arcs have strict margins.

The unique solution has nonsingular derivative. With $p_0,p_3$ fixed on their transverse anchor edges, completion fixes $p_1,p_2$ uniquely as above. The two affine equations for the active lines are one of $(u=b,v=a)$, $(u+v=W,v=a)$, or $(u=b,u+v=H)$. Each pair has independent normals in $(u,v)$. The remaining two edge parameters are determined at nonzero speed. Thus the four-parameter linearized residual has trivial kernel. The nearby regular polygon contributes one zero modulo two. Uniform continuity of the residual under the registered deformation preserves local degree, also when the original family slides or when all three lengths tie.

Finally take any cutoff below the minimum side on the family and shrink the neighborhood so every trial quadruple there has its first two points farther apart than the cutoff. A sufficiently close arbitrary registered redraw keeps the residual nonzero on the neighborhood boundary, preserves its odd degree, and still clears the cutoff there. This argument follows the whole family, so the surviving square need not approach any one preselected member.

**Comparison with the disk.** If $t=m$, the squared side is at least $m^2/2$. Every disk in $L$ has diameter at most $m$. If $t<m$, the side is $\sqrt{a^2+b^2}\ge\max(a,b)$. Any axis-aligned square contained in $L$ has side at most $\max(a,b)$, since a larger one would have a top-right corner strictly inside the missing quadrant. A disk of radius $r$ contains an axis-aligned square of side $\sqrt2 r$. Both cases give the claimed bound.

### When one arm spans the shorter dimension

If $b\ge m$, then $W=m<H$. When $b>W$ the two end walls are taller than $W$ and [v340's tall-wall theorem](v339-v340-end-wall-selection.md) gives robust side exactly $W$. When $b=W$, reflect left to right: the initial horizontal roof at height $W$ has length $W-a$ and rises to $H$. The v340 end-plateau proof applies directly to its normalized real-parameter model. It gives robust side at least $\sqrt{d^2+(W-d)^2}$ with $d=\min(W-a,W/2)$. The minimum possible value is $W/\sqrt2$. In both cases $2r\le W$. When $a\ge m$ interchange the two axes. This exhausts the remaining cases, including the transition at equality.

For $W=2,H=7,a=1,b=2$, the disk radius is 1 and the plateau bound is $\sqrt2$. The [v337 sharp example](v337-terrain-component-degree.md), reflected horizontally, supplies the matching upper bound on the robust side. Thus the constant is sharp even within this nonconvex class.

### Computing the inradius independently

For a disk of proposed radius $x$, its center coordinates must lie at least $x$ from the bottom and left edges. Of these possible centers, $(x,x)$ is farthest from the missing northeast quadrant, because distance to that quadrant is nonincreasing in each coordinate. Thus the disk fits exactly when
$$
2x\le m,\qquad (a-x)_+^2+(b-x)_+^2\ge x^2 .
$$
For $A=\min(a,b)$ and $B=\max(a,b)$ this yields
$$
r=\min(m/2,r_\infty),\qquad
r_\infty=\begin{cases}
B/2 & B\ge 2A,\\
A+B-\sqrt{2AB} & B<2A .
\end{cases}
$$
In the first case the quadrant constraint becomes $B-x\ge x$ once $x\ge A$. In the second the threshold lies below $A$; solving the quadratic gives the smaller root. For instance, with $W=9,H=7,a=b=1$, the largest disk reaches around the elbow and has radius $2-\sqrt2$. The selected square has side $\sqrt2$.

## Verification and next question

The [rational implementation](../code/convex_l_disk_v344_v345.py) constructs the marked interval and computes its minimum squared side. The [tests](../tests/test_convex_l_disk_v344_v345.py) compare the floor with the exact algebraic inradius for all 441 integer L shapes with $2\le W,H\le7$, check reflection and rational rescaling, and use the independent v341 polygon chart tracer to recover the minimum side and absolute degree one on ten selected components. These include a continuous sliding family, a three-way tie and fractional lengths. The proofs above cover real lengths and arbitrary registered redraws.

```sh
python -m unittest discover -s tests -p 'test_convex_l_disk_v344_v345.py' -v
```

A next test is to tilt the inner wall of the L while keeping the other five edges in their original directions. The two anchor corners still determine the square. The contact condition $v=a$ becomes a linear equation coupling $v$ and $t$. The open question is whether a marked family still has odd degree **and** reaches the disk target through every transition for that larger class. Arbitrary Jordan boundaries remain open.
