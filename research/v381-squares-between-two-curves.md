# v381. Squares between two curves

[v380](v380-comb-laws-for-any-arms.md) found, numerically, that aligned squares, with two corners on each of two neighbouring turns, sit exactly where the gap between the turns is stationary. This note proves it, for any two curves, smooth or polygonal.

The idea is to stand the square on the lower curve. Its bottom-left corner is a single unknown $x_1$. The top-left corner touching the upper curve then fixes the side. Whether the top-right corner also lands on the upper curve becomes one real function $H(x_1)$, and the squares are exactly its zeros. $H$ is the change of the gap across the square, plus an error controlled by how steep the curves are and how much the upper one bends. That gives both halves at once. Where the gap rises and then falls by more than the error, $H$ changes sign and a square must exist. And wherever a square exists, the gap barely changes across it.

**Results.**
1. **Reduction (Theorem 381.1).** For graphs $f_1<f_2$ with slopes at most $\sigma<1$, the squares standing on $f_1$ with their top on $f_2$ are exactly the zeros of a continuous function
$$H(x_1)=-\big[g(x_1+u)-g(x_1)\big]+R,\qquad |R|\le|s|\cdot{\rm TV}(f_2'),\quad |s|\le\sigma u,$$
where $g=f_2-f_1$, $u=u(x_1)$ is the horizontal extent of the side, $s$ its rise, and ${\rm TV}$ is the total variation of $f_2'$ over the square's footprint.
2. **Existence (Corollary 381.2).** If the gap across one side-length rises by more than the error bound at one position and falls by more than it at another, there is a square in between.
3. **Localization (Corollary 381.3).** Every such square has $|g(x_1+u)-g(x_1)|\le|s|\,{\rm TV}(f_2')$. For smooth curves this means $|g'|\le\sigma\,\|f_2''\|\,u$ somewhere under the square.
4. **Kinks (Corollary 381.4).** For polygons, ${\rm TV}(f_2')$ is the sum of the upper curve's turning angles under the square. So kink squares sit where a vertex makes the gap's slope change sign. This is the mechanism of v379's kink squares.
5. **Checks:**
   - wobbling spirals: in 60 windows at stationary gaps, $H$ has a zero in 58, and all 58 polish to exact inscribed squares. The explicit test of Corollary 381.2 certifies 43 of them (every gap minimum, and the maxima of the strongest wobble). Without stationary gaps there are no zeros;
   - localization holds on all 508 aligned squares from v380, including squares straddling turns up to 12 apart. The bound is tight: the gap change equals the bound to $10^{-8}$;
   - polygon drawings: the zeros of $H$ are exactly the aligned squares, 178 of 178, one-to-one, with no extras.

## The reduction

Let $f_1<f_2$ be graphs on an interval $I$: $C^1$, or piecewise linear, with $|f_i'|\le\sigma<1$. A square stands on $f_1$ if its corners are $P_1=(x_1,f_1(x_1))$, $P_2=(x_1+u,f_1(x_1+u))$ with $u>0$, and $P_4=P_1+J(P_2-P_1)$, $P_3=P_2+J(P_2-P_1)$, where $J$ is the quarter turn upward. Write $s=f_1(x_1+u)-f_1(x_1)$, so $P_2-P_1=(u,s)$ and $J(u,s)=(-s,u)$.

**Theorem 381.1.** For each $x_1$ there is exactly one $u(x_1)\in(0,g(x_1)/(1-\sigma^2)]$ with $P_4$ on $f_2$, and $u$ is continuous. The square is inscribed (with $P_3$ on $f_2$ too) iff $H(x_1)=0$, where
$$H(x_1)=s-\big[f_2(x_1+u-s)-f_2(x_1-s)\big]=-\big[g(x_1+u)-g(x_1)\big]+R,$$
$$R=\int_0^{s}\big[f_2'(x_1+u-\tau)-f_2'(x_1-\tau)\big]\,d\tau,\qquad |R|\le|s|\cdot{\rm TV}\big(f_2';W\big),$$
and $W$ is $[x_1,x_1+u]$ widened by $|s|$ on the side $s$ points to.

*Proof.* $P_4=(x_1-s,f_1(x_1)+u)$ lies on $f_2$ iff $F(u)=f_1(x_1)+u-f_2(x_1-s(u))=0$. Here $F(0)=-g(x_1)<0$, and $F$ increases with slope at least $1-\sigma^2$ (one-sided slopes for polygons). So there is a unique root, at most $g/(1-\sigma^2)$. Then $P_3=(x_1+u-s,f_1(x_1+u)+u)$ lies on $f_2$ iff $f_1(x_1+u)+u=f_2(x_1+u-s)$. Subtracting the $P_4$ equation gives $H=0$.

Adding and subtracting $f_2(x_1+u)-f_2(x_1)$ and $f_1(x_1+u)-f_1(x_1)=s$ gives the second form. Then $R$ is the stated integral, and $|f_2'(a)-f_2'(b)|$ is at most the total variation between $a$ and $b$. ∎

**Corollary 381.2 (existence).** Write $\Delta(x)=g(x+u(x))-g(x)$ and $B(x)=|s(x)|\,{\rm TV}(f_2';W(x))\le\sigma u\,{\rm TV}$. If $\Delta(x_-)>B(x_-)$ and $\Delta(x_+)<-B(x_+)$ (or the reverse), with all corners for $x_1\in[x_-,x_+]$ inside $I$, then a square stands on $f_1$ with bottom-left corner in $(x_-,x_+)$. This follows from the intermediate value theorem applied to $H$.

**Corollary 381.3 (localization).** If a square stands on $f_1$ at $x_1$, then $|\Delta(x_1)|\le B(x_1)$. For $f_2\in C^2$ this gives $|g(x_1+u)-g(x_1)|\le\sigma\|f_2''\|u^2$, so $|g'(\xi)|\le\sigma\|f_2''\|u$ for some $\xi$ under the square. The gap is stationary to within the pitch-curvature product.

**Corollary 381.4 (kinks).** If $f_2$ is piecewise linear, ${\rm TV}(f_2';W)$ is the sum of its turning angles in $W$, and $R=0$ when $W$ contains no vertex of $f_2$. A square then needs $g(x_1+u)=g(x_1)$ exactly for the piecewise-linear $g$. That happens where a vertex of either curve makes $g'$ change sign. Every vertex pair across neighbouring turns can do this, which is why every small square of the drawings sits at a vertex (v379).

So Principle 380.2 is now a theorem: aligned squares exist where the gap rises and falls by more than the curvature error, and only where the gap is stationary up to that error. It holds for any pair of graphs, not just spirals, and for squares straddling non-neighbouring turns as well.

## Checks

**Existence at stationary gaps.** The wobbling family of v380 is $r=\Theta^{-1/2}(1+\frac\kappa\Theta\sin\omega\Theta)$. In a frame at each stationary point of the radial gap ($\Theta\in[8,40]$), the two turns are sampled as graphs. $H$ is scanned over bottom-left positions within three gaps of the stationary point, and each zero is polished on the true spiral by Newton's method.

| $\kappa,\omega$ | stationary gaps (min / max) | windows with a zero | polished to exact squares | certified by Cor. 381.2 |
|---|---|---|---|---|
| 0.3, 1 | 5 / 5 | 10 | 10 | 5 minima, 0 maxima |
| 0.7, 2 | 10 / 10 | 20 | 20 | 10 minima, 0 maxima |
| 0.5, 3 | 15 / 15 | 28 | 28 | 15 minima, 13 maxima |
| 0, – and 0.2, 2 | none | 0 | 0 | – |

The explicit test is conservative at gap maxima of the gentler wobbles: $H$ changes sign there, but the worst-case bound $B$ is too large to certify it. In the strongest wobble the two uncertified maxima also have no zero inside the window. Those squares are too large for it.

**Localization.** For each aligned square from v380's search, put the frame on the two turns it straddles. These are neighbours for 148 of the 508, and 2 to 12 turns apart for the rest. Then compare the gap change across the square with $B$. All 508 satisfy $|\Delta|\le B+|H|$, where $|H|\le2\times10^{-8}$ is the residual left by sampling. The largest ratio $(|\Delta|-|H|)/B$ is $1+4\times10^{-11}$. So the bound is attained: on these curves $f_2'$ is monotone across each square, and the variation bound is exact.

**Polygons.** For the drawings ($\lambda=0,1$ at $s=1/20$ depth 4, and $\lambda=0$ at $s=1/30$ depth 5), $H$ is computed exactly on the piecewise-linear turns. Frames are placed at 96 angles per turn, for turn gaps 1–3. The zeros with side below $s$ coincide with the exact aligned squares from the square counter: 56/56, 62/62 and 60/60, with nothing extra.

## What this changes

- **The comb picture is now proved except for one piece.** The general law (v380) bounds where squares can be. The two-curve theorem says exactly where aligned squares are. v379's selection rule governs the tilted ones. The open piece is the matching existence statement for tilted squares at their resonant radii. It needs a similar reduction for squares spread over several turns.
- **Polygon kink squares are now explained exactly.** That matters for the asymmetric drawing data of v374 and v378, all of which involve them.
- **A counting tool.** $H$ turns aligned squares into zeros of one real function per pair of turns, which can be counted and tracked along deformations. That is the bookkeeping the surgery rule (v378) and the path condition (v374) need.

[Code](../code/two_curve_v381.py) · [Tests](../tests/test_two_curve_v381.py) · [Data](../results/v381/two_curve_check.json)
