# v374. Asymmetric wild points: the gap path and the scale mismatch

v370 and v372 handle a wild point with a symmetry (a half-turn or a mirror), v371 extends them to curves joined to such a point without passing a square of side $s$, and v373 shows what such a path must sweep through. A generic wild point has no symmetry. The obvious plan is to deform it into a symmetric one and carry the protection back along the deformation.

This note shows that half of the plan works and pins down exactly where the other half fails. The deformation always exists, and if it is made small enough it cannot change anything at a fixed cutoff $s$. But Theorem 370.2 only protects cutoffs whose small squares sit deep inside the symmetric ball. Near a spiral point the small squares of side $s$ do not sit near the point. They fill a comb out to radius about $s^{a/(1+a)}$, far larger than $s$. So the plan closes exactly when a curve's "invisibility radius" at cutoff $s$ beats that comb radius, and nothing we know forces this.

**Results.**
1. **Gap path (Lemma 374.1).** Near a point where the curve is two radially monotone arms, the angular gap between the arms can be moved to $\pi$ through Jordan curves. The end of the path is half-turn symmetric near the point, and the path can be confined to any small ball.
2. **Local invisibility (Lemma 374.2).** If a Jordan curve has no square of side $s$, then every curve that differs from it only inside a small enough ball around $p$ has none either. So a small deformation at $p$ keeps the large degrees at cutoff $s$.
3. **Reduction (Theorem 374.3).** An asymmetric two-arm point inscribes a square as soon as, for one cutoff $s$, the curve symmetrized inside a ball of radius $\varepsilon$ has all its squares of side below $s$ deep inside that ball while $\varepsilon$ is still below the invisibility radius at $s$.
4. **Comb law (Claim 374.4, heuristic, checked exactly).** For two spiral arms $r=\theta^{-a}$ with gap $\pi\pm\varphi$, small squares of side below $s$ reach out to radius $r_{\rm comb}(s)\approx\big(s/(a\,g)\big)^{a/(1+a)}$ with $g=\pi-|\varphi|$. So the reduction needs an invisibility radius of order $s^{a/(1+a)}$, which is much larger than $s$.
5. **Exact checks:**
   - on eight spiral drawings, symmetric and asymmetric, the outermost small square sits within 2% of the predicted comb radius;
   - along an exact gap path at cutoff $1/20$, protection is lost and regained; both events are aligned squares at radius $6.4\,s$, which is the comb radius, not the wild point.

## Setting

Put $p$ at the origin. Say $p$ is a **two-arm point** of a Jordan curve $\gamma$ if, for some $R>0$:
- $\gamma\cap\bar B(p,R)$ is one arc $\alpha$ through $p$, and $\gamma$ is polygonal away from $p$ (hypothesis (H1) of Theorem 370.2);
- each of the two arms of $\alpha$ meets every circle $|z|=r$, $0<r\le R$, exactly once.

Then each arm is a graph over the radius: arm $j$ is $r\mapsto re^{i\theta_j(r)}$ with $\theta_j$ continuous on $(0,R]$. The arms meet only at $p$, so the gap $\delta(r)=\theta_2(r)-\theta_1(r)$ is never a multiple of $2\pi$, and we lift it to $(0,2\pi)$. The angles may wind without bound as $r\to0$. The spirals of v370 are two-arm points with $\delta\equiv\pi$.

## The gap path

**Lemma 374.1.** Let $p$ be a two-arm point and $0<\varepsilon<\varepsilon'\le R$. There is a uniformly continuous family $\gamma_\lambda$, $\lambda\in[0,1]$, of Jordan curves with $\gamma_0=\gamma$, all equal to $\gamma$ outside $B(p,\varepsilon')$, such that $\gamma_1$ is polygonal away from $p$ and $\gamma_1\cap\bar B(p,\varepsilon)$ is half-turn symmetric.

*Proof.* Take a continuous $\chi:(0,R]\to[0,1]$ equal to $1$ on $(0,\varepsilon]$ and to $0$ on $[\varepsilon',R]$. Move arm two to the angle
$$\theta_2^\lambda(r)=\theta_2(r)+\lambda\chi(r)\big(\pi-\delta(r)\big).$$
The new gap $\delta+\lambda\chi(\pi-\delta)$ is a convex combination of $\delta(r)$ and $\pi$, so it stays in $(0,2\pi)$. Each point of arm two moves along its own circle, so the new arm still meets every circle once and never meets arm one. Nothing changes for $r\ge\varepsilon'$, and the rest of $\gamma$ stays outside $\bar B(p,R)$. So every $\gamma_\lambda$ is a Jordan curve, and the family is uniformly continuous in $\lambda$.

At $\lambda=1$, arm two is $-(\text{arm one})$ on $\bar B(p,\varepsilon)$. That part is polygonal because arm one is. On the annulus $\varepsilon\le r\le\varepsilon'$ the moved arm is piecewise smooth, compact, and at positive distance from arm one. Replace it by a fine inscribed polygon with the same ends. With a fine enough mesh this stays in the annulus up to an arbitrarily small error, avoids arm one, and the straight homotopy to it passes only through Jordan curves. Append that homotopy to the path. ∎

The convex-combination step is the whole point of radial monotonicity: it lets the gap be interpolated pointwise without the arms ever touching. A point whose arms double back radially would need a different construction.

## Local invisibility

**Lemma 374.2.** Let $\gamma$ be a Jordan curve with no inscribed square of side $s$, and let $p\in\gamma$. There is $\varepsilon_{\rm inv}(s)>0$ such that no closed curve $\gamma'\subset\gamma\cup B(p,\varepsilon_{\rm inv}(s))$ has an inscribed square of side $s$.

*Proof.* Otherwise there are $\varepsilon_n\to0$ and squares $Q_n$ of side $s$ with corners in $\gamma\cup B(p,\varepsilon_n)$. Pass to a subsequence where the corners converge. Corners on $\gamma$ converge to points of $\gamma$, since $\gamma$ is closed. Corners in $B(p,\varepsilon_n)$ converge to $p\in\gamma$. The limit is a square of side $s$ inscribed in $\gamma$. ∎

Nothing is assumed about $\gamma'$ inside the ball, so the lemma covers every deformation supported there, including the gap path of Lemma 374.1 with $\varepsilon'<\varepsilon_{\rm inv}(s)$. By the path principle (Corollary 371.2), such a deformation keeps both large degrees at cutoff $s$.

The radius $\varepsilon_{\rm inv}(s)$ is not quantitative. It is at most the distance from $p$ to the free corner of any near-square of side $s$ with three corners on $\gamma$, and it can be tiny.

## The reduction

**Theorem 374.3.** Let $p$ be a two-arm point of $\gamma$. Suppose that for some $s>0$ and some $\varepsilon'<\min(R,\varepsilon_{\rm inv}(s))$, with $\varepsilon<\varepsilon'$, the curve $\gamma_1$ of Lemma 374.1 satisfies $s\le s_0(\gamma_1)$, where $s_0$ is the threshold of Theorem 370.2 at radius $\varepsilon$. (In particular, every square of side below $s$ on $\gamma_1$ or on its drawings has a corner in $B(p,\varepsilon/8)$.) Then $\gamma$ inscribes a square of side at least $s$.

*Proof.* If $\gamma$ has a square of side $s$ we are done. Otherwise $\varepsilon_{\rm inv}(s)$ is defined. Every curve on the gap path lies in $\gamma\cup B(p,\varepsilon')$, so by Lemma 374.2 none has a square of side $s$. The curve $\gamma_1$ satisfies Theorem 370.2 at radius $\varepsilon$: (H1) because it is polygonal away from $p$, (H2) by construction, and (H3) is automatic (v372, item 4). Since $s\le s_0(\gamma_1)$, Corollary 371.3 says that every Jordan curve joined to $\gamma_1$ without passing a square of side $s$ inscribes a square of side greater than $s$. The gap path, run backwards, is such a join to $\gamma$. ∎

So the asymmetric case reduces to one inequality between two radii at a single cutoff:
- the **confinement radius**: how far out the squares of side below $s$ reach on the symmetrized curve; Theorem 370.2 needs this below $\varepsilon/8$;
- the **invisibility radius** $\varepsilon_{\rm inv}(s)$ of the original curve.

If $8\,r_{\rm conf}(s)<\varepsilon_{\rm inv}(s)$ for one $s$, the point inscribes a square.

## The comb law

**Claim 374.4 (heuristic).** Near two spiral arms $r=\theta^{-a}$ with gaps $\pi\pm\varphi(r)$, the squares of side below $s$ reach out to
$$r_{\rm comb}(s)\approx\Big(\frac{s}{a\,g}\Big)^{a/(1+a)},\qquad g=\pi-|\varphi(r_{\rm comb})|.$$

*Derivation.* One turn of an arm changes the radius by about $2\pi a\,r^{1+1/a}$. The other arm splits that turn in the ratio of the two gaps, so the closest pair of neighbouring passes is about $a\,g\,r^{1+1/a}$ apart. At radius $r\gg s$ the passes are nearly straight on the scale $s$. Four corners on nearly parallel lines form a square only if two neighbouring lines are at most about $s$ apart: two corners on each of two lines force the separation to equal the side, and spreading the corners over more lines needs closer lines. So small squares stop where the closest spacing reaches $s$. This is not a proof: it ignores curvature corrections of relative size $s/r$, and it does not rule out squares that use the outer closure.

**Exact check.** Drawings use $a=1/2$, $\theta_0=3$, 24 vertices per turn, and $\varphi(r)=0.4\sin(5\log(1/r))$ scaled by $\lambda$ ($\lambda=0$ symmetric, $\lambda=1$ asymmetric). Every drawing reaches inside the comb radius. The predicted radius uses the fixed point of the formula above.

| $s$ | depth | $\lambda$ | small squares | aligned / reversed | outermost centre | predicted | centre / $s$ |
|---|---|---|---|---|---|---|---|
| 1/20 | 4 | 0 | 58 | 29 / 29 | 0.3127 | 0.3169 | 6.3 |
| 1/20 | 4 | 1 | 64 | 35 / 29 | 0.3286 | 0.3258 | 6.6 |
| 1/25 | 4 | 0 | 34 | 17 / 17 | 0.2906 | 0.2942 | 7.3 |
| 1/25 | 4 | 1 | 42 | 22 / 20 | 0.3009 | 0.2969 | 7.5 |
| 1/30 | 5 | 0 | 62 | 31 / 31 | 0.2755 | 0.2769 | 8.3 |
| 1/30 | 5 | 1 | 68 | 33 / 35 | 0.2790 | 0.2782 | 8.4 |
| 1/40 | 6 | 0 | 66 | 33 / 33 | 0.2492 | 0.2515 | 10.0 |
| 1/40 | 6 | 1 | 70 | 29 / 41 | 0.2556 | 0.2570 | 10.2 |

Every outermost centre is within 2% of the prediction. The ratio centre/side grows as $s$ shrinks, as $s^{-1/(1+a)}=s^{-2/3}$ predicts. No square lies in the crossed sector.

**What this means for the reduction.** For a spiral point, $r_{\rm conf}(s)\approx r_{\rm comb}(s)$, so Theorem 374.3 needs
$$\varepsilon_{\rm inv}(s)\gtrsim 8\Big(\frac{s}{a\,g}\Big)^{a/(1+a)},$$
that is, an invisibility radius larger than $s$ by a factor growing like $s^{-1/(1+a)}$. A curve can miss squares of side $s$ by much less than that. So the deformation that symmetrizes the point is visible at exactly the cutoffs that Theorem 370.2 would protect.

## Where the events happen

The full gap path (support the whole spiral, out to radius $0.58$) on the depth-4 drawing at $s=1/20$:

| $\lambda$ | 0 | 1/4 | 1/2 | 3/4 | 1 |
|---|---|---|---|---|---|
| small aligned / reversed | 29 / 29 | 31 / 27 | 32 / 27 | 33 / 27 | 35 / 29 |
| large graceful parity | 1 | 1 | 0 | 1 | 1 |

Bisection finds two parity events, at $\lambda\in(2567/8192,\,5135/16384)$ and $\lambda\in(4975/8192,\,9951/16384)$. At each, one aligned square shrinks through side $1/20$ (from $0.0500003$ to $0.0499998$ or $0.0499999$). Its centre is at $6.4\,s\approx0.32$, the comb radius, on the outer turns. None of the events is near the wild point.

These are depth-truncated polygons, so they are not counterexamples: every curve here inscribes squares, and v373 says each parity change is a square of side $1/20$ on an intermediate curve. What the data shows is where those squares are. Symmetrizing a spiral point rotates its comb, and the squares the rotation creates and destroys live at the comb radius.

## What would close asymmetric points

Theorem 374.3 isolates two ways forward:
- **A quantitative invisibility radius.** Show that if $\gamma$ has no square of side $s$, then near-squares of side $s$ with a corner near $p$ stay away from $p$ by at least $C\,r_{\rm comb}(s)$. This is a "near-square forces a square" statement at the comb scale.
- **A deformation that leaves the comb alone.** Instead of rotating whole turns, change the curve only below the pass spacing, so the confinement radius of the deformed region is set by the deformation, not by the comb.

Either would make two-arm points as tractable as the symmetric ones. Points whose arms are not radially monotone need a different gap path altogether.

[Code](../code/asymmetric_spiral_v374.py) · [Tests](../tests/test_asymmetric_spiral_v374.py) · [Data](../results/v374/comb_scan.json)
