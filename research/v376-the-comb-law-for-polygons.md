# v376. The comb law for polygons

[v375](v375-the-comb-law.md) proved the comb law for smooth spiral arcs. The theorems that use it, 370.2 and 374.3, need curves that are polygonal away from the wild point. This note carries the law over to polygons inscribed in the arcs. The cost is one error term: how far an edge sags away from the arc it spans.

The same error term answers a question raised in review. A fixed number of vertices per turn eventually makes the drawing cross itself, because the edges sag by more than the gap between turns. The comb law with sag gives a simple rule: refine the mesh so the sag stays below the turn spacing. Then the polygon is simple and obeys the comb law at every depth.

**Results.**
1. **Sag bound (Lemma 376.0).** Every point of an edge spanning $\Delta\Theta$ in lifted angle lies within $h=\frac{\Delta\Theta^2}8\sup|\gamma''|+\nu$ of the arc, where $\nu$ is the vertex rounding.
2. **Polygonal comb law (Theorem 376.1).** For a square of side $t$ on the polygon, $t\ge d_g(\rho^*)-2h$, where $\rho^*$ is the smallest radius of the corresponding points on the arc. The smallness conditions are those of v375, with $t$ replaced by $t+\sqrt2h$ and an extra $2\mu$ (twice the mesh) in the straightness condition.
3. **Confinement and simplicity (Corollaries 376.2, 376.3).** Every small square has a corner within $d_g^{-1}(t+2h)+h$. An inscribed polygon whose edges all satisfy $2h<d_g$ is simple.
4. **Mesh rule (Corollary 376.4).** Both hold at every depth if the mesh shrinks like $r^{1/(2a)}$. At 24 vertices per turn the criterion first fails at depth 28 (symmetric) and 24 (asymmetric). The asymmetric drawing actually crosses itself somewhere between depths 70 and 75.
5. **Exact checks:**
   - on the depth-4 drawings, every corner lies within $h$ of its arc point, and $h$ is within 7–10% of the true largest sag;
   - on the eight v374 drawings, 414 of the 464 small squares meet the hypotheses, and every square with corners on the arms (448) satisfies the bound;
   - with the refined mesh, the asymmetric drawing is simple at depths 75 and 100 (10 827 and 16 597 vertices), and $2h/d_g\le0.15$ all the way in.

## Setting

Keep the notation of v375. The arc $\alpha$ has arms $\gamma_1(\Theta)=\Theta^{-a}e^{i\Theta}$ and $\gamma_2(\Theta)=\Theta^{-a}e^{i(\Theta+\pi+\varphi(\Theta^{-a}))}$ for $\Theta\ge\theta_0$, with $|\varphi|\le\bar\varphi<\pi$, $|r\varphi'|\le L$, $g=\pi-\bar\varphi$, and $d_g$, $K$, $\varepsilon_1$, $\tau$, $\eta$ as there. Now also assume $|r^2\varphi''|\le L_2$.

An **inscribed polygon** $P$ has, on each arm $j$, vertices $\hat z_{j,k}$ with $|\hat z_{j,k}-\gamma_j(\Theta_k)|\le\nu$, where $\theta_0=\Theta_0<\Theta_1<\cdots$ and $\Theta_{k+1}-\Theta_k\le\mu$. Consecutive vertices on an arm are joined by edges. (The drawings also have an inner connector and an outer closure. A square with a corner on those is not covered here.)

## The sag bound

**Lemma 376.0.** Let $e$ be the edge from $\hat z_{j,k}$ to $\hat z_{j,k+1}$, and let $x=(1-\sigma)\hat z_{j,k}+\sigma\hat z_{j,k+1}$. Put $x^*=\gamma_j(\Theta_k+\sigma\Delta\Theta)$ with $\Delta\Theta=\Theta_{k+1}-\Theta_k$. Then
$$|x-x^*|\le h_e:=\frac{\Delta\Theta^2}8\sup_{[\Theta_k,\Theta_{k+1}]}|\gamma_j''|+\nu,$$
and on arm $j$
$$\frac{|\gamma_j''(\Theta)|}{r}\le\Big(1+\frac1a\Big)K^2+(1+KL)^2+2K(1+KL)+K^2\Big(L_2+\big(1+\tfrac1a\big)L\Big),\qquad K=K(r)$$
(with $L=L_2=0$ on arm 1).

*Proof.* For $f(u)=\gamma_j(\Theta_k+u\Delta\Theta)$, the linear interpolation error is $(1-\sigma)f(0)+\sigma f(1)-f(\sigma)=-\int_0^1G(\sigma,u)f''(u)\,du$. Here $G\ge0$ is the Green's function of $-d^2/du^2$ on $[0,1]$, with $\int G=\sigma(1-\sigma)/2\le\frac18$. So the error is at most $\frac18\Delta\Theta^2\sup|\gamma_j''|$, and rounding adds at most $\nu$.

For the curvature, write $\gamma_j=re^{i\Psi}$. Then $\gamma_j''=e^{i\Psi}\big(r''-r\Psi'^2+i(2r'\Psi'+r\Psi'')\big)$, with $r'=-Kr$, $r''=(1+\frac1a)K^2r$, $\Psi'=1+d\varphi/d\Theta\in[1-KL,1+KL]$, and $|\Psi''|=|\varphi''r'^2+\varphi'r''|\le K^2(L_2+(1+\frac1a)L)$. ∎

For the drawings, $\mu=\pi/12$, $\nu<10^{-7}$, and the sag is about $\mu^2r/8$, roughly $0.003$ at the radii where the small squares sit.

## The polygonal comb law

**Theorem 376.1.** Let $Q$ be a square of side $t$ with corners $x_i$ on edges $e_i$ of $P$, and let $x_i^*$ be their arc points from Lemma 376.0. Put
$$h=\max_ih_{e_i},\qquad \rho^*=\min_i|x_i^*|,\qquad t'=t+\sqrt2h,$$
and let $\omega$ bound the angle by which rounding turns any edge, so $\omega\le\arcsin(2\nu/\text{edge length})$. Suppose the v375 conditions hold with $t'$ in place of $t$ and the straightness condition strengthened to
$$\big(\eta(t',\rho^*)+2\mu\big)(1+\varepsilon_1)+\tau+\omega<\frac\pi2.$$
Then $t\ge d_g(\rho^*)-2h$.

*Proof.* Run the proof of Theorem 375.1 on the arc points, with neighbours defined through them. The arc points are at most $\sqrt2t+2h=\sqrt2t'$ apart, so Lemma 375.A applies with $t'$. Two corners whose arc points are not neighbours therefore satisfy $|x_i-x_j|\ge|x_i^*-x_j^*|-2h\ge d_g(\rho^*)-2h$.

For Lemma 375.B, three pairwise-neighbour corners have arc points in a $\Theta$-window of length at most $\eta(t',\rho^*)$. Their edges span $\Theta$-ranges of length at most $\mu$ that contain those points, so all edges between two of the corners lie within a window of length $\eta+2\mu$. The vector from one corner to another along $P$ is a sum of positive multiples of edge vectors. Each edge vector is a chord of the arc over its range, turned by at most $\omega$. So it points into the tangent cone of the arc over the window, widened by $\omega$. As in v375, that cone has opening less than $\pi/2$, so no two such vectors are perpendicular. The rest of the proof is unchanged. ∎

**Localization.** Nothing in either proof uses the arc beyond the largest radius met by the corners' edges. So $K_R$, $\varepsilon_1$ and $\tau$ may be evaluated at that radius $R'$ instead of $R$. The checks below do this, and the same holds for Theorem 375.1.

**Corollary 376.2 (confinement).** Under the conditions of Theorem 376.1, every square on $P$ of side $t$ has a corner within $d_g^{-1}(t+2h)+h$ of $p$.

*Proof.* $|x_i|\ge|x_i^*|-h$, and $t\ge d_g(\rho^*)-2h$ gives $\rho^*\le d_g^{-1}(t+2h)$. ∎

**Corollary 376.3 (simplicity).** Suppose every pair of edges $e,e'$ on the two arm paths satisfies $2\max(h_e,h_{e'})<d_g(r)$, where $r$ is the smaller radius of their arc points, together with the straightness condition with $t'=2h$. Then the two arm paths of $P$ are simple and disjoint.

*Proof.* Suppose non-adjacent edges meet at a point $x=y$. Their arc points are at most $2h$ apart. If they are not neighbours, Lemma 375.A (with $\sqrt2t=2h$) gives $|x^*-y^*|\ge d_g(r)>2h$, which is impossible. If they are neighbours, both edges are on one arm, inside a short window. The path along $P$ from $x$ back to $y$ contains at least one whole edge, and it is a sum of positive multiples of vectors in an open cone of opening less than $\pi/2$. Such a sum cannot vanish. ∎

## How fine the mesh must be

**Corollary 376.4 (mesh rule).** To leading order $h\approx\mu^2r/8$ and $d_g\approx agr^{1+1/a}$. So $2h\le\theta\,d_g$ holds when
$$\mu\lesssim2\sqrt{\theta\,g\,a}\;r^{1/(2a)}.$$
A mesh that shrinks like $r^{1/(2a)}$ ($\propto r$ for $a=1/2$) keeps both corollaries at every depth, with confinement radius $d_g^{-1}(t)(1+O(\theta))$. A fixed mesh does not. Near $p$ the turn spacing $\propto r^{1+1/a}$ falls below the sag $\propto r$, and then the edges of neighbouring turns can cross.

**At 24 vertices per turn** ($\mu=\pi/12$, $a=1/2$), the criterion $2h<d_g$ first fails at depth 28 for $\varphi=0$ and at depth 24 for the asymmetric $\varphi$. The criterion is only sufficient. The asymmetric drawing is still simple at depth 70 by bisection, and it crosses itself by depth 75, as found in review: edges 830 and 959 cross. The symmetric drawing is still simple at depth 75. Its two arms have vertices on common rays, which this argument does not use.

**With the refined mesh** $\mu(\Theta)=\frac\pi{12}\sqrt{3/\Theta}$, which is proportional to $r$ and agrees with the drawings at $\Theta_0=3$, the ratio $2h/d_g$ stays below $0.066$ ($\varphi=0$) and $0.146$ (asymmetric) out to depth 100. The asymmetric drawing is simple at depths 75 and 100.

## What this settles

- **v375, Corollary 375.3, items 2 and 3.** These now hold for inscribed polygons whose mesh obeys the rule above, which can always be arranged when drawings are built for a proof. One correction to item 2: Lemma 369.6's threshold is limited by the turn spacing as well as by the edge length, since edges on neighbouring turns are non-adjacent. What the comb law removes is the dependence on edge length. It also places the threshold at radius $R/8$ rather than at the $R/16$ where Lemma 369.6 is applied.
- **The v374 drawings.** At depth at most 6 they are well inside the safe range, and their parity data stand.
- **Still open.** The lower half of the comb law, and the path condition (b) of Theorem 374.3. Both are unchanged by this note.

## Exact checks

**Arc points and sag.** On the depth-4 drawings at $s=1/20$, every corner of every small square on the arms lies within $h$ of its arc point, at most $0.90h$ (224 corners for $\varphi=0$, 248 asymmetric). The largest true sag of any edge, sampled at 41 points, is at most $0.93$ of its bound $h_e$.

**Theorem 376.1 on the v374 drawings.** Here $L=2\lambda$ and $L_2=\sqrt{104}\,\lambda$, since $|r^2\varphi''|=|{-10}\sin u+2\cos u|\le\sqrt{104}$ for $\varphi=0.4\sin u$, $u=5\log(1/r)$. The hypotheses use the localized $R'$.

| $s$ | depth | $\lambda$ | squares | on the arms | meet hypotheses | bound holds | largest $h$ | largest $\eta$ |
|---|---|---|---|---|---|---|---|---|
| 1/20 | 4 | 0 | 58 | 56 | 56 | all | 0.0033 | 0.41 |
| 1/20 | 4 | 1 | 64 | 62 | 34 | all | 0.0046 | 0.98 |
| 1/25 | 4 | 0 | 34 | 32 | 32 | all | 0.0030 | 0.35 |
| 1/25 | 4 | 1 | 42 | 40 | 34 | all | 0.0039 | 0.83 |
| 1/30 | 5 | 0 | 62 | 60 | 60 | all | 0.0028 | 0.32 |
| 1/30 | 5 | 1 | 68 | 66 | 66 | all | 0.0034 | 0.73 |
| 1/40 | 6 | 0 | 66 | 64 | 64 | all | 0.0024 | 0.26 |
| 1/40 | 6 | 1 | 70 | 68 | 68 | all | 0.0029 | 0.62 |

The two off-arm squares in each drawing touch the inner connector. Where the hypotheses hold, the side exceeds the bound $d_g(\rho^*)-2h$ by a factor of at least 1.59. The asymmetric rows at $s=1/20$ and $1/25$ fail the straightness condition for some squares, because the extra $2\mu=\pi/6$ is large at 24 vertices per turn. A finer mesh would remove this.

[Code](../code/polygonal_comb_v376.py) · [Tests](../tests/test_polygonal_comb_v376.py) · [Data](../results/v376/polygonal_comb_check.json)
