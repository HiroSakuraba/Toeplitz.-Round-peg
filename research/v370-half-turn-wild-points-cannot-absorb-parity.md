# v370. A half-turn-symmetric wild point cannot absorb the odd square parity

[v369](v369-relabeling-dichotomy-and-collapse-localization.md) showed that a noncollapse proof can fail only when the odd parity of approximating polygons collapses into one point. This note tests that at a single wild point, on a curve that is polygonal everywhere else. The wild points used here lie outside the sufficient conditions of the two recent general results. The result is a theorem for one geometric kind of wild point, plus exact data on what happens without that structure.

**Results.**
1. **Sector degrees of polygons (Lemma 370.1).** Every Jordan polygon has mod-two degree 1 for squares whose corners occur along the curve in the square's own orientation (*aligned*), and 0 for the opposite order (*reversed*).
2. **Half-turn wild points (Theorem 370.2).** Suppose a Jordan curve is polygonal away from one point $p$, and near $p$ it is a single arc invariant under the half-turn about $p$. Then it inscribes a square of side at least an explicit $s_0>0$, however wild the arc is at $p$.
3. **The mechanism.** The half-turn reverses the arc while preserving square orientation, so it trades small aligned squares for small reversed ones. Polygons have sector degrees $(1,0)$, so small squares cannot make the combined large count even. The exact counts match this pairing at every depth and cutoff tested.
4. **Without symmetry the parity is not local.** Changing the inner connecting arc does change the small-square parity, through squares that straddle the innermost disk: 77 versus 78 aligned squares at one depth. And even with symmetry, the small *aligned* count alone is odd at many depths. So a lemma stated for one sector only would fail here, while the two-sector version holds.

## Which wild points are new

Asano and Ike's theorem requires smooth approximations whose primitives of the Liouville form, $\int y\,dx$ along the parametrization, converge uniformly on compact sets (their Theorem 1.1). It covers every rectifiable curve and every locally monotone curve. A follow-up preprint covers curves of finite $p$-variation for some $p<2$.

Take a radially monotone double spiral at $p$ with turn radii $r_k$. Each turn adds roughly the area of its disk, $\mp\pi r_k^2$, to the primitive, with the same sign every turn. So the primitive converges only if $\sum r_k^2<\infty$. Finite $p$-variation for some $p<2$ needs $\sum r_k^p<\infty$. For radii $r\approx\theta^{-a}$, both conditions fail exactly when $a\le1/2$. In the repo's power spirals $t\,e^{i/t^p}$ that is $p\ge2$. The spiral centre is not locally monotone either.

These are sufficient conditions, and I have not checked every other criterion in the literature. The computation of the primitive is an estimate of the turn contributions, not a proof that their condition fails. Any single curve built here has explicit squares away from $p$, so its existence is not in question. The content is a theorem for the whole class of such points.

## Lemma 370.1: sector degrees of polygons

For a Jordan polygon $C$, let $X_\pm/T$ be the quotient row spaces of v369 for the two orientations. The *aligned* copy uses the curve's own orientation. Remove the $T$-invariant set $N$ of rows whose corners have diameter less than $\tau/2$. Lemma 369.6 makes $\partial N$ zero-free, and zeros on $\partial X$ are collapsed, so they lie in $N$. The degrees $\pi_{\rm al}(C)$ and $\pi_{\rm rev}(C)$ over the two quotients are therefore defined and independent of the choice of $N$.

**Lemma 370.1.** $\pi_{\rm al}(C)=1$ and $\pi_{\rm rev}(C)=0$ for every Jordan polygon $C$.

*Proof.* A simple polygon with at least four vertices has an ear: a vertex $v$ whose triangle with its neighbours contains no other part of the polygon. Move $v$ straight to the midpoint of its neighbours. The two moving edges stay inside the ear, so every intermediate polygon is simple, and the final one has a straight angle at $v$, which can be deleted. Lemma 369.6 applies throughout, including at straight angles, since adjacent collinear edges carry no square. The smallest non-adjacent edge distance stays positive along each motion, so one fixed $N$ works for the whole motion and the degrees do not change.

Repeating this reaches a triangle. For a generic acute triangle the exact census gives three transverse aligned squares and no reversed ones. ∎

The tests compare the exact counter with the complete census on two triangles, a convex quadrilateral and a dart, and check three spiral drawings with 21, 29 and 41 vertices. Every case has odd aligned and even reversed counts: for example, 21 aligned and 10 reversed squares on the 21-vertex drawing. Classical results give odd square counts for generic smooth curves. This is the polygon version, split by sector, with an elementary proof.

## Theorem 370.2: half-turn wild points

Let $g$ be the half-turn about $p$.

**Theorem 370.2.** Let $\gamma$ be a Jordan curve and $p\in\gamma$ such that:
- **(H1)** $\gamma\setminus B(p,r)$ consists of finitely many segments for every $r>0$;
- **(H2)** for some $R>0$, $\gamma\cap\bar B(p,R)$ is a single arc $\alpha$ through $p$ with $g(\alpha)=\alpha$;
- **(H3)** there are points $a_n\to p$ on $\alpha$ such that replacing the subarc of $\alpha$ between $a_n$ and $g(a_n)$ by the segment $[a_n,g(a_n)]$ gives a Jordan curve $D_n$, and the removed subarcs shrink to $p$.

Choose $s_0>0$ with two properties: every square of side below $s_0$ on $\gamma$ or on any $D_n$ has a corner in $B(p,R/8)$, and $R/8+s_0\sqrt2<R$. (By (H1), the polygonal part outside $B(p,R/16)$ has a positive non-adjacent edge distance, so Lemma 369.6 provides such an $s_0$.) Then $\gamma$ inscribes a square of side at least $s_0$.

*Proof.* Each $D_n$ is a Jordan polygon, by (H1) and (H3). If $D_n$ has a square of side exactly $s_0$, keep it. Otherwise each sector degree splits into the part with side below $s_0$ and the part above it, because the cutoff and the collapse tube are zero-free.

Every small square of $D_n$ lies in $B(p,R)$, where $D_n$ is the $g$-invariant arc obtained from $\alpha$. Reversing an arc reverses the order in which a square's corners occur along it, while $g$ preserves the square's orientation. So $g$ maps small aligned rows homeomorphically onto small reversed rows, and the residual transforms by a linear isomorphism. Hence the small aligned and small reversed degrees are equal mod 2. A square fixed by $g$ would have to lie in both sectors, which is impossible, so the pairing has no fixed points.

By Lemma 370.1, the large aligned plus large reversed degree is $1+0-2\cdot(\text{small})\equiv1$. So $D_n$ has a square of side at least $s_0$. As $n\to\infty$, $D_n\subset\gamma\cup B(p,\varepsilon_n)$ with $\varepsilon_n\to0$. A subsequence of these squares converges to a square of side at least $s_0$ with corners on $\gamma$. ∎

The argument never controls the small squares themselves. They may number in the hundreds and carry odd parity in each sector. It only uses that they come in pairs across sectors.

## The exact computation

The drawings are double spirals: arm one at radius $\theta^{-1/2}$ from $\theta=3$ inward, with 24 rational vertices per turn, arm two its exact half-turn image, a connector through $p$, and a polygonal outer closure. The counter finds every square with each corner inside an edge and a nonsingular contact system. numpy prunes, and every reported square is recomputed exactly. It reproduces the full census on the polygons above. No degenerate candidates occurred.

| depth (half-turns) | vertices | side $<1/10$: aligned, reversed | side $<1/20$ | side $<1/40$ |
|---|---|---|---|---|
| 2 | 53 | 1, 1 | 0, 0 | 0, 0 |
| 3 | 77 | 9, 9 | 5, 5 | 0, 0 |
| 4 | 101 | 34, 34 | 29, 29 | 0, 0 |
| 5 | 125 | 60, 60 | 53, 53 | 9, 9 |
| 6 | 149 | 87, 87 | 78, 78 | 33, 33 |
| 7 | 173 | 132, 132 | 104, 104 | 57, 57 |
| 8 | 197 | 209, 209 | 159, 159 | 81, 81 |

The tests also check that the half-turn maps each small aligned square exactly onto a small reversed one, and each crossed square onto a crossed one. On the 21-, 29- and 41-vertex drawings the full census confirms that the large graceful count is odd.

**What symmetry is doing.** At depth 6 with cutoff $1/20$, a bent connector through a point off the diameter gives 77 aligned and 78 reversed small squares. Their total is odd, so on that drawing the large graceful count is *even*. Moving the inner arc changes the parity through squares with a corner on it. Only a symmetric connector guarantees the pairing. And the small aligned count is odd at depths 3, 4, 5 and 8 for cutoff $1/20$. So "small squares are even in the odd sector" is false here, while the version summed over both graceful sectors always holds.

## What this suggests

The cancellation comes from one property: an involution of the local arc that reverses its direction and preserves square orientation. Among isometries fixing $p$, only the half-turn does both. The natural next question is whether approximate half-turn symmetry suffices. Suppose the arc near $p$ is within $\delta(r)$ of its half-turn image at scale $r$. Then the v368 margins bound how far the pairing can fail. A rate $\delta(r)=o(r)$ might keep the two small-sector degrees equal at every depth. That would replace an exact symmetry by an asymptotic one, and would include wild points with no symmetry at all.

[Code](../code/wild_point_parity_v370.py) · [Tests](../tests/test_wild_point_parity_v370.py) · [Scan](../results/v370/spiral_small_square_scan.json)
