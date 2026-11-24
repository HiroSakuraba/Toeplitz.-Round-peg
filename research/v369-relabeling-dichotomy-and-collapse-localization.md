# v369. The relabeling dichotomy and where collapse can happen

[v368](v368-certificate-margin-and-one-stage-criterion.md) left a target: on every Jordan curve, either some three-stretch certificate is isolated, or some square family comes back to itself under a change of starting corner. This note proves a version of that dichotomy for the exact zero set of the curve itself, with no approximations, and uses it to localize any failure of noncollapse.

**Results.**
1. **Invariant components give squares (Theorem 369.2).** Suppose a connected component of the exact zero set is carried to itself by a nontrivial relabeling. Then the curve has a square with two corners at least $1/3$ apart in parameter, and side at least $m_\gamma/\sqrt2$. [v362](v362-relabeling-path-noncollapse.md) proved this along approximating paths. Here it holds directly on the curve.
2. **The residual is equivariant (Lemma 369.3).** Relabeling acts freely on rows and acts on the square residual by an explicit linear map $L$ with $L^4=I$. So the residual is a section of a vector bundle over the quotient, and each square in a sector is one zero there.
3. **Otherwise the zero set splits (Proposition 369.4).** If no component is invariant, the four circles of collapsed rows lie in four distinct components. The zero set then splits into clopen pieces: one piece $V$ around one collapse circle, its three relabelings, and a compact set $R$ of genuine squares.
4. **Collapse is localized (Theorem 369.5).** The mod-two degree $d_R$ of the residual near $R$ is the same for every curve close to $\gamma$. If $d_R=1$, every nearby curve has a square with separated corners, so noncollapse holds. If $d_R=0$, any odd square parity of nearby polygons is carried by squares near one collapse piece. In that case either the curve has squares of every size up to some $s_0>0$, or those squares genuinely collapse to a point.
5. **Polygons have no small squares (Lemma 369.6).** A polygon with $\tau$ the smallest distance between non-adjacent edges has no inscribed square of side below $\tau/\sqrt2$.

The checks are exact: the identity for $L$, freeness, the one-third bound on every balanced row with denominator 24, the invariant family on the symmetric octagon, the single square of the v247 decagon, and the small-square bound on random star polygons.

## Rows and the zero set

Fix a registration $\gamma:\mathbb R/\mathbb Z\to\mathbb R^2$ of a Jordan curve. A **row** is a lifted quadruple $t_0\le t_1\le t_2\le t_3\le t_0+1$ modulo a common integer shift. The rows form $X\cong S^1\times\Delta^3$, a compact 4-manifold with boundary; the boundary consists of the rows with a zero gap. Write $p_i=\gamma(t_i)$, $v=p_1-p_0$, $J$ for the counterclockwise quarter turn, and
$$
G(x)=(r_1,r_2)=(p_2-p_1-Jv,\ p_3-p_0-Jv)\in\mathbb R^4 .
$$
A zero is a square whose corners appear along the curve in their square order, positively oriented. Reversing $J$ treats the reflected sector the same way. Squares whose corners appear in crossed order are outside this note.

Let $\hat Z=G^{-1}(0)$. A zero with two equal consecutive times has two equal adjacent corners, so all four corners are equal. Hence the zeros on $\partial X$ are exactly the **collapsed rows**. They form four circles $\Delta_0,\dots,\Delta_3$, according to which gap equals one. Every other zero is a genuine square.

The relabeling is $T(t_0,t_1,t_2,t_3)=(t_1,t_2,t_3,t_0+1)$. It preserves $\hat Z$ and permutes the $\Delta_i$ cyclically.

## Invariant components

**Lemma 369.1.** Let $K\subset X$ be connected, and suppose $x$ and $T^kx$ lie in $K$ for some $k\in\{1,2,3\}$. Then some row of $K$ has two times at circular distance at least $1/3$.

*Proof.* If a gap exceeds $1/2$ it is unique. On the open set of rows that have such a gap, its index is locally constant. If every row of $K$ had one, the index would be constant on $K$. But $T$ lowers it by one modulo four, so it differs at $x$ and $T^kx$. Hence some row of $K$ has every gap at most $1/2$. Such a row has a pair at distance at least $1/3$, by the counting argument of v362. ∎

The tests check the last step on every balanced row with denominator 24, along with v362's sharp example.

**Theorem 369.2.** If a component $K$ of $\hat Z$ satisfies $T^kK=K$ for some $k\in\{1,2,3\}$, then $\gamma$ has an inscribed square with two corners at parameter distance at least $1/3$. Its side is at least $m_\gamma/\sqrt2$, where $m_\gamma=\min_{d(s,t)\ge1/3}|\gamma(s)-\gamma(t)|$.

*Proof.* Apply Lemma 369.1 to $K$. The resulting row is an exact zero with two distinct corners, so it is a genuine square. Those two corners are adjacent or opposite. ∎

On the quarter-turn-symmetric octagon with equal times, $\gamma(t+1/4)=J\gamma(t)$ exactly. So the rows $(t,t+1/4,t+1/2,t+3/4)$ form a connected family of exact squares, and $T$ maps the family to itself. This is the circle's mechanism, and Theorem 369.2 applies to it.

## Equivariance

**Lemma 369.3.** $T$ acts freely on $X$. For any four points,
$$
G(Tx)=L\,G(x),\qquad L(r_1,r_2)=\big(r_2-(I+J)r_1,\ -Jr_1\big),\qquad L^4=I .
$$

*Proof.* If $Tx$ equalled $x$ up to an integer shift, then $t_1=t_0$, $t_2=t_1$, $t_3=t_2$ and $t_0+1=t_3$, which is impossible. For the identity, write $p_2=p_1+Jv+r_1$ and $p_3=p_0+Jv+r_2$. The relabeled row has first side $v'=Jv+r_1$. Its residuals are $p_3-p_2-Jv'=r_2-r_1-Jr_1$ and $p_0-p_1-Jv'=-Jr_1$. The relation $L^4=I$ follows from $T^4=\mathrm{id}$. The tests check both identities exactly on random rational quadruples, in both orientations. ∎

So $G$ is a section of the bundle $(X\times\mathbb R^4)/T$ over the compact manifold $X/T$, with $T$ acting on the fibre by $L$. One geometric square in the sector is exactly one zero in $X/T$. Its mod-two degree over any region of $X/T$ with zero-free boundary is well defined and does not change under small perturbations. In $X$ itself every count is a multiple of four, which is why mod-two degrees must be taken on the quotient.

## When no component is invariant

**Proposition 369.4.** Suppose no component of $\hat Z$ is invariant under $T$, $T^2$ or $T^3$. Let $K_i$ be the component containing $\Delta_i$. Then:
- the $K_i$ are pairwise distinct, and $TK_i=K_{i+1}$;
- for every $\eta>0$ there is a clopen $V\subset\hat Z$ with $K_0\subset V$, $V$ within $\eta$ of $K_0$, and $V,TV,T^2V,T^3V$ pairwise disjoint;
- the set $R=\hat Z\setminus\bigcup_jT^jV$ is compact, $T$-invariant, and consists of genuine squares.

*Proof.* If $K_i=K_j$ with $i\ne j$, then $T^{j-i}K_i=K_j=K_i$ would be invariant. In a compact Hausdorff space, components equal quasi-components. So $K_0$ is the intersection of the clopen sets containing it. By compactness there is a clopen $U$ with $K_0\subset U$, $U$ within $\eta$ of $K_0$, and $U$ disjoint from the compact set $\bigcup_{j=1}^3T^jK_0$. Put $V=U\setminus\bigcup_{j=1}^3T^{-j}U$. It is clopen, it contains $K_0$, and $V\cap T^jV=\varnothing$: a point $x=T^jy$ with $y\in V$ has $x=T^jy\notin U$. Collapsed rows lie in $\bigcup_jT^jK_0$, so $R$ contains none. ∎

**Corollary.** If $K_0\ne\Delta_0$, then $\gamma$ has inscribed squares of every side in $[0,s_0]$, for some $s_0>0$. Indeed, the side $|p_1-p_0|$ is continuous on the connected set $K_0$, vanishes on $\Delta_0$ and is positive somewhere.

## Localizing collapse

Assume the situation of Proposition 369.4. Choose $T$-invariant open sets in $X$:
- $O=\bigsqcup_jT^jO_V$, containing $\bigcup_jT^jV$, with each $T^jO_V$ near $T^jV$;
- $P$, with closure a compact PL 4-manifold inside the interior of $X$, containing $R$;
- with $\bar O\cap\bar P=\varnothing$ and $\partial O\cup\partial P$ disjoint from $\hat Z$.

These exist because the pieces are clopen and pairwise at positive distance, and $R$ avoids $\partial X$. Let $d_R$ be the mod-two degree of $G$ over $\bar P/T$.

**Theorem 369.5.** There is $\varepsilon>0$ such that for every registered Jordan curve $\gamma'$ with $\sup_t|\gamma'(t)-\gamma(t)|<\varepsilon$:
1. every zero of $G_{\gamma'}$ lies in $O\cup P$, and the degree of $G_{\gamma'}$ over $\bar P/T$ equals $d_R$;
2. if $d_R=1$, then $\gamma'$, and $\gamma$ itself, have a square whose row lies in $\bar P$. Since $\bar P$ avoids collapsed rows, the corners of such a square are separated by a fixed parameter distance $\kappa>0$, uniformly in $\gamma'$;
3. if $\gamma'$ is a polygon with $\tau'$ as in Lemma 369.6, and $N$ is the $T$-invariant set of rows whose four corners have diameter less than $\tau'/2$, then $\deg(G_{\gamma'},O_V\setminus N)=\pi(\gamma')+d_R \pmod 2$. Here $\pi(\gamma')$ is the degree of $G_{\gamma'}$ over $(X\setminus N)/T$.

*Proof.* Each residual coordinate changes by at most $4\delta$ when every point moves by at most $\delta$. On the compact set $X\setminus(O\cup P)$, $|G|$ is bounded below, so (1) holds for small $\varepsilon$, by homotopy. For (2), an odd degree forces a zero in $\bar P$. Apply this to $\gamma'$, and to $\gamma$ itself at $\delta=0$. Rows in $\bar P$ are uniformly away from collapse, so they have two times at circular distance at least some $\kappa>0$.

For (3), the zeros of $G_{\gamma'}$ in $X\setminus N$ lie in $P$ and in $O\setminus N$. The quotient of $O$ is $O_V$, because $T$ permutes the pieces freely. Additivity of degree gives the formula. ∎

**Lemma 369.6.** Let $\gamma$ be a polygon with $n\ge4$ edges, and let $\tau$ be the smallest distance between non-adjacent edges. Then $\gamma$ has no inscribed square of side $s<\tau/\sqrt2$.

*Proof.* The corners lie in a set of diameter $s\sqrt2<\tau$. Such a set meets no two non-adjacent edges, so it meets at most two edges, and two only if they are adjacent. The corners therefore lie on one segment or on two segments leaving a common vertex.

Four points on one line are not a square, and neither are three collinear points plus one more. With two corners on each of two lines, either the two pairs are opposite sides, which would make the lines parallel; or they are the diagonals, whose common midpoint is the vertex where the lines meet. Then the two corners of one diagonal lie on opposite sides of the vertex, but each edge lies on one side. ∎

The tests check the lemma against complete censuses of random star polygons. On the v247 decagon, $\tau\approx0.309$. Its only inscribed square has side about $1.410$, and its corners are in square order.

So for polygons the collapse circles are isolated components. A zero on $\partial N$ would be a square of side $\tau'/(2\sqrt2)$, so $\partial N$ is zero-free and the degree $\pi$ in (3) is defined. Whether $\pi$ is odd for every polygon is the classical parity question. v369 does not use it. If $\pi(\gamma_n)$ is odd along an approximating sequence and $d_R=0$, then (3) forces odd degree near $K_0$ at every stage. Those squares either stay in $K_0\setminus\Delta_0$, which by the corollary means the target has squares of every small size, or they collapse to points on $\Delta_0$.

## What this says

Every Jordan curve is in one of three regimes:
- some square family is invariant under relabeling, as for the circle; the curve has a square with corners separated by a third of the curve;
- no family is invariant and the degree near the genuine squares is odd; then squares persist with separated corners on every nearby curve;
- no family is invariant, that degree is even, and any odd parity of approximating polygons sits in squares near a single collapse circle.

Only the third regime can defeat a noncollapse argument. There, the obstruction is concentrated at one collapse piece and its three relabelings, not spread over the zero set. The v368 certificates are product-box tools for finding odd degree inside $R$. They can miss it when the time projections of an odd component overlap, and no box sees an isolated square of even local degree. Theorem 369.5 counts the degree of all of $R$ at once, so isolated squares of even degree simply contribute zero. Whether the third regime with genuine collapse ($K_0=\Delta_0$) can occur for a Jordan curve whose approximations have odd parity is the question this reduction leaves.

[Code](../code/relabel_dichotomy_v369.py) · [Tests](../tests/test_relabel_dichotomy_v369.py)
