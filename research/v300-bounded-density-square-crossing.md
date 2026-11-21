# v300 — Bounded boundary mass and a square crossing with moving contacts

Benjamin John Schulz · 20 November 2026

A bounded sweep of two boundary measures on any parametrized Jordan curve admits a single compact connected family of balanced, positive-size rectangles reaching both measures. For a class of curves described by two arcs at paired levels, two weights chosen before solving for any contacts force a square on that family. The same weights work through arbitrarily many polygonal refinements and on an infinite-length limit, with square side at least $1/9$ in the example's units. The general Jordan case still requires a geometric argument forcing the aspect ratio to cross one.

## 1. A connected family for arbitrary Jordan curves

Give the source circle length one and let $\gamma$ be a Jordan embedding. Label rectangle corners in geometric cyclic order and require the boundary to encounter them in that cyclic order (orient the source counterclockwise). Write $I_{ij}$ for the forward source arc from corner $i$ to its successor $j$, and color $I_{01},I_{23}$ red. For two probability measures $\mu_0,\mu_1$ on the source, each with density at most $K$ relative to normalized length, set $\mu_u=(1-u)\mu_0+u\mu_1$. A rectangle is balanced at $u$ when its total red mass is $1/2$. Cyclic relabeling exchanges the colors and preserves balance.

**Proposition.** There is a compact connected subset of $[0,1]\times(S^1)^4$ of balanced labeled rectangles on $\gamma$ whose projection is $[0,1]$. On the whole set each side length is at least
$$
m_\gamma(1/(4K)),\qquad
m_\gamma(a)=\min_{d_{S^1}(s,t)\ge a}|\gamma(s)-\gamma(t)|>0 .
$$
This is a continuum statement. There can be many rectangles at a given weight.

**Size.** One red arc has mass at least $1/4$ and at most $1/2$. Its length and the length of its complement are therefore at least $1/(4K)$ and $1/(2K)$ respectively. Its chord has the indicated lower bound. Apply the same argument to the blue arcs. Opposite sides are equal. On an approximating embedding within Euclidean sup distance $e$ of $\gamma$, the corresponding floor is at least $m_\gamma(1/(4K))-2e$.

**Polygonal continuation.** Schwartz's Theorems 1.5–1.7 give every polygon in his open dense generic set a graceful global component of labeled rectangles. Such a component is either a proper arc with aspect running from zero to infinity at its two ends, or a circle preserved by cyclic relabeling. His Corollary 1.4 gives single-measure balance on arbitrary Jordan curves. We use the polygon component structure to retain connectedness while the measure varies. [Schwartz, *A Trichotomy for Rectangles Inscribed in Jordan Loops*, §§1–2](https://www.math.brown.edu/~res/Papers/sweepout.pdf)

On the arc let $B_i=\mu_i(\mathrm{red})-1/2$. Both functions are continuous. At a flat limit the two adjacent short sides approach two distinct points. Their source arcs have length tending to zero because the complementary route between either short pair contains the other pair near a different source point. Thus the two $B_i$ extend with values $-1/2$ at one end and $+1/2$ at the other, up to exchanging the labels. The scalar $(1-u)B_0+uB_1$ has opposite strict signs on the lateral edges of the compactified arc times $[0,1]$. Its zero set contains a continuum joining the two weight edges.

Here is the elementary continuation fact used both here and below. For a continuous scalar on a rectangle with opposite strict signs on its vertical sides, approximate uniformly by piecewise affine functions on finer triangulations. Perturb finitely many vertex values away from zero while retaining both side signs. Zeros in each triangle are disjoint line segments and glue to arcs or loops. The lower horizontal side has an odd number of sign changes. Arcs returning to it account for an even number of ends, so at least one arc reaches the upper side. A Hausdorff-convergent subsequence of such compact arcs in the parameter rectangle yields a connected zero set meeting both sides, by uniform convergence. This proves the fact even for a nonisolated zero set.

On an elliptic component let $T$ cyclically relabel the corners. Then $B_i(Tp)=-B_i(p)$. A common zero of $B_0,B_1$ gives a single rectangle for the entire sweep. Otherwise normalize $(B_0,B_1)$ to a map into the unit circle. The free order-four action on the component circle can be parametrized by quarter shifts, reversing the parameter if necessary: a finite-order orientation-reversing circle map would have a fixed point, and the translates of a fundamental arc give this parametrization. For an angle lift $b(t)$ of the normalized pair, $b(t+1/4)-b(t)$ is a constant odd multiple of $\pi$ by continuity. Hence the winding is nonzero. The lifted angle traverses a full quadrant with $B_0\le0, B_1\ge0$, from one bounding axis to the other, in either order. Take the first hit of the far boundary after a traversal starts and the last preceding hit of the near boundary. Along that subarc,
$$
u=\frac{-B_0}{B_1-B_0}
$$
is well defined and continuous, equaling zero and one at its ends. Its graph gives the desired connected family in consistent labels.

Finally approximate the given parameterized Jordan embedding uniformly by generic Jordan polygons, retaining the same source circle, orientation and two source measures. Use the compatible polygonal approximations in Schwartz’s §2.3, followed if needed by arbitrarily small generic vertex perturbations in each polygon space. Take a Hausdorff-convergent subsequence of the polygon continua in the compact product $[0,1]\times(S^1)^4$. The limit is connected and meets both weight faces. Rectangle equations pass to the limit. The common side floor above makes every limiting rectangle genuine, so its source corners are distinct and retain cyclic order. The source arc masses converge because the two measures are nonatomic and the weights converge. This proves the proposition.

**What remains.** Cyclic relabeling reverses the sign of the first side minus the second. To deduce a square one needs opposite aspect signs on the *same consistently labeled continuum*. One cannot assign a sign to an unlabeled endpoint by arbitrarily choosing its first corner. The following construction supplies both signs in a geometric labeling.

## 2. How to prescribe the weights without locating a square

Let $h:[0,1]\to[0,\infty)$ be continuous, vanish at the two ends and be positive inside. Trace the lower arc from left to right and the upper arc back. The basic example is
$$
\gamma_\pm(x)=(x,g(h(x))\pm h(x))
$$
for any continuous $g$. These graphs enclose a Jordan domain. More generally take $\gamma_\pm(x)=T_{h(x)}(x,\pm h(x))$ for a continuous family of orientation-preserving plane isometries $T_y$, provided the two resulting arcs form a Jordan embedding with the stated cyclic order. Two positions $s<t$ with the same level $y=h(s)=h(t)$ give an actual rectangle by applying the one isometry $T_y$ to the four corners of $[s,t]\times[-y,y]$. Its consistently labeled width and height are $t-s$ and $2y$. This shared motion of equal-level points is the structural hypothesis.

Give each copy of the $x$-interval half the probability measure $\rho(x)\,dx$ and use the corresponding source parametrization (each copy occupies half a turn, linearly in $x$). Its normalized source density is $\rho$. The red arcs of this rectangle have total mass $H(t)-H(s)$, where $H$ is the cumulative distribution on $[0,1]$.

**Weight selection for an isolated high region.** A sufficient profile condition is a threshold $b>0$, an interval $J$ of positive length on which $h\ge b$, and $\operatorname{diam}\{x:h(x)\ge b\}<2b$. Every profile with a unique maximum satisfies this condition. As $b$ increases to its positive maximum the nested closed superlevel sets shrink to the maximizer, and for each $b$ below that maximum a neighborhood has height above $b$.

Choose a positive baseline density and three disjoint packets near the left end, the middle and the right end. Give the baseline mass $1/4$ and each packet mass $1/4$, uniformly on its interval. For sufficiently small packets, the gap between consecutive packets exceeds twice the supremum of $h$ on the two end packets. An interval of mass $1/2$ must meet at least two packets in positive mass because a proper interval gets strictly less than $1/4$ from the uniform baseline. Thus it reaches an end packet, has one endpoint near the corresponding end of the domain, and spans at least the packet gap. More precisely take the end packets abutting 0 and 1. Equal endpoint heights then force width greater than height. For the second density assign half the mass uniformly to the domain and half uniformly to $J$. A half-mass interval can neither contain $J$ nor miss it, since $J$ itself has more than half the mass. Thus some endpoint is in $J$. Equal endpoint heights put both endpoints in the indicated superlevel set, forcing height greater than width. All choices use profile bounds, before computing any rectangle.

Mix the two densities linearly in $u$ and write $H_u$ for the cumulative distribution. Strict lower bounds on the densities make $(u,p)\mapsto H_u^{-1}(p)$ continuous. For $0\le p\le1/2$ set
$$
s=H_u^{-1}(p),\quad t=H_u^{-1}(p+1/2),\quad F(u,p)=h(s)-h(t).
$$
At $p=0$ the sign is negative for every $u$, and at $p=1/2$ it is positive. The scalar continuation argument gives a compact connected zero set reaching $u=0,1$. On it $s<t$ and $h(s)>0$ (a zero-height pair would use both ends and have mass one). These are balanced genuine rectangles. The width minus height takes both signs on this one set, so it contains a square. If both densities have upper bound $K$, every width and in particular every square side on it is at least $1/(2K)$. The two source coordinates are solved over the full half-mass interval domain; no four corner windows are inputs.

## 3. A single dyadic sweep for every refinement and its limit

Use the three wide packets $[0,1/16]$, $[15/32,17/32]$, $[15/16,1]$ and the tall packet $J=[7/16,9/16]$. The actual densities are
$$
\rho_W=1/4+4\,\mathbf1_{\text{three packets}},\qquad
\rho_T=1/2+4\,\mathbf1_J .
$$
They integrate to one and their mixtures have $K=9/2$. The following coarse inequalities suffice at all depths.

* $h\le1/8$ on the end packets.
* $h\ge3/8$ on $J$, and $\{h\ge3/8\}\subset[1/4,3/4]$.

The wide gate has width minus height at least $13/32-1/4=5/32$. The tall gate has width minus height at most $1/2-3/4=-1/4$. The selected square has side at least $1/9$.

For a concrete refining family put $T(x)=\min(x,1-x)$, $v(z)=2\operatorname{dist}(z,\mathbb Z)$ and $a(x)=1$ on the left half and $2$ on the right half. Set
$$
h_N(x)=\tfrac98 T(x)+a(x)\sum_{n=1}^N 2^{-n-6}v(4^n x),
\qquad g(y)=y/3 .
$$
All added terms vanish at the midpoint, so the profiles are continuous. Old vertices are retained, and the perturbations continue on every open interval. The graphs are rational polygons at finite depth. The uniform tail is at most $2^{-N-5}$, so their parametrized boundaries converge to the corresponding Jordan curve. Every profile lies between $9T/8$ and $9T/8+1/32$. This implies all the displayed gate inequalities. The midline varies with thickness, and the two half-profiles differ from depth one onward. The construction requires neither horizontal nor vertical reflection symmetry.

**Infinite variation.** On each cell of length $4^{-n}$, all terms before level $n$ are affine. At its endpoints and midpoint all later terms vanish. The midpoint of the limiting profile exceeds the mean of its two endpoint values by $a\,2^{-n-6}$ with the cell's interior multiplier. The variation on the two half-cells is at least twice this excess. Summing over $4^n$ cells gives variation at least $2^{n-5}$, which diverges. Each graph has infinite length since its height is a nonzero constant multiple of the profile. The same estimates for truncated contours and Hausdorff compactness in $(u,p)$ also give a subsequential limit of spanning witnesses. There is no need to identify one particular square across refinements.

The exact replay checks depths zero through four. The unique *equal-level* square moves from $(s,t)=(4/17,13/17)$ at depth zero to $(63/268,821/1072)$ at depth one. Both coordinates are free to change. At depth four they are $(863/3664,2831/3664)$ and the mixture weight is $53/176$. The computation exhausts the equal-level square equation at these depths. Oblique squares outside this construction are outside that census.

## 4. Controls and limits

* On the complete v274 concave atlas give each of the four edges one quarter of the parameter interval. If $q$ is a rectangle's red source fraction, the least balancing cap is $1/(2\min(q,1-q))$, achieved by constant densities on the two colors. The vertex square has $q=23/48$, cap $24/23$. The persistent square has $q=409/528$, cap $264/119$. Hence uniform weighting misses both squares, and the cheapest square at this time is on the folded component that disappears for negative time. Substituting the four v288 source coordinates along the persistent band gives $q(x)=x(169+71x)/(48(1+3x+x^2))$ for $x=w/(w+h)$. Its derivative has numerator proportional to $169+142x+44x^2>0$. On $[2/5,3/5]$ the greatest required cap is $316/103$. These retrospective caps benchmark the computation. The weights in §3 are prescribed in advance. [Atlas](v274-full-rectangle-binormal-atlas.md) · [Moving fold](v280-v283-moving-relative-class.md) · [Source coordinates](v288-v291-rough-arc-continuation.md)
* Two exact checks exercise nonuniqueness and cyclic relabeling. A folded rational profile in the tests has five isolated balanced rectangles of the selected form at mixture $1/4$. On the rotating square loop of a circle, take two source densities alternating $3/2$ and $1/2$ on quarters, shifting the second density by $1/8$ of a turn. The signed red masses change sign under a quarter turn, have no common zero, and on turns $1/8\le t\le1/4$ their balancing mixture runs from zero to one.
* Even an arbitrarily small drift can remove every rectangle using the selected vertical pairs. Replace the midline $g(h(x))$ by $g(h(x))+\delta x$, $\delta\ne0$. A rectangle with two vertical sides at $s<t$ requires equal half-thicknesses and equal midline heights. The residual midline difference at equal level is $\delta(t-s)\ne0$. The resulting graphs still form a Jordan curve and converge uniformly to the unperturbed example as $\delta\to0$. The obstruction concerns that particular pairing of supports. Other rectangle orientations remain possible, and §1 still applies.
* The cap is relative to a specified source measure. Preserving the parametrization across refinements matters. The folded profile demonstrates why uniqueness cannot be assumed. The concave atlas with both weights uniform gives a constant sweep with no balanced square. For a general embedding the open issue is to choose weights and prove an aspect crossing on a single consistently labeled family.

## Replay

Run `python -m unittest discover -s tests -p 'test_balanced_sweep_v300.py' -v` and `python code/balanced_sweep_v300.py --check` from the repository root. All calculations use rational arithmetic in the Python standard library. The [certificate](../results/v300/balanced-sweep.json) contains the finite checks. Sections 1–3 provide the topological and limiting proofs.
