# v302 — A connected rectangle crossing on the full rough family

Benjamin John Schulz · November 2026

**Theorem.** Take any sheared v301 rough refinement, including depth zero and the infinite-length limit, with shear of either sign, $|d|\le 1/100$. Independently displace its parametrized points by at most $1/1000$ in Euclidean distance, with a continuous Jordan embedding as the result. On each resulting curve there is a compact connected set of actual inscribed rectangles on which the length of the more vertical side divided by the length of the more horizontal side runs through the whole interval $[99/100,100/99]$. Every side is longer than $407/1000$ in units of the original tip separation. The source encounters the four corners in geometric order. A continuous choice of a balancing packet of source length $1/12$ exists on the whole set.

This gives wide rectangles, a square and tall rectangles on the same set even at zero shear, where the seed has a continuum of squares. Connectedness permits branches and folds. The theorem applies to each fixed curve and allows independent choices of perturbation at different depths. [v301 family and square census](v301-tilted-squares-and-boundary-weights.md)

## 1. A position envelope closed under interpolation

Let $P$ be the convex seed with successive vertices $(0,0),(1/2,-3/8),(1,0),(1/2,3/4)$. Write $T(u)=\min(u,1-u)$, extended by the same formula to all real $u$. The lower and upper reference points above $u$ have heights $-3T(u)/4$ and $3T(u)/2$. On the unit source circle the lower branch goes from $u=0$ to $1$ in half a turn and the upper returns in the other half. Set
$$
O=(1/2,1/8),\qquad r=29/100,\quad R=69/100,\quad E=79/1500 .
$$

The proof works under the following broader hypothesis. In these same source coordinates, the Euclidean distance to $P(u)$ is at most $E$, $|x-u|\le1/1000$, and the signed height deviation $y-P_y(u)$ belongs to the interval corresponding to its branch.

| Source branch | Lower deviation | Upper deviation |
|---|---:|---:|
| Lower | $-1/48-11/1000$ | $11/1000$ |
| Upper | $-11/1000$ | $1/24+11/1000$ |

Indeed, the v300 profiles have $h=9T/8+e$ with $0\le e\le1/32$. Their graph heights are $-2h/3$ and $4h/3$. Shear adds $du$ and the additional displacement contributes at most $1/1000$ to either coordinate. The triangle inequality gives $E=1/24+1/100+1/1000$. These pointwise constraints are convex, so linear interpolation from the seed with shear $+1/100$ to any target stays inside them. The intermediate parameter maps are allowed to self-intersect. The seed edge lines have distances $2/5$ and $5/(4\sqrt{13})$ from $O$, and the largest vertex radius is $5/8$. Since $r+E$ is below both line distances and $5/8+E<R$, every intermediate curve stays in the annulus $r<|p-O|<R$. [v300 profiles and uniform tail](v300-bounded-density-square-crossing.md)

## 2. Keep the center in the protected region

For a rectangle let $C$ be its center, $\rho$ its half-diagonal, and $D=|C-O|$. Set $t=99/100$. If the side ratio belongs to $[t,1/t]$, either side divided by the diagonal is at least $p_0=t/\sqrt{1+t^2}$. For $D>0$ let $p$ be the largest absolute projection of a unit center-to-corner direction onto the line $OC$. If $a,b$ are the half-side lengths and the line makes angle $\theta$ with the first side, then
$$
p=(a|\cos\theta|+b|\sin\theta|)/\sqrt{a^2+b^2}\ \ge p_0 .
$$
The minimum and maximum squared vertex radii about $O$ equal $D^2+\rho^2\mp2D\rho p$. At $D=\rho>0$ their ratio would require $R/r\ge\sqrt{(1+p_0)/(1-p_0)}=\sqrt{1+t^2}+t$. Our annulus violates this condition, since $r^2+2trR-R^2=2099/500000>0$. A collapsed rectangle at $D=\rho=0$ would lie at $O$. Thus a zero cannot cross the boundary $D=\rho$.

Inside $D<\rho$, the minimum squared vertex radius is at most $\rho^2$, since $p_0>1/2$ (also immediate if $D=0$). Opposite corners are at most $2R$ apart, so $r\le\rho\le R$. Each side has squared length at least
$$
4r^2t^2/(1+t^2)=8242641/49502500>(407/1000)^2 .
$$

## 3. Exclude the orientation at which the side labels would exchange

Suppose a rectangle in the protected region has a side at $45^\circ$ to the horizontal, including slope $-1$. Write its center as $(X,Y)$ and let $a,b$ now denote its **full** side lengths. Its vertical tips have heights $Y\pm A$ at horizontal offsets of magnitude $|B|$ from $X$. Its horizontal tips have abscissas $X\pm A$ at vertical offsets of magnitude $|B|$ from $Y$, where
$$
A=(a+b)/(2\sqrt2),\qquad |B|=|a-b|/(2\sqrt2).
$$
Since $|B|/A\le1/199$ and $A^2+B^2=\rho^2$, we have $A>0.289$ and $|B|<0.004$. All decimals below are exact rationals.

At a vertical tip the source coordinate differs from $X$ by less than $0.005$. By the 1-Lipschitz property of $T$ and the table in §1, an upper-branch tip has height in $1.5T(X)+[-0.019,0.061]$. A lower-branch tip has height in $-0.75T(X)+[-0.036,0.015]$. The separation $2A>0.578$ exceeds the width of either interval. An upside-down assignment is also impossible, since all lower points have $y\le0.011$ and all upper points have $y\ge-0.011$. The top tip is upper and the bottom tip lower. Subtracting and adding their height bounds gives
$$
9T(X)/8-0.017\le A\le9T(X)/8+0.0485,\qquad
Y\ge3T(X)/8-0.0275 .
$$
Consequently $T(X)>0.21$ and $Y-|B|>0.011$. Both horizontal tips must belong to the upper branch.

At any upper point the deviation from $1.5T(x)$ at its **actual** abscissa lies in $[-0.0125,1/24+0.0125]$, an interval of width $1/15$. The two horizontal tips differ in height by at most $2|B|$. It follows that
$$
|T(X-A)-T(X+A)|\le (0.008+1/15)/1.5=56/1125<1/20 .
$$
If both tips were on one side of the tent's crest, the difference would be $2A$. They therefore straddle it, where the difference in absolute value equals $|2X-1|$. Hence $T(X)>19/40$ and $A>4139/8000$. But all actual abscissas lie between $-0.001$ and $1.001$, so $A\le501/1000$. The bounds conflict by $131/8000$. This excludes the orientation for every allowed ratio on the entire interpolation.

## 4. One degree selects one spanning continuum

On ordered source parameters $(s_0,s_1,s_2,s_3)$ let $p_i=\gamma(s_i)$, $v=p_1-p_0$ and $J$ be counterclockwise quarter turn. Set
$$
F_k=(p_2-p_1-kJv,\ p_3-p_0-kJv).
$$
For arbitrary tuples define $C=(p_0+p_2)/2$ and $\rho=|p_2-p_0|/2$. Select the open domain $|C-O|<\rho$ and $v_x>|v_y|$. At a noncollapsed zero this singles out the counterclockwise labeling starting with the right-going, more horizontal side. That side is unique up to translation to the opposite side, on which the counterclockwise direction is reversed. The ratio $k$ has the stated vertical-over-horizontal meaning. The first inequality excludes collapsed zeros even in the closure. Sections 2 and 3 exclude every zero on either lateral boundary.

At shear $d=1/100$, $k=1$, the v301 census gives exactly one geometric square, with all four parameters strictly interior to their reference edges. In cyclic order starting on its right-going side those edges are $(3,1,2,3)$. The exact substitution gives $\rho^2-|C-O|^2=14458021697071/73277024040000$ and $v_x-|v_y|=2550/14267$. In the four oriented edge parameters the Jacobian has determinant $385209/640000\ne0$. The square contributes degree one modulo two and is the only zero in the selected domain. [Exact replay record](../results/v302/rough-crossing.json)

For an explicit Euclidean degree, replace each $s_i$ by the angle of $z_i\in\mathbb R^2$ with $1/2<|z_i|<3/2$, and append the four equations $|z_i|-1=0$ to $F_k$. The new zeros have all radii one. Radial variables give four independent derivatives at the seed zero, preserving its regularity. The product domain is bounded. Along the linear curve interpolation with $k=1$, and then as $k$ varies over $[t,1/t]$ on the target, the zero set in the closed selected domain is compact and misses its lateral boundary. Homotopy invariance on these open slices therefore gives odd degree at every ratio.

For the stronger connected conclusion, apply the continuation principle on $[t,1/t]$ with the target fixed. If no component met both ends, Whyburn's separation lemma would split the compact zero set into two disjoint compact sets containing the respective end fibers. An open neighborhood isolating the first set would have odd degree at the initial end and zero degree at the final end, contradicting homotopy invariance. Thus a compact component meets both ends and projects onto the full interval. Removing the radial variables gives the claimed source configurations. This proof applies directly to the limiting continuous curve. With a uniformly convergent sequence of registered targets, a subsequential Hausdorff limit of their continua also spans both ends: compactness preserves connectedness, the equations pass to the limit, and the two boundary exclusions keep the limit in the selected domain. [Smets, *Brouwer Degree and Applications*, chapter 7](https://www.ljll.fr/smets/ULM/Brouwer_Degree_and_applications.pdf)

## 5. Source order and one balancing choice on the entire set

Three corners of a rectangle of full side lengths $a,b$ form a right triangle of minimum strip width $ab/\sqrt{a^2+b^2}$ (its smallest altitude). This is at least $2rt/(1+t^2)=5742/19801>2E$. If three vertices became collinear while moving back along straight segments to their unsheared reference points, the original three vertices would be in a strip of width $2E$ about that line. Thus every oriented triple keeps its sign. On the convex counterclockwise seed these four reference points occur in geometric order, so the original source parameters do too. The seed has Lipschitz constant $\sqrt{13}$ for circular source distance. Each consecutive source gap exceeds $(407/1000-2E)/\sqrt{13}>1/12$.

Write $q$ for the total source length of the red forward arcs $I_{01}$ and $I_{23}$ and $\ell=1/12$. Prescribe the unit-mass densities $1/2+6\mathbf 1_{[w,w+\ell]}$ on the unit source circle. Choose the phase
$$
w=s_1-\ell(1-q)\pmod 1 .
$$
The packet lies across the red-to-blue junction $s_1$, with red overlap $\ell(1-q)$ and blue overlap $\ell q$. Both adjacent gaps exceed $\ell$. Thus the red mass is $q/2+6\ell(1-q)=1/2$. Source gaps and $w$ vary continuously on the selected configurations, including across the circle seam. Adjoining $w$ to the continuum preserves compactness and connectedness. Balance here uses the specified boundary parameter, which remains meaningful on the infinite-length curve.

## 6. A directly computed crossing on a rough polygon

Take the sixteen-edge depth-one polygon with shear $+1/100$. Number its vertices from the left tip in source order. Move vertex $i$ by $10^{-4}((3i\bmod 7)-3,(2i\bmod 5)-2)$ and interpolate the errors on each edge. Its Euclidean error is at most $\sqrt{13}/10000$ and the resulting polygon passes the exact simplicity check. On this input the ordered edge tuple $(15,4,8,13)$ supports an explicit continuous branch over the **whole** ratio interval. All four within-edge parameters stay strictly between $1/100$ and $99/100$.

To verify this, form the affine contact matrix of $F_k$ on those four edges. Its determinant and the four Cramer numerators have degree at most four in $k$, since every matrix and right-hand-side entry is affine in $k$. Interpolate at five rational equally spaced ratios, including both endpoints. The replay recovers the degree-four Bernstein coefficients of the determinant and both cleared inequalities for each tag. All 45 coefficients have the determinant's sign. The determinant stays nonzero, the four rational functions of $k$ remain inside the edges, and substitution at $k=1$ places the square in the protected region. Sections 2 and 3 then keep that branch inside throughout. This example supplies actual continuously moving contact positions as an additional check on the general topological argument.

## Replay and scope

Run `python code/rough_crossing_v302.py --check` and `python -m unittest discover -s tests -p 'test_rough_crossing_v302.py' -v`. The arithmetic replay checks the bounds and seed input. The tests independently solve for actual endpoint rectangles on the positive, negative and zero shear seeds. The continuum and all-depth results use the proof above. For curves outside the position envelope, a new argument is needed to separate or track the two side labels. [v303 direction search and its obstruction](v303-finding-a-direction.md)
