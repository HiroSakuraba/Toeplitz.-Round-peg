# v301 — What becomes of the square when the vertical rectangles disappear?

Benjamin John Schulz · November 2026

An arbitrarily small nonzero shear of the v300 seed leaves exactly one geometric square. It tilts, moves to another set of edges, and escapes the entire original interpolation of boundary weights. Positive and negative shears select different ends of a continuum of squares on the unsheared seed. Replacing one boundary weight recovers a connected band of tilted rectangles through the square. A second construction gives a macroscopic balanced square on every sheared rough refinement and its infinite-length limit.

## 1. The lost square and the surviving one

Number the directed edges of the convex quadrilateral with successive vertices $(0,0),(1/2,-3/8),(1,0),(1/2,3/4)$ by $0,1,2,3$. Apply the shear $S_d(x,y)=(x,y+dx)$. Use the same source as [v300](v300-bounded-density-square-crossing.md), with the bottom graph traversed in increasing $x$ and the top in decreasing $x$. The source densities stay attached to these parameters.

**At zero drift.** Write $w=9/17$, $m=4/17$, $s=m+z$ and $t=m-z/2$, with $-4/17\le z\le4/17$. Take the bottom two corners at $s,s+w$ and the top two at $t+w,t$. The first two successive side vectors are $(w,3z/2)$ and $(-3z/2,w)$. Thus the v300 vertical-pair square at $z=0$ lies inside a continuum of genuine squares. The v300 counts explicitly concerned equal-level pairs. In particular, the full four-edge square chart is singular.

**For $0<d\le1/100$.** There is exactly one square up to cyclic and reflected labelings. Its contacts in increasing boundary order belong to edges $(1,2,3,3)$. Put $D=4d^2+12d+17$. Its affine parameters within those edges are
$$
(a,b,c,e)=D^{-1}\big(4d^2+12d+1,\ 4,\ 5-4d^2,\ 17-4d^2\big).
$$
All are strictly interior and $c<e$. The alternative chart using all four edges has tags $((2d+3)/(2d),3/(2d),(4d+3)/(4d),3/(4d))$. Its affine solution escapes the edges.

**Exhaustion.** The seed polygon is convex with strict vertex turns, and the four vertices of a positive square are extreme in their own convex hull. Sorting them along the convex boundary gives the geometric cyclic order even when a support edge contains two of them. Three square corners cannot lie on a line. If two polygon edges alone carried all four corners, they would either carry parallel opposite sides or contain intersecting diagonals intersecting in the interiors of both segments. No two seed edge lines are parallel and the boundary is simple. These properties persist under the shear. Thus the sorted edge assignments consist of the all-distinct chart and twelve charts using three edges.

For each of the thirteen charts the counterclockwise square equations are a $4\times4$ affine system in its edge parameters. Its entries and right side are affine in $d$. The determinant and four Cramer numerators therefore have degree at most four. The replay evaluates each at $100d=0,1/4,1/2,3/4,1$, then recovers its exact degree-four Bernstein coefficients. Nonnegative coefficients and a positive last coefficient give strict positivity on $(0,1]$. The reversed signs prove negativity. The determinants have fixed nonzero signs throughout the punctured interval. For each excluded chart the certificate records one tag strictly outside the unit interval for every drift. For the winner it checks all eight strict bounds and the order of the two repeated-edge contacts. Substitution of the displayed rational functions in the original equations verifies their identities (after clearing $D$, each residual has degree at most three). [Per-chart certificate](../results/v301/tilted-sweep.json)

Reflecting in $x=1/2$ and translating vertically by $-d$ sends the $d$ curve to the $-d$ curve. The unique positive-drift square tends to the displayed seed square at $z=-4/17$ (the left tip). The negative-drift square tends to the one at $z=4/17$ (the right tip). These are distinct even as unordered point sets. Consequently there is **no continuous choice of a single geometric square through zero drift** on this family. The square sizes stay positive on both sides.

## 2. Every original mixture misses the tilted square

For positive drift, cyclically relabel the winning contacts starting with the later contact on edge 3, so the ordered tags are $(e,a,b,c)$. In that order their $x$ coordinates $(l,m,n,p)$ are
$$
D^{-1}\big(2d(2d+3),\ (2d+3)^2,\ 4d^2+12d+15,\ 2(2d^2+3d+3)\big).
$$
For any density used below, write $H$ for its cumulative distribution on $[0,1]$ and give each boundary graph half the probability measure $dH$. The two red boundary arcs run from the upper left contact through the left tip to the bottom contact, and back along the upper graph from $n$ to $p$. Their total mass is
$$
q_H=\tfrac12\big(H(l)+H(m)+H(n)-H(p)\big).
$$
One cyclic relabel exchanges red and blue, so the existence of a balancing mixture is label-invariant.

The exact coordinate enclosures are $[0,1/256]$, $[9/17,8/15]$, $[15/17,8/9]$ and $[7/20,71/200]$, respectively. Positive-denominator polynomial signs verify them on the whole punctured interval. Evaluating the piecewise-linear CDFs at the enclosure endpoints, with the fourth sign reversed, yields the following sufficient bounds. This also handles the bottom contact crossing a density breakpoint.

| Measure | Lower red fraction | Upper red fraction |
|---|---:|---:|
| Original wide | $13693/27200$ | $47677/92160$ |
| Original tall | $9493/13600$ | $32749/46080$ |
| New end-packet measure | $5293/13600$ | $18517/46080$ |

The first two lower bounds exceed $1/2$. Since mass is affine in the mixture and the geometric square is unique, **no original v300 mixture balances any square** at these shears. Reflection preserves the two symmetric source densities and balance, proving the same for negative drift. For $d=1/100$, the two actual fractions are $58807/114136$ and $242015/342408$. In particular, the general v300 continuum of positive-size balanced rectangles can miss the square angle after an arbitrarily small perturbation.

## 3. A repaired weight pair and a connected tilted crossing

Retain the tall measure. Replace the wide measure by the symmetric density
$$
\rho_E(x)=1/2+4\,\mathbf1_{[0,1/16]\,\cup\,[15/16,1]}(x).
$$
Both densities integrate to one and are bounded by $9/2$. On the unique square we have $q_E<1/2<q_T$. The mixture weight on the tall measure is exactly $(1/2-q_E)/(q_T-q_E)$. The same two measures work for every drift in the interval. At $d=1/100$ the tall measure receives weight $34783/105594$.

There is also a connected family crossing the square angle under independent positional changes. Fix $0<d\le1/100$ and the edge labeling $(3,1,2,3)$ from §2. Allow any compatible continuous Jordan embedding at Euclidean sup distance at most $d/100000$ from the sheared seed. In the contact box of halfwidth $d/8$ around its four *reference* tags, prescribe side ratio
$$
1-d/20000\ \le k\le\ 1+d/20000.
$$
Use the residual with two vector entries $p_2-p_1-kJ(p_1-p_0)$ and $p_3-p_0-kJ(p_1-p_0)$, where $J$ is counterclockwise quarter turn. At a zero the second side is $k$ times the first side rotated through a right angle.

Here are uniform bounds for the degree calculation.

1. The center tags stay strictly inside their edges after the box enlargement. In particular $1-e=(8d^2+12d)/D>2d/3$, and $a,b,c$ and their distances from 1 exceed $1/20$. Also $e>9/10$ and $e-c>1/2$. This also keeps the repeated-edge contacts in order.
2. At $k=1$ the reference residual has constant matrix $A$ with determinant $9D/256$. Its entries have absolute value below 2, since the two individual edge coordinate speeds are at most $1/2$ and $3/4+d/2$. Each cofactor has absolute value at most 48. Since $9D/256\ge153/256$, $\|A^{-1}\|_\infty<400$.
3. Throughout the box each component of the reference first side is at most 1 in absolute value. The horizontal span is at most 1. Its vertical difference uses one upper point at $l\le1/256+d/16$ and one lower-right point with $m\ge9/17-d/16$, so it is bounded by $3/8+2(1/256+d/16)+d<1$. Moving $k$ changes each row by at most $d/20000$. An independent Euclidean positional error $\epsilon$ at the four tagged points contributes at most $(2+2k)\epsilon\le6\epsilon$. Preconditioning therefore changes the translated identity by at most $400(d/20000+6d/100000)=11d/250<d/8$.
4. The resulting Brouwer degree is $+1$ after preconditioning, with a zero-free lateral boundary for the entire ratio interval. The connected continuation lemma from [v284](v284-v287-rectangle-family-stability.md) gives a compact connected set meeting both ratio endpoints. It contains a square. Its first side stays tilted: the reference upper contact has nonnegative height and the lower contact has height at most $-3(1-8/15-d/16)/4+d<-3/10$, leaving ample margin for the two positional errors. The first side’s horizontal projection on the perturbed curve is at least $9/17-1/256-d/8-2d/100000>1/2$.
5. Mass is computed in the transported source coordinates, even when the perturbed contacts move away from their reference edges. Each tag changes its reference $x$ by at most $d/16$. The formula of §2 changes by at most $(9/2)d/8$. The strict bounds $q_E<1/2<q_T$ still hold. The unique balancing weight is continuous on the entire crossing. Thus the lifted set in ratio, contacts and weight is connected and contains wide, square and tall rectangles under the same two preassigned densities. This statement concerns the weights realized on the ratio band. It makes no claim that they cover the entire mixture interval $[0,1]$.

Reflection gives identical allowances in $|d|$ for negative drift, with cyclic labels adjusted. Both the allowed positional error and the ratio bandwidth shrink with the shear. The next result instead covers the full rough amplitudes.

## 4. The full sheared rough family and a weight sweep fixed in advance

Let $h_N$ and $h_\infty$ be the exact profiles of v300, and shear the two graphs $-2h_N/3$, $4h_N/3$ by $|d|\le1/100$. Allow a further independent Euclidean displacement of at most $1/1000$ in the fixed parametrization, assuming the resulting loop is Jordan. Every such curve has a square of side at least $29\sqrt2/100$. Every square of this size on these curves can be balanced by a member of the same sliding-packet family, with density cap $13/2$.

**Existence.** Take $O=(1/2,1/8)$. Its distances to the lower and upper supporting lines of the unsheared polygon are $2/5$ and $5/(4\sqrt{13})$. Its farthest polygon vertex is at distance $5/8$. The sup distance from the seed to any of the stated curves is at most $E=1/24+1/100+1/1000=79/1500$, using the v300 profile tail bound from depth zero and the graph height multipliers. The entire curve lies in the annulus about $O$ of inner radius $29/100$ and outer radius $69/100$. Straight interpolation from the seed stays away from $O$, so the winding is one. The ratio $69/29$ is less than $1+\sqrt2$. Matschke's annulus theorem therefore gives the claimed side floor. For the sheared two graphs themselves at nonzero drift, any axis-aligned rectangle would require equal half-thickness and equal midline at two distinct $x$ values. Their midline difference at equal thickness is $d(t-s)$, so every square on these graphs is tilted. Additional independent perturbations preserve the size conclusion. [Matschke, *A Survey on the Square Peg Problem*, Theorem 5](https://www.math.u-bordeaux.fr/~bmatschke/SurveyOnSquarePeg_extension14.pdf)

**Order and separation of the contacts.** Let a square on the perturbed curve have side $L\ge29\sqrt2/100$. Each triple of its vertices is a right isosceles triangle of minimum strip width $L/\sqrt2$ (the shortest altitude). Since $2E<29/100$, moving its vertices back to their seed positions cannot make any triple collinear along the straight interpolation: collinearity would put the original triple in a strip of width $2E$. All four turn signs are preserved. The seed points therefore have the same cyclic order as the square, by convexity. Following the source gives the geometric cyclic order. In the normalized source the seed is $\sqrt{13}$-Lipschitz. A pair of consecutive square corners therefore has circular source distance at least $(L-2E)/\sqrt{13}>1/12$. In particular each of the four boundary arcs contains a full source interval of length $1/12$.

**Balancing.** Before choosing a curve or a square, prescribe the source-circle weights
$$
\rho_v(t)=1/2+6\,\mathbf1_{[v,v+1/12]}(t),\qquad v\in\mathbb R/\mathbb Z .
$$
For any of these squares with red source fraction $q\in(0,1)$, putting the packet in a red arc gives red mass $q/2+1/2$. Putting it in a blue arc gives $q/2$. The red mass is continuous as the packet travels, so it takes the value $1/2$. The argument obtains the square from the annulus geometry and balances it with a prescribed family. A connected near-square crossing for the full rough amplitudes is still open. In particular the packet argument does not imply that an arbitrary Jordan curve admits a macroscopic square or that one given pair of packets forces it.

## Replay

Run `python -m unittest discover -s tests -p 'test_tilted_sweep_v301.py' -v` and `python code/tilted_sweep_v301.py --check`. The rational replay exhausts all thirteen charts over the continuous punctured interval, and the tests compare against the independent v223 full-label census. The annulus theorem and the connected degree argument are mathematical inputs. [Verification records](../results/v301/tilted-sweep.json)
