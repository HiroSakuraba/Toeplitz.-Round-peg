# v362. A change of starting corner prevents collapse

On a circle a square can rotate until the next corner takes the place of the first. V361 shows why that motion defeats a local three-window certificate. Here the motion itself gives a lower bound on square size, even on approximations for which the moving figures are only approximately square.

## The one-third lemma

Use registered times on $\mathbb R/\mathbb Z$, with circular distance $d$. Represent four boundary-ordered labels by lifted times $t_0\le t_1\le t_2\le t_3\le t_0+1$ modulo a common integer shift. Allow equal times. The four forward gaps, including the wraparound gap, sum to one. A cyclic change $T$ of the starting label shifts the gaps and replaces the time row by $(t_1,t_2,t_3,t_0+1)$.

**Lemma.** Any continuous path from a row $x$ to $T^k x$ for $k=1,2,3$ has a row at which two parameters are at circular distance at least $1/3$.

If a gap exceeds $1/2$, its label is unique. If every row had such a gap, the label would be locally constant along the path, whereas the endpoint permutes it without a fixed label. Hence at some row all gaps are at most $1/2$. Suppose all pair distances there were less than $1/3$. With one time as origin, represent all times in $(-1/3,1/3)$. The span between the extreme representatives is less than $2/3$. Since those two representatives have circular distance less than $1/3$, their span must be less than $1/3$. This leaves a gap greater than $2/3$, a contradiction. The constant for this balanced-row implication is sharp, as four times can approach three equally spaced positions with one position repeated.

## A sufficient condition on approximations to any Jordan curve

Let $\gamma$ be a continuous injective loop and $P_n$ uniformly registered polygonal approximations with error $\delta_n\to0$. On each path as above, evaluate $p_i=P_n(t_i)$ and put $v=p_1-p_0$. In this criterion the cyclic boundary labels are also in positive geometric square order at a noncollapsed zero. Require, everywhere along the path,
$$
 \|p_2-p_1-Jv\|,\quad \|p_3-p_0-Jv\|\ \le\varepsilon_n,
 \qquad \varepsilon_n\to0 .
$$
A subsequence of approximations suffices. Reversing the convention gives the reflected sector. Crossing boundary-order sectors are outside this criterion.

Define $m_\gamma=\min_{d(s,t)\ge1/3}\|\gamma(s)-\gamma(t)\|$. Compactness and injectivity give $m_\gamma>0$. At the balanced row the separated physical pair is at least $m_\gamma-2\delta_n$ apart. The ideal square on the first side puts every measured corner within $\varepsilon_n$ of its ideal position. Therefore
$$
 \sqrt2\,\|v\|\ \ge\ m_\gamma-2\delta_n-2\varepsilon_n .
$$
Take a convergent subsequence of the four circle parameters. The residuals vanish in the limit and the bound keeps the square nondegenerate. Its side is at least $m_\gamma/\sqrt2$. This curve-dependent quantity uses a fixed registration and needs no inradius estimate. Recording only square-like endpoints would be insufficient, because a path can leave the square-like configurations between them.

## What would it take to construct the paths?

There is a precise restriction. On the compact space of lifted cyclic rows, a circle times the closed simplex of forward gaps, let $Z$ be the exact zero set on the target, including collapsed rows. Paths of vanishing uniform residual with a nontrivial relabeling exist along uniformly converging approximations **if and only if** some connected component of $Z$ is invariant under a nonidentity $T^k$.

For necessity, fix $k$ on a subsequence and take a Hausdorff limit of the compact connected path images and a limit $x$ of starting rows. The limit is connected, lies in $Z$ by uniform convergence and contains $x,T^kx$. The component of $Z$ containing $x$ intersects its translate, hence equals it. For sufficiency, in the ambient row space take the open set of rows with both target residuals below any given tolerance. Its component containing the invariant connected set is path-connected by local path connectivity. It contains $x,T^kx$, so join them there. Uniformly small redraws perturb both residuals uniformly. Taking tolerances to zero gives the required paths. This argument asserts existence and does not produce paths for general input.

A noncircular ellipse shows the restriction matters. In centered principal coordinates its equation is $p^THp=1$ with distinct positive diagonal entries. Opposite vertices of an inscribed square are $c\pm v,c\pm Jv$. Subtraction gives $c^THv=c^THJv=0$, so $c=0$. Equality of the two remaining quadratic values gives $v_x^2=v_y^2$, which yields exactly one geometric square. Its four choices of first corner are isolated. The collapsed rows form four disjoint circles, one for each choice of gap equal to one. Every nonidentity relabeling permutes these components without fixing one. Consequently there are **no arbitrarily accurate relabeling paths** for this ellipse, although it contains its square. V268's obstruction treats paths of genuine rectangles; the present statement concerns approximately square configurations and their possible collapse.

## Exact path replay

The verifier takes rational polygon vertices, their registration knots and a piecewise-linear path in lifted parameter coordinates. It checks cyclic order and endpoint relabeling. On each path segment it splits at every polygon knot crossed by any corner, including across the parameter seam. Between breaks both residual vectors are affine, so their maximum squared Euclidean norms occur at endpoints. Intersecting the four linear gap inequalities finds a balanced row and an explicit separated pair.

The blind search also builds a finite graph directly from the ordered polygon vertices, including rows with repeated vertices. It keeps rows below a chosen squared-error ceiling and joins them by single-corner or synchronized four-corner steps of one edge, respecting cyclic order in the lift. Along each such step the physical corner coordinates interpolate affinely, so checking the error at the two ends certifies the connecting edge. A connected component containing a row and a nontrivial relabeling supplies an explicit path for the independent replay. This exhausts the chosen graph only. The enumeration takes order $n^4$ vertex rows.

On a rational quarter-turn-symmetric octagon all three nonidentity relabelings replay with zero residual. Moving just one vertex horizontally by $h$ gives maximum squared residual exactly $2h^2$ on the tested quarter-turn path for $h=1/16,1/64,1/256$. The residual vanishes at both path endpoints on these tests. The blind search also finds paths on the reference octagon and the $h=1/64$ redraw, without a seed or a supplied path. A path beginning with all four parameters in $[0,3/100]$ and cycling them incurs squared error greater than one in the interior on the reference polygon. This is why comparing only the endpoints can falsely suggest a useful path.

**Next research target.** On approximations for which the labeled square components are isolated, group square configurations above a fixed parameter-separation cutoff and track their mod-two degree through refinement. V359 provides the fixed-arc compactness argument, and v362 handles a component that reconnects to another labeling. The missing implication is that at least one positive-separation mechanism survives on some sequence for each target. The circle and ellipse rule out making either of these particular mechanisms a universal requirement by itself.

[Exact path verifier](../code/square_relabel_path_v362.py) · [Tests](../tests/test_square_relabel_path_v362.py) · [v342 cyclic components](v342-unmarked-components-and-square-frames.md) · [v268 rectangle paths](v268-rectangle-monodromy-and-protected-square-criterion.md)
