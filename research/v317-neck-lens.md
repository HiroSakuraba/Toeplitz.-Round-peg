# v317. Selecting two graph arcs and excluding artificial contacts

Benjamin John Schulz · November 2026

Take two disjoint portions of a Jordan boundary that are graphs over the same interval. They can bound a wide middle region even if their ends never meet and the rest of the curve folds back across the interval. We give four height comparisons that allow a large-square theorem for *closed* pairs of graphs to be used on these *open* portions. A square-grid selector obtains the arcs and cuts from the vertices of the polygon. A second test computes odd local degree and a sharper side bound when the selected graphs are close to four affine segments.

## Local theorem

Use orthogonal coordinates with equal length units. Let $f,g$ be 1-Lipschitz on $[L,R]$, with $f<g$ in the interior and $f\le g$ at the ends. Suppose their graphs are portions of the target. Put $H=\max(g-f)>0$ and $c=9/500$. Choose $0<d\le cH/\sqrt2$ with $d\le R-L$. At each end with a positive gap require the following inequalities. Omit the test at an end where the graphs already meet.

| End | Two height comparisons |
| --- | --- |
| Left | $g(L)-d<f(L+d)$ and $g(L+d)<f(L)+d$ |
| Right | $g(R)-d<f(R-d)$ and $g(R-d)<f(R)+d$ |

**Conclusion.** The two original graphs contain a square of side at least $cH$, with an adjacent pair of vertices on each graph. The remaining boundary of the target is unrestricted. The verifier uses $d=H/80$, which satisfies $d<cH/\sqrt2$. A certified disk of radius $r$ between the graphs whose center abscissa lies in the interval gives $H\ge 2r$ by taking its vertical diameter.

**External input.** Rifford's Theorem 1.1 supplies the constant $0.018$ for two 1-Lipschitz graphs meeting at both ends and bounding a Jordan curve. We use that theorem on a temporary closed curve. The open-arc conclusion requires the additional argument below. [L. Rifford, *A quantitative version of Tao's result on the Toeplitz Square Peg Problem*, Theorem 1.1](https://arxiv.org/abs/2106.01914)

### Proof

At the left end let $e=g(L)-f(L)$. For $e>0$, extend each graph to their common midpoint height at $x=L-e/2$ with straight segments of slopes $-1$ on the lower graph and $+1$ on the upper graph. Cap the right end in the same way, reversing the slopes. Both graphs remain 1-Lipschitz and the maximum gap stays $H$.

First establish that Rifford supplies a square that can be labeled with lower vertices $(t,f(t))$ and $(t+a,f(t)+b)$, $a>0$, and corresponding upper vertices obtained by shifting $(-b,a)$. In this paragraph $f,g$ denote the extended graphs. Multiply their ordinates by $0<\lambda<1$. On either rescaled graph every secant has absolute slope less than one, so three corners of a square cannot lie on it: two of their joining lines are perpendicular. Two diagonally opposite corners on each graph are likewise impossible. Hence each graph carries an adjacent pair, and a common tip cannot be a vertex (the union would then contain at most three of the four). Orient the lower pair from left to right. A downward completion would put an upper vertex at $(t+b,f_\lambda(t)-a)$ below the lower graph there, since $|b|<a$. The upward completion follows. Now let $\lambda\uparrow1$. Compactness retains this labeling and gives a square of side at least $cH$ on the completed pair, including at the critical Lipschitz constant one. Its horizontal run $a$ satisfies $a\ge cH/\sqrt2\ge d$ and $|b|\le a$.

Let $W(x)=g(x)-f(x)$ on the completed interval. At the abscissa of **each** corner one has $a-|b|\le W(x)$. For instance the upper partner of the lower-left corner has ordinate $f(t)+a$ and abscissa $t-b$, so $g(t)\ge f(t)+a-|b|$. Apply the same bound to $f$ at $t-b$ and to the other two corners.

Suppose some vertex lies in the added left cap. Take the leftmost corner at $x=L-z$, $0<z\le e/2$. For $b\ge0$ it starts the upper side on $g$, and for $b\le0$ it starts the lower side on $f$. Ties cause no problem. On that graph $F$ take the signed slope direction $\sigma=+1$ or $-1$ respectively. Define the deficit of a chord over $[u,v]$ by
$$
D(u,v)=(v-u)-\sigma(F(v)-F(u)).
$$
The deficit is nonnegative, additive on consecutive intervals, and at most twice the interval length. For the whole horizontal run it is $a-|b|\le W(L-z)=e-2z$. On the real arc, the left comparisons give $D(L,L+d)>e$ for the applicable graph and sign. Also $d>e/2$, since the two comparisons contain both graph values in an open interval of length $2d-e$. Therefore the first length-$d$ piece of the offending chord contains $[L,L+d-z]$; its deficit exceeds $e-2z$ after removing a terminal interval of length $z$. This contradicts the whole-chord bound. At the right end the rightmost corner ends a chord on the lower graph if $b\ge0$ and on the upper graph if $b\le0$. Use the last length-$d$ piece and remove the initial length $z$ from $[R-d,R]$. The two right comparisons give exactly the same contradiction. All corners are genuine.

**Why cap exclusion matters.** Take the diamond with tips $(0,0),(2,0)$ and other corners $(1,\pm1)$. Trim its graphs to $[\epsilon,2-\epsilon]$. The endpoint gaps tend to zero but the caps restore the diamond's two outer tips and its full square of side $\sqrt2$. That particular square uses two artificial vertices at every trimming. The height tests reject these ends. Other squares may still survive trimming.

## Automatic selection on a square grid

The implementation takes a simple polygon in cyclic order, with no supplied cuts or marked arcs. In frame $u=(x+y)/2$, $v=(y-x)/2$, every nonzero horizontal or vertical grid edge becomes a graph segment of slope $+1$ or $-1$. Split at reversals of $u$ and orient each maximal run by increasing $u$. The other diagonal frame interchanges $u,v$ up to sign and can be selected separately. Lengths in either frame must be multiplied by $\sqrt2$ to recover physical lengths.

For each pair of runs, take their overlap and order the two graphs by height. Compute the piecewise linear gap at their merged abscissas. Try the full overlap and the closures of connected components where the gap exceeds a dyadic level. This picks out a large middle lobe even when another lobe lies beyond either narrow opening. Clip the ends by exact interpolation, find the maximum gap, and evaluate the height comparisons. Keep the largest certified $H$ in the tested frame. The finite level search can miss a certificate between its tested levels. A miss carries no inference about square existence. A closed axis-aligned intersection sweep first checks the input polygon for overlaps and self-contacts.

In particular, on a grid loop consisting of two complete coordinatewise northeast paths with common southwest and northeast ends, the full pair meets at both tips. Rifford already yields $M(P)\ge(9/500)\sqrt2\,i(P)$, where $i(P)$ is the largest side of an axis-aligned square in the interior and $M(P)$ is the largest inscribed-square side. The physical diagonal slice through the center of the interior square has length at least $\sqrt2\,i(P)$ in that direction (use a supremum if necessary). The localized test permits other boundary portions which violate this two-path description.

## An all-refinement example and a larger local bound

In the $(u,v)$ frame use two broken lines with abscissas $0,4,8$. The lower ordinates are $0,-1,1$ and the upper ones are $e,2+e,1+e$ with $e=1/256$. On each interval of horizontal length $1/n$, replace a line of slope $m$ by two segments of slopes $+1,-1$ on the upper graph or $-1,+1$ on the lower graph, choosing durations to match the old endpoints. The vertical departure is at most $(1-m^2)/(2n)\le1/(2n)$. The teeth point outward, so the gap exceeds $e$ throughout $(0,8)$. Its maximum $H_n$ is at least $3+e$.

Attach a folded detour at each end. At the right end use local coordinates based at the lower endpoint and follow $(0,0),(2,-2),(1,-3),(2,-4),(6,0)$. Connect to its reflection about $v=e/2$ by a tip at $(6+e/2,e/2)$, traversing the reflected half in reverse. At the left end use half the detour's displacements before reflection and point it to the left, keeping the seam gap $e$. These are simple disjoint attachments. Each lower half stays below the seam level except at its ends, the upper half stays above, and each detour lies on the exterior side of its seam. The two center arcs have disjoint interiors by their gap. Mapping back by $(x,y)=(u-v,u+v)$ gives a simple square-grid polygon up to the common lattice unit. The first diagonal frame has ten maximal runs.

For all $n\ge64$, the comparisons on $[0,8]$ have margin at least $d/2-e-1/(2n)>0$, using $d=H_n/80$, the nominal slope bound $1/2$ and exact endpoint values. Both seams have gap $e$, and on a short extension of the central runs beyond either seam the gap increases. Thus $(0,8)$ is a component at the tested dyadic level $e$. The local theorem supplies a physical square of side at least $\sqrt2(9/500)(3+e)$ independently of refinement.

**Upgrading the witness to degree one.** On a selected interval the verifier constructs secants on each graph from the two ends and the midpoint. Try all 16 assignments of the four square vertices to these two halves on their respective graphs. Use physical square variables in the *working frame*, with vertices
$$
(X,Y),\quad(X+a,Y+b),\quad(X+a-b,Y+a+b),\quad(X-b,Y+a).
$$
For each assigned line $v=m_j u+c_j$ let $z_j=v_j-m_j u_j-c_j$. If this affine map is invertible, use the box $|z_j|\le p$, with $p=(R-L)/512$. Check its 16 inverse images of box vertices. Every assigned abscissa must lie strictly inside its own interval on the entire box and $a$ must stay positive. Require each true graph to deviate from its line by strictly less than $p$ there (checking the whole assigned segment suffices). On each opposite pair of faces its true residual has the sign of the corresponding $z_j$. Straight interpolation therefore gives mod-two degree one, even for a continuum of squares. The lower-left to lower-right ordering and the distinct upper graph select one labeling. This degree computation is independent of the external size theorem once the four supports are chosen.

For the family the successful supports are lower-left half, lower-right half, upper-right half, upper-left half. The four abscissa ranges on the full box are
$$
[8333/2688,8557/2688],\quad[28201/5376,28649/5376],\quad
[12947/2688,13171/2688],\quad[14359/5376,14807/5376].
$$
The secants are identical for every $n$, and $1/(2n)<p=1/64$ for every $n\ge64$. Exact linear optimization on the box gives $a\ge3805/1792$ and $b\ge681/1792$. Consequently the **physical squared side is at least $7470893/802816>9$** throughout the marked region. Its odd degree therefore protects a square above the fixed physical cutoff 3. The same calculation works for any continuous pair of graphs strictly inside the four vertical error bands, including the limiting broken lines.

| Refinement | Polygon vertices | Selected interval | Local theorem, squared physical floor | Degree box, squared physical floor |
| ---: | ---: | --- | ---: | ---: |
| 64 | 2,068 | $[0,8]$ | $48149721/8192000000$ | $7470893/802816$ |
| 128 | 4,116 | $[0,8]$ | $480249/81920000$ | $7470893/802816$ |

The [exact output](../results/v317/neck-lens.json) gives the endpoint gaps, signed margins and support intervals. The tests also cover cyclic relabeling, reversing traversal, a false cap witness, a singular refinement and invalid polygon contacts.

## What survives when both diagonal frames fold on every edge?

Once one has a verified affine box, the target arcs themselves can fold. In the working frame let $[A_j,B_j]$ be the full abscissa range of corner $j$ in the box. For each support find a target traversal strictly inside $|v-m_j u-c_j|<p$, with one endpoint abscissa below $A_j$ and the other above $B_j$. Each traversal may cross its vertical slices many times. This suffices for the same side floor and an odd **marked** degree.

To prove it, use the four box coordinates and four independent arc parameters. Set each of the four corner positions equal to its corresponding arc position, giving eight scalar equations. Deform the traversals to straight portions of their support lines, using linear interpolation of parametrizations and retaining the endpoint abscissas. The open strip is convex. A box face cannot contain a zero because one corner there has line residual exactly $p$ in absolute value, while its entire traversal stays strictly inside the strip. An arc parameter endpoint cannot contain a zero by the verified abscissa bounds. At the straight end the unique square coordinates solve the invertible four-line system, and each remaining parameter is uniquely determined by its nonconstant linear abscissa. That zero is regular by elimination. Degree one transfers to the folded target. Only existence is needed if other target portions also cross the same strips.

The exact scan takes an unmarked target and searches consecutive vertices in each open strip for a traversal passing both ends of the required abscissa range. Convexity certifies the intervening segments. This supplies parameters from the actual target independently of its registration with the reference. Subdividing long edges may expose further traversals missed by a vertex scan.

For a stress test take $n$ a multiple of 64. The preceding grid fixture has lattice unit $G=1/(4n)$ in the physical frame. Subdivide its straight edges into lengths between $G$ and $4G$ on that lattice. On the middle third of each short edge insert a rectangular detour to its left of depth $\tau=G/24$, returning to the same edge. The resulting polygon is Jordan. Nonincident lattice segments are at least $G$ apart, and detours on incident perpendicular edges stay at least $G/3$ from their shared endpoint in the longitudinal coordinate, beyond the depth of the other detour. The collinear consecutive case has disjoint middle thirds.

Each replacement contains a reversal of projection in **each** diagonal frame. A monotone run in either frame can therefore use parts of at most two successive blocks, with abscissa span at most $4G+2\tau$. This tends to zero. Nevertheless, in the working frame a detour changes each coordinate from a point of its old segment by at most $\tau/2$. The nominal support-line residual along each of the four central half-arcs stays bounded by $1/(2n)+(1+1/2)\tau/2<1/64$. Its ends retain abscissas 0 and 4, or 4 and 8, beyond the box ranges. The traversal theorem therefore keeps the same physical side bound **above 3 for every such refinement**. The detoured polygons converge to the two nominal broken lines with the exterior attachments in registered distance: the original teeth and every replacement block have vanishing diameter and corresponding block endpoints are retained.

The replay transfers the box from $n=128$ to a separately folded $n=64$ target with **30,740 vertices**. The longest monotone span in either diagonal frame is only $65/6144$ in working units. All four traversals are recovered directly on that target. This use of a coarse reference succeeds where insisting on a single macroscopic raw monotone run would fail.

## Remaining obstacle

This gives a geometric selector and a size floor for curves having a suitable local pair of graph arcs. General Jordan curves may contain folds at every scale, and even a grid approximation may have monotone runs that all become short in both frames. The traversal test handles that failure when a coarse reference supplies four suitable support lines. Finding such a reference, or a different certificate of positive size, from an arbitrary Jordan input remains open. The [v315 cutoff lemma](v315-cutoff-payment.md) would suffice if the failures could be paired or converted to simultaneous near contacts at one positive size.

[Checker](../code/neck_lens_v317.py) · [Regression tests](../tests/test_neck_lens_v317.py) · [Noncollapse goal](general-noncollapse-goal.md)
