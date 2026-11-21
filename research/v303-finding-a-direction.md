# v303 — Finding the direction from an unmarked polygon

Benjamin John Schulz · November 2026

The v302 proof needs a direction in which the two sides of a near-square rectangle cannot exchange their horizontal and vertical roles. Here is a finite exact test for finding such a direction on a rational polygon supplied in cyclic vertex order. It also recovers a source parametrization, which is needed to interpret the density. A smooth fourfold-symmetric example shows that some curves admit a tie in **every** fixed direction.

## A certificate from two vertices

Orient the input polygon counterclockwise by its signed area. Try every ordered pair of distinct vertices $A,B$. With $d=B-A$ use coordinates
$$
x=\frac{(p-A)\cdot d}{d\cdot d},\qquad y=\frac{\det(d,p-A)}{d\cdot d}.
$$
This is an orientation-preserving similarity that takes the candidate endpoints to $(0,0)$ and $(1,0)$. No square contacts are inputs. Traverse the counterclockwise arc from $A$ to $B$ as the lower graph and the complementary arc backwards as the upper graph. Require $x$ to increase strictly on both traversals. If it does, set the branch parameter $u=x$ and give each branch half the source circle as in v302.

For each graph segment check its endpoints and any interior point with $x=1/2$. At each such point form $\delta=y-P_y(x)$ using the lower or upper v302 seed height and the corresponding allowed interval $[L,H]$. Return the minimum over all checked points of
$$
\min\{1/1000,\ \delta-L,\ H-\delta,\ E-|\delta|\},\qquad E=79/1500 ,
$$
provided it is nonnegative. Otherwise reject this pair. On each side of the crest, $\delta$ is affine, so these tests certify the entire curve. If the returned number is $b$, any registered continuous perturbation of Euclidean size at most $b$ in normalized coordinates has $|x-u|\le 1/1000$, height within $[L,H]$ relative to the original source, and Euclidean distance at most $E$ from the seed. These are precisely the inputs of [v302, §1](v302-connected-rough-rectangle-crossing.md). The perturbed graph may double back after the source has been assigned. The target must still be a Jordan embedding. In the original coordinates the allowed displacement is $b|B-A|$ and every side in the spanning continuum has length greater than $0.407|B-A|$.

The search costs $O(n^3)$ rational operations on $n$ vertices, including at most $n(n-1)$ candidate frames. A missing certificate leaves the fixed-direction question undecided for that input. Reversal or cyclic reindexing of the vertex list preserves the set of geometric answers. Orientation-preserving similarities take every successful pair to the corresponding successful pair and preserve the dimensionless budget.

## Positive example with independent arcs

The exact seven-vertex example has the following cyclic vertices.

$$
(0,0),\ (1/4,-19/100),\ (1/2,-3/8),\ (1,0),\
(3/4,39/100),\ (1/2,77/100),\ (1/4,3/8).
$$

Exactly one of the 42 ordered pairs passes. Its chord is $(0,0)$ to $(1,0)$ and its normalized perturbation budget is $1/1000$. After applying $(x,y)\mapsto(7+3x-4y,-2+4x+3y)$, the search returns $(7,-2)$ to $(10,2)$, yielding an actual displacement budget $1/200$ and side floor greater than $2.035$. The rotated input supplies its vertices alone. The upper and lower graphs have independent deviations: for instance, at $x=1/2$ the lower point equals its reference height while the upper one is $1/50$ higher. [Exact outputs](../results/v303/polygon-directions.json)

## A fixed direction can fail on every axis

Take the smooth polar Jordan curve of radius $f(\theta)=(1+\cos(4\theta)/100)/2$ around $O=0$. It lies in the v302 annulus $0.29<|p|<0.69$. Since $f(\theta+\pi/2)=f(\theta)>0$, the four points at $\theta+j\pi/2$ are the corners of a square centered at $O$. Its first side makes angle $\theta+3\pi/4$. For any proposed frame angle $\alpha$, choosing $\theta=\alpha-\pi/2$ puts that side at $45^\circ$ relative to the frame. This square is in the protected region $|C-O|<\rho$. The tie persists for every choice of frame, despite the arbitrarily thin annulus obtainable by decreasing the amplitude. The same construction applies to the boundary of a square about its center, whose radial function also has period $\pi/2$. Here the obstruction follows from the whole family of squares, independently of the search returning no certificate.

## The same curve still has a crossing

For $0<\theta<\pi/2$, put $a=f(\theta)\cos\theta$ and $b=f(\theta)\sin\theta$. The points $(\pm a,\pm b)$ are all on the curve by its reflection symmetries. They give one continuous branch of actual axis-aligned rectangles with side ratio $b/a=\tan\theta$. It spans every positive ratio and passes through a square at $\theta=\pi/4$. Thus excluding ties from the entire protected region is stronger than what existence on a single component requires.

One can locate the branch intersections exactly in the space of **centered** rectangles. Let $\theta$ and $\theta+\phi$, $0<\phi<\pi$, be two successive corner directions. The two antipodal points complete a rectangle precisely when the two radii agree. The equation factors as
$$
f(\theta+\phi)-f(\theta)
=-\tfrac{1}{100}\sin(4\theta+2\phi)\sin(2\phi)=0 .
$$
The square circle $\phi=\pi/2$ meets the other branches at $\theta\in(\pi/4)\mathbb Z$. These eight labeled junctions give two geometric squares modulo a quarter turn. Along the square circle, advancing $\theta$ by $\pi/2$ returns to the same geometric square with the corners cyclically permuted. The ratio of second to first side is $\cot(\phi/2)$; cyclic relabeling exchanges it with its reciprocal. This describes the full centered locus, while off-center configurations may also occur.

**Next question.** Can one transfer the degree between overlapping choices of frame through a junction where the side labels exchange? The factored centered example gives explicit branches and a known permutation against which to check such a rule. A subsequent test would perturb the fourfold symmetry and examine whether a spanning component survives after the square circle splits. For a general Jordan curve, collapse outside the present annulus class is a separate open issue.

## Replay

Use `python code/direction_search_v303.py --polygon polygon.json` with an array of cyclic coordinate pairs. Integers, decimal numbers and quoted fractions are accepted. The result lists all successful ordered tip pairs, their squared separations and their certified normalized budgets. For the examples run `python code/direction_search_v303.py --check` and `python -m unittest discover -s tests -p 'test_direction_search_v303.py' -v`.
