# v307 — A square-size gap in a convex boundary layer

Benjamin John Schulz · November 2026

Let $K$ be a compact convex body in the plane and fix $B(O,r)\subset K\subset B(O,R)$, using closed disks with $r>0$. The layer $S_\delta$ consists of the points of $K$ whose distance from $\partial K$ is at most $\delta$. We use $\rho$ for a square's half-diagonal; its side is $\sqrt2\rho$.

**Convex-layer lemma.** If all four vertices of any genuine square belong to $S_\delta$, then
$$
G(\rho)=\rho^2-(r-\delta)\rho+\delta R\ \ge 0.
$$
**Large-square theorem.** Let $\Gamma\subset S_\delta$ be a Jordan curve enclosing $O$. If $0<q<r$ and
$$
\delta < \frac{q(r-q)}{R+q},
$$
then $\Gamma$ inscribes a square with half-diagonal strictly greater than $q$. In fact its half-diagonal is at least
$$
\rho_+=\frac{r-\delta+\sqrt{(r-\delta)^2-4\delta R}}2 .
$$
There are no differentiability, graph, parametrization or isolated-square hypotheses on the input. The inequalities concern one enclosing convex body; they need not hold for every Jordan curve. For an intrinsic sufficient test take $K=\operatorname{conv}\Gamma$ and measure the greatest depth of the curve below its hull boundary.

## 1. Proof of the gap

Write $C$ for the square center and use axes along its sides. Reverse their signs as necessary so that both coordinates of $O-C$ are nonnegative. Take the square corner $p=C+(s/2,s/2)$ in these coordinates, where $s=\sqrt2\rho$. There is a unit vector $n$ with $0\le h_K(n)-n\cdot p\le\delta$, where $h_K(n)=\max_{x\in K}n\cdot x$. Indeed, take a supporting line at a nearest boundary point. The ball about $p$ of radius $\operatorname{dist}(p,\partial K)$ lies inside $K$ and touches there, giving equality between this radius and the supporting-line slack. This also covers a boundary point by taking any support there.

The adjacent vertices $p-se_i$ lie in $K$. Hence each $n_i\ge-\delta/s$. By convexity $C\in K$, so $|O-C|\le R$. Using both nonnegative coordinates of $O-C$ gives
$$
n\cdot(C-O)\le (\delta/s)\|O-C\|_1\le\delta R/\rho .
$$
Also $n\cdot(p-C)\le\rho$. The inball gives $h_K(n)\ge n\cdot O+r$. Combining the three estimates yields $r-\delta\le\rho+\delta R/\rho$, as claimed. Only the corner facing $O$ needs small support slack; membership in $K$ suffices for its two neighbors. For $\delta=0$, every inscribed square on the boundary of a convex body containing the disk has half-diagonal at least $r$. A disk realizes equality.

## 2. Why a square occurs on the large side

The gap alone permits arbitrarily small squares. Here is the existence step, with the labeling convention explicit.

Let $E$ be the ellipse centered at $O$ of semiaxes $2r,r$. Deform it to $K$ through the Minkowski sums $K_t=(1-t)E+tK$. They all contain the same inball and fit in $B(O,M)$, $M=\max(2r,R)$. Use defining functions $f_t(x)=\max_{|n|=1}(n\cdot x-h_{K_t}(n))$, which vary continuously since support functions interpolate linearly. They are negative inside, zero on the boundary and positive outside. By the zero-depth lemma, none of these boundaries has a square with $\rho=q$.

Next interpolate linearly from $f_K$ to signed distance from $\Gamma$, negative on its bounded interior. Both functions are positive outside $K$ and negative on the open core $\{x\in K:\operatorname{dist}(x,\partial K)>\delta\}$. Indeed this core is convex (an intersection of inward support half-planes, equivalently the interior of the eroded body); it contains $O$ because $\delta<r$, and avoids $\Gamma$. The exterior of $K$ lies in the unbounded component of the complement of $\Gamma$. Every zero of an intermediate function is consequently in $S_\delta$. The intermediate zero sets need not be curves. Since $G(q)<0$, none can contain all four vertices of a square at the cutoff.

For precision, parametrize a geometric square by $C,\rho,\theta\bmod \pi/2$, its ordered trial corners on an angle lift being $C+\rho e_{\theta+j\pi/2}$ for $j=0,1,2,3$. Use the compact four-dimensional domain $|C-O|\le L$, $q\le\rho\le L$, with $L>M$. Evaluate the defining function at all four corners. At the angle seam the four entries are cyclically permuted: they form a section of the rank-four real vector bundle obtained by this same permutation on the fiber. The section has no boundary zero. The inner face is excluded above; any four zeros have $C\in B(O,M)$ and $\rho\le M$ by averaging and taking an opposite pair. Counting on the angle quotient counts one geometric square, avoiding the fourfold overcount.

The parity of interior zeros of a transverse section is invariant under a homotopy nonzero on the boundary. One proof smooths and perturbs the section homotopy transversely, using local bundle charts and perturbations supported away from the zero-free boundary, after arbitrarily close smooth approximation there if necessary. The compact one-dimensional zero manifold has an even total number of boundary points. This works with mod-two counts regardless of the seam permutation. In particular, if the final continuous section were everywhere nonzero it would have a positive minimum norm and could be approximated with zero count zero.

The ellipse starts with one regular geometric square. To see regularity, translate $O$ to zero and use the defining function $x^TAx-1$ with diagonal positive $A$ having unequal entries. The two opposite-corner differences are $4\rho C^TAu$ and $4\rho C^TAv$, where $u=e_\theta,v=Ju$. They force $C=0$ with invertible center derivative. The two pair averages then give $\rho^2u^TAu=\rho^2v^TAv=1$. Equality of the two quadratic forms has a simple angle root at $\pi/4$ modulo $\pi/2$, and the sum has nonzero radius derivative. Passing from the quadratic defining function to $f_E$ by convex interpolation leaves the same ellipse as its zero set, including on the boundary cutoff. Thus the parity starts odd and forces a final zero at $\rho>q$. The gap implies $\rho\ge\rho_+$ since $q$ lies strictly between the two real roots. This proves the theorem, including for nonisolated final squares.

## 3. Exact coarse example and controls

Use $K=[-4,4]\times[-1,1]$, $O=0$, $r=1$, $R=5$, $\delta=1/32$, $q=3/4$. The cutoff polynomial equals $-1/128$. Every Jordan loop in this inner boundary layer enclosing $O$ therefore has square side greater than $3\sqrt2/4$. The stronger bound is $\rho\ge(31+\sqrt{321})/64$. The same lower bound holds on every eligible refinement and limit. For a simple polygon fixture, traverse the box counterclockwise from $(-4,-1)$, replacing the upper side by the path $(4,1),(1,1),(1/2,31/32),(-1,1),(-4,1)$. All vertices of its hull are retained and every modified point is in the layer. The geometry allows arbitrary folds within the layer so long as the loop is Jordan and separates the core. This fixture's diameter is at least eight and strip width is two, excluding the $1+\sqrt2$ circular-annulus criterion at every interior center.

The unmodified long rectangle illustrates why a global argument matters. It has no isolated geometric square. If four square corners met four distinct sides without sharing an edge, the square would span both the full width 8 and the full height 2, whereas its axis-aligned bounding box has equal width and height. Otherwise two corners share a supporting edge, and must be adjacent because a diagonal cannot support a square. The square is then axis-aligned. To touch the boundary at all four corners it must span the height (spanning the width would require height 8). The side-two squares translate continuously across the box, including up to both ends. The layer theorem applies before choosing any member of this family.

Enclosure is essential. In the unit disk take the circle of radius $\delta/4$ with center $(1-\delta/2,0)$ for $\delta=1/32$. Its entire boundary lies in the layer; its squares have half-diagonal $\delta/4$. It does not enclose $O$. Positive layer thickness also allows small squares on loops that *do* enclose $O$: cut a rectangular notch of width and depth $\delta$ into the top side of the long box. All four notch corners lie on the loop and form a square of side $\delta$. The theorem forces a second, large square there. These two tests distinguish a geometric gap from existence on its large side.

The [rational replay](../code/convex_shell_v307.py) checks the cutoff, radical and two controls; the section argument supplies existence. The executable also accepts a rational Jordan polygon and computes its convex hull. For each curve segment it proves that the union of the hull support strips covers the entire parameter interval. It uses a rational lower bound on the Euclidean length of each hull edge, so a rejected strip test is inconclusive. Run `python code/convex_shell_v307.py` and `python -m unittest discover -s tests -p 'test_convex_shell_v307.py' -v`. With `--input input.json`, supply keys `polygon` (cyclic coordinate pairs), `center` (one coordinate pair), and rational `r`, `R`, `depth`, `q`; fractions may be quoted. The proof gives squares; a consistently labeled near-square rectangle continuum would require an additional argument, since cyclic relabeling exchanges reciprocal ratios.

## Context and remaining question

Matschke identifies a scale bound in terms of an interior disk as a possible route to the general conjecture, and proves an existence and size result in a circular annulus. Chambers proves persistence near a smooth reference using a curvature and winding condition. A public discussion of tangential convex polygons also uses the square corner facing the disk. Here the support slack allows a nonconvex target, and a quotient parity argument selects the large side of the gap. No priority or optimal-thickness claim is made. [Matschke, survey](https://www.math.u-bordeaux.fr/~bmatschke/SurveyOnSquarePeg_extension14.pdf) · [Chambers, *On the square peg problem*](https://arxiv.org/abs/2203.02613) · [Convex polygon discussion](https://math.stackexchange.com/questions/5140338/does-every-convex-polygon-have-only-large-inscribed-squares)

The proof suggests a narrower next question for curves with deep indentations. At a single candidate scale, can one exclude a square on the deformation using just the supporting-line slack of its corner facing $O$? Section 1 permits large slack elsewhere. A deformation-independent way to certify that directional condition for general curves is still missing.
