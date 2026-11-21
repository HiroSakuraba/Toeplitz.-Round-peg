# v314. Local square degree after an unpaired crossing

Benjamin John Schulz · November 2026

At the asymmetric endpoint of [v313](v313-half-turn-degree.md), the degree counting all geometric squares of side greater than $2/5$ is even. Changing the reflection center preserves that global degree. Restricting the **centers of the squares counted** gives a different outcome. The central strip and the upper protrusion have odd degree separately. Both force a square, even under further nonsymmetric Jordan perturbations. The central square has side at least $31\sqrt2/32$.

We give a general formula for localizing the degree at a half-turn, then two ways to keep the resulting center and size boundaries free of squares. One uses exact contact equations and yields three disjoint size intervals. The other uses elementary span inequalities and allows arbitrary rough changes in two prescribed parts of a rectangle.

## Local half-turn formula

Write $R$ for counterclockwise rotation by $\pi/2$. Coordinates for a geometric square are its center $C$, half-diagonal $\rho>0$, and angle $\theta$ modulo $\pi/2$. Its vertices are $x_j=C+\rho R^ju_\theta$, $0\le j<4$, where $u_\theta=(\cos\theta,\sin\theta)$. Its side is $\sqrt2\rho$. Evaluating a continuous real field $f$ at the four vertices gives a section of the rank-four real bundle with cyclic permutation at the angle seam. For a Jordan curve use signed distance, positive outside. A zero represents one geometric square regardless of the choice of first vertex. The mod-two degree on a compact region with zero-free boundary is the parity of the zeros of a sufficiently close transverse section. It applies even if the original section has a continuum of zeros.

Let $U$ be a bounded half-turn-invariant center region about $O$, with piecewise smooth boundary and $O$ in its interior. Suppose $0<a<b$, $f(2O-x)=f(x)$ and the section is zero-free on **every boundary face** of $\overline U\times[a,b]\times(\mathbb R/(\pi/2)\mathbb Z)$. Define on each end circle
$$
w_r(\theta)=(f(O+r u_\theta), f(O+r R u_\theta)),\qquad 0\le\theta\le\pi/2.
$$
The endpoint equals $S w_r(0)$, where $S(x,y)=(y,x)$. Join $e=(1,1)$ to $w_r(0)$ by a path $p$ in the punctured plane. Follow $p$, then $w_r$, then the reverse of $Sp$. Let $\tau(r)$ be the winding of that closed walk modulo two. Replacing $p$ changes the integer winding by twice a winding, since $S$ reverses orientation. Rotating the starting angle gives a continuous family of twisted loops and preserves $\tau$. A strictly exterior-positive circle has $\tau=0$. A strictly interior-negative circle has $\tau=1$, as one sees by contracting to $(-1,-1)$ and connecting around one side of the origin.

**Theorem.** The local geometric-square degree is $\tau(a)+\tau(b)\pmod2$.

**Proof.** On configurations the half-turn reflects $C$ about $O$ and acts trivially on the size and angle class. Identify the fibers by shifting the four entries by two. The fixed base is the cylinder $C=O$. Over it the fixed fibers have coordinates $(x,y,x,y)$ and the complementary fibers have coordinates $(x,y,-x,-y)$. Both the normal plane to the cylinder and the complementary fiber carry minus the identity.

Approximate the section smoothly and equivariantly (average a close smooth approximation with its transform), keeping the positive boundary margin. Perturb its restriction to the fixed cylinder transversely within the fixed subbundle, extending local bump perturbations equivariantly. The fixed zeros are now finite and interior. At each, a local perturbation linear in $C-O$ in the complementary fiber can make the normal derivative invertible. It vanishes on the fixed cylinder. Equivariance makes the other two derivative blocks between unlike eigenspaces zero, so the full four-dimensional zero is regular. Retain small disjoint neighborhoods of these regular zeros. The rest of the fixed locus is zero-free away from these neighborhoods, and the remaining zeros are in the free part of the action. Perturb there downstairs in the quotient and lift. Each resulting off-center zero has a distinct mate. Thus the degree equals the parity of the zeros on the fixed cylinder. Only the proof perturbations of sections need be smooth.

Cut the fixed cylinder along an angular seam missing the transverse zeros. In this rectangle the restricted equation takes values in an ordinary plane. Its zero parity is the boundary winding modulo two. Write $h$ for the radial trace from the inner to the outer radius on the starting seam. The other seam has trace $Sh$. Use connector $p$ on the inner boundary and $p$ followed by $h$ on the outer boundary. Subtracting the two closed end walks cancels the common connector terms and gives the rectangle boundary winding up to sign. Reduction modulo two proves the formula. The end windings agree with those before perturbation by the boundary margin. An invariant center region disjoint from $O$ has even degree by the free-pair argument alone.

**Transfer to a nonsymmetric curve.** Choose a trial point $O$ and an invariant center region as above. At every configuration on its size and center faces set $A_j=f(x_j)$ and $B_j=f(2O-x_j)$. If the closed segment $[A,(A+B)/2]$ misses zero on all these faces, the straight interpolation to $g(x)=(f(x)+f(2O-x))/2$ preserves the local degree. Compute the two twists for $g$. A difference forces a square of the original curve in the open region. The midpoint field need not define a Jordan curve. For an already even field a difference of twists alone forces a centered square by the cylinder argument. The center walls are needed for the full local degree and its transfer.

A second transfer test uses a reference curve directly. Say that curves are registered within $\epsilon$ if they have Jordan parametrizations with coordinate sup distance at most $\epsilon$ at each common parameter. Outside the closed coordinate-sup $\epsilon$-tube of the reference, the straight homotopy of parametrizations avoids the test point, so the two winding numbers, hence the inside/outside signs, agree. Interpolation of signed distances has all its point zeros in that tube even if intermediate parametrizations self-intersect. Excluding boundary configurations whose four vertices all lie in the tube therefore suffices.

## Three disjoint size intervals beyond the centered inradius

For a reproducible symmetric example, take rays $(3,0),(3,1),(2,2),(1,3),(0,3)$ in that order. Multiply them respectively by $3,1,1,3,1$ in the first quadrant. In the second quadrant rotate a copy by $R$ and multiply its five original rays respectively by $1,2,3,2,3$. Use the common endpoint only once. Complete the polygon by a half-turn, again omitting duplicate endpoints. The resulting 16-vertex Jordan polygon has $(0,3)$ on its boundary. All four scan radii below exceed its centered inradius.

| Half-diagonal $r$ | $\tau(r)$ | Certified quarter-circle subarcs |
| ---: | ---: | ---: |
| 4 | 1 | 6 |
| 5 | 0 | 6 |
| $33/5$ | 1 | 12 |
| 7 | 0 | 8 |

Parameterize the unit quarter-circle by $((1-t^2)/(1+t^2),2t/(1+t^2))$, $0\le t\le1$. Each subarc and its rotated copy lie in their corresponding endpoint boxes by coordinate monotonicity. Subdivide until at least one box misses every polygon edge, certified by coordinate or edge-normal separation. Thus one evaluation coordinate has constant strict sign on each accepted subarc. Replacing its evaluation path by the segment of endpoint sign vectors preserves winding: interpolate each endpoint to its sign vector, then straighten within the corresponding open half-plane. Zero signs in the other coordinate are allowed. If the depth budget is exhausted the verifier reports an unresolved cutoff.

**Center and radius clearance.** Take $|C|_\infty\le h=1/8$ and $\epsilon=1/128$. For each tuple of four edges $E_j$ supporting candidate tube contacts the common variable $v=\rho u_\theta$ must satisfy
$$
v\in\bigcap_{j=0}^3\big(R^{-j}E_j+[-h-\epsilon,h+\epsilon]^2\big).
$$
Intersect the six support half-planes of each enlarged segment (the coordinate and segment normals with both signs). Exact clipping leaves 12 tuples, forming three orbits of four by cyclic relabeling. For each surviving tuple choose an edge endpoint $p_j$ and normal $n_j$. The four nominal equations are
$$
n_j\cdot C+(R^{-j}n_j)\cdot v=n_j\cdot p_j.
$$
Their matrices are invertible. At the exact roots the center is zero, the four contacts are strictly inside their respective edges, and the three squared radii are $457/25,2088/49,2169/49$. On the tube each row has error at most $\epsilon\|n_j\|_1$. Multiplication by the entrywise absolute inverse matrix bounds both coordinates of $C$ and $v$. The largest center bound is $15/128<1/8$. The three squared-radius enclosures are contained respectively in $(16,25)$, $(25,(33/5)^2)$ and $((33/5)^2,49)$. The [exact endpoints and displaced sample](../results/v314/local-twist.json) allow direct comparison.

Every registered Jordan perturbation of size at most $1/128$ therefore contains at least three distinct squares with centers in the open box, one with half-diagonal in **each** of $(4,5)$, $(5,33/5)$ and $(33/5,7)$. All center and radius faces remain clear throughout signed-distance interpolation. Each interval has odd degree by its two end twists. In particular the largest of these three has side greater than $33\sqrt2/5$. This certifies existence without following an individual square through the perturbation. As a cross-check the program solves the four-contact equations again on an explicit vertex-displaced polygon and finds exactly three squares in this box, all with displaced centers.

## A protected central square under arbitrary cap changes

Start with $[-4,4]\times[-1,1]$. Allow a Jordan boundary with arbitrary changes in the two closed cap boxes
$$
|x|\le 3/2,\qquad 1\le |y|\le H=25/16,
$$
whose boundary and inside/outside status agree with the rectangle off the boxes. Then allow a further registered perturbation of at most $\eta=1/32$. The final curve has a square with $|C_y|<1/2$ and **side at least $31\sqrt2/32$**. The shape of the boundary inside a cap and the relation between the two caps are unrestricted subject to these assumptions.

For the first homotopy interpolate the signed distances of the rectangle and the cap-modified curve. For the second interpolate to the registered target. Throughout, every point zero lies within $\eta$ of the rectangle boundary or cap boxes. Call this union of thickened sets $T$. We use the center slab $|C_y|<1/2$ and side cutoff $s=2/5$, truncating $|C_x|$ and the upper size far beyond any possible zeros. Write $w$ for the square's horizontal span, which equals its vertical span. For side $\ell$ one has $\ell\le w\le\sqrt2\ell$.

* At $\ell=s$ in the closed slab, all vertices have $|y|\le1/2+s<1-\eta$. Contacts with $T$ must therefore be in the two vertical side strips. One strip gives $w\le2\eta<s$ and touching both requires $w\ge8-2\eta>\sqrt2s$.
* On $|C_y|=1/2$, the overall vertical bound gives $w\le2(H+\eta-1/2)$. A thick cap and a vertical side strip would require $w\ge4-3/2-2\eta$, strictly larger. With no vertical-side contact, all contacts are in the positive or negative height bands $[1-\eta,H+\eta]$ and its reflection. In one band $|C_y|\ge1-\eta$. In both bands the midpoint of the extrema gives $|C_y|\le(H-1)/2+\eta<1/2$.
* It remains to consider contacts supported only by the thickened rectangle boundary. On the center wall one cannot reach both vertical sides by the overall span bound, or both horizontal rows since now $w\le2(1+\eta-1/2)<2-2\eta$. At most two adjacent strips suffice. Use the outer lines of those strips as coordinate axes, reflecting signs so that all vertices lie in the positive quadrant. Each vertex is within $2\eta$ of at least one axis. Each center coordinate is at least $w/2$. Among four quarter-turn-related displacements one is northeast of the center, so that vertex has both coordinates at least $w/2$. Hence $w\le4\eta<s$. The same bound holds for a single strip.

These estimates exclude the size and center walls along both homotopies. For the reference rectangle the local twist theorem gives odd degree. The circle of radius $s/\sqrt2$ lies strictly inside and a sufficiently large outer circle lies outside. Thus some square above the cutoff survives in the open slab.

To improve its side bound, consider a square there with $\ell>s$. If it visits both height bands, $w\ge2-2\eta$. If it visits a cap and a vertical strip, $w\ge4-3/2-2\eta>2-2\eta$. A single height band cannot contain a square centered in the slab. For rectangle-only contacts, two opposite strips likewise give $w\ge2-2\eta$, while the adjacent-strip estimate would contradict $\ell>s$. Hence $\ell\ge(2-2\eta)/\sqrt2=31\sqrt2/32$.

## The asymmetric endpoint with even total degree

Take the two protrusions of v313 with upper scale 1 and lower scale $9/10$. At side cutoff $2/5$ the global degree is zero by its single-crossing calculation. The cap theorem gives local degree one in the central slab. In the exterior slabs $|C_y|>1/2$, the bound $w\le2(H-|C_y|)$ prevents reaching a cap and a vertical side. Opposite height bands would put the center inside the central slab. Squares supported by the rectangle alone have side two and center height zero, as shown in [v307](v307-convex-shell-square-gap.md). The remaining exterior squares are each supported on one protrusion and its horizontal baseline. Their span is at most the height of that protrusion. Any such square must touch the protruding part, so all contacts lie within the normalized profile used for the v313 census (the maximal absolute horizontal offset is at most $7/16+17/32<1$ before scaling). The upper protrusion has one regular square above cutoff, with squared side $533/2890$. The lower has none.

| Center height | Local degree above side $2/5$ |
| --- | ---: |
| $C_y<-1/2$ | 0 |
| $-1/2<C_y<1/2$ | 1 |
| $C_y>1/2$ | 1 |

Re-running the v313 closed-tube test **at this asymmetric endpoint** excludes side $2/5$ at coordinate-sup halfwidth $1/8192$. It closes all 2,164 leaves, reaching depth 28. The center walls stay clear by the cap estimates. Thus every Jordan curve registered within $1/8192$ of this endpoint has two separately forced squares: one in the central slab with side at least $31\sqrt2/32$, and another above it with side greater than $2/5$. Additional squares are possible and global parity stays even.

The open problem for this approach is to obtain a center region and two positive size boundaries with odd relative degree, or a square on a tested boundary, for an arbitrary Jordan curve. The exact examples show what information survives after a global cancellation. The center-wall exclusions above still use particular geometry.

[Exact verifier](../code/local_twist_v314.py) · [Regression tests](../tests/test_local_twist_v314.py) · [Rational results](../results/v314/local-twist.json)
