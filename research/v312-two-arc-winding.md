# v312. A boundary winding test for squares on rough arcs

Benjamin John Schulz · November 2026

Two marked arcs suffice to test a candidate region of square configurations on an arbitrary Jordan curve. The test only evaluates trials on the boundary of a two-parameter rectangle. A second version prescribes all four contact arcs and admits a sharp corridor criterion. On the accumulating-notch family, both formerly straight contact arcs can fold backwards throughout strips of halfwidth $1/80$; a square of side at least $39/40$ survives at every stage and on the limit.

## 1. Complete an oriented side and take a winding number

Let $Q=J(s)$ and $B=K(t)$ range over disjoint compact subarcs of a Jordan curve $\Gamma$. Write $R$ for a counterclockwise quarter turn. The remaining vertices for this orientation are given uniquely by
$$
P=Q+R(B-Q),\qquad D=B+R(B-Q).
$$
Let $\phi$ be signed distance to $\Gamma$, positive inside. On the closed parameter rectangle define $F=(\phi(P),\phi(D))$. If $F$ avoids zero on the perimeter, the integer winding of the perimeter image about zero is defined. **Nonzero winding forces an inscribed square of side at least $\operatorname{dist}(J,K)$.** Indeed winding equals planar Brouwer degree. A zero gives four genuine vertices, and the positive separation gives a nonzero side. No regularity beyond an embedding is needed. To obtain the square on the opposite side reverse the two roles. For finitely many interior regular zeros, the winding modulo two counts the squares in this oriented side window. An individual isolated transverse square in a finite census can always be detected by sufficiently small windows about an oriented adjacent pair: that pair uniquely determines this completion, so no other census zero survives in the shrinking rectangle. The widths need not stay bounded below across different curves.

Only four one-parameter sweeps are required. On each one, one side endpoint stays fixed and the other traverses its arc. For polygons in general position, count crossings of $P$ through $\Gamma$ for which $D$ is inside, with the sign of crossing of the positive vertical ray in the $(\phi(P),\phi(D))$ plane. This computes winding without enumerating interior squares. The supplied verifier uses exact segment geometry and rejects ambiguous ray events. The continuous winding theorem also covers zero-free perimeters for which this particular ray is unsuitable.

**Quantitative stability.** On the perimeter put
$$
m=\min\max\big(\operatorname{dist}(P,\Gamma),\operatorname{dist}(D,\Gamma)\big).
$$
A registered Jordan perturbation by Euclidean sup distance $\delta$ preserves winding if $(2+\sqrt2)\delta<m$. Its side floor is $\operatorname{dist}(J,K)-2\delta$ when positive. To prove this, interpolate the positions of the registered side endpoints and independently interpolate the two signed-distance functions. Every point zero of the interpolated distance function lies in the closed $\delta$-tube of the reference curve: outside the tube the registered loops have the same inside/outside status by winding invariance. Both completed vertices move by at most $(1+\sqrt2)\delta$, since the coefficient maps are $R$ and $1\pm R$. The vertex with reference distance at least $m$ therefore stays outside the tube on each perimeter trial. The interpolated two-component map avoids zero there.

**Exact triangle example.** On $(0,0),(4,0),(1,3)$ choose $Q=(4/7+a,0)$ and $B=(16/7+b,0)$ with $|a|,|b|\le1/16$. At the two completions the signed line equations for the left and right sloping edges are $4a-b$ and $a-2b$. Their matrix has determinant $-7$; the unique zero has winding $-1$. The inverse matrix has sup norm $5/7$. On the perimeter at least one absolute residual is at least $7/80$. The corresponding Euclidean distance to its supporting line exceeds $1/40$ after normalization by at most $\sqrt{10}$. The other edges cannot be closer than $1/40$ to that vertex: the heights are at least $12/7-1/8$, the right-edge residual at $P$ is at least $12/7-1/16$, and the left-edge residual at $D$ is at least $36/7-3/16$. Thus $m>1/40$. Any registered Euclidean perturbation by at most $1/160$ has a square anchored on the two marked arcs with side at least $883/560$. This calculation uses only the two side windows, without specifying contacts on either sloping edge in the perturbed curve.

## 2. Four rough arcs and a sharp corridor criterion

To prescribe $P\in I$ and $D\in L$, use two separating curves instead of the whole boundary. Extend $I$ and $L$ independently to proper embedded lines and give each a continuous signed distance. It suffices that their added tails introduce no simultaneous false zero. In the construction below neither completion can reach its corresponding tails. The same perimeter-winding proof then gives a square meeting all four marked arcs. This is a local degree. Other squares elsewhere on the Jordan curve may contribute to different tests.

Fix $0\le e<h$, $2h+e<a<1-h$, and $a<1-e$. Place four arbitrary simple arcs in the following boxes, with each endpoint on its designated transverse end face. The middle of an arc can turn backwards. For the vertical and horizontal tail extensions of $I,L$ require that only the endpoints meet their longitudinal end faces.

| Arc | Horizontal range | Vertical range | End faces |
| --- | --- | --- | --- |
| $I$ | $[-e,e]$ | $[-a,a]$ | bottom and top |
| $J$ | $[1-e,1+e]$ | $[-h,h]$ | bottom and top |
| $K$ | $[1-h,1+h]$ | $[1-e,1+e]$ | left and right |
| $L$ | $[-a,a]$ | $[1-e,1+e]$ | left and right |

Assume they are subarcs of one Jordan curve. Set the two residuals positive on the right of $I$ and above $L$. The vertical coordinate of $P$ and horizontal coordinate of $D$ are always in $[-2h-e,2h+e]$. The completed vertices cannot reach artificial tails.

**Theorem.** If $h>3e$, the four-arc problem has degree one with $J$ oriented up and $K$ right. It forces a square of side at least $1-2e$. The same existence bound holds on closed arcs at $h=3e$.

In fact $P_x=Q_x+Q_y-B_y$. On the bottom face of $J$ it is at most $2e-h<-e$; on its top face it exceeds $e$. Likewise $D_y=B_x+B_y-Q_x$ lies below $1-e$ on the left face of $K$ and above $1+e$ on its right face. Straight interpolation of the residual map to the two centered parameter coordinates avoids zero on the perimeter. This gives winding $+1$. The horizontal separation of $P$ from $Q$ supplies the side floor. At equality, shrink each arc's transverse deviation from its centerline by a positive factor tending to one, keeping longitudinal coordinates. Each resulting four-arc problem satisfies the strict test. Compactness in the parameter rectangle gives a limiting square on the original closed arcs; continuation at intermediate stages uses only the four arcs.

**Sharpness.** For $h<3e$ take vertical $I$ at $x=e$, vertical $J$ at $x=1-e$, horizontal $K$ at $y=1+e$, and horizontal $L$ at $y=1+e$. The equation for $P_x$ would require $Q_y=3e$, outside $J$. At equality the square lies on the end of $J$; the strict boundary degree is then undefined. These four segments can be joined in cyclic order outside the boxes. The side bound is also attained: put both horizontal arcs at $1-e$ instead. Then $Q=(1-e,e)$ and $B=(1-e,1-e)$ give the unique marked square of side $1-2e$ in the strict regime.

**A winding certificate beyond the sign test.** For $e=1/16,h=1/4,a=5/8$, use straight $I,L$ on their centerlines and $J$ from $(1,-1/32)$ to $(1,1/4)$. Take $K$ through $(3/4,1),(13/16,15/16),(15/16,1),(5/4,1)$. On the first endpoint of $J$, $P_x$ takes both signs as $B$ moves. The perimeter nevertheless winds once; the independent four-edge calculation gives one marked square. Two further fixtures have winding zero, with respectively two and zero marked squares. This shows exactly what an inconclusive boundary test permits.

## 3. Four folded arcs, including an infinite-length limit

Use $e=1/16,h=1/4,a=5/8$. The [fixture and verifier](../code/rough_arc_winding_v312.py) specify rational six-vertex paths on each of the four arcs and six additional connecting vertices. Every marked arc reverses its longitudinal direction on an interior segment. All four contain transverse displacements greater than $1/512$ from their centerlines. The 30-vertex loop passes an exact Jordan simplicity check. The four boundary sweeps give three positive-ray hits on one face, winding $+1$. An independent exhaustive solution of the $5^4=625$ marked edge products finds **31** transverse squares, each with one vertex on each designated arc. Using the entire loop as both separators gives winding $-1$ on the same two side arcs. An independent search allowing both remaining vertices anywhere on the 30 edges also finds exactly 31 squares in this oriented side window. Neither count includes squares whose oriented side misses the window. The theorem protects at least one square with side at least $7/8$.

The same boxes also carry simple arcs of infinite length on **every open subarc**. On each marked segment $A+t v$, $0\le t\le1$, replace it by
$$
A+t v+\alpha Rv\,t(1-t)W(t),\quad
W(t)=\sum_{n\ge1}2^{-n}\operatorname{dist}(2^{n^2}t,\mathbb Z),
\qquad \alpha=1/4096 .
$$
The series converges uniformly and $0\le W\le1/2$. Projection along $v$ is strictly increasing on that segment. Every perturbed point lies in the thin diamond with vertices $A,A+v,A+v/2\pm\alpha Rv/2$. Exact convex-polygon tests verify that the 20 diamonds and the remaining connector segments meet only at intended common endpoints, and that each diamond remains in its marked corridor. Endpoints and the longitudinal reversals of the original paths persist.

For infinite variation, sample any interior interval on the mesh $2^{-n^2-1}$. All terms with index greater than $n$ vanish at mesh points. The $n$th term alternates between zero and $2^{-n-1}$. Lower terms have total Lipschitz constant at most $\sum_{k<n}2^{k^2-k}=o(2^{n^2-n})$, so eventually each mesh-step change has magnitude at least $2^{-n-2}$. The number of steps on the interval grows on the order of $2^{n^2}$. Multiplication by $t(1-t)$ preserves infinite variation on compact interior intervals, since its reciprocal there is bounded and Lipschitz. Thus each replacement has infinite length on every nontrivial subarc. Finite series truncations and the infinite Jordan limit all obey the same $7/8$ bound; no straight subarc is needed at any of the four contacts.

## 4. Return to the accumulating notches

Keep $I,L$ and the complement of $J,K$ from [v311](v311-separated-arc-parity.md) on any nonempty notch stage or its infinite limit. Replace the bottom and top marked segments by arbitrary simple arcs with the original endpoints, inside the same horizontal ranges and the height strips $|Q_y|\le e,\ |B_y-1|\le e$, where $e=1/80$. Assume the result is Jordan. These arcs may backtrack in $x$, with no parametrization or slope bound inside the strips.

For $P=(x,y)$, put $Q=(z,j)$, $B=(u,1+k)$. At a zero of the four-arc square equations,
$$
z=1+x+k-j,\quad u=1+x+y+k-2j,\quad D=(x+y-j,1+y+k-j).
$$
Along straight interpolation to any fixed linear traversal with the same endpoints, both arcs remain in their respective boxes. The two longitudinal predictions move by at most $2e,3e$ and the fourth-corner prediction by at most $2e$ in sup norm. The original face margins for $I,J,K,L$ are respectively $1/32,1/8,1/8,1/16$. The smallest remaining margin is $1/160$. Hence the degree of v311 stays one, and the vertical $Q,B$ separation is at least **$39/40$**. This permits transverse excursions $6.4$ times the earlier $1/512$ registration allowance on these two arcs.

For completeness, this degree still gives the combined v311 four-cut bit on finite regular censuses with admissible cuts. Write $v=z-x\ge3/8$. If $P,Q$ were opposite, a lower completion would have height at most $(1/8+e-3/8)/2<-1/16$; the downward adjacent completion has a corner at height at most $e-3/8$. Both miss the whole curve. In the upward completion the upper right corner $(z+y-j,j+v)$ is strictly between the vertical sides and above the bottom region. It is on the top region; set $k=0$ on its unchanged portion. Its predicted abscissa $u$ lies strictly inside $K$ by the $3e$ bound, so it belongs to $K$. The fourth corner is high, and its abscissa is at most $5/8+e<3/4$; it must be on the unchanged top or left. On the top, $|y|\le2e$, so its abscissa is at most $1/2+3e<5/8$. On the left its height exceeds $3/4$. Thus it belongs to $L$, and every square hitting $I,J$ contributes once to the marked degree.

## 5. Why universal nonzero winding fails, even on the circle

On the unit circle the two coordinates of $F$ agree identically. In fact, for any center $C$, translating $Q,B$ together by $R(B-Q)$ gives
$$
|P-C|^2-|D-C|^2=|Q-C|^2-|B-C|^2.
$$
The cross term vanishes by perpendicularity. Equal input radii therefore imply equal output radii. Signed distance to the circle depends only on radius. Every zero-free perimeter image lies on the diagonal in the residual plane and has winding zero. Yet every $Q$ on the circle gives a square by taking $B=RQ$. Conversely the nondegenerate square completions on the circle obey precisely this relation, since opposite corners must have the origin as midpoint. The resulting one-parameter family passes through the boundary of any proper parameter rectangle that meets the family in its interior. Thus the condition “every Jordan curve has two windows with defined, nonzero winding” is false.

**A limiting alternative.** Fix two disjoint closed parameter windows with image distance $d>0$ on a target Jordan curve. Suppose arbitrarily close registered approximants have nonzero whole-curve winding on these windows. Either the target has an exact oriented square on their perimeter, or its perimeter is zero-free and its winding equals that of all sufficiently close approximants. In the latter case the target has an interior square. Both alternatives give side at least $d$. Indeed, if there is no perimeter square, compactness gives $m>0$ and the quantitative stability lemma applies with the target as reference. This leaves open how to force one alternative at positive separation for arbitrary curves; the circle rules out demanding the second alternative alone.

## Scope and replay

The whole-curve test applies to any specified pair of disjoint parameter windows with nonzero boundary winding. The outstanding problem is to force either a boundary square or nonzero winding on separated windows of suitably chosen approximations. The exact finite census detects isolated transverse squares and the folded examples give sufficient geometric conditions. Neither supplies universal windows.

The rotating-side construction also features in [Tao's integration approach](https://doi.org/10.1017/fms.2017.23); [Greene and Lobb](https://doi.org/10.4171/cmh/619) treat curves formed from two Lipschitz graphs. Here the certificate is a winding on two designated parameter intervals, and permits folds on all marked arcs. Run the Python verifier and its unittest file; [recorded results](../results/v312/rough-arc-winding.json) contain the rational output.
