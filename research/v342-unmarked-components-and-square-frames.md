# v342. Counting a square once when its four labels belong to one family

Benjamin John Schulz · November 2026

On a square boundary of side $S$, the connected labeled family from [v341](v341-ledge-and-component-certificates.md) has signed degree zero. Counting each geometric square once gives **local parity one**. Its robust largest-square side is exactly $S/\sqrt2$. There is also a sharp finite-error guarantee. Among Jordan boundaries admitting a registration with coordinate sup error at most $\varepsilon$ from the reference, the infimum of the largest inscribed square side is
$$
\frac{S-2\varepsilon}{\sqrt2},
\qquad 0<\varepsilon\le S/4.
$$
At zero error the reference itself has largest side $S$. The finite-error lower bound uses Matschke's existing square-annulus theorem. An explicit family of oblique rhombi gives the upper bound.

## Relative parity on one component of geometric squares

Use the [v341 tracer](../code/polygon_component_v341.py) on the four marked boundary parameters to find an isolated positive-size zero component $C$. Write $g$ for cyclic relabeling and $H=\{g^j:g^jC=C\}$. The action on individual positive-size squares is free even when the whole component is invariant. Take the union of the translates of $C$, and choose a small invariant open neighborhood with compact closure, separated from parameter collisions and from every other zero component. Its quotient by the four relabelings counts each geometric counterclockwise square once. This includes squares of any curve-order sector represented in the selected components.

The residual $F=(p_2-p_1-J(p_1-p_0),p_3-p_0-J(p_1-p_0))$ transforms by the order-four target matrix of [v262, section 3](v262-cyclic-quotient-relative-parity.md). Thus it defines a section of a rank-four real vector bundle on that quotient. Define the **relative parity** by the number modulo two of zeros of a transverse section sufficiently close to it. Independence follows by joining two perturbations through sections nonzero on the boundary, smoothing relative to the ends and making the homotopy transverse. The compact one-dimensional zero manifold has an even number of endpoints. The same reasoning applies to the continuous geometric section of a sufficiently close registered curve. If the reference parity is odd, every strict cutoff below the component's minimum side survives. This is local to the isolating neighborhood. The [v219 vanishing absolute obstruction](v219-equivariant-quotient-collision-link.md) is compatible with this construction.

**Covering rule.** The parity of the labeled degree on $C$ equals $|H|$ times the geometric parity modulo two. Choose the invariant neighborhood small enough that translates of the neighborhood of $C$ are identical or disjoint according to $H$. Lift a regular quotient section. Each geometric zero lifts to exactly $|H|$ zeros in the neighborhood of $C$, since the group is free there. When the labeled degree is odd, $H$ is trivial and the geometric parity is odd. When $|H|$ is even, one needs a quotient computation. In particular, dividing v341's signed count by four would be invalid: its constant target shift need not respect relabeling.

All **476** recorded v341 certificates have trivial component stabilizer, verified directly from their edge charts. Their existing odd labeled degrees therefore transfer to geometric parity without changing any size floor. In v262 the selected boxes had disjoint cyclic images. The ordinary square here is the complementary case: one connected family contains all four choices of first vertex.

### Exact replay on a connected family

For a rational polygon with vertices $P_i$, deform the actual vertices to $P_i+tV_i$, $t\downarrow0$. The resulting geometric residual respects the cyclic action. There are finitely many charts, so any nearby root limiting to the selected component union belongs to its saturated chart list. Conversely any feasible root in a listed closed chart limits to that chart's reference slice, which lies entirely in the selected component union by convexity. A sufficiently small deformation remains embedded for the fixtures used here.

Test one cyclic representative of each edge tuple in the saturation. The chart equations are $(A+tA')u=b+tb'$ because the residual is linear in the evaluated points. The [replay](../code/quotient_component_v342.py) uses the [v337 rational-function box solver](../code/terrain_component_degree_v337.py) over the one-sided germ. It checks singular systems for feasibility as well. It accepts a surviving chart only when its root is unique, interior and transverse, and its edge tuple has four distinct cyclic rotations. Otherwise it declines the deformation. With these checks, different surviving representatives cannot encode the same geometric square. The count modulo two is the quotient parity. A collection of other deformation directions need not be searched for the conclusion on an accepted replay.

| Reference component | Stabilizer order | Charts in one labeled component | Chart orbits in its saturation | Surviving orbit count for the chosen deformation |
|---|---:|---:|---:|---:|
| Whole large family on the $2\times2$ boundary | 4 | 60 | 15 | **1** |
| Translating unit squares in the tower on (2,7) | 1 | 5 | 5 | **0** |

For the first row the map is $(x,y)\mapsto(x+ty,y+tx)$; negating $t$ gives the same count. Only the four-different-edges chart 0123 survives. The reference floor of squared side is 2. The tower row uses the exact vertex velocities in the tests. It illustrates that quotient parity can itself vanish. The branch geometry and registered neighborhood in the two rows are different.

## Why the square bound is sharp

For an oblique rhombus let $A,B$ denote its unequal half-diagonals. Its unique inscribed square has side $2AB/(A+B)$ and four strictly interior edge contacts.

To see uniqueness, suppose two square vertices lie on one closed side of the rhombus. They must be adjacent: the other two vertices of a square straddle the line of a diagonal. Take the rhombus side as a horizontal base of length $L$, with the opposite side shifted by $d$ and at height $h$, so $L^2=d^2+h^2$ and $d\ne0$. A square on this base with side $s<h$ must have both upper vertices on the two lateral sides, hence $s=L>h$. If $s=h$, fitting the upper and lower pairs requires $L-|d|\ge h$, whereas $|d|+h>L$. Thus each rhombus side contains exactly one contact, strictly in its interior. Convexity gives their cyclic order. The square center, as the midpoint of opposite contacts, belongs to both rhombus midlines and is its center.

Use the diagonals as axes. Two consecutive centered vertices $(x,y),(-y,x)$ on $|x|/A+|y|/B=1$ satisfy $|x|/A+|y|/B=|y|/A+|x|/B=1$. Thus $|x|=|y|=AB/(A+B)$. This also proves transversality: opposite-pair equations fix both center coordinates by independent normals; the remaining equations in the fixed open sign quadrant have nonzero determinant $1/A^2-1/B^2$ up to sign. Such rhombi converging to a square of side $S$ have largest square side tending to $S/\sqrt2$. For the side-2 replay the unique squared side is $2(1-t^2)^2$, independently checked by exhaustive rational polygon censuses at three finite $t$.

All inscribed squares on the square boundary span its full width. Indeed, if two vertices use one side, they form a square edge, and the other pair either lies on the opposite side or spans the two lateral sides, forcing side $S$. Otherwise each boundary side has a contact. A Euclidean square of side $s$ has axis-aligned bounding width at most $\sqrt2s$. Thus this family's minimum is $S/\sqrt2$, attained by the four side midpoints. Alternatively its full interval up to cyclic relabeling has consecutive edge vector $(u,S-u)$ with $0\le u\le S$. The odd geometric parity above gives the robust lower bound, and the rhombi give equality.

## A quantitative guarantee for redrawing a square

Matschke's Theorem 2.1 states that a continuous closed loop of winding one in $1\le\|x\|_\infty\le3$ inscribes a square of side at least $\sqrt2$. Scale and translate it: an inner square of half-side $r>0$ and a concentric parallel outer square of half-side $R\le3r$ give a side floor $\sqrt2r$ for any loop in that annulus winding once. Reversing the traversal handles winding minus one. Theorem 1.3 of the same paper is the distinct circular-annulus version. The square-annulus theorem and the parity mechanism in its Lemma 2.3 are prior results. [Matschke, *On the square peg problem and some relatives*, §2, Theorem 2.1 and Lemma 2.3](https://www.math.u-bordeaux.fr/~bmatschke/SquarePegAndRelatives2_05.pdf)

If a curve is registered within $\varepsilon<S/2$ of the side-$S$ square boundary, its sup radius about the center stays between $S/2-\varepsilon$ and $S/2+\varepsilon$. Straight interpolation between corresponding points avoids the center, so its winding is preserved. The ratio is at most three exactly when $\varepsilon\le S/4$. This proves the stated finite-error lower bound. It covers arbitrarily fine oscillations.

For the matching upper infimum fix positive $\varepsilon$ in this interval. First shrink the centered reference square to side $S-2\varepsilon+\eta$ with $0<\eta<2\varepsilon$. This costs $\varepsilon-\eta/2$ in registered max norm. Apply a nonzero centered rhombus deformation of displacement less than $\eta/2$. Choose its strain tending to zero with $\eta$. As $\eta\downarrow0$, the unique square on those curves has side tending to $(S-2\varepsilon)/\sqrt2$. No attainment assertion is needed.

## Finding the square frame without marking contacts

For fixed axes and a polygon made from filled consecutive columns, the [frame search](../code/terrain_square_annulus_v342.py) finds the largest admissible $r$ over **all centers**. For a block of columns $[a,b)$ let $B$ be the largest bottom height and $T$ the smallest top. A contained axis-aligned open square over these columns fits in $[a,b]\times[B,T]$. Every contained square is covered by the block of columns its open horizontal interval meets. For a candidate $r$, admissible centers for the inner square form the box
$$
[a+r,b-r]\times[B+r,T-r].
$$
To enclose the entire polygon in outer half-side $3r$, intersect with
$$
[x_{\max}-3r,x_{\min}+3r]\times[y_{\max}-3r,y_{\min}+3r].
$$
Within any one block it suffices to try $r=\min(b-a,T-B)/2$. If a smaller $r$ worked, the inner box's faces move inward at speed one up to this maximum while the outer constraints expand at speed three. Equivalently the two cross-inequalities in each axis get weaker and the individual intervals stay nonempty. Clipping the outer bounding-box midpoint to the inner center box minimizes the required outer radius. Among equal inner radii, choose the least outer radius. The center is interior to the Jordan polygon, so the winding condition holds. The procedure searches a fixed orientation; rotations can change its answer.

| v336 sample | Frame at target with $R<3r$ | Frame at target with $R=3r$ | No admissible fixed-axis frame |
|---|---:|---:|---:|
| 282 single terrains | 67 | 65 | 150 |
| 194 glued terrains | 7 | 6 | 181 |

Here the target is $r=i/2$. No sample row had a best admissible frame strictly below the target. The [ledger](../results/v342/square-frames.csv) records the chosen center and radii. A strict frame gives an explicit registered neighborhood: for coordinate error $\varepsilon\le(3r-R)/4$, the winding persists and the perturbed side floor is $\sqrt2(r-\varepsilon)$. At equality the frame alone proves existence on the reference. For all sample rows, v341 separately gives persistence, including those where this frame search fails.

The outstanding selection problem is to find a geometrically large region of odd unmarked parity along approximations to an arbitrary loop. The ordinary square resolves a failure caused purely by relabeling. The tower control and the failed frame searches show why dropping labels alone does not complete the general proof.

## Replay

```sh
python code/quotient_component_v342.py
python code/terrain_square_annulus_v342.py
python -m unittest discover -s tests -p 'test_quotient_component_v342.py' -v
```
