# v339–v340. Selecting large squares on end-wall terrains

Benjamin John Schulz · November 2026

A grid terrain has a horizontal base and positive unit-width columns of heights $H_0,\ldots,H_{w-1}$. Let $i$ denote its largest contained open axis-aligned square side. A certificate with lower bound $s$ means that for each $s'<s$, every sufficiently close registered Jordan boundary has an inscribed square of side greater than $s'$. Registration means uniform convergence of boundary parametrizations. The perturbed boundary can have arbitrarily many additional vertices or be nonpolygonal.

**Result.** On the 252 end-wall rows of [v338](v338-robust-terrain-proof.md), 131 have a strict four-contact certificate at or above $i/\sqrt2$. The component and corner lemmas below certify another 97. The resulting 228 certificates use exact boundary geometry. The count is for this fixed sample, which has 237 distinct height sequences. The general end-wall conjecture remains open.

| First certificate applied | Rows | Smallest newly relevant example |
|---|---:|---|
| v339 strict crossing | 131 | Includes all three v338 global-parity exceptions |
| v340 sliding family between tall end walls | 42 | (3,3), a rectangle of width 2 and height 3 |
| v340 plateau followed by a rise | 21 | (2,5) and its reflection are members |
| v340 one wall corner | 34 | (2,5,2) |
| Unresolved by these tests | 24 | (1,1,5,2) |

Smallest means least area (sum of column heights), then width, lexicographic heights and seed. The [252-row certificate file](../results/v340/end-wall.csv) retains repeated seeds' shapes, all misses, the squared side of the best strict crossing, and the support coordinates used for each certificate. The [selector](../code/terrain_selector_v339.py) regenerates seeds 0 through 299 with v338's generator, interior-side filter and excursion classifier. It uses no randomly perturbed coordinates. The first nonrectangular strict miss is (2,5). Both (1,1,5,2) and (1,4,2,2) remain at area 9.

## v339. Enumerating strict contacts

For integer $a<k<b$ put $u=k-a$, $v=b-k$. The candidate vertices in geometric order are
$$
(a+b-k,0),\quad (b,v),\quad (k,u+v),\quad (a,u).
$$
They form a square of squared side $u^2+v^2$. An interior or exterior wall at $x$ spans the heights of the columns before and after $x$, using exterior height zero. Accept precisely when the three indicated ordinates lie strictly between those heights at their respective walls. The base contact is automatically strictly interior. Every strict crossing of the [v337 four-contact form](v337-terrain-component-degree.md) has these coordinates. The cutoff test is $2(u^2+v^2)\ge i^2$. This enumerates a particular contact pattern. Other patterns can also have nonzero degree.

| Heights | $(a,k,b)$ | Squared side |
|---|---|---:|
| (8,4,3,1,1,6) | (0,1,6) | 26 |
| (5,1,7,4,6,4) | (0,3,6) | 18 |
| (4,3,1,6,3,8) | (2,4,6) | 8 |

For all degree arguments below we use one local labeling in geometric counterclockwise order. For four labeled boundary parameters the continuous square residual is $(p_2-p_1-J(p_1-p_0),p_3-p_0-J(p_1-p_0))$, where $J$ rotates by a positive quarter turn. An isolated compact labeled zero set has a neighborhood with no boundary zero. A sufficiently small regular deformation computes its local degree. A registered uniform perturbation changes the residual uniformly and preserves that degree. By shrinking the neighborhood, every strict cutoff below the minimum reference side of the component persists. These are local statements, so other square components may coexist.

## v340. Tall end walls

**Theorem.** If both exterior walls are strictly taller than the terrain width $w$, its robust largest-square side is exactly $w$. The intervening roof may dip below $w$.

**Proof.** Normalize the width to 1 and let $m>1$ be the shorter end height. The axis-aligned unit square slides vertically, with bottom height $s$ ranging from 0 to $m-1$. Label its bottom-left and bottom-right vertices first. Their first edge points to the right throughout the interval, giving a continuous lift disjoint from its other three cyclic labelings.

This interval is isolated. Near any member with $s>0$ both lower contacts are strictly inside the end walls. Near $s=0$ both upper contacts are strictly inside them. If either pair is on its respective walls, the horizontal displacement between the pair is 1. The other pair has the same horizontal displacement and must stay inside the unit strip. Consequently its endpoints also lie on the two walls. The edge is horizontal, the side is 1, and the zero belongs to the displayed interval.

Apply the taper $(x,y)\mapsto(ty+(1-2ty)x,y)$, for $t>0$ small enough to preserve the order of the walls at every height. This is an embedding of the terrain strip up to its maximum height. Every vertex stays between the new walls $x=ty$ and $x=1-ty$. Write the lower edge vector of a square near the family as $(d,e)$ with $d>0$. The upper vertices are obtained by adding $(-e,d)$. If both lower vertices are on their walls, the left upper vertex would require $-e\ge td$ and the right one $-e\le-td$. Thus no such square survives. Compactness leaves only the bottom end as a possible limit.

Near that end, both upper vertices stay on their wall interiors. There are three remaining possibilities for the lower pair. When both are on the base they give an axis-aligned square of side $1/(1+2t)$. When the lower left is on its wall at height $y$ and the lower right is on the base, $e=-y$. The two upper-wall equations give $y=td$ and $(1+2t+t^2)d=1$. These equations place all four points strictly inside their specified segments for small $t$. Reflection gives the third square. Each of these affine contact systems has a unique solution by nonsingular linear elimination and therefore contributes one modulo two. The three squares give odd local degree. The limiting lower bound is 1. Conversely the horizontal extent of any square is at least its side, so the reference strip itself gives the matching upper bound. Rescaling proves the theorem.

## v340. A plateau at one end

**Theorem.** Suppose the left end height equals $w$, the initial horizontal roof extends a distance $A$ before a step *up*, and the right end height exceeds $w$. Then the robust size is at least $\sqrt{d^2+(w-d)^2}$ with $d=\min(A,w/2)$. The reflected statement also holds. This includes arbitrary changes to the roof away from the contact traces.

Normalize $w=1$ and put $a=A/w\in(0,1)$. The component from v337 extends to the whole interval
$$
(1-u,0),\ (1,1-u),\ (u,1),\ (0,u),\qquad 0\le u\le a .
$$
The minimum is attained at $u=\min(a,1/2)$. The orientations of its first edges span strictly less than a quarter turn, so one can isolate a single labeling. To check the incident edge charts, use the model polygon $(0,0),(1,0),(1,H),(a,H),(a,1),(0,1)$ with any $H>1$. In the original terrain the contact traces have identical local geometry up to shortening unused portions of their incident edges. At the right wall only heights up to 1 and a neighborhood above 1 are needed. At the roof step only a neighborhood immediately above height 1 is needed. Compactness permits restriction to these common subarcs.

In the model, number edges from zero. At $u=0$ there are two choices of incident edge at each of three corners. At $u=a$ there are two choices at the roof endpoint. Up to cyclic relabeling the resulting chart list is the nine rows below. At zero shear chart 0145 gives the interval, 1155 can only fit a unit axis-aligned square on the left wall of height 1, and the seven remaining charts have nonsingular linear equations and give just the indicated endpoint when feasible. This isolates the family for every $0<a<1$, $H>1$.

Use the shear $(x,y)\mapsto(x+ty,y)$ for $t\downarrow0$. The table gives one decisive equation or infeasibility per chart. Each parameter runs from 0 to 1 in the edge orientation listed above. It follows by substituting the four edge parametrizations in the square residual.

| Chart | Feasibility for small positive $t$ |
|---|---|
| 0014, 0015 | The second base parameter is $1+t$. |
| 0114 | The first parameter is $-t^2$. |
| 0115 | The second parameter is $-t/[H(1+t^2)]$. |
| 0155 | The third parameter is $t(t-1)/(1+t^2)$. |
| 1145 | The first parameter is $t(t-1)/[H(1+t^2)]$. |
| 1155 | If the fourth parameter is $z$, nonnegativity of the third requires $(1+t^2)z\ge1$. The first is $[1+t^2-t-(1+t^2)z]/[H(1+t^2)]<0$. |
| 0145 | The top point's pre-shear abscissa is $(1-t)/2$. The solution lies in the open plateau exactly when $a>(1-t)/2$. The other three parameters are $(1+t)/2,1/(2H),1/2$, in their contact order with the top one omitted. |
| 0135 | The top point's height is $[1+(1-2a)t]/(1+t^2)$. It lies above the step foot exactly when $a<(1-t)/2$. The other contacts limit to the interiors of the base and walls. |

The two feasible alternatives solve nonsingular affine systems (their determinant magnitudes are $2aHt$ and $H(H-1)(1+t^2)$, respectively, in the normalized model). Exactly one is interior for all sufficiently small positive $t$; at $a=1/2$ it is 0145. It approaches the minimum-size member of the interval. This gives odd degree and proves the stated lower bound. Other components may affect the largest square on the full terrain.

## v340. One wall contact at a corner

**Lemma.** Keep the base, roof step and one side wall as strict interior contacts in the v339 construction. If the other side contact is a corner, a transverse crossing of the following scalar function along its two incident boundary rays gives local degree of absolute value one.

| Free contact $(X,Y)$ | Scalar crossing |
|---|---|
| Left | $X+Y-k$ changes sign |
| Right | $X-Y-k$ changes sign |

For example, use a base point $(s,0)$ and an opposite roof point $(k,h)$. If the right wall stays fixed at $x=b$, its equation gives $s=2b-k-h$ and the free left point becomes $(b-h,h+k-b)$. Its trace is the line $X+Y=k$. Its horizontal coordinate fixes $h$ uniquely, and the other wall parameter is fixed by its ordinate. Elimination leaves precisely the scalar crossing equation on the marked corner arc, multiplied by nonsingular linear factors in the full residual. A sign change gives one-dimensional local degree $\pm1$. Fixing the left wall instead gives the free right point $(a+h,a+h-k)$ and the second test. A corner touching the line from one side has no such certificate.

On (2,5,2), the square $(1,0),(3,1),(2,3),(0,2)$ has side $\sqrt5$. Along the left corner the quantity $X+Y-2$ goes from positive to negative. Its reflection gives a second local square. The implementation checks the signs on the incident horizontal and vertical rays. The independent test rebuilds the actual polygon edges for every reported crossing on every end-wall row, computes the four-coordinate Jacobian in each incident chart, and verifies both determinants have the same nonzero sign at accepted corners.

## Remaining work and replay

Of the 24 remaining rows, the first two by area have heights (1,1,5,2) and (1,4,2,2). Their rows have no strict base/three-wall witness. A useful next test is to allow the marked horizontal contact to move to a roof ledge or rotate the whole contact pattern. For a cutoff reached only along a family, enumerate the incident charts of that family with the v337 census. A miss in the CSV is an absence of these certificates, not a bound on the terrain's inscribed squares.

```sh
python code/terrain_selector_v339.py
python -m unittest discover -s tests -p 'test_terrain_selector_v339_v340.py' -v
```
