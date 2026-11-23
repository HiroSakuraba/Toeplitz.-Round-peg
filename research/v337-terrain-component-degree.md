# v337. A sharp robust square bound for two terrains in the archive

Benjamin John Schulz · November 2026

Let $M(C)$ be the largest inscribed square side on a closed Jordan curve (zero if there is none). Fix a boundary traversal of a polygon $P$. Its **robust size** is the limit, as $\delta\downarrow0$, of the infimum of $M(C)$ over Jordan boundaries admitting a registration with $\|C-P\|_\infty<\delta$. Registration compares corresponding boundary points and permits arbitrarily many additional vertices or a nonpolygonal curve. The reference size $i(P)$ is the side of the largest contained open axis-aligned square.

**Theorem.** The robust size equals $i(P)/\sqrt2$ for both of the following polygons.

| Reference polygon | $i(P)$ | $M(P)$ | Robust size | Large component's mod-2 degree |
|---|---:|---:|---:|---:|
| Terrain with column heights (2,7) | 2 | 2 | $\sqrt2$ | 1 |
| [v334 glued witness A](v334-two-terrain-gluing.md) | 7 | 5 | $7/\sqrt2$ | 1 |

Every strict cutoff below the indicated robust size therefore holds throughout some registered neighborhood. The neighborhood may depend on the cutoff. A shear gives the matching upper bound. The theorem settles two fixtures of [Conjecture 336.1](v336-robust-terrain-squares.md). The uniform assertion over all terrains and their unions remains open.

## 1. The two components to isolate

Number edges from zero in the cyclic orders printed below. In both examples let $w$ denote the span of the two walls used by the large family. Four points of the boundary in counterclockwise geometric order are
$$
B=(w-u,0),\quad R=(w,w-u),\quad T=(u,w),\quad L=(0,u).
$$
They form a square of squared side $u^2+(w-u)^2=w^2/2+2(u-w/2)^2$.

**Single terrain.** The boundary is (0,0), (2,0), (2,7), (1,7), (1,2), (0,2). Set $w=2$ and $0\le u\le1$. The corner square and diamond are the two ends of one closed interval of square orbits. There is precisely one other component of positive-size squares, the vertically translating unit square on the two tower walls, with bottom height from 2 to 6. A squared-side cutoff of 3/2 separates them.

The diamond alone can disappear. Move the right wall left to $x=2-\varepsilon$. A square converging to the diamond would still use the base and both walls in their interiors. Perpendicular equal diagonals would put its top point at height $2-\varepsilon$, while all nearby roof points have height at least 2. The corner end still admits a square of side $2-\varepsilon$.

**Glued witness.** The boundary is (0,0), (6,0), (6,−1), (8,−1), (8,3), (7,3), (7,7), (3,7), (3,8), (0,8). Here $w=7$ and $3\le u\le4$. Both ends have side 5, and the middle has side $7/\sqrt2$. All charts outside this interval and its endpoint contacts have squared sides strictly below 16. We use cutoff 24. The original upper terrain's side-7 square lost a base contact in the open seam. The replacement uses the exposed base between abscissas 3 and 4.

## 2. Exhaustion, including chart boundaries

The finite algebra below is part of the proof. [Exact replay](../code/terrain_component_degree_v337.py) uses the rational edge equations from [v223](../code/polygon_orbit_census_v223.py). For every ordered choice of four edges solve the four affine square equations within the closed unit parameter box. We fix counterclockwise *square* labels, retaining all curve orders, and take one cyclic representative per edge tuple. Three vertices on one straight edge are impossible at positive size. Rank-deficient slices are bounded convex polytopes, fully determined by their enumerated vertices. Squared side is convex on each slice, so a chart with all vertex sizes below cutoff cannot reach it. The explicit interval formula above gives the **minimum** on the surviving segment.

| Polygon | Charts meeting cutoff, in cyclic representatives | Chart containing the nonconstant segment |
|---|---|---|
| (2,7) | 0014, 0015, 0114, 0115, 0135, 0145, 0155, 1145, 1155 | 0145 |
| Glued A | 0469, 0569, 0579 | 0569 |

Each digit is an edge number. Every other chart in these rows contains just an endpoint at the reference polygon. For (2,7), there are 14 positive-size charts in total. The other nonconstant segment is 1133. The test compares the rational endpoint orbits and checks the size separation. Gluing endpoint charts is essential to obtaining the component's degree.

## 3. A single labeling and an isolating neighborhood

Work on the four-dimensional torus of four marked boundary parameters, each increasing along the given traversal. Write $p_j=P(s_j)$, $e=p_1-p_0$, and $J(x,y)=(-y,x)$. The continuous square residual is
$$
F_P(s)=(p_2-p_1-Je,\ p_3-p_0-Je).
$$
Take the open set $U$ given by $|e|^2>c$, $2e_x+e_y>0$, and $2e_y-e_x>0$, with the relevant squared cutoff $c$. The last two inequalities describe the open right angle between the rays (2,1) and (−1,2). The whole large family has its edge from $B$ to $R$ strictly in this sector. The other three counterclockwise labelings rotate that edge by successive quarter-turns. Thus each orbit contributes exactly one lift. The census and strict inequalities show that $F_P$ is zero-free on $\partial U$. This construction uses the **reference** polygon to define $U$ throughout the deformation.

Compactness gives a positive residual margin on $\partial U$. Under any registered uniform perturbation of the curve, the residual changes uniformly, since it is linear in the four evaluated points. The straight homotopy between sufficiently close residuals keeps the boundary zero-free. This defines an ordinary local degree, including for perturbations with no tangent or edge decomposition.

## 4. Resolve the family with a shear

Use the registered linear shear $(x,y)\mapsto(x+\sigma t y,y)$ for $t>0$, with $\sigma=+1$ for (2,7) and $-1$ for glued A. If arbitrarily small shears had an additional square above the cutoff, a subsequence on one fixed edge tuple would limit to a square at or above cutoff in that closed chart. Thus only the table in section 2 needs checking. This compactness argument also covers changes of edge label near a polygon vertex.

| Polygon | Unique feasible chart near $0^+$ | Its four edge parameters in the displayed edge order | Jacobian determinant | Squared side |
|---|---|---|---:|---:|
| (2,7) | 0145 | $((1+t)/2,1/7,t,1/2)$ | $-56t$ | 2 |
| Glued A | 0569 | $(7(1-t)/12,1/8,7(1-t)/8,9/16)$ | $1536t$ | $49/2$ |

Here feasibility is checked for **all sufficiently small positive** $t$. The replay works over rational functions of $t$. It forms a basis of independent equations and adds active unit-box faces in every way that completes a square system. Each resulting vertex is checked against every original equation and box inequality. The sign at $0^+$ is the sign of the ratio of the lowest nonzero polynomial coefficients. This also rejects the singular sheared chart 1155 on (2,7). A nonempty closed slice of the box has a vertex. On the two surviving charts all four parameters are strictly interior and the displayed determinants are nonzero.

The unique sheared zeros converge to the marked diamond at $u=w/2$, in $U$. With the coordinate and residual orientations above the local degrees are respectively −1 and +1. In particular both are odd. Homotopy invariance transfers the degree to the reference polygon and to sufficiently close registered curves.

For any such sequence of curves approaching the reference, take zeros in $U$. Every subsequential parameter limit is a reference zero in $\overline U$, hence belongs to the large family. Its squared side is at least $w^2/2$. This proves the robust lower bound. For the upper bound, the shear has exactly one square above the chosen cutoff for every sufficiently small positive $t$, and its side is exactly $w/\sqrt2$. Hence the limiting infimum is sharp. No sampling of perturbations is used.

## 5. A strict-crossing lemma for other boundaries

Suppose a boundary has a horizontal base contact at $(w-a,0)$, two vertical wall contacts at $(0,a)$ and $(w,w-a)$, and a vertical roof step crossing $(a,w)$, where $0<a<w$. Assume each contact is in the relative interior of its segment. The four marked local arcs force a square of limiting side at least $w/\sqrt2$ under registered perturbations.

Indeed put the base and roof points at $(s,0)$ and $(a,h)$. The two completed vertices have abscissas $(s+a-h)/2$ and $(s+a+h)/2$. Equating them to 0 and $w$ gives a nonsingular system in $(s,h)$, with solution $(w-a,w)$. The wall ordinates then fix their parameters uniquely and transversely. Thus the four-contact map has degree of absolute value one locally. Its reference side is $\sqrt{a^2+(w-a)^2}$. For a terrain of width $w$, $i(P)\le w$. For example the terrain (2,5,5) gives a regular square of side $\sqrt5$. The same proof applies after gluing whenever the four contacts remain on exposed open segments.

More generally, the component certificate survives a change elsewhere on the boundary if the boundary agrees on marked open subarcs covering all four contact traces of the component, with their registrations preserved. Restrict to a smaller isolating neighborhood of its compact labeled lift contained in the product of these subarcs, using a union of parameter boxes if needed. Excision preserves the degree, and the residual on this neighborhood is identical. For (2,7) these neighborhoods must include both sides of each visited corner. Establishing that an arbitrary gluing retains such a component, or supplies another one when a trace is swallowed, is the remaining geometric selection problem.

## Replay

Install SymPy 1.14, then run from the repository root.

```sh
python code/terrain_component_degree_v337.py
python -m unittest discover -s tests -p 'test_terrain_component_degree_v337.py' -v
```
