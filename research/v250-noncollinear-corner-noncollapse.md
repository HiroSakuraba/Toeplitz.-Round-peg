# v250: A genuine-corner square handoff and a branch-independent noncollapse bound

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact-rational theorems, an explicitly changing strictly convex polygon family with a genuinely non-collinear corner, full real-parameter square continuation, and a separate no-collapse bound applying to **every** inscribed square of a fixed polygon. The result is **restricted** and does not establish a universal Big Peg constant or solve the Toeplitz square-peg conjecture.

## 1. The question

Our v248 local Brouwer result protects a square while its corners lie in the interiors of four fixed edges. v249 verifies a complete two-chart continuation across a *collinear subdivision vertex*. An actual polygon corner is harder because its two incident supporting lines are different, so the square solutions on either side of the transition need not describe the same algebraic function.

v250 separates two mechanisms:

1. A particular positive-size square branch can cross a **genuine corner** and continue in a different, independently certified edge chart.
2. Irrespective of branches, **all** squares on a fixed polygon are bounded away from zero by a positive geometric feature scale, provided the polygon has at least four sides. This is not a uniform theorem across arbitrarily intricate polygons, because the feature scale can shrink while the inradius stays large.

The second mechanism makes precise what it would mean for extreme corners to prevent square collapse, while demonstrating why corner angles alone do not yield a universal result.

## 2. Exact genuine-corner construction

Start with the rational octagon from v248:

\[
A_0=(39,21)/20,\ A_1=(18,40)/20,\ A_2=(-17,39)/20,\
A_3=(-38,20)/20,
\]
\[
A_4=(-43,-19)/20,\ A_5=(-20,-38)/20,\
A_6=(21,-38)/20,\ A_7=(39,-23)/20.
\]

For the four original supporting edges \((0,2,4,6)\), the exact Cramer square root \(z=(p_x,p_y,u,v)\) is

\[
z_*=\left(
\frac{21561}{15136},\frac{808131}{529760},
-\frac{44719}{15136},-\frac{20205}{105952}
\right).
\]

The vector

\[
K=\left(\frac{3241}{8989},-\frac{5455}{8989},
\frac{665}{8989},1\right)
\]

is in the **exact nullspace of the three unchanged incidence equations for edges (2,4,6)**. This assertion is verified directly in code, including all three resulting square-corner segment conditions.

Define the desired corner at the event by displacing the square along that kernel:

\[
V_*=(z_{*,x},z_{*,y})-\tfrac12(K_x,K_y).
\]

Let \(E=A_1-A_0\), and move this polygon vertex as

\[
\boxed{V(\lambda)=V_*+\frac{\lambda-\tfrac12}{4}E,\qquad 0\le\lambda\le1.}
\]

Use the cyclic nine-vertex polygon

\[
P(\lambda)=(A_0,V(\lambda),A_1,A_2,A_3,A_4,A_5,A_6,A_7).
\]

Only the new vertex moves, but **the geometric boundary changes**, not merely its edge parametrization. The two incident edges at \(V\) are never collinear: their exact signed cross product is

\[
\boxed{(V-A_0)\times(A_1-V)=\frac{6622}{44945}>0}.
\]

In contrast, v249's vertex had an identically zero turn because it lay *on* the original edge.

### The all-time Jordan certificate

For each pair of *nonadjacent* edges, take a global axis-aligned bounding rectangle containing the edge throughout all real times. Because every vertex is affine in \(\lambda\), it suffices to consider its coordinates at \(\lambda=0,1\). The **minimum coordinate separation** of these global edge rectangles is the exact rational

\[
\boxed{g=\frac{1810915}{8503594}>0}.
\]

Therefore no nonadjacent edges can intersect at any time. Additionally, exact Bernstein positivity checks show that all nine successive turns remain *strictly positive* and all edge lengths stay nonzero for every \(\lambda\). Thus no adjacent backtracking occurs and \(P(\lambda)\) is a strictly convex simple closed polygon for every real time.

Each vertex at \(\lambda=0,1\) lies strictly in the radius-three disk centered at the origin. Because the disk is convex and every vertex moves affinely, the entire polygon and its interior lie in that disk for every \(\lambda\). Accordingly, its **true inradius is at most three**.

## 3. The genuine two-chart square handoff

Reuse the independent exact Cramer/Bernstein finite-atlas verifier from v249. It solves the four supporting-line equations for the square parameters as rational functions of \(\lambda\), then proves the rational **closed-segment** inequalities for every parameter time, with strict positive square-side lower bound.

The two required ordered square edge charts are:

| Time | Square-supporting edges | Location of first square corner |
|---|---|---|
| \(0\le\lambda\le1/2\) | \((1,3,5,7)\) | Segment from kink \(V\) to \(A_1\) |
| \(1/2\le\lambda\le1\) | \((0,3,5,7)\) | Segment from \(A_0\) to kink \(V\) |

At the exact handoff time \(\lambda=1/2\), the two Cramer roots **agree exactly** and their first corner equals \(V_*\). Its segment parameter equals zero in the first chart and one in the second, so it is a valid boundary corner in **both** closed-edge descriptions.

Crucially, **the two square branches are algebraically distinct away from the handoff**. Cross-multiplying the Cramer numerators and determinants shows they are not the same rational function. This is a genuine change of supporting lines, not v249's one branch rewritten with a collinear edge subdivision.

For both chart intervals, exact rational Bernstein inequalities certify:

\[
u(\lambda)^2+v(\lambda)^2 > 9
\quad\text{for **every real** }\lambda\in[0,1].
\]

The square has side **strictly greater than three** throughout. The Jordan polygon stays within a radius-three disk, hence its true inradius is at most three. Therefore

\[
\boxed{J(P(\lambda))=
\frac{s_{\max}(P(\lambda))}{r_{\rm in}(P(\lambda))}
>1\qquad(\forall\,\lambda\in[0,1]).}
\]

This is a new restricted, quantitative all-time Big Peg result across a genuine polygon corner. It is **not** a result for arbitrary polygon families.

## 4. A no-collapse bound for **every** square on a fixed polygon

Let \(P\) be any simple Jordan polygon with at least four edges, and define its exact geometric feature scale

\[
d(P)=\min_{\substack{e_i,e_j\\ \text{nonadjacent edges}}}
\mathrm{dist}(e_i,e_j)>0.
\]

**Elementary lemma:** Every nondegenerate square whose four distinct corners lie on \(\partial P\), regardless of which edge charts or square branches describe it, has side

\[
\boxed{s\ge d(P)/\sqrt2}.
\]

**Proof:** Assign each square corner to an incident polygon edge (if it is itself a polygon vertex, choose either incident edge). A nondegenerate square's four corners cannot all lie on one line, nor on the union of two *adjacent* closed polygon edges sharing a vertex. Indeed, if two corners belong to each of two adjacent supporting rays, those two pairs must be pairs of square sides (impossible because their lines are not parallel), or diagonals (impossible because their common midpoint would have to be the common ray endpoint, whereas both corners on each ray are on the same side of that endpoint). Three collinear square corners are impossible. Thus at least two of the four corners lie on **nonadjacent** polygon edges. Their Euclidean distance is at least \(d(P)\) and at most the square's diagonal \(\sqrt2\,s\). Therefore \(s\ge d(P)/\sqrt2\). QED.

Because two nonadjacent closed polygon edges are disjoint and there are finitely many of them, \(d(P)>0\) for every fixed simple polygon. The exact verifier computes \(d(P)^2\) as the minimum of all four endpoint-to-opposite-segment distances for every nonadjacent edge pair; disjoint segment distance always occurs at such a pair of points.

### A uniform family-wise bound, independent of a square branch

For our nine-vertex moving-kink family, the **global edge-box separation** \(g>0\) above implies \(d(P(\lambda))\ge g\) for every real \(\lambda\). Therefore **every** inscribed square (not just our two-chart branch) satisfies

\[
\boxed{s^2\ge \frac{g^2}{2}
=\frac12\left(\frac{1810915}{8503594}\right)^2>0.}
\]

This is the branch-independent noncollapse principle. The stronger explicit branch theorem gives \(s_{\max}>3\), but this separate feature-gap bound would remain valid even if the square branches underwent births, deaths, or chart exchanges. It requires no tracking of individual square roots, only a fixed feature gap and the known existence of squares in this specific family.

## 5. Why extreme corner angles alone are not the answer

There is a clean, exact counterexample to **using the feature-gap lemma alone as a universal argument**.

Start with the square \([-1,1]^2\). For \(0<\delta<1/4\), add a small *outward isosceles triangular tooth* along its upper boundary. Its boundary vertices in cyclic order are

\[
(-1,-1),(1,-1),(1,1),(\delta,1),(0,1+\delta),
(-\delta,1),(-1,1).
\]

The tooth's apex angle is **exactly 90 degrees**, and its two foot interior angles are **exactly 225 degrees**, independent of \(\delta\). The four original corner angles remain exactly 90 degrees. Thus the polygon's full turning-angle pattern does not vary when the tooth is made arbitrarily small.

Yet the two nonadjacent tooth-related edges have endpoints at distance at most \(\sqrt2\delta\), so

\[
0<d(P_\delta)^2\le2\delta^2\longrightarrow0.
\]

The polygon still **contains the closed unit disk** centered at the origin, since it contains the original square. Therefore a positive inradius and fixed extreme angles do **not** prevent the *feature-gap* \(d(P)\) from tending to zero.

**Crucial logical distinction:** This does **not** show that the largest inscribed square on these toothed polygons gets small. In fact, their large square boundary remains essentially intact. It only shows that the elementary feature-gap lower bound cannot produce a universal positive Big Peg constant from corner-angle data alone.

## 6. Reproducibility and negative controls

All proof decisions are made with standard-library Python exact Fractions:

~~~sh
python code/big_peg_corner_v250.py
python -m unittest -v tests.test_big_peg_corner_v250
~~~

The new CI workflow replays both commands, the predecessor v249 two-chart certificate, and the v248 rough-curve persistence certificate. Tests reject a flattened (collinear) kink, wrong supporting chart, a missing/invalid feature gap, false oversize square thresholds, nonaffine vertex motions, self-intersecting polygons, and malformed or floating-point tooth scales.

### Research status

- **Proved:** A genuine noncollinear corner does not necessarily destroy a large square, even when the chart's supporting line changes. A finite exact atlas can prove all-time continuation and a uniform lower side for this specific family.
- **Proved:** With a uniform separation between nonadjacent edges, **all** squares in a family are bounded away from zero, independently of which branch realizes them.
- **Disproved as a sufficient strategy:** A uniform positive separation between nonadjacent edges follows from fixed vertex angles and fixed positive inradius. The fixed-angle tooth family refutes this implication.
- **Still open:** A geometric/topological replacement for the feature-gap assumption that survives refinement to infinitely many corners. In particular, no theorem is established for arbitrary singular square-orbit bifurcations or all Jordan curves. The unrestricted Big Peg positive-c conjecture and Toeplitz remain unsolved.

## Next research direction

The real task is to build a quantitative **multi-branch** invariant of the square configuration space: if an isolated large square disappears, can total parity, degree, or another conserved quantity force at least one *other large square* to remain? This invariant must exclude the accumulation of all large square branches into features of vanishing size. A feature-gap lower bound alone is inadequate, as the tooth family explicitly demonstrates.
