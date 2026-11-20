# v254: A genuine positive-size square-pair death at a noncollinear polygon corner

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Result:** An explicit, continuously varying family of **strictly convex Jordan pentagons** in which **two distinct positive-size inscribed squares collide into a single geometric square at a polygon vertex, then both disappear as valid boundary squares**. Their four-edge square-incidence matrices remain nonsingular at the collision. The singularity occurs in the *glued geometric square-orbit configuration space*, not as a rank drop of the individual Cramer charts. A third square persists.

**Scope:** The algebraic local event and its one-sided branches are certified at **every real parameter**, with a complete independent geometric-square census at four **rational parameter slices**. We do **not** claim an all-time census of every other possible square, a generic signed cobordism, or a solution of Toeplitz.

## 1. Why this is not another cutoff crossing

v252–v253 tracked two square configurations that both continued to exist while their sizes crossed a selected positive cutoff. Those events modified the *filtered* square count but did not extinguish any geometric square orbit.

v254 exhibits a fundamentally different phenomenon:

> Two different geometric squares exist on one side of a critical time. Their four corners converge to the **same nonzero square**, one corner reaches an actual **noncollinear polygon vertex**, and the two distinct supporting-edge solutions become invalid on the opposite side.

The positive square itself exists at the collision instant. There are two regular algebraic branches before, one common geometric square exactly at the event, and no feasible continuation of either of those two particular branches after. This is an actual pair annihilation of geometric square configurations, not a choice of size cutoff.

## 2. An exact rational moving pentagon

Let \(0.6\le t\le0.65\). The counterclockwise polygon vertices are

\[
\begin{aligned}
V_0(t)&=(5355/10000-t/5,\;6134/10000),\\
V_1&=(-5770/10000,\;10693/10000),\\
V_2&=(-12571/10000,\;2233/10000),\\
V_3&=(-8305/10000,\;-3562/10000),\\
V_4&=(285/10000,\;-7637/10000).
\end{aligned}
\]

Only \(V_0\) moves. All coordinates remain exact rational-affine functions of the real parameter \(t\).

The new [v254 verifier](../code/big_peg_pair_death_v254.py) proves the polygon is **strictly convex and Jordan for every real** \(t\in[3/5,13/20]\): every successive signed turn is strictly positive by exact Bernstein signs, every edge has strictly positive length, and axis-aligned enclosures of every pair of nonadjacent moving edges remain disjoint with positive rational separation. Endpoint Jordan checks independently verify orientation.

### The three tracked square configurations

Represent an ordered counterclockwise square by \(z=(p_x,p_y,u,v)\), with corners

\[
C_0=p,\quad C_1=p+(u,v),\quad
C_2=p+(u-v,u+v),\quad C_3=p+(-v,u).
\]

The three edge choices are

\[
\begin{array}{ll}
C=(0,1,2,4)&\text{persistent square},\\
A=(0,1,3,4)&\text{first disappearing square},\\
B=(0,2,3,4)&\text{second disappearing square}.
\end{array}
\]

Each has a rational Cramer root \(z_E(t)=N_E(t)/D_E(t)\), with polynomial \(N_E,D_E\). The verifier **independently recomputes** all three matrices and proves:

1. Every \(D_E(t)\) has constant **strictly nonzero** sign throughout the entire time interval, even at the collision.
2. All three supporting-line square candidates have squared side **strictly greater than 1** throughout the interval, with exact Bernstein signs.
3. The persistent square \(C\) has **all eight** strict segment-membership inequalities positive at every real time.
4. Squares \(A,B\) each satisfy **seven** strict inequalities everywhere. Their eighth inequalities, which determine whether \(C_1\) lies on the correct side of vertex \(V_2\), have an exact common algebraic factor.

## 3. The genuine pair-death event

That common factor is the explicit **rational quadratic**

\[
\boxed{
Q(t)=72468080000000t^2
-496738335230000t
+291143052416293.
}
\]

Exact rational Bernstein signs establish \(Q'(t)<0\) throughout \([3/5,13/20]\). The endpoint signs are opposite; therefore \(Q\) has a **unique** real zero \(\alpha\) in this interval.

The existing v229 algebraic-time verifier independently checks the strict rational isolating bracket

\[
\boxed{\frac{64722}{100000}<\alpha<
\frac{64723}{100000}.}
\]

**The two exceptional edge-contact conditions factor exactly as**

\[
\begin{aligned}
\text{A's distance from the end of edge 1, signed}&=Q(t)\,H_A(t),\\
\text{B's distance from the start of edge 2, signed}&=Q(t)\,H_B(t),
\end{aligned}
\]

where *both* \(H_A\) and \(H_B\) are proved **strictly positive everywhere in the full interval** by exact rational Bernstein coefficients.

Accordingly:

- For \(t<\alpha\), both squares \(A\) and \(B\) lie strictly on all their four selected edges, are positive-size and represent **two distinct geometric squares**.
- At \(t=\alpha\), the exceptional inequalities are equalities and both squares touch \(V_2\).
- For \(t>\alpha\), both become strictly infeasible as squares of the polygon, because their respective second corners fall *off* their selected boundary segments.

This proves the one-sided appearance/death by exact inequalities **for every real parameter**, not sampled numerical interpolation.

### They actually collide geometrically, not merely in side length

For all four Cramer coordinates \(j=0,1,2,3\), the verifier checks the polynomial identity

\[
Q(t)\ \bigm|\ N_{A,j}(t)D_B(t)-N_{B,j}(t)D_A(t).
\]

Because \(D_A(\alpha),D_B(\alpha)\ne0\), we have

\[
\boxed{z_A(\alpha)=z_B(\alpha).}
\]

The square's **second corner** at that root equals the original pentagon vertex \(V_2\), proved by an additional exact polynomial divisibility test in *both* Cartesian coordinates.

Conversely, the **first-coordinate difference divided by \(Q(t)\)** has a certified nonzero sign on the whole interval. Thus the two square configurations are *genuinely different* for every \(t\ne\alpha\), including all \(t<\alpha\), rather than duplicate edge labels for one orbit.

The third square \(C\) is **not** the same as this collision square: its Cramer first-corner difference has **nonzero algebraic sign at \(\alpha\)**, checked with the existing v229 real algebraic root-sign method.

All squares at the event have nonzero side, with \(s^2>1\). Nothing is shrinking toward a point.

## 4. Independent COMPLETE unfiltered square censuses at rational times

The pre-existing [v223 exact geometric-orbit census](../code/polygon_orbit_census_v223.py) does not use our selected Cramer branches. It searches all edge assignments (including repeated edges), all square orientations, polygon vertex hits and singular systems, deduplicating by the four geometric corners. The v254 verifier reruns it with **zero size cutoff** at four rational times.

| Parameter \(t\) | Position relative to \(\alpha\) | ALL positive geometric square orbits |
|---|---|---:|
| \(3/5=0.6\) | Before | **3** |
| \(64722/100000=0.64722\) | Immediately before | **3** |
| \(64723/100000=0.64723\) | Immediately after | **1** |
| \(13/20=0.65\) | After | **1** |

Every square counted by the independent rational-slice censuses has \(s^2>1\). The sample counts confirm the local algebraic pair death is accompanied by a **3→1 decrease** in the total genuine geometric square count across the isolated event.

**Important distinction:** These four complete *slices* do not, by themselves, prove that exactly three and exactly one squares exist at *every* real time on the respective sides. We are not claiming that additional configurations cannot have their own unrelated events elsewhere in the window. The all-real-time proof covers the three explicitly identified branches; the independent v223 all-square counts are complete at the four stated rational times.

A stronger global pentagon event ledger, across *all* repeated and singular edge charts for every algebraic parameter time, is still future work.

## 5. Mathematical significance for Toeplitz

We now have two distinctly different types of certified event:

| Event type | Established milestone | Effect on geometric square orbits |
|---|---|---|
| Size-boundary crossing | v252–v253 | A square remains present, but enters/leaves the **filtered** count |
| Actual geometric pair death at polygon corner | **v254** | **Two distinct squares merge into one** at the vertex event, then both disappear as valid square configurations |

The regular-time orbit count changes by **two**, so its parity at rational times on opposite sides agrees, as expected from common generic polygonal parity phenomena. The single shared geometric square at the *singular event instant* need not obey ordinary generic-time parity rules.

Equally importantly, the two branch Cramer determinants do **not** vanish. A purported square-orbit birth/death classifier that looks **only for rank-deficient 4×4 incidence matrices would MISS this event.** Topological gluing across polygon vertices is an independent source of singular behavior and must be included in any complete event theory.

This demonstrates a genuine local mechanism for appearance/disappearance of positive-size squares without disproving their total existence: the third square persists throughout the deformation.

## 6. Reproduce and limitations

Pure Python standard library, exact Fraction arithmetic and prior v223/v224/v229/v249/v252 modules:

~~~sh
python code/big_peg_pair_death_v254.py
python -m unittest -v tests.test_big_peg_pair_death_v254
~~~

The [19 exact tests](../tests/test_big_peg_pair_death_v254.py) cover all-time Jordan geometry, uniqueness/isolation of the algebraic event, all eight edge-contact inequalities for each branch, exact factor/divisibility checks, distinctness except at collision, positivity, persistent third square, all four independent complete rational-slice censuses, and fail-closed negative controls.

The dedicated CI also replays v253's global square-size bottleneck theorem and v229's singular rational-function germ audit.

**Not proved:** A full all-real-time census of *all* square orbits for this pentagon homotopy; a general polygonal square-orbit cobordism including all vertex and rank-degenerate events; a universal positive Big Peg constant; or the unrestricted Toeplitz square-peg conjecture.

### Future research

Construct a **finite event graph** with both rank-singular and vertex-gluing strata, classify local orbit multiplicities, and develop a degree/parity accounting rule with two kinds of boundary:

1. square-size cutoff crossings; and
2. the actual geometric orbit-branch bifurcations established here.

The remaining obstruction to Toeplitz is uniformity across arbitrary polygon complexity, where vertex-gluing events can accumulate at ever smaller scales.
