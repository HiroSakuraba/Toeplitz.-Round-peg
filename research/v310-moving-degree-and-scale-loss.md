# v310 — What a change of square-count parity forces on the limiting curve

Benjamin John Schulz · November 2026

Two sequences of polygons can converge rapidly to the same convex quadrilateral while disagreeing, at every step, about the parity of the square count above an entire interval of sizes. This obstruction comes with a positive conclusion for arbitrary Jordan limits. If the disagreement persists at one positive size, the limit realizes that exact size. If a positive measure of cutoffs keeps changing, the limit realizes at least that measure of distinct squared side lengths. A different estimate below certifies the square affected by a local edit even if the new arc oscillates infinitely often.

## 1. Transfer across a wall

Use physical squares, each counted once. Write the four vertices as the cyclic quarter turns of a radius vector about the center. The parameters are center, positive half-diagonal and angle modulo a quarter turn. Write $s$ for squared **side** length, twice the squared half-diagonal. Evaluating a continuous signed-distance function at the four vertices gives a section of a real rank-four bundle over this four-dimensional space. A change of first vertex cyclically permutes the four coordinates. Work modulo two. At isolated transverse zeros the degree counts geometric squares modulo two. Local degree still exists on compact regions with zero-free boundary when the zero set is singular. See [the bundle construction](v307-convex-shell-square-gap.md).

For a motion of these sections, restrict to a compact region bounded by center and size faces. Require the side-face intersections, including their seams with the two end faces, to be zero-free. The difference of the two end degrees equals the sum, modulo two, of the degrees on the side faces times the time interval. To see this, make one sufficiently small transverse perturbation on the whole product, compatible on faces and disjoint from the seams at its zeros. The zero set is a compact one-manifold and has an even number of boundary points. Its counts give the stated degrees by homotopy invariance. This argument allows actual, even nonisolated, zeros in the interior of a moving wall. Adjacent regions use the same restricted section, so a shared wall contributes equally to both. A signed-distance interpolation can have auxiliary zeros off an actual curve. When the homotopy is itself a path of Jordan boundaries, a wall zero is an actual square on an intermediate curve.

## 2. The limiting curve pays for repeated parity changes

Fix $0<a<b$. Define $S(\Gamma)$ to be the set of $s\in[a,b]$ attained by geometric squares on a compact Jordan curve $\Gamma$. Let $E_\delta$ be the corresponding set for the closed Euclidean $\delta$-neighborhood of the curve, allowing the four vertices anywhere in the neighborhood. The centers of the candidate squares stay bounded. The size window and angle circle are compact. Consequently the $E_\delta$ are compact, nested as $\delta$ decreases, and their intersection is exactly $S(\Gamma)$. In particular their lengths decrease to $|S(\Gamma)|$, where length means one-dimensional Lebesgue measure.

For a polygon with finite transverse squares write $P(s)$ for the parity above $s$. At regular cutoffs this also makes sense as the degree above $s$, for any Jordan curve with that cutoff free of squares. Truncate at a sufficiently large radius and center box, which introduce no zeros. Registration means a common parametrization with pointwise Euclidean displacement bounded by $\delta$.

**Theorem.** Suppose two Jordan curves are each within registered distance $\delta$ of $\Gamma$. For every $s\in[a,b]\setminus E_\delta$ their degrees above $s$ agree. Hence for regular polygon endpoints the length of their parity-disagreement set is at most $|E_\delta|$.

**Proof.** Outside the closed tube, both endpoint winding numbers agree with that of $\Gamma$ up to orientation, by straight interpolation of their parametrizations with the reference parametrization and invariance of inside/outside. Thus the endpoint signed distances have the same strict sign there. Every point zero of their linear interpolation belongs to the closed tube. A four-corner zero at scale $s$ would put $s$ in $E_\delta$. Compact degree continuation above $s$ now gives equality. No embedding assumption on an interpolation of the parametrizations is required.

**Consequences.** Let pairs of approximants converge in registered distance to the same Jordan curve.

- Every positive cutoff with infinitely many parity disagreements is realized by a square of precisely that size on the limiting curve. Apply the theorem in every fixed tube.
- If the limsup of the lengths of the disagreement sets in $[a,b]$ is $m$, then $|S(\Gamma)|\ge m$. The cutoffs can move at every step. Positive $m$ therefore forces uncountably many distinct sizes. This conclusion alone need not give an interval.
- If $S(\Gamma)$ has measure zero in the window, all finite-transverse approximations converge in $L^1$ there to the degree profile of $\Gamma$ (define it arbitrarily on its null spectrum). Use the same argument with the second endpoint equal to $\Gamma$ and let the tube shrink.

All three statements concern a **fixed positive window**. They allow every sufficiently close approximation to have zero degree above every cutoff in that window. Ruling out escape through sizes tending to zero remains essential.

## 3. A convex discontinuity and its explicit limit

Start with $K=((0,0),(1/2,-3/8),(1,0),(1/2,3/4))$, in cyclic order. Apply either $A_d(x,y)=(x,y+dx)$ or $B_d(x,y)=(x,(1+d)y)$, with $0<d\le 1/100$. Each resulting quadrilateral has precisely one square. The first has $s>39/100$ and center $x<47/100$. The second has $s<3/10$ and center $x=1/2$.

Here are details of the certification. Convexity forces cyclic contact order. Two distinct supporting edges alone cannot contain all four vertices unless they are parallel and carry opposite sides: a supporting line through two square vertices carries a side. These quadrilaterals have no parallel edges. After cyclic relabeling the 13 nondecreasing four-edge tuples with at least three distinct edges exhaust the possibilities. A repeated edge split by the starting label can be put first or last by restarting at its first contact. On each chart the four square equations are affine in the four edge tags. Their determinants and Cramer numerators are polynomials of degree at most four in $d$. The replay checks exact Bernstein signs on the whole half-open interval. The shear winner is $(1,2,3,3)$ with strictly ordered repeated contacts, as in [v301](v301-tilted-squares-and-boundary-weights.md); the stretch winner is $(0,1,2,3)$. The other charts have a tag strictly outside $[0,1]$. Both winners have interior tags and nonzero determinant, hence local degree one. Transversality can equivalently be checked by intersecting the linear subspace of labeled squares in eight corner coordinates with the product of the four edge lines.

For the stretch put $\lambda=1+d$. Its winning square is axis-aligned with left coordinate $4/(8+9\lambda)$ and side $9\lambda/(8+9\lambda)$. For the shear, substituting the v301 tags in the side from the final edge to edge 1 gives $v_x>21/40$ and $-v_y>43/125$; the sum of their squared lower bounds exceeds $39/100$. The signs, center inequality and stretch bound are checked after clearing positive denominators.

Along the straight interpolation of these affine transformations, each quadrilateral is convex and moves at most $d$ from $K$. The disk centered at $(1/2,1/8)$ of radius $1/4$ stays inside: the minimum supporting-line distance on $K$ is $5/(4\sqrt{13})$ and convex support functions change by at most $d$. The [convex inball bound](v307-convex-shell-square-gap.md) gives $s\ge 1/8$ throughout. The region $s>1/10,\ x<47/100$ changes degree from one to zero. Only its center wall can account for this change. With all centers allowed, the degree above any $q\in[3/10,39/100]$ changes from one to zero. Here the only contributing wall is $s=q$. Both transfers concern actual intermediate squares.

Take the shear for even $n$ and the stretch for odd $n$, with $d_n=2^{-n}/100$. Registered displacements between successive polygons sum to at most $\sum_n(d_n+d_{n+1})<\infty$. The profiles disagree on a window of length $9/100$ at **every** step and lose that length on each shear-to-stretch step. This rules out any modulus tending to zero that would bound the $L^1$ change by registered displacement alone, even for this convex class with a common inball. These are inscribed-square approximations on moving polygons. They do not require the new polygon to retain every earlier vertex.

The limit shows exactly where the sizes went. On successive edges of $K$ use the tags
$$
(2t-10/17,\;2t-9/17,\;t-1/17,\;t),\qquad 5/17\le t\le 9/17.
$$
They lie in $[0,1]$. Substitution gives a square whose first side vector is $(9/17,\tfrac32(t-9/17))$. The resulting squares realize **every** squared side length from $81/289$ to $117/289$. The shear and stretch select opposite ends of this subfamily as $d\downarrow0$.

## 4. A geometric sufficient condition for retaining positive size

For finite-transverse polygons converging uniformly to a Jordan curve, write $M_i=\int_a^b P_i(s)\,ds$ and $L_i=|\{s:P_i(s)=1, P_{i+1}(s)=0\}|$ in that window. If $\sum L_i\le L<M_0$, the limit contains a square of squared side at least $a+M_0-L$. Indeed $M_n\ge M_0-L$ by telescoping. The support of $P_n$ lies below the largest squared side, so the latter is at least $a+M_n$. Compactness gives four limiting vertices on the curve, with the positive lower bound intact. This is a sufficient criterion. Section 3 shows why a finite length of boundary motion cannot supply its hypothesis.

**Bounding a stage before knowing its squares.** Replace part of a boundary inside a closed square box $D$ with center $p_0$ and halfwidth $e$. Assume the two boundaries and their inside/outside designations agree off $D$. Assume the common boundary there lies on a supplied finite set of lines. Under signed-distance interpolation, point zeros off $D$ lie on the unchanged arcs. At a cutoff $s>8e^2$, at most one corner fits in $D$. If all corners avoid $D$, their square already belonged to the old curve. For finite source spectra, these unchanged cutoffs have measure zero. For a square meeting $D$, write its corners as $p,p+v,p+Jv,p+(I+J)v$, with $p\in D$ and $J$ counterclockwise quarter turn. Assign the other three corners to an ordered triple of supporting lines, allowing repeated lines.

For a line $n_i\cdot x=h_i$ and its corner matrix $T_i=I,J,I+J$, the equation is $(n_i^T T_i)v=h_i-n_i\cdot p$. Each invertible pair gives $v=c+B(p-p_0)$, so component $k$ differs from $c_k$ by at most $e\sum_j|B_{kj}|$. Interval squaring and summation enclose $s=|v|^2$. Intersect the enclosures from all invertible pairs and union over triples. Retain the whole window for a triple having no invertible pair, and retain the part below $8e^2$. This cover uses rational arithmetic for rational inputs. Outside it, apart from the finite unchanged spectrum, degree cannot change. Its intersection with the original odd set bounds $L_i$. Repeating with the actual fixed arcs at each stage permits varying boxes and different edge directions.

**A square forced to touch a rough replacement.** Take the triangle $(0,0),(4,0),(1,3)$, the box centered at $(4/7,12/7)$ with $e=1/4096$, and the size window $[1,4]$. Exact census gives three transverse squares, of squared sides $72/25,144/49,360/121$. The middle square has the box center as its upper left corner. The other two lie outside the narrower band $(2.9,2.95)$ in size. The union over all 27 line triples, after clipping to the large window, consists of just one interval
$$
U=[9662300209/3288334336,\;4832526581/1644167168]
\subset (2.93835,2.93920).
$$
Any Jordan replacement inside the box with the same boundary and sides off the box has a square touching the box with squared side in $U$. In fact, the narrower band's boundary has no homotopy zeros, and its initial degree is one. Every square in the band avoiding the box would be an original triangle square, yet the only original square in the band touches the box. The single-corner estimate therefore puts the guaranteed square in $U$. This conclusion applies directly to an infinitely oscillating replacement, including a Jordan limit of infinitely many edits in the box. If the replacement removes the old contact, the forced square is geometrically different. The unchanged triangle squares cannot witness the assertion. For the separate integral estimate the initial mass is $284063/148225$ and one-step loss is at most $98339/234881024$.

## 5. All four original contacts can be removed at once

For a finite union of edit boxes use the one-box line estimate for configurations with just one corner in the union. For two corners in the same box retain sizes up to its squared diameter. For two distinct boxes centered at $p,q$ with halfwidths $e,f$, put $r=e+f$. Their squared distance is in
$$
[A,B]=\left[\sum_{k=1}^2 (|p_k-q_k|-r)_+^2,\;\sum_{k=1}^2 (|p_k-q_k|+r)^2\right].
$$
Add both $[A,B]$ and $[A/2,B/2]$ to the size cover, since these two corners can be adjacent or opposite. Take the union over unordered box pairs. This covers three or four moving corners as well. The proof uses no derivatives on the replacing arcs. As before, cutoffs from squares avoiding every box are unchanged and must be accounted for when testing a particular wall.

Put a box of halfwidth $1/4096$ around **each** vertex of the triangle's middle square, at $x=4/7,16/7$ and $y=0,12/7$. Compute the combined cover in the isolating band $(2.9,2.95)$. It is a single closed interval
$$
[603635761/205520896,\;302161969/102760448].
$$
It lies strictly inside the band. The same degree-one argument forces a square in this closed interval on **every** Jordan replacement with boundary and inside/outside unchanged off the union of boxes. The four arcs can all be rough and all four original contacts can be removed. Any square in the band must touch at least one box because the unique triangle square in the band has all its contacts in the boxes. This result does not specify which boxes the individual vertices of a new square occupy.

The replay also constructs a 19-vertex Jordan polygon by cutting out the four original contacts. Each replacement has a middle segment parallel to its original edge, shifted up by a different multiple of $1/(128\cdot4096)$. Exact segment tests exclude all four old contacts. Solving on the four new middle segments gives an actual square of squared side $39582494097493/13469017440256$, strictly different from $144/49$ and inside the certified interval. This last computation verifies one witness without making a uniqueness claim.

## What remains open

The disagreement theorem recovers a positive-size square whenever parity keeps changing at a positive cutoff. The remaining scenario is stabilization to zero on every fixed positive window while all polygon squares concentrate near zero. The multi-box estimate permits all four corners to change, but its union of pair-distance intervals can fill the whole window as the number of boxes grows. The next target is to control that union for the successive refinements of one arbitrary Jordan curve, or to force disagreement at a fixed positive cutoff. Neither follows from the rate of uniform approximation. For the classical compactness obstruction see [Matschke's survey](https://www.ams.org/notices/201404/rnoti-p346.pdf).

Replay with `python code/scale_transfer_v310.py` and `python -m unittest discover -s tests -p 'test_scale_transfer_v310.py' -v`. [Exact output](../results/v310/scale-transfer.json).
