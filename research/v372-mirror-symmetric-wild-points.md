# v372. Mirror-symmetric wild points

[v370](v370-half-turn-wild-points-cannot-absorb-parity.md) used the half-turn: it reverses the local arc and keeps square orientation, so it trades small aligned squares for small reversed ones. Among isometries fixing a point, the only other kind that can reverse an arc through that point is a reflection. This note works out that case.

A reflection reverses both the arc and the square's orientation, so it keeps every small square in its own sector. Pairs then cancel within each sector. The only unpaired squares are those symmetric about the mirror, and they correspond to chords of a single arm. If each arm stays inside a 90° cone around the perpendicular to the mirror, no such chord exists, and the point cannot absorb the odd aligned parity.

**Results.**
1. **Reduction (Lemma 372.1).** At a mirror-symmetric wild point, the small aligned degree has the parity of the number of small mirror-symmetric aligned squares, and likewise for reversed. Mirror-symmetric squares are exactly the chords $P\mapsto P+(0,2x_P)$ of one arm.
2. **Cone points (Theorem 372.2).** Suppose a Jordan curve is polygonal away from $p$, and near $p$ it is an arc symmetric under a reflection through $p$, with each arm inside the open 45° cone about the perpendicular to the mirror. Then it inscribes a square of explicit size. Inside the cone the arms may curl arbitrarily.
3. **Several wild points (Theorem 372.3).** Finitely many wild points, each of half-turn or mirror-cone type, can be combined.
4. **(H3) is automatic.** In Theorem 370.2 the drawing can be cut at the last exit of an arm from a small ball, so hypothesis (H3) is not needed.
5. **Exact checks:**
   - on curled cone arms, the small aligned and reversed counts are both even, and the reflection maps each sector onto itself with no fixed square;
   - on a drawing small enough for a full census, the large aligned count is odd;
   - when an arm leaves the cone and creates exactly one small chord, exactly one small sector becomes odd, and its extra square is the mirror-symmetric one.

## Setting and reduction

Put $p$ at the origin and the mirror on the $y$-axis, with reflection $\rho(x,y)=(-x,y)$. Suppose $\gamma\cap\bar B(p,R)$ is one arc $\alpha$ through $p$ with $\rho(\alpha)=\alpha$, and $\rho$ swaps its two arms.

A point of an arm on the axis would be fixed by $\rho$, so it would lie on both arms. Hence $\alpha$ meets the axis only at $p$, and arm one lies in $x>0$.

**Drawings without (H3).** For small $r$, let $a$ be the last point of arm one in $\bar B(p,r)$, going outward. Replace the subarc of $\alpha$ between $\rho(a)$ and $a$ by $[\rho(a),p]\cup[p,a]$. The outer remainder of the arm lies outside $\bar B(p,r)$ except at $a$, and the new piece lies inside. So the result $D_r$ is a Jordan polygon, mirror-symmetric in $B(p,R)$, and $D_r\to\gamma$ uniformly as $r\to0$. The same cut, with the segment $[a,-a]$, gives the drawings of Theorem 370.2 without assuming (H3).

**Lemma 372.1.** Let $s$ be such that every square of side below $s$ on $D_r$ lies in $B(p,R)$, and suppose no square has side exactly $s$. Then:
- the small aligned degree of $D_r$ is congruent to the number of small aligned squares $Q$ with $\rho(Q)=Q$, and likewise for reversed;
- such squares are exactly those with corners $(\pm u,y)$ and $(\pm u,y+2u)$, where $P=(u,y)$ and $S(P)=(u,y+2u)$ both lie on the right half of $D_r\cap B(p,R)$, and $S(x,y)=(x,y+2x)$.

*Proof.* $\rho$ maps the local arc to itself and reverses its direction, so it reverses the order in which a square's corners occur. It also reverses the square's orientation. The two reversals cancel, so $\rho$ maps small aligned squares to small aligned squares, and reversed to reversed. Away from symmetric squares the involution is free, and an equivariant perturbation near the zero set pairs the zeros. So each small sector degree has the parity of its symmetric squares.

A square symmetric about the $y$-axis has the axis either as a diagonal or as a midline. As a diagonal, two corners would lie on the axis, which meets the local arc only at $p$. As a midline, the square is axis-aligned, with the corners described. ∎

This turns a four-corner question into a two-point question about one arm. On polygons the chords are found exactly, by solving a $2\times2$ system for each pair of edges.

## Cone points

**Theorem 372.2.** Let $\gamma$ be a Jordan curve with the following properties:
- $\gamma\setminus B(p,r)$ consists of finitely many segments for every $r>0$;
- $\gamma\cap\bar B(p,R)$ is an arc symmetric under the reflection $\rho$ in a line $L$ through $p$, with $\rho$ swapping the arms;
- each arm, apart from $p$, lies in the open cone of half-angle $45°$ about the perpendicular to $L$ on its side.

Then $\gamma$ inscribes a square of side at least $s_0$, where $s_0$ is chosen as in Theorem 370.2 from the polygonal part outside $B(p,R/16)$.

*Proof.* In coordinates, the right half $\beta$ of the local arc of $D_r$ lies in $\{|y|<x\}\cup\{p\}$, since the segment $[p,a]$ stays in the convex cone. Then $S(\beta)\subset\{y>x\}\cup\{p\}$ and $S^{-1}(\beta)\subset\{y<-x\}\cup\{p\}$. So a chord $P$, $S(P)$ inside $\beta$ would have $P=p$ and side zero. By Lemma 372.1, both small sector degrees of $D_r$ are even.

Lemma 370.1 gives aligned degree 1 for $D_r$. So the large aligned degree is odd, and $D_r$ has a square of side at least $s_0$ unless one has side exactly $s_0$. A subsequence converges as $r\to0$ to a square of side at least $s_0$ on $\gamma$. ∎

Here the large *aligned* count is odd by itself. In the half-turn case only the sum over the two graceful sectors was.

**Which points are new.** A cone arm cannot wind around $p$, so a slow spiral is impossible. It can still curl. Take finite double spirals of size $c_k\to0$ along the arm, each with about $c_k^{-2}$ turns. Each curl swings the primitive $\int y\,dx$ up by roughly a fixed amount and back down. So the primitive has non-vanishing oscillation at $p$ and no continuous extension, which is outside Asano and Ike's condition as I read it. The curls give infinite $p$-variation for every $p<2$, and the point is not locally monotone. As in v370, these are estimates against sufficient conditions. Any one such curve may have obvious squares elsewhere, so the content is the theorem for the class.

**Several points.**

**Theorem 372.3.** Let $\gamma$ be polygonal away from finitely many points $p_i$, each of the type in Theorem 370.2 (half-turn) or Theorem 372.2 (mirror cone), in disjoint balls. Then $\gamma$ inscribes a square.

*Proof.* Choose $s$ so small that every small square has a corner near exactly one $p_i$, and lies in that point's ball. Cut every point at once. At a half-turn point the small squares contribute equally to both graceful sectors. At a mirror-cone point they contribute an even number to each. So the small aligned plus small reversed degree is even, and the large graceful degree is $1+0\equiv1$. ∎

## Exact checks

Arm one runs from the origin along the $x$-axis through curls at $x=2^{-k}$, each a polygonal double spiral of size $2^{-k}/4$. Arm two is its mirror image. The outer closure is generic and not symmetric. Every drawing is checked to be Jordan, and the counter is the exact transverse counter of v370.

| drawing | vertices | side $<1/16$: aligned, reversed, crossed | fixed squares |
|---|---|---|---|
| depth 1, 2 turns per curl | 72 | 60, 60, 150 | none |
| depth 2, 2 turns per curl | 138 | 210, 204, 482 | none |

Both sector counts are even at each depth. The reflection maps each sector onto itself, and a full census of a small curled drawing confirms that its large aligned count is odd.

**Controls outside the cone.** Two zigzag arms leave the cone, and each has exactly one chord with $u<1/8$. In one, the small squares below side $1/4$ number 1 aligned and 0 reversed. In the other, 0 aligned and 1 reversed. In both, the odd square is the mirror-symmetric square built on that chord. So the cone is doing the work. Without it, the parity in a sector is the parity of the chord count, and nothing forces that to be even.

## Where this leaves the local question

The isometric involutions fixing a point and reversing an arc through it are the half-turn and the reflections. Both kinds now give theorems:

| symmetry | effect on sectors | condition | what stays odd |
|---|---|---|---|
| half-turn | swaps them | none | large aligned plus reversed |
| reflection | keeps them | no chords of the arm (the 45° cone suffices) | large aligned |

v371 extends both to curves joined to such a point without passing a square of side $s$.

A wild point with no isometric symmetry needs a different handle. The reflection case suggests one. There, the obstruction reduced to the parity of a one-arm chord count. For a general point, the analogue would be a reduction of the small-square parity to a lower-dimensional count that can be controlled at every depth.

[Code](../code/mirror_wild_point_v372.py) · [Tests](../tests/test_mirror_wild_point_v372.py)
