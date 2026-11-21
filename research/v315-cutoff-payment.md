# v315. Which cutoff crossings force a square on the target curve?

Benjamin John Schulz · November 2026

A reference curve can have odd degree above a side length $s$ while every chosen approximation to a different target has even degree there. A crossing somewhere along a path between them only gives a square for the intermediate defining field. We give a sufficient condition that places a side-$s$ square on the **target**. It allows some intermediate squares to stay a fixed positive distance away forever. We also show how to turn a near-target cutoff configuration into a square on an arbitrarily close polygonal approximation with nonzero local degree.

Throughout, $s$ denotes **side length**. The parameter called $s$ in [v310](v310-moving-degree-and-scale-loss.md) was squared side length. Degree counts each physical square once modulo two. The geometric-square bundle and the half-turn seed are described in [v313](v313-half-turn-degree.md).

## Relative degree with only part of the cutoff wall paired

Let $P_n$ converge in registered uniform distance to a Jordan curve $\Gamma$. Choose cutoffs $s_n\to s>0$. On each stage start with a section of odd degree above $s_n$ and end with the signed-distance section for $P_n$. Join the two by a continuous section homotopy. Truncate centers and upper sizes to a compact domain containing all relevant zeros. Assume the outer faces stay zero-free. For the moment require both ends to have no zeros at the cutoff. The **cutoff wall** $W_n$ has coordinates (center, angle modulo a quarter turn, homotopy time) at the single side length $s_n$. Its dimension and the section rank are both four.

For a square configuration $Q$ write $d_\Gamma(Q)=\max_{x\text{ vertex of }Q}\operatorname{dist}_\infty(x,\Gamma)$. Suppose that on $W_n$ there is an isolating open part $V_n$ with compact closure away from the wall boundary, zero-free frontier and mod-two degree zero. One may use a finite disjoint union of regions here. Require that **every wall zero outside $V_n$** satisfy $d_\Gamma(Q)\le\eta_n$, where $\eta_n\to0$. Either part may be empty. All data may change with $n$.

**Theorem.** The limit curve has a square of side exactly $s$, or its degree above $s$ is defined and odd. In particular it has a square of side at least $s$. If the endpoint degrees are even infinitely often, the limit realizes $s$ exactly. If an endpoint $P_n$ already realizes $s_n$ on infinitely many stages, the same exact-size conclusion follows directly. A cutoff zero only at the starting reference has no such implication.

**Proof.** At a stage with even endpoint degree, the cobordism boundary identity from v310 gives total wall degree one. Excision removes $V_n$ with contribution zero. The remainder has a zero. Its four corners lie within $\eta_n$ of $\Gamma$. Along any infinite sequence of such stages, centers are bounded, and angles and side lengths have a convergent subsequence. The limit is a square on $\Gamma$ of side $s$. The same compactness argument applies to any endpoint squares of side $s_n$.

If $\Gamma$ misses size $s$, compactness gives a closed neighborhood of $s$ that is free of squares even in a sufficiently thin closed tube around $\Gamma$. Otherwise square configurations in tubes of widths tending to zero with sizes tending to $s$ would converge to an exact square there. The preceding argument now excludes even endpoints for all sufficiently large $n$. Interpolation of signed distances between $P_n$ and $\Gamma$ has point zeros only in a shrinking tube by registration and winding. The gap therefore identifies the odd endpoint degree above $s_n$ with the degree of $\Gamma$ above $s$. This also explains why neither a common reference nor tracking branches is required.

**How to verify degree zero away from the target.** Suppose a free involution of $\overline V_n$ lifts linearly to the bundle, and the wall section $F$ is equivariant there. On the quotient, perturb the descended section transversely away from the zero-free boundary, then lift. Its zeros occur in pairs upstairs. This proves the required parity for singular as well as isolated zeros. There is also a boundary-only test. Let $TF$ denote the transported section under the involution and $G=(F+TF)/2$. If the straight segment from $F$ to $G$ avoids zero on the frontier, their relative degrees agree. Since $G$ is equivariant its degree is zero. With an invariant norm it suffices that $\sup_{\partial V_n}|F-TF|<2\inf_{\partial V_n}|F|$. This test allows asymmetric sections inside $V_n$ and a different pairing region on each stage. Two disjoint boxes with equal local degrees are another sufficient certificate.

The useful search condition is therefore to cover the cutoff wall by excluded cells, isolated collections of even degree, and cells in which **all four** vertices are close to the target. Only the last error must tend to zero. A vertex close to an edited arc while its three companions stay uncontrolled gives no limit square. A symmetric seed can be quite far from $\Gamma$ so long as the wall accounting is valid. Exact symmetry of arbitrarily close reference *curves* about bounded centers would itself force symmetry on their limit.

## Turning a tube witness into an odd-degree approximation

Let $Q=(q_0,q_1,q_2,q_3)$ have side $\ell$ and let $P$ be a registered polygonal approximation. Suppose each $q_j$ has a supplied contact $p_j\in P$ within coordinate sup distance $e>0$. If $\ell>16e$, one can move $P$ by at most $e$ in registered coordinate sup distance to a simple polygon containing **these four exact corners**. The corners may originally miss $P$ and the contacts may be vertices. When a wall witness is within $\eta$ of $\Gamma$ and $P$ is within $\delta$ of $\Gamma$, one may use $e=\eta+\delta$ (or an arbitrarily small strict enlargement).

For a single contact use the explicit map
$$
H_j(x)=x+(q_j-p_j)\max\!\left(0,\min\!\left(1,2-\frac{\|x-p_j\|_\infty}{2e}\right)\right).
$$
Its perturbation of the identity is $1/2$-Lipschitz in the sup norm. The map is injective with Lipschitz inverse on its image. For any desired output $y$, the inverse equation $x=y-(H_j(x)-x)$ is a contraction, so it is onto. It is a translation within radius $2e$, and the identity outside radius $4e$. Distinct $q_j$ are at sup distance at least $\ell/\sqrt2$, so distinct $p_j$ are separated by at least $\ell/\sqrt2-2e>8e$. The four open supports are disjoint and each homeomorphism preserves its support box. The composite moves every point at most $e$, sends each $p_j$ to $q_j$, and is piecewise linear. Evaluating it on $P$ with the original edge parametrization gives the claimed registered polygon. The construction uses rational arithmetic for rational data.

If the resulting square has singular contacts, the [marked polygon lemma](general-noncollapse-goal.md) regularizes it with arbitrarily small further displacement. Thus a sequence of positive-scale tube witnesses supplies explicitly repairable approximations with odd **local** degree above every fixed cutoff below $s$. On stages where $P_n$ already has odd degree above $s_n$, its degree region itself works above any smaller fixed size floor. This connects the exact-cutoff outcome to the project's noncollapse goal even when the unmodified approximations have even global degree.

## Exact example with distant paired crossings and even endpoints

Use the five-point protrusion profile in v313. An upper protrusion at $u$ has coordinates $(u,1)+\lambda p$ for profile points $p$, and the lower one at $b$ has $(b,-1)-\lambda p$. Horizontal boundary continues on either side of each mouth. Take the base rectangle $[-4,d]\times[-1,1]$, with $4\le d\le6$ and $|u|,|b|\le1$. All stages enclose the open unit disk. Fix $s=2/5$.

The exact normalized squared sides for a square on one profile and its baseline are $1386334197/9194892100$ and $533/2890$. Only the latter can reach $s$ for $9/10\le\lambda\le1$, at $\alpha^2=2312/2665$. Both contact systems are regular with strictly interior tags. The larger square crosses transversely in the scale variable because its squared side varies as $\lambda^2(533/2890)$. Completeness here follows as in v313. A cutoff square cannot be supported by the unmodified rectangle, which has only side-two squares. A protrusion contact rules out the opposite row and the vertical walls by distance. On the remaining baseline and protrusion the horizontal span equals the vertical span and is at most $17\lambda/32$. Hence the normalized contacts all lie within horizontal distance $7/16+17/32<1$ of the attachment. The existing finite profile census covers them throughout every motion below.

Write $e_n=1/(10\cdot2^n)$ and bracket $\alpha$ by rationals $l_n<\alpha<h_n$ of distance $e_n$ by bisection from $[9/10,1]$. Begin with $d=4$, $u=-1$, $b=1$ and both scales $9/10$. This seed is centrally symmetric with zero-free cutoff and has odd degree by v313.

| Motion | Cutoff events | Degree at end |
| --- | --- | ---: |
| Extend the right wall to $d=6$ | None | 1 |
| Increase both protrusion scales to 1 | Two regular crossings at reflected squares | 1 |
| Translate attachments to $u=1$, $b=-1+e_n$, keeping both scales 1 | None | 1 |
| Decrease the lower scale to $l_n$, keeping the upper scale 1 | One regular crossing at lower scale $\alpha$ | 0 |

The extension breaks global half-turn symmetry before the paired crossings occur. The two contact charts in the growing protrusions still correspond under the half-turn near their zeros. Each has degree one, so their union is a valid $V_n$. The limiting target has upper attachment $1$ at scale 1 and lower attachment $-1$ at scale $\alpha$. The linear edge registrations of the endpoints have error at most $2e_n$ by the vertex bounds. Each distant square has a vertex of absolute height at least $1+s$ since its vertical span is at least its side and its contacts are all on or beyond the original baseline. At that height the corresponding target protrusion is horizontally separated by at least $2-(7/16+17/32)-7/16=19/32>s$. The flat row is at least $s$ away and the side walls and opposite row are farther. Thus these wall squares have $d_\Gamma(Q)\ge2/5$ for every $n$. The theorem discards their paired contribution.

At the remaining crossing, the lower attachment differs from the target's by exactly $e_n$ and the scale equals $\alpha$. All four corners lie within $e_n$ of $\Gamma$. The theorem gives a square of side **exactly $2/5$** on this asymmetric target even though the endpoints all have degree zero above that cutoff. In this example, the explicit limit of the lower profile squares verifies the conclusion independently.

The rational replay goes one step further. On the lower profile of scale $l_n$ mark the four contacts corresponding to the normalized larger square. Choose as desired corners those of the same profile square at scale $h_n$ and the same attachment. Each displacement is at most $e_n$ because all normalized contacts lie within the unit sup box. The explicit PL repair preserves the supporting lines' directions near the contacts and moves the four marks to the desired square. The verified four-line matrix is invertible, so the repaired square has odd local degree and side **strictly greater** than $2/5$. At depth eight the displacement budget is $1/2560$ and the exact repaired squared side is $121272957/757596160>4/25$. This calculation does not assert a change of the repaired polygon's global degree.

## What must still be proved

For each arbitrary Jordan input, select a positive limiting cutoff and construct comparisons in which every unpaired wall configuration has all four corners in a shrinking target tube. The example establishes that shrinking errors are needed only after removing degree-zero regions. It supplies no general selection rule or universal lower bound in terms of inradius. A direct test for the next constructions is to look for unpaired configurations that stay a fixed positive distance from the curve, even after a proposed local pairing.

[Exact checker](../code/cutoff_payment_v315.py) · [Regression tests](../tests/test_cutoff_payment_v315.py) · [Rational output](../results/v315/cutoff-payment.json)
