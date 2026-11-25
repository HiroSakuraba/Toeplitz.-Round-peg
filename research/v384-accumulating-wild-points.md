# v384. Wild points that pile up

So far every wild curve in this series had one wild point, or a few. This note takes curves with infinitely many spiral points $p_k$ piling up at a limit point $p_\infty$. Two questions are asked.
- At a fixed cutoff $s$, does the large graceful parity $L(\gamma,s)$ even make sense, and how many of the points does it see?
- As $s\to0$, does the parity settle, or flip forever?

**The short answer.** At each cutoff only finitely many points matter: those bigger than about $s$. When those are well separated, they act independently, and the parity is a sum of single-point contributions. So, in the sparse regime, piling up adds nothing new: the question reduces to one wild point at a time, read at its own relative cutoff. In the self-similar regime, the few points at the cutoff scale are too close for the separation argument, so they could interact. In the model they don't: exact counts at five scales show no interaction term.

**Results.**
1. **The invisibility radius, pinned (Lemma 384.1).** At any point $p$ of a curve with no square of side $s$, the largest ball inside which arbitrary changes leave $L(\cdot,s)$ alone has radius between $\min(s/2,d_N)$ and $\min(s/\sqrt2,d_N)$. Here $d_N$ is the distance from $p$ to the nearest free corner of a near-square (side $s$, three corners on the curve). v374 had only the upper bound.
2. **Flattening in a strip (Lemma 384.2).** Wild points whose disks lie in a strip of half-width $\eta<s/(2\sqrt2)$ around a line can all be flattened at once, infinitely many included, provided no other part of the curve outside the strip comes within $\sqrt2\,s$ of them. The reason: three corners of a side-$s$ square never fit in such a strip.
3. **Finite range (Theorem 384.3).** For any curve with wild points $p_k\to p_\infty$, and any cutoff $s$ at which it has no square of side $s$, flattening all points of size below some $\delta(s)>0$ leaves $L(\gamma,s)$ unchanged. So the parity sees finitely many points. On the test models, Lemma 384.2 gives $\delta(s)$ explicitly, about $0.35\,s$.
4. **Additivity (Theorem 384.4).** If the points that remain are more than $\sqrt2\,s$ apart, then
$$L(\gamma,s)-1\equiv\sum_k\big(L(\gamma^{(k)},s)-1\big)\pmod 2,$$
   where $\gamma^{(k)}$ keeps point $k$ alone. A cluster of interacting points can be kept together as one unit.
5. **Symmetric points stay protected (Corollary 384.5).** If every point is half-turn symmetric and the separation holds at a cutoff $s$, then $L(\gamma,s)$ is odd. On the sparse model this holds for every cutoff below $0.0304$. So that curve inscribes a square of side at least $0.0304$, although it is polygonal away from infinitely many wild points. This extends Theorem 370.2, which allowed one.
6. **A recursion for the sparse model (Corollary 384.6).** With sizes $\rho_k=4^{-k}$,
$$L(\gamma,s/4)\equiv L(\gamma,s)+\delta\big(s/(4\rho_{k_0})\big),$$
   where $\delta(r)$ is the parity defect of one spiral at relative cutoff $r$. So the parity flips forever exactly when the largest spiral's defect does, along a geometric sequence of cutoffs.
7. **Checks** (exact counts on rational Jordan polygons with up to six spirals):
   - **Sparse, twisted and symmetric, four cutoffs each:**
     - the full drawing, the drawing with the strip-flattenable spirals removed, and the additive prediction agree every time;
     - the per-spiral counts of small aligned and reversed squares equal those of each spiral alone, exactly;
     - the twisted parity flips (0, 0, 1, 0 at $s=1/60,1/120,1/240,1/480$), as the recursion predicts;
     - the symmetric parity is always 1.
   - **Self-similar, five cutoffs covering one period in $\log_2 s$, twisted and symmetric:** the strip lemma cannot clear the one or two spirals at the cutoff scale. Yet exact counts find no interaction. No square reaches across two spirals, the cluster's parity equals the sum of its members', and each member alone leaves $L$ odd.

**A bug fix in the square counter.** The exact counter of v370 skipped candidate squares whose $4\times4$ contact matrix had absolute determinant below $10^{-14}$. That matrix is built from edge vectors, so on edges shorter than about $10^{-4}$ every candidate was skipped. Tiny spirals therefore looked as if they flipped the parity. The test is now relative to the edge lengths. All 129 tests of the earlier notes that use the counter pass unchanged, so no earlier result was affected; their drawings never had edges that short.

## Setting

$L(\gamma,s)$ is the mod-2 degree of the large aligned plus large reversed squares (side above $s$), defined whenever $\gamma$ has no square of side exactly $s$ (v371, v373). It is constant along any uniformly continuous path of closed curves, simple or not, none of which has a square of side exactly $s$. This is the path principle, Corollary 371.2. For a Jordan polygon it equals $1+(\text{small aligned})+(\text{small reversed})$, by Lemma 370.1.

## The invisibility radius

**Lemma 384.1.** Let $\gamma$ have no square of side $s$, let $p\in\gamma$, and let $N_s$ be the set of free corners $A$ of squares of side $s$ with the other three corners on $\gamma$. Let $d_N={\rm dist}(p,N_s)$. Then:
- every closed curve inside $\gamma\cup B(p,\rho)$ with $\rho<\min(s/2,d_N)$ has no square of side $s$;
- for $\rho>\min(s/\sqrt2,d_N)$ some closed curve inside $\gamma\cup B(p,\rho)$ has one.

*Proof.* Take $\rho<\min(s/2,d_N)$ and a square of side $s$ on such a curve. Two of its corners are at least $s>2\rho$ apart, so at most one lies in the ball off $\gamma$. If none does, it is a square of $\gamma$, which is excluded. If one does, it is a free corner within $\rho<d_N$ of $p$, which is also excluded. The upper bound is v374's: a detour through a nearby free corner, or the circle of radius $s/\sqrt2$. ∎

$N_s$ is closed, by compactness, and does not contain $p$, since that would be a square of $\gamma$. So $d_N>0$. In the regime $\rho<s/2$ the lemma says exactly which local changes are invisible: those that stay away from free corners of near-squares.

## Flattening in a strip

**Lemma 384.2.** Let $\gamma$ have no square of side $s$. Let $S$ be the strip of half-width $\eta<s/(2\sqrt2)$ around a line, and let $B_j$ be disks (finitely or infinitely many) contained in $S$, such that every point of $\gamma\setminus S$ is more than $\sqrt2\,s$ from every $B_j$. Then no closed curve inside $\gamma\cup\bigcup_jB_j$ has a square of side $s$. Hence, if each $\gamma\cap B_j$ is an arc whose diameters tend to 0, replacing every such arc by its chord leaves $L(\cdot,s)$ unchanged.

*Proof.* Any three corners of a side-$s$ square form a right isosceles triangle with legs $s$. Its smallest width is $s/\sqrt2>2\eta$, so at most two corners lie in $S$. Take a square with a corner off $\gamma$. That corner lies in some $B_j\subset S$. Every other corner is within $\sqrt2\,s$ of it, so it lies in another disk (inside $S$) or on $\gamma$ within $\sqrt2\,s$ of $B_j$, which is also inside $S$. So all four corners lie in $S$, which is impossible. A square with every corner on $\gamma$ is excluded by hypothesis. The straight homotopy of each arc to its chord stays in its disk. Because the diameters tend to 0, the homotopy is uniformly continuous. The path principle finishes the proof. ∎

## Finite range

**Theorem 384.3.** Let $\gamma$ be a Jordan curve with no square of side $s$. Let $p_k\to p_\infty$ be points of $\gamma$, each with a disk $B_k=B(p_k,\rho_k)$ meeting $\gamma$ in one arc, with $\rho_k\to0$. Then there is $\delta(s)>0$ such that flattening every $B_k$ with $\rho_k<\delta(s)$ leaves $L(\gamma,s)$ unchanged. Only finitely many points have $\rho_k\ge\delta(s)$.

*Proof.* Suppose not. Then there are $\delta_n\to0$ and squares of side $s$ with corners in $\gamma\cup\bigcup_{\rho_k<\delta_n}B_k$. The disks involved have radius below $\delta_n$ and centres on $\gamma$. Pass to a subsequence of convergent corners. The limit corners lie on $\gamma$, which is closed. So the limit is a square of side $s$ on $\gamma$, a contradiction. Once no closed curve in that set has a square of side $s$, the straight flattening homotopy keeps $L$, by the path principle. ∎

So at each cutoff the parity is a statement about a curve with finitely many wild points. Compactness does not say how small $\delta(s)$ is. On the models below, Lemma 384.2 makes it explicit.

## Additivity

**Theorem 384.4.** Suppose the curve is polygonal away from finitely many wild points $p_1,\dots,p_m$ (after flattening), plus straight pieces along a line. Suppose the disks $B_k$ are more than $\sqrt2\,s$ apart and more than $\sqrt2\,s$ from every non-straight part of the polygonal remainder. Suppose also that the remainder has no square of side below $s$. Then, with $\gamma^{(k)}$ the curve with only $p_k$ kept,
$$L(\gamma,s)-1\equiv\sum_{k=1}^m\big(L(\gamma^{(k)},s)-1\big)\pmod2,$$
whenever $\gamma$ and all $\gamma^{(k)}$ have no square of side $s$.

*Proof.*
- Truncate each spiral at depth $n$ by a chord through its centre, giving polygons $D_n$ and $D_n^{(k)}$. For large $n$ these have the same $L$ as $\gamma$ and $\gamma^{(k)}$, by Lemma 384.1 applied at the finitely many centres. (Equivalently, by compactness, as in Lemma 374.2.)
- In $D_n$, a square of side below $s$ has diameter below $\sqrt2\,s$, and it cannot have all four corners on straight pieces. So it has a corner in exactly one disk $B_k$, and all its corners lie within $\sqrt2\,s$ of $B_k$. There $D_n$ and $D_n^{(k)}$ coincide and are traversed in the same order. So the small squares of $D_n$ split into the small squares of the $D_n^{(k)}$, sector by sector.
- Lemma 370.1 turns this into the identity. ∎

A group of points that are closer than $\sqrt2\,s$ can be kept together as one unit $\gamma^{(K)}$. The proof is the same.

**Corollary 384.5 (symmetric points).** If, in addition, each $\gamma\cap B(p_k,R_k)$ is invariant under the half-turn about $p_k$, with $R_k$ beyond $\rho_k+\sqrt2\,s$, then $L(\gamma,s)$ is odd.

*Proof.* The small squares of $D_n$ near $p_k$ lie inside $B(p_k,R_k)$, where the curve is half-turn invariant. As in Theorem 370.2, the half-turn pairs small aligned squares with small reversed ones there. So the small graceful count is even, and $L(D_n,s)=1$. ∎

## The test models

The base is the rectangle $[-1,1]\times[0,1]$, traversed counterclockwise. Its bottom edge runs through $p_\infty=0$. At $x_k=2^{-k}$, for $k\ge2$, a double spiral of radius $\rho_k$ is spliced into the edge:
- the edge enters the disk on the left;
- arm one, $r=\rho_k(\theta/3\pi)^{-1/2}$, winds in from angle $\pi$;
- a chord crosses the centre;
- arm two winds out to the right.

When arm two is the half-turn image of arm one, the whole local arc, edge pieces included, is invariant under the half-turn about $p_k$. The twisted version turns arm two by $\pi+\phi(r/\rho_k)$ instead, with $\phi(u)=0.4\sin(5\log(1/u))$ as in v374, so it is asymmetric at every scale. Twist and depth are scale-free, so spiral $k+1$ is an exact scaled copy of spiral $k$, up to rounding.

- **Sparse:** $\rho_k=4^{-k}$. Each point is isolated at its own scale.
- **Self-similar:** $\rho_k=2^{-k}/8$. The picture is the same at every scale.

The sparse drawings keep spirals $2\le k\le7$, kept ones to depth 3 half-turns and flattenable ones to depth 2. The self-similar drawings keep $2\le k\le8$, all to depth 2. Every drawing has 12 vertices per turn, and all are exact rational Jordan polygons.

**Sparse hypotheses for every small cutoff.** Write $\eta=0.999\,s/(2\sqrt2)$. Spirals with $\rho_k>\eta$ are kept, and the rest are flattened by Lemma 384.2. Let $b$ be the smallest kept spiral.
- If $b\ge4$, the gap from spiral $b$ to the next disk in is $2^{-b-1}-\tfrac54 4^{-b}$. Since $s<4^{-b}/(0.999/(2\sqrt2))$, so that $\sqrt2\,s<4.004\cdot4^{-b}$, this exceeds $\sqrt2\,s$ because $2^b\ge16>10.51$.
- Kept neighbours $k<k+1\le b$ with $k\ge3$ are separated by the same comparison: their gap $2^{-k-1}-\tfrac54 4^{-k}$ beats $\sqrt2\,s<1.001\cdot4^{-k}$ once $2^k>4.5$.
- If $b=3$, the gaps are $0.0430$ (spiral 3 to spiral 4) and $0.0469$ (spiral 2 to spiral 3). These exceed $\sqrt2\,s$ for $s<0.0304$.

So Theorem 384.4 applies at every cutoff below $0.0304$ (also checked on 3 000 cutoffs from $10^{-7}$). With symmetric spirals, Corollary 384.5 gives $L(\gamma,s)=1$ for every such $s$ at which the curve has no square of side $s$. A limit of squares of side above $s_n\uparrow0.0304$ then gives a square of side at least $0.0304$.

**The recursion (Corollary 384.6).** In the sparse model let $\delta(r)=1+L(\gamma_1,r\rho)$ be the defect of a single spiral of radius $\rho$ on a straight edge, at relative cutoff $r$. This is scale-free. Additivity gives $L(\gamma,s)\equiv1+\sum_k\delta(s/\rho_k)$. Only finitely many terms are nonzero: $\delta(r)=0$ once the spiral fits in the strip, and the list of spirals starts at $k_0$. Replacing $s$ by $s/4$ shifts every relative cutoff to the next spiral and adds one term, the largest spiral at its new, smaller relative cutoff. That is the recursion in the summary. Iterating it,
$$L(\gamma,4^{-m}s)\equiv L(\gamma,s)+\sum_{i=1}^m\delta\big(4^{-i}s/\rho_{k_0}\big).$$
So the parity of the whole accumulating curve, along a geometric sequence of cutoffs, is a running sum of one spiral's defects read deeper and deeper. For symmetric spirals every defect is 0 (v370). For asymmetric ones the defect sequence is exactly the single-point question of v374 and v377.

## Checks

**Sparse model.** Spirals $k=2,\dots,7$. $L$ is computed exactly on the full drawing, on the drawing with only the kept spirals, and on each kept spiral alone. "Prediction" is $1+\sum_k(L_k-1)$.

| cutoff | kept | flattened (Lemma 384.2) | twisted: full / kept only / prediction | singles $L_k$ | symmetric: full / prediction |
|---|---|---|---|---|---|
| 1/60 | 2, 3 | 4–7 | 0 / 0 / 0 | 1, 0 | 1 / 1 |
| 1/120 | 2, 3, 4 | 5–7 | 0 / 0 / 0 | 1, 0, 1 | 1 / 1 |
| 1/240 | 2, 3, 4 | 5–7 | 1 / 1 / 1 | 0, 1, 0 | 1 / 1 |
| 1/480 | 2–5 | 6, 7 | 0 / 0 / 0 | 1, 1, 0, 1 | 1 / 1 |

- **Same squares in the full drawing and alone.** In all eight runs, every kept spiral's small aligned and reversed counts in the full drawing equal its counts alone. No square has corners near two spirals.
- **Flattened spirals.** These are drawn to depth 2. Each carries 30 aligned and 36 reversed squares (twisted), or 16 and 16 (symmetric), all smaller than $s$. That is an even graceful number, as Lemma 384.2 requires.
- **Scale invariance.** Spiral $k$ at cutoff $s$ has exactly the counts of spiral $k+1$ at $s/4$. For example, 44 aligned and 41 reversed for spiral 3 at $1/60$ and for spiral 4 at $1/240$.
- **Recursion.** $L(1/240)=L(1/60)+\delta=0+1$ and $L(1/480)=L(1/120)+\delta=0+0$, with $\delta$ the defect of spiral 2 alone at the smaller cutoff.
- **Symmetric drawings.** Aligned and reversed counts are equal spiral by spiral.

**Self-similar model.** Spirals $k=2,\dots,8$, depth 2. Cutoffs $1/64,1/76,1/91,1/108,1/128$ span one period of $\log_2 s$. "Big" spirals are separated by more than $\sqrt2\,s$ from everything smaller, so Theorem 384.4 applies to them. The rest form the cluster at $p_\infty$. In the cluster, Lemma 384.2 clears the smallest spirals, and the "hard" ones remain.

| cutoff | big | cluster | cleared by Lemma 384.2 | hard | hard together: $L$ | prediction from singles | squares across two spirals |
|---|---|---|---|---|---|---|---|
| 1/64 | 2, 3 | 4–8 | 6–8 | 4, 5 | 1 (both twists) | 1 | 0 |
| 1/76 | 2–4 | 5–8 | 5–8 | none | 1 | 1 | – |
| 1/91 | 2–4 | 5–8 | 7, 8 | 5, 6 | 1 | 1 | 0 |
| 1/108 | 2–4 | 5–8 | 7, 8 | 5, 6 | 1 | 1 | 0 |
| 1/128 | 2–4 | 5–8 | 7, 8 | 5, 6 | 1 | 1 | 0 |

- At $1/64$ the full symmetric drawing (all seven spirals) was also counted: $L=1$, equal to the prediction, with no square across two spirals.
- The hard spirals have radius between $0.35\,s$ and $0.7\,s$. All their squares are smaller than $s$, with an even graceful count (30+36 twisted, 16+16 symmetric). So they act as if flattened, although Lemma 384.2 cannot show it.
- Lemma 384.1 would show it if no near-square of side $s$ has its free corner in their disks, which is what the counts suggest.

## What this changes

- **Accumulation is not a new obstacle in the sparse regime.** At each cutoff the parity is a finite sum of single-point defects. Protected points stay protected when infinitely many of them pile up, as long as they separate at every scale.
- **The self-similar regime is where interaction could enter.** The separation argument fails at the cutoff scale in every period, so a log-periodic interaction term is possible in principle. In the model it is zero at every scale tested, for symmetric and twisted spirals alike.
- **Open:** the single-point defect $\delta(r)$ of an asymmetric spiral as $r\to0$, which is v374's path condition in another guise; and the interaction term in general, beyond the model.

[Code](../code/accumulating_v384.py) · [Tests](../tests/test_accumulating_v384.py) · [Data](../results/v384/trials.json)
