# v375. The comb law

[v374](v374-asymmetric-wild-points.md) found that near a spiral wild point, small squares do not crowd around the point. They spread over a comb of nearly parallel turns, out to a radius much larger than their side. That was a heuristic. This note proves its upper half: every small square has a corner inside an explicit radius, and that radius has the predicted size.

The proof uses just two facts about the comb. First, two points that are not close neighbours along the same arm are at least one turn spacing apart. Second, a short stretch of arm is too straight to hold three corners of a square. So some side of the square must jump from one turn to another, and its length is at least the turn spacing there. The spacing grows with the radius, so a small square must reach in to where the spacing is small.

**Results.**
1. **Comb law (Theorem 375.1).** For two spiral arms $r=\theta^{-a}$ with angular gaps $\pi\pm\varphi(r)$, $|\varphi|\le\bar\varphi<\pi$, a square of side $t$ inscribed in the arc with all corners at radius at least $\rho$ has $t\ge d_g(\rho)$, where
$$d_g(\rho)=\frac{g\,a\,\rho^{1+1/a}}{1+\tfrac\pi2 a\rho^{1/a}},\qquad g=\pi-\bar\varphi,$$
under explicit smallness conditions.
2. **Confinement (Corollary 375.2).** Every inscribed square of side below $d_g(r)$ has a corner in $B(p,r)$. Equivalently, a square of side $t$ has a corner within $\rho(t)=d_g^{-1}(t)=(t/(ag))^{a/(1+a)}(1+o(1))$. This is the upper half of v374's Claim 374.4, now proved.
3. **Consequences (Corollary 375.3).** On the smooth arc, the threshold of Theorem 370.2 is set by the turn spacing, not by a polygon's edge length, and condition (a) of Theorem 374.3 holds as soon as $s\le d_g(\varepsilon/8)$. So the open path condition of v374 is pinned to an explicit radius, about $8\,(s/(ag))^{a/(1+a)}$. Carrying this over to inscribed polygons, which those theorems need, takes an error term that is sketched but not written out.
4. **Exact checks:**
   - on the continuous arcs, about 2.8 million sampled pairs of non-neighbours closer than $0.08$ are all at least $d_g(r_{\min})$ apart, with margin 7.5% (symmetric) and 12.7% (asymmetric);
   - on the eight v374 drawings, all 464 small squares have a corner inside $\rho(\text{side})$. The outermost square's nearest corner reaches about 89–90% of the bound.

## Setting

Put $p$ at the origin. Fix $a>0$, $\theta_0>0$, $R=\theta_0^{-a}$, and a $C^1$ function $\varphi$ on $(0,R]$ with
$$|\varphi|\le\bar\varphi<\pi,\qquad |r\varphi'(r)|\le L.$$
The arc $\alpha$ consists of $p$ and two arms, parametrized by the lifted angle $\Theta\ge\theta_0$:
$$\text{arm 1: }\Theta^{-a}e^{i\Theta},\qquad \text{arm 2: }\Theta^{-a}e^{i(\Theta+\pi+\varphi(\Theta^{-a}))}.$$
Arm 2 turns arm 1 by $\pi+\varphi$ at each radius, as in v374; $\varphi=0$ is the half-turn spiral of v370. Every point $z\neq p$ of $\alpha$ has a lifted angle $\Theta_z=|z|^{-1/a}$, determined by its radius alone.

Write $K(r)=a\,r^{1/a}$, the slope of a turn relative to the circle, and set
$$K_R=K(R),\quad \varepsilon_1=K_RL,\quad \tau=\arctan\frac{K_R}{1-\varepsilon_1},\quad \eta(t,\rho)=\Big(\frac{\pi}{\sqrt2}+\sqrt2L\Big)\frac t\rho.$$
The **smallness conditions at $(t,\rho)$** are
$$K_R\le\tfrac14,\qquad \varepsilon_1<1,\qquad \eta\le\tfrac\pi4,\qquad \eta(1+\varepsilon_1)+\tau<\tfrac\pi2.$$
The first two are conditions on $R$ only. The last two hold whenever $t/\rho$ is small. Note $|d\varphi/d\Theta|=K(r)\,|r\varphi'(r)|\le\varepsilon_1$.

Call two points $z,w\neq p$ **neighbours**, $z\sim w$, if they lie on the same arm and $|\Theta_z-\Theta_w|\le\pi/2$.

## Two lemmas

Throughout, $z,w\neq p$ are points of $\alpha$, $r=\min(|z|,|w|)$, and $\beta\in(-\pi,\pi]$ is the angle from $z$ to $w$ seen from $p$. Two elementary bounds are used repeatedly:
- **radial:** $|w-z|\ge\big||w|-|z|\big|=|\Theta_w^{-a}-\Theta_z^{-a}|\ge K(r)\,r\,|\Theta_w-\Theta_z|$, by the mean value theorem;
- **angular:** $|w-z|\ge2\sqrt{|w||z|}\,|\sin(\beta/2)|\ge\frac2\pi r|\beta|$.

**Lemma 375.A (spacing).** Suppose $|w-z|\le\sqrt2\,t$ and the smallness conditions hold at $(t,r)$. If $z\not\sim w$, then $|w-z|\ge d_g(r)$.

*Proof.* Write $\Delta=\Theta_w-\Theta_z$ and $K=K(r)$. Comparing polar angles,
$$\Delta\equiv\beta-(c_w-c_z)-(\varphi_w-\varphi_z)\pmod{2\pi},$$
where $c=0$, $\varphi_\cdot=0$ on arm 1 and $c=\pi$, $\varphi_\cdot=\varphi(|\cdot|)$ on arm 2.

*Different arms.* The constant difference is $\pm\pi$ and only one value of $\varphi$ appears. So $|\Delta|\ge\pi-\bar\varphi-|\beta|=g-|\beta|$.

*Same arm, $|\Delta|>\pi/2$.* The $\varphi$ difference $\delta$ is $0$ on arm 1. On arm 2, $|\delta|\le L|\log(|w|/|z|)|\le L|w-z|/r\le\sqrt2Lt/r\le\eta$. Write $\Delta=\beta-\delta+2\pi n$. If $n=0$, then $|\beta|>\pi/2-\eta\ge\pi/4$, so $|w-z|\ge2r\sin(\pi/8)>0.76\,r$. If $n\ne0$, then $|\Delta|\ge2\pi-|\beta|-\eta\ge g-|\beta|$.

In the cases with $|\Delta|\ge g-|\beta|$, put $x=|\beta|$. The two bounds give
$$|w-z|\ge\max\Big(Kr(g-x),\ \tfrac2\pi rx\Big)\ge\frac{gKr}{1+\pi K/2}=d_g(r),$$
since a decreasing and an increasing linear function have a maximum of at least their crossing value. In the remaining case, $0.76\,r\ge d_g(r)$ because $K\le K_R\le\frac14$ gives $d_g(r)\le\frac{\pi/4}{1+\pi/8}r<0.57\,r$. ∎

**Lemma 375.B (short stretches are straight).** Suppose the smallness conditions hold at $(t,\rho)$. Then no three corners of a square of side $t$, all at radius at least $\rho$, are pairwise neighbours.

*Proof.* For neighbouring corners, $|\Delta|\le\pi/2$, so $n=0$ above and $\beta=\Delta+\delta$. The angular bound gives $|\beta|\le\frac\pi2|w-z|/r\le\frac{\pi}{\sqrt2}t/\rho$, so $|\Delta|\le|\beta|+|\delta|\le\eta$. Three pairwise neighbours therefore lie on one arm, inside a $\Theta$-interval of length at most $\eta$.

On arm $j$ the point is $re^{i\Psi}$ with $r=\Theta^{-a}$ and $\Psi=\Theta+c_j+\varphi_j$. Its derivative is $e^{i\Psi}(-K r+ir\Psi')$, so the tangent direction is
$$\Psi+\frac\pi2+\arctan\frac{K(r)}{\Psi'},\qquad \Psi'=1+\frac{d\varphi_j}{d\Theta}\in[1-\varepsilon_1,1+\varepsilon_1].$$
Over the interval, $\Psi$ varies by at most $\eta(1+\varepsilon_1)$, and the arctangent stays in $(0,\tau]$. So all tangent directions on that stretch lie in an interval of length less than $\eta(1+\varepsilon_1)+\tau<\pi/2$. A chord is the integral of the unit tangent, so every chord of the stretch points into the same interval, and no two chords are perpendicular. But any three corners of a square include a vertex together with both of its neighbours, and those two sides are perpendicular. ∎

## The comb law

**Theorem 375.1.** Let $Q$ be a square of side $t$ with all four corners on $\alpha\setminus\{p\}$, all at radius at least $\rho$, and suppose the smallness conditions hold at $(t,\rho)$. Then $t\ge d_g(\rho)$.

*Proof.* On the corners of $Q$, $\sim$ is an equivalence relation. For transitivity, neighbouring corners have $|\Delta|\le\eta$ (proof of Lemma 375.B), so $A\sim B\sim C$ gives the same arm and $|\Theta_A-\Theta_C|\le2\eta\le\pi/2$. By Lemma 375.B each class has at most two corners. The four sides of $Q$ form a 4-cycle, and a class of two corners contains at most one side. So at least two sides of $Q$ join corners that are not neighbours. Corners are at most $\sqrt2t$ apart, and the conditions at $(t,\rho)$ imply them at $(t,r)$ for every $r\ge\rho$. So Lemma 375.A applies to such a side and gives $t\ge d_g(r_{\min})\ge d_g(\rho)$, since $d_g$ is increasing. ∎

**Corollary 375.2 (confinement).** Let $0<r\le R$, and suppose the smallness conditions hold at $(d_g(r),r)$. Then every square inscribed in $\alpha$ with side below $d_g(r)$ has a corner in $B(p,r)$. (A corner at $p$ counts.)

*Proof.* Otherwise all corners are at radius at least $r$. The conditions at $(t,r)$ follow from those at $(d_g(r),r)$, because $\eta$ increases with $t$. Theorem 375.1 then gives $t\ge d_g(r)$. ∎

Since $d_g(r)=ag\,r^{1+1/a}(1+O(r^{1/a}))$, the confinement radius for side $t$ is
$$\rho(t)=d_g^{-1}(t)=\Big(\frac t{ag}\Big)^{a/(1+a)}\big(1+O(t^{1/(1+a)})\big),$$
and $t/\rho(t)\to0$, so the conditions hold for every small $r$ once $K_R\le\frac14$ and $\varepsilon_1<1$. These two can always be arranged by shrinking $R$ (starting the arc at a larger $\theta_0$). For the half-turn spiral ($\varphi=0$, $L=0$), $g=\pi$ and $\eta=\frac\pi{\sqrt2}t/\rho$.

## What this settles

**Corollary 375.3.**
1. **Claim 374.4, upper half.** Near a two-arm spiral point, squares of side below $s$ reach out at most to $d_g^{-1}(s)$, which has the predicted size $(s/(ag))^{a/(1+a)}$. Here $g$ uses the global bound $\bar\varphi$; shrinking $R$ makes it local.
2. **Theorem 370.2 on spirals.** Theorem 370.2 needs the curve to be polygonal away from $p$, so it applies to polygons inscribed in $\alpha$, not to $\alpha$ itself. There, Lemma 369.6 gives a threshold of the order of the polygon's edge length, which shrinks as the polygon is refined. The comb law says that the right threshold is set by the turn spacing instead: on the smooth arc, squares of side below $d_g(R/8)$ have a corner in $B(p,R/8)$, however fine the approximation.
3. **Theorem 374.3, condition (a), arc part.** The locally symmetrized curve $\gamma_1$ is again a two-arm spiral arc, with $\varphi_1=(1-\chi)\varphi$. If the cutoff $\chi$ is chosen as a function of $\log r$, then $|r\varphi_1'|$ stays bounded. So on the smooth arc, (a) holds as soon as $s\le d_g(\varepsilon/8)$. The open condition (b) therefore concerns a path supported out to the explicit radius $\varepsilon\approx8\,(s/(ag))^{a/(1+a)}$.

**Polygonal arcs.** Items 2 and 3 use the theorem for polygons inscribed in $\alpha$, and the proof above is written for the smooth arc. It should carry over with an error term. Lemma 375.B only uses chord directions, and every chord of an inscribed polygon points into the tangent cone of the sub-arc it spans. Lemma 375.A loses at most twice the distance $h$ between the polygon and $\alpha$, so it needs $h$ small compared with the turn spacing. We have not written that version out. The exact checks below are on inscribed polygons with 24 vertices per turn, and they satisfy the smooth bound with no error term at all.

**What is not proved.** The matching lower bound, that squares of side near $s$ really do reach out to about $\rho(s)$, remains observational. On the depth-4 drawings at $s=1/20$, the outermost square's nearest corner is at $0.289$ against the bound $0.322$ (symmetric) and at $0.303$ against $0.341$ (asymmetric). Its centre is within 4% of v374's heuristic radius. That lower half is the part that makes the path condition hard.

## Exact checks

**Lemma 375.A on the continuous arcs.** Sample 7 turns of each arm at 6000 points. Among pairs closer than $0.08$ that are not neighbours (about 2.77 million for each arc), the smallest ratio $|w-z|/d_g(r_{\min})$ is $1.0754$ for $\varphi=0$ and $1.1268$ for $\varphi=0.4\sin(5\log(1/r))$. The 7.5% margin is roughly the factor $1+\frac\pi2K$ given up in the proof.

**Theorem 375.1 on the v374 drawings.** These are polygons inscribed in the arcs, so this is a consistency check, not a verification. Here $g=\pi-0.4\lambda$ and $L=2\lambda$, and the smallness conditions hold at $(s,\rho(s))$ for every row.

| $s$ | depth | $\lambda$ | squares | outside bound | worst min-corner / $\rho$ | $\rho(s)$ | $\eta$ at $s$ |
|---|---|---|---|---|---|---|---|
| 1/20 | 4 | 0 | 58 | 0 | 0.922 | 0.3255 | 0.34 |
| 1/20 | 4 | 1 | 64 | 0 | 0.905 | 0.3415 | 0.74 |
| 1/25 | 4 | 0 | 34 | 0 | 0.922 | 0.3010 | 0.30 |
| 1/25 | 4 | 1 | 42 | 0 | 0.905 | 0.3157 | 0.64 |
| 1/30 | 5 | 0 | 62 | 0 | 0.936 | 0.2825 | 0.26 |
| 1/30 | 5 | 1 | 68 | 0 | 0.944 | 0.2962 | 0.57 |
| 1/40 | 6 | 0 | 66 | 0 | 0.945 | 0.2558 | 0.22 |
| 1/40 | 6 | 1 | 70 | 0 | 0.972 | 0.2681 | 0.47 |

The worst ratios come from small squares at the truncated inner end of each drawing, where the polygon turns around. The theorem does not cover those squares, yet they still satisfy the bound.

## Where this leaves asymmetric points

v374 reduced a two-arm point to a path that symmetrizes it out to radius $\varepsilon$ without passing a square of side $s$. This note fixes $\varepsilon$: about $8\,(s/(ag))^{a/(1+a)}$, with a proof. On that scale the path rotates whole turns of the comb. The next step is the lower half of the comb law together with a count of the side-$s$ squares that such a rotation creates and destroys. If they come in pairs, or can be avoided by choosing $s$ between the side lengths the comb produces, the path condition follows.

[Code](../code/comb_law_v375.py) · [Tests](../tests/test_comb_law_v375.py) · [Data](../results/v375/comb_bound_check.json)
