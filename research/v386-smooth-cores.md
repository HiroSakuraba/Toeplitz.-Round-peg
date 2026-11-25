# v386. Smooth cores: what regularity buys, and what the twist takes back

[v385](v385-asymmetry-below-the-invisibility-radius.md) left the hard case of the path condition: points twisted at every scale. There the symmetrization must reach the confinement radius, which for polygonal arms is the comb radius $\sim s^{a/(1+a)}$. [v379](v379-squares-on-a-smooth-spiral.md) observed that smooth arms confine small squares much more tightly, to $r_{\rm smooth}(s)\sim s^{2a/(2a+1)}$. It suggested that extending Theorem 370.2 to smooth curves would shrink the radius the path condition must reach. This note does that extension and makes the smooth confinement rigorous. It then finds the catch: a twist brings comb-scale squares back.

The tool is again the British flag identity of v383. Written as a mixed second difference, it becomes an exact version of v379's selection rule, with no error terms and with the twist included. From it:
- a symmetric smooth core confines squares as tightly as v379 predicted;
- a twisted core does too for graceful squares, unless the twist varies quickly;
- squares that escape to the comb scale through the twist's size are always crossed, so they cannot change the parity.

**Results.**
1. **Exact selection identity (Theorem 386.1).** For any square on the half-turn spiral, with corners on passes $k_A,\dots,k_D$,
$$\pi j+\Delta\theta=(2a+1)\,\frac{\zeta^{2a+1}}{\xi^{2a+2}}\,\delta_B\delta_D\qquad\text{exactly},$$
   where:
   - $j=k_A+k_C-k_B-k_D$;
   - $\Delta\theta$ is the "angular flag" $\theta_A+\theta_C-\theta_B-\theta_D$, with $|\Delta\theta|\le2\arcsin(t^2/2|Z|^2)$;
   - $\delta_X=\Theta_X-\Theta_A$ are lift differences;
   - $\xi,\zeta$ lie in an explicit interval around the lifts.

   With a twist, $\pi j$ becomes $\pi j-\Sigma^*$, where $\Sigma^*$ is the alternating sum of the twist over the corners on arm two.
2. **Smooth confinement, rigorous (Theorem 386.2).** On the plain half-turn spiral, in v383's range, every square of side $t$ with $t\le0.2\,|Z|$ satisfies one of two explicit inequalities. Both force
$$t\ \ge\ c_2\,r^{1+1/(2a)}\,(1-o(1)),\qquad c_2=a\sqrt{2\pi/(2a+1)},$$
   where $r$ is its smallest corner radius. That is, every corner set reaches within $r_{\rm smooth}(t)=(t/c_2)^{2a/(2a+1)}$ of the point. This is v379's Corollary 379.4, now a theorem. Squares with $t>0.2|Z|$ have a corner within $6t\ll r_{\rm smooth}(t)$.
3. **Smooth cores in Theorem 370.2 (Lemma 386.3, Theorem 386.4).** Every $C^2$ Jordan curve has sector degrees $(1,0)$ below its smallest square, and $L=1+$(small graceful) as for polygons. So Theorems 370.2 and 374.3 hold with "polygonal away from $p$" replaced by "$C^2$ away from $p$". Their confinement condition (a) needs only the **graceful** small squares inside the symmetric ball. Crossed ones never enter the parity.
4. **What the twist does (Proposition 386.5).** On twisted smooth arms there are two cases:
   - **Crossed squares can escape through the twist's size.** These have one diagonal pair of corners on each arm. Their effective level is $2\varphi/\pi$, so near a zero of $\varphi$ they come out to the comb scale.
   - **Graceful squares can escape only through the twist's gradient.** If arm two holds one or three corners, $j$ is odd and $|j_{\rm eff}|\ge1-\Phi/\pi$. Otherwise the twist enters only through differences, bounded by $|r\varphi'|\,\log(R/r)$.

   So graceful squares are confined to about $\max\big(r_{\rm smooth}(t),\,(L\,t/c_L)^{a/(1+a)}\big)$, with $L=\sup|r\varphi'|$.
5. **Checks:**
   - all 807 labelled squares of v382 satisfy the identity to $4\times10^{-17}$. All 428 with $t\le0.2|Z|$ have the implied ratio inside its interval and satisfy the confinement inequality. The worst has $t/(c_2r^2)=1.061$, so the constant is nearly sharp.
   - on twisted arms ($\varphi=0.4\sin(5\log(1/r))$), Newton finds 9 squares far below the plain bound, down to $t/(c_2r^2)=0.20$, at side $1.46$ turn spacings near the twist zero $r=e^{-4\pi/5}$. All 9 are of the predicted kind ($j=0$, diagonal pair on arm two) and all are crossed.

## The exact identity

Use the setting of v383: passes $\rho=(\theta+\pi k)^{-a}$, lifts $\Theta=\theta+\pi k$, and $F(\Theta)=\Theta^{-2a}$, the squared radius. Choose a corner $A$ whose two neighbours lie on other passes. In v383's range every small square has one, since its pattern is 1-1-1-1 or 2-1-1 diagonal. Then $|\delta_B|,|\delta_D|\ge\pi-W$.

**Theorem 386.1.** Put $e=\pi j+\Delta\theta$, so that $\Theta_C=\Theta_A+\delta_B+\delta_D+e$. The flag identity $F(\Theta_A)+F(\Theta_C)=F(\Theta_B)+F(\Theta_D)$ reads
$$\underbrace{F(\Theta_A+\delta_B+\delta_D)-F(\Theta_A+\delta_B)-F(\Theta_A+\delta_D)+F(\Theta_A)}_{=\ \delta_B\delta_D\,F''(\xi)}+\underbrace{F(\Theta_C)-F(\Theta_C-e)}_{=\ e\,F'(\zeta)}=0.$$
With $F''/(-F')=(2a+1)/\Theta$ at the respective points, this is the identity of the summary. Here $\xi$ lies in the hull of $\Theta_A,\Theta_B,\Theta_D$ and $\Theta_B+\Theta_D-\Theta_A$, and $\zeta$ between $\Theta_C$ and $\Theta_B+\Theta_D-\Theta_A$. All of these lie in $[2\Theta_{\rm lo}-\Theta_{\rm hi},\,2\Theta_{\rm hi}-\Theta_{\rm lo}]$, where $\Theta_{\rm lo},\Theta_{\rm hi}$ are the smallest and largest corner lifts.

For the angular flag, write the corners as $Z+v_X$ with $v_C=-v_A$, $v_D=-v_B$, $v_B=iv_A$, and $u=v_A/Z$. Then $\theta_A+\theta_C-2\arg Z=\arg(1-u^2)$ and $\theta_B+\theta_D-2\arg Z=\arg(1+u^2)$. So $|\Delta\theta|\le2\arcsin|u|^2=2\arcsin\frac{t^2}{2|Z|^2}$. ∎

v379's Theorem 379.2 is the leading-order form of this: $\xi\approx\zeta\approx\Theta$, $\delta\approx\pi n$, $\Delta\theta\approx0$. Here nothing is approximated.

**With a twist.** A corner on arm two (odd pass) has lift $\Theta_X=\theta_X+\pi k_X-\varphi(\rho_X)$. The same computation gives $\pi j-\Sigma^*+\Delta\theta=(2a+1)\zeta^{2a+1}\xi^{-2a-2}\delta_B\delta_D$, with
$$\Sigma^*=\varphi^*_A+\varphi^*_C-\varphi^*_B-\varphi^*_D,\qquad\varphi^*_X=[k_X\text{ odd}]\,\varphi(\rho_X).$$
Call $j_{\rm eff}=j-\Sigma^*/\pi$ the effective level.

## Smooth confinement

**Theorem 386.2.** Let a square of side $t$ lie on the plain half-turn spiral in v383's range, with $t\le\tau|Z|$ and $2\Theta_{\rm lo}>\Theta_{\rm hi}$. (For $a=1/2$ the second holds when $t<0.24|Z|$; take $\tau=0.2$.) Write $r$, $R$ for its smallest and largest corner radius, and
$$\Lambda_1=2\Theta_{\rm lo}-\Theta_{\rm hi},\qquad\Lambda_2=2\Theta_{\rm hi}-\Theta_{\rm lo}.$$
Then either:
- **$j=0$, and** $2\arcsin\frac{t^2}{2|Z|^2}\ \ge\ (2a+1)(\pi-W)^2\,\Lambda_1^{2a+1}/\Lambda_2^{2a+2}$; or
- **$j\ne0$, and**
$$\pi-2\arcsin\frac{t^2}{2|Z|^2}\ \le\ (2a+1)\frac{\Lambda_2^{2a+1}}{\Lambda_1^{2a+2}}\cdot\frac{t^2}{2a^2}\,r^{-2/a-2}\Big(1+\frac t{\sqrt2r}\Big)^2 .$$

*Proof.* The case $j=0$ is the identity with $|\delta_B\delta_D|\ge(\pi-W)^2$.

For $j\ne0$, $|e|\ge\pi-|\Delta\theta|$. Bound $|\delta_X|\le\frac1ar^{-1/a-1}|\rho_X-\rho_A|$ by the mean value theorem. Then use
$$|\rho_B-\rho_A|\,|\rho_D-\rho_A|\le\big(|p\cdot\hat A|+\tfrac{t^2}{2r}\big)\big(|q\cdot\hat A|+\tfrac{t^2}{2r}\big)\le\tfrac{t^2}2\Big(1+\tfrac t{\sqrt2r}\Big)^2,$$
which holds because $p\perp q$, $|p|=|q|=t$, and $(p\cdot\hat A)^2+(q\cdot\hat A)^2=t^2$. ∎

To leading order, the first case needs $t\gtrsim\pi\sqrt{2a+1}\,r^{1+1/(2a)}$, and the second $t\gtrsim c_2\,r^{1+1/(2a)}$. Since $c_2<\pi\sqrt{2a+1}$, every square in range has $t\ge c_2r^{1+1/(2a)}(1-o(1))$.

Squares outside the range are of two kinds:
- those with $t>0.2|Z|$ have $|Z|<5t$, so a corner lies within $6t$;
- those with larger windows or a smaller outer lift are covered the same way, or are not small.

For $a=1/2$, $c_2=\sqrt\pi/2\approx0.886$ and $r_{\rm smooth}(t)=(t/c_2)^{1/2}$, against the comb radius $(2t/\pi)^{1/3}$.

## Smooth cores in the parity theorems

**Lemma 386.3.** Let $C$ be a $C^2$ Jordan curve with curvature at most $\kappa$ and reach $r_0>0$. Then:
- $C$ has no inscribed square of side below some $t_0(\kappa,r_0)>0$;
- for every $s$ at which $C$ has no square of side $s$, $L(C,s)=1+$(number of graceful squares of side below $s$) mod 2, and the aligned and reversed large degrees are those of a polygon.

*Proof.* Take a square of side $t<t_0$.
- Its corners are within $\sqrt2t$ of each other. By the reach, they lie on one subarc of length $O(t)$, whose total turning is $O(\kappa t)<\pi/4$ for small $t$.
- Cut the curve away from that subarc to get an arc. By v380's general comb law, an inscribed square on an arc has two sides whose connecting subarcs each turn by at least $\pi/4$. That is impossible here.

The same holds uniformly for inscribed polygons $P_m$ of $C$ with mesh tending to 0:
- the polygonal subarc between two vertices turns by at most the smooth turning over a slightly longer arc;
- the reach transfers to $P_m$ for large $m$.

So $P_m$ has no squares of side below $t_0/2$, and Lemma 370.1 gives $L(P_m,t)=1$ for $t<t_0/2$. Since $C$ has no squares of side $t$, Proposition 371.1 transfers this to $C$ for large $m$. The cutoff version of the path principle then moves $t$ up to $s$, flipping $L$ at each graceful square's side. ∎

**Theorem 386.4.** Theorem 370.2 holds for a curve that is $C^2$ away from $p$, if:
- the truncations $D_n$ are made with $C^2$ half-turn-symmetric caps, so each $D_n$ is a $C^2$ Jordan curve;
- the threshold condition is required of graceful squares only: every graceful square of side below $s_0$ on $\gamma$ or on any $D_n$ has a corner in $B(p,R/8)$.

The same goes for Theorem 374.3 and its condition (a).

*Proof.* As in Theorem 370.2, with Lemma 386.3 in place of Lemma 370.1. The half-turn pairs small aligned squares with small reversed ones. That needs only the graceful small squares to lie in the symmetric ball, and crossed squares do not enter $L$. ∎

For a symmetric smooth core, Theorem 386.2 makes condition (a) hold at radius $\varepsilon\approx8\,r_{\rm smooth}(s)$, instead of v375's $8\,\rho(s)$. That improvement is the one v379 suggested.

## What the twist takes back

**Proposition 386.5.** On twisted smooth arms, $\Sigma^*$ depends only on which corners lie on arm two:
- **none:** $\Sigma^*=0$, as on the plain spiral;
- **a diagonal pair** ($\{A,C\}$ or $\{B,D\}$): $\Sigma^*\approx\pm2\varphi$, and $j$ is even. With $j=0$, the identity allows $|\delta_B\delta_D|$ as small as $2\Theta|\varphi|/(2a+1)$. So such squares exist with side $t\approx\sqrt{2|\varphi|/\pi}\;c_2r^{1+1/(2a)}$, and near zeros of $\varphi$ down to $\sqrt2$ turn spacings. **These squares are crossed:** along the curve, arm one comes first and arm two second, so a diagonal pair on one arm is adjacent in curve order but not in the square.
- **one or three corners:** $j$ is odd, and $|\pi j_{\rm eff}|\ge\pi-\Phi-|\text{difference}|$. Smooth confinement holds with $c_2$ replaced by $c_2\sqrt{1-\Phi/\pi}$ to leading order.
- **an adjacent pair, or all four:** $\Sigma^*$ is a difference (or second difference) of $\varphi$ over corners at most $t\sqrt2$ apart, so $|\Sigma^*|\le2L\log(R/r)$ with $L=\sup|r\varphi'|$. If $j=0$, the identity gives
$$(2a+1)(\pi-W)^2\Lambda_1^{2a+1}\Lambda_2^{-2a-2}\le2L\log\tfrac Rr+2\arcsin\tfrac{t^2}{2|Z|^2},$$
   hence $t\gtrsim\min\big(c_L\,r^{1+1/a}/L,\ \pi\sqrt{2a+1}\,r^{1+1/(2a)}\big)$, with $c_L=(2a+1)\pi^2/(2\sqrt2)$.

Graceful squares have arm-one corners contiguous in the square's cyclic order. So they fall in the first, third or fourth case, and their confinement radius is about $\max\big(r_{\rm smooth}(t)(1-\Phi/\pi)^{-a/(2a+1)},\ (Lt/c_L)^{a/(1+a)}\big)$. ∎

So for the parity, the twist's size is harmless: everything it lets escape is crossed. Its gradient is not. A twist varying by $L$ per unit of $\log r$ lets graceful squares reach a fraction $(L\,ag/c_L)^{a/(1+a)}$ of the comb radius. For v374's profile $L=2$, which gives about $0.77$ of the comb radius for $a=1/2$, so essentially no gain. A slowly varying twist does gain: $L=0.1$ gives $0.28$.

## Checks

**Plain spiral.** v382's 807 labelled squares, rebuilt from their three-corner branches, and the 40 squares of a v380 random search:

| set | squares | with $t\le0.2|Z|$ | ratio inside interval | inequality holds | min $t/(c_2r^2)$ |
|---|---|---|---|---|---|
| v382 branches | 807 | 428 | 428 | 428 (all $j\ne0$) | 1.061 |
| v380 search | 40 | 7 | 7 | 7 | 1.231 |

- The flag identity holds to $4\times10^{-17}$.
- The 33 search squares with larger $t/|Z|$, up to 0.50, are near the point, as the proof requires. No square with $t\le0.2|Z|$ has $j=0$, as v383 and the first alternative predict.

**Twisted arms.** These use $\varphi(r)=0.4\sin(5\log(1/r))$, the default twist of v374, with Newton's method on corners held to prescribed passes. Seeds are the straight-comb guesses for labels with $n_1n_3\ne0$:
- at 120 radii in $[0.06,0.3]$;
- densely around the zeros of $\varphi$ at $r=e^{-2\pi/5},e^{-3\pi/5},e^{-4\pi/5}$.

Nine distinct squares were found. All have $j=0$ with a diagonal pair on arm two, all are crossed, all satisfy the twisted identity, and all violate the plain bound:

| $t/(c_2r^2)$ | $j_{\rm eff}$ | $r_{\min}$ | side / turn spacing |
|---|---|---|---|
| 0.201 | 0.039 | 0.0778 | 1.46 |
| 0.213 | −0.044 | 0.0829 | 1.45 |
| 0.371 | −0.122 | 0.1332 | 1.57 |
| 0.513 | 0.212 | 0.1748 | 1.66 |
| 0.545 | −0.157 | 0.0902 | 3.41 |

The two smallest sit beside the zero at $r=e^{-4\pi/5}=0.081$, with side $\sqrt2$ turn spacings: they are comb-scale squares. In every row $t/(c_2r^2)\ge\sqrt{|j_{\rm eff}|}$ roughly, as the identity requires. This search was aimed at the $j\approx0$ positions, so it says nothing about graceful twisted squares.

## What this changes

- **v379's smooth comb law is proved,** in v383's range: $t\ge c_2r^{1+1/(2a)}(1-o(1))$, with the constant attained to within 6%.
- **Theorems 370.2 and 374.3 now allow smooth arms,** and their confinement condition concerns graceful squares only. For a symmetric smooth core the path condition need only reach $8\,r_{\rm smooth}(s)\sim s^{2a/(2a+1)}$, not $8\rho(s)\sim s^{a/(1+a)}$.
- **For twisted points, the twist's gradient is what matters.** Its size only creates crossed squares. Its gradient $L$ lets graceful squares out to about $(Ls)^{a/(1+a)}$. So the path condition's reach is governed by how fast the asymmetry varies across scales. v374's events at the comb radius are what a steep twist predicts.
- **Open:** whether the graceful events of a steep twist cancel. A census of graceful squares on twisted smooth arms, aimed at resonances rather than at $j\approx0$, is the next measurement.

[Code](../code/smooth_cores_v386.py) · [Tests](../tests/test_smooth_cores_v386.py) · [Data](../results/v386/smooth_check.json)
