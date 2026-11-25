# v385. Asymmetry below the invisibility radius is harmless

The path condition is the open step in v374's route to squares at an asymmetric two-arm point. Theorem 374.3 asks for a cutoff $s$ at which symmetrizing the point passes no square of side exactly $s$. v374 noted that blanket invisibility cannot supply this, because the symmetrization has to reach out to about $8\,r_{\rm conf}(s)$, far beyond $s$.

This note settles the path condition in one regime and pins down where the obstacle can live. The idea is to symmetrize only what is actually asymmetric. If the twist vanishes everywhere between some tiny radius $\rho$ and the confinement radius, the symmetrizing path moves the curve only inside $B(p,\rho)$. Such a move is invisible once $\rho$ is below the distance from $p$ to the nearest near-square free corner. That distance turns out to have a sharp, explicit lower bound.

**Results.**
1. **Nearest free corner (Lemma 385.1).** Near a two-arm spiral point, every near-square of side $s$ has its free corner at distance at least
$$\eta^*(s)=\frac{\pi/2-\Phi^*-2\nu}{\sqrt2}\;a\,s^{1+1/a}\,\big(1+O(s^{1/a})\big)$$
   from $p$, given by an explicit inequality. Here $\Phi^*$ is the size of the twist at radius $s$ ($\Phi^*=0$ for a symmetric point), and $\nu$ is the phase defect of an inscribed polygon ($\nu=0$ for smooth arms). The reason: the free corner's two neighbours sit at radius $s$, a quarter-turn apart. Two points at the same radius a quarter-turn apart are never on the comb at once. The free corner has to move by $\eta^*$ to fix that.
2. **Invisibility at a spiral point (Corollary 385.2).** By Lemma 384.1 the invisibility radius at $p$ is at least $\min(s/2,\eta^*(s))$. For $a=1/2$ this is $\eta^*\approx0.555\,s^3$ when the point is symmetric at radius $s$. The bound is sharp: on the near-square lattice of v378 the nearest free corner comes within $0.008\%$ of $\eta^*$ (symmetric) and $0.06\%$ (twisted).
3. **Symmetric window (Theorem 385.3).** Suppose the twist vanishes on an annulus $\rho\le r\le\varepsilon$, condition (a) of Theorem 374.3 holds at radius $\varepsilon$, and $\rho<\min(s/2,\eta^*(s))$. Then path invisibility holds automatically, and the point inscribes a square of side greater than $s$.
4. **Scale gaps (Corollary 385.4).** On smooth arms, symmetry on one annulus $[\rho,\varepsilon]$ with $\rho<\rho_{\max}(\varepsilon)$ is enough, where $\rho_{\max}(\varepsilon)$ is computed in closed form. It behaves like $\kappa_a\,\varepsilon^{((1+a)/a)^2}$, that is $\varepsilon^9$ for $a=1/2$.
5. **Where a counterexample would have to be asymmetric (Corollary 385.5).** A two-arm point with no inscribed square must be asymmetric in every window $[\rho_{\max}(\varepsilon),\varepsilon]$. In log scale, its asymmetric scales cannot leave a gap from $\varepsilon$ down to about $\varepsilon^{9}$.
6. **Checks:**
   - v378's near-square lattice at 80 cutoffs between 0.003 and 0.03, symmetric and twisted: the nearest free corner is always at least $\eta^*(s)$. The ratio has minimum $1.00008$ (symmetric) and $1.0006$ (twisted), median about 1.4.
   - v377's independent solver at three cutoffs agrees: minimum ratio 1.076 (symmetric) and 1.002 (twisted).

## Setting

Use the notation of v374–v376. The point $p=0$ has arms
$$\Theta^{-a}e^{i\Theta}\quad\text{and}\quad\Theta^{-a}e^{i(\Theta+\pi+\varphi(\Theta^{-a}))},\qquad\Theta\ge\theta_0,$$
with $|r\varphi'(r)|\le L$. A point $z$ lies on arm $j$ iff its comb phase
$$h_j(z)=|z|^{-1/a}-\arg z+[j=2]\,\big(\pi+\varphi(|z|)\big)$$
is in $2\pi\mathbb Z$. For a polygon inscribed in the arms (v376), every point lies within the sag $h$ of an arm point. Its phase then misses $2\pi\mathbb Z$ by at most
$$\nu\le\tfrac1a(r-h)^{-1/a-1}h+\arcsin\tfrac h{r-h}+L\log\tfrac r{r-h}.$$
Under v376's mesh rule $h\lesssim d_g(r)$, so $\nu$ is a fixed fraction that can be made small.

## The nearest free corner

**Lemma 385.1.** Take $s$ so small that the curve inside $B(p,2s)$ consists of the two arms. Let a square of side $s$ have three corners on the curve and its fourth corner $A$ with $|A|=\eta<s/2$. Put
$$\Phi^*=\max\Big(\sup_{[s-\eta,\,s+\eta]}|\varphi|,\ \operatorname{osc}_{[s-\eta,\,s+\eta]}\varphi\Big).$$
Then
$$2\arcsin\frac\eta s+\frac{2\sin\!\big(\frac\pi4+\arcsin\frac\eta s\big)}{a}\,\eta\,(s-\eta)^{-1/a-1}+2\nu\ \ge\ \frac\pi2-\Phi^* .$$
So every free corner is at distance at least $\eta^*(s)$ from $p$, where $\eta^*(s)$ is the largest $\eta$ violating this inequality.

*Proof.* Write the neighbours of $A$ as $B=A+se^{i\theta}$ and $D=A+ise^{i\theta}$. They are curve points within $s+\eta<2s$ of $p$, so they lie on arms $j_B$ and $j_D$.

**Angles.** $B=se^{i\theta}(1+w)$ with $|w|=\eta/s$. So $|\arg B-\theta|\le\arcsin(\eta/s)$, and likewise $|\arg D-\theta-\pi/2|\le\arcsin(\eta/s)$. Hence $\arg D-\arg B=\pi/2+e$ with $|e|\le2\arcsin(\eta/s)$.

**Radii.** $f(A)=|A+se^{i\theta}|-|A+ise^{i\theta}|$ vanishes at $A=0$. Its gradient is the difference of two unit vectors at angle $\pi/2+e$, so it has norm at most $2\sin(\pi/4+\arcsin(\eta/s))$. Hence $\big||B|-|D|\big|\le2\eta\sin(\pi/4+\arcsin(\eta/s))$. Both radii lie in $[s-\eta,s+\eta]$, so the mean value theorem bounds $\big||B|^{-1/a}-|D|^{-1/a}\big|$ by the second term of the inequality.

**Phases.** Subtracting the arm conditions for $B$ and $D$,
$$\tfrac\pi2+e+\big(|B|^{-1/a}-|D|^{-1/a}\big)+o_B-o_D+(\text{defects})\in2\pi\mathbb Z,$$
where $o$ is the arm offset and the defects total at most $2\nu$. The offset difference $o_B-o_D$ is one of four things:
- $0$, if both are on arm one;
- $\varphi(|B|)-\varphi(|D|)$, if both are on arm two, at most the oscillation;
- $\pi+\varphi(|B|)$, if $B$ is on arm two and $D$ on arm one;
- $-\pi-\varphi(|D|)$, if $B$ is on arm one and $D$ on arm two.

So $\frac\pi2+(o_B-o_D)$ lies at distance at least $\frac\pi2-\Phi^*$ from $2\pi\mathbb Z$. The remaining terms must close that distance, which is the inequality. ∎

The left side increases with $\eta$, and it is dominated by the radial term $\sqrt2\,\eta\,s^{-1/a-1}/a$. That gives the leading order in the summary. Only $B$ and $D$ are used, so the bound holds for all eight arm triples at once. It also holds whatever $C$ does.

**Corollary 385.2.** If the curve has no square of side $s$, every closed curve inside $\gamma\cup B(p,\rho)$ with $\rho<\min(s/2,\eta^*(s))$ has no square of side $s$.

*Proof.* Lemma 384.1, with $d_N\ge\eta^*(s)$ by Lemma 385.1. ∎

v374 showed the invisibility radius is at most $d_N$. Together with the lattice data below, the invisibility radius at a spiral point is $\eta^*(s)$ up to a factor that comes arbitrarily close to 1 at the worst cutoffs. This also turns v377's observation, nearest free corners at $0.32$–$1.16$ turn spacings, into a theorem with the sharp constant. For a symmetric point the bound is $\eta^*/d(s)=1/(2\sqrt2)\approx0.354$ turn spacings.

## The symmetric window

**Theorem 385.3.** Let $p$ be a two-arm point in the setting of Theorem 374.3. Let $0<\rho<\varepsilon$ and $s>0$ satisfy:
- the twist vanishes on $[\rho,\varepsilon]$;
- condition (a) of Theorem 374.3 holds at radius $\varepsilon$ for the curve $\gamma_1$ that is half-turn symmetric on $B(p,\varepsilon)$ (by Corollary 375.3 this holds for smooth arms when $s\le d_g(\varepsilon/8)$, and by v376 with a sag correction for polygons);
- $\rho<\min(s/2,\eta^*(s))$, with $\eta^*$ computed for $\Phi^*=0$.

Then $\gamma$ inscribes a square of side greater than $s$.

*Proof.*
- If $\gamma$ has a square of side $s$, we are done.
- Otherwise, run the gap path of Lemma 374.1 with its cutoff equal to 1 on $[0,\rho]$ and 0 beyond some $\rho'$ with $\rho<\rho'<\min(s/2,\eta^*(s))$. On $[\rho,\rho']$ the twist is already 0, so the interpolation changes nothing there. Every curve on the path lies inside $\gamma\cup B(p,\rho')$, and the end curve is $\gamma_1$, symmetric on all of $B(p,\varepsilon)$.
- Since $\rho'<s/2$, and condition (a) forces $s\sqrt2<\varepsilon$, the window $[s-\eta^*,s+\eta^*]$ lies inside $[\rho,\varepsilon]$. So $\Phi^*=0$ there, and $\eta^*$ is the symmetric value. Corollary 385.2 gives condition (b) of Theorem 374.3.
- Condition (a) is assumed. Theorem 374.3 finishes the proof. ∎

**Corollary 385.4 (smooth arms).** Let $s_\varepsilon=d_g(\varepsilon/8)$. If the twist vanishes on $[\rho,\varepsilon]$ with
$$\rho<\rho_{\max}(\varepsilon)=\min\big(s_\varepsilon/2,\ \eta^*(s_\varepsilon)\big),$$
the point inscribes a square. The best cutoff is the largest one allowed, because $\eta^*$ grows with $s$. Since $s_\varepsilon\asymp\varepsilon^{1+1/a}$ and $\eta^*(s)\asymp s^{1+1/a}$, we get $\rho_{\max}(\varepsilon)\asymp\varepsilon^{((1+a)/a)^2}$. For $a=1/2$:

| $\varepsilon$ | best cutoff $s_\varepsilon$ | $\rho_{\max}$ | $\log\rho_{\max}/\log\varepsilon$ |
|---|---|---|---|
| 0.3 | $8.3\times10^{-5}$ | $3.1\times10^{-13}$ | 23.9 |
| 0.1 | $3.1\times10^{-6}$ | $1.6\times10^{-17}$ | 16.8 |
| 0.01 | $3.1\times10^{-9}$ | $1.6\times10^{-26}$ | 12.9 |
| $10^{-8}$ | $3.1\times10^{-27}$ | $1.6\times10^{-80}$ | 9.97 |

The ratio tends to 9 slowly, because of constant factors such as $8^9$. The numbers are tiny, but they are explicit.

**Corollary 385.5.** If a two-arm point inscribes no square, then for every small $\varepsilon$ its twist is nonzero somewhere in $[\rho_{\max}(\varepsilon),\varepsilon]$. So a counterexample must be asymmetric at scales recurring at least once in every such window. The plain asymmetric profile of v374, $\varphi(r)=0.4\sin(5\log(1/r))$, is asymmetric at almost every scale, so the corollary says nothing about it. Points whose asymmetry thins out, with symmetric stretches spanning a factor of $\varepsilon^{-8}$ in radius, are now covered.

## Checks

**The lattice.** For each cutoff, v378's lattice enumeration lists every near-square whose free corner lies within $4K(s)s$ of $p$, across all eight arm triples, on the smooth arcs. The table compares the nearest free corner with $\eta^*(s)$, computed from the exact inequality with $\Phi^*$ evaluated on its own window.

| twist | cutoffs | min nearest / $\eta^*$ | median |
|---|---|---|---|
| symmetric ($\lambda=0$) | 80, $s\in[0.003,0.03]$ | 1.00008 | 1.44 |
| twisted ($\lambda=1$) | 80, same | 1.0006 | 1.37 |

**An independent solver.** v377's solver places the free corner in the comb cell nearest $p$ and turns the square until the opposite corner closes. At $s=0.02,0.01,0.005$ the minimum ratio is 1.076 (symmetric) and 1.002 (twisted).

**Leading order.** The exact $\eta^*$ and the leading-order formula agree to $0.7\%$ at $s=0.05$ and to $0.007\%$ at $s=0.005$.

## What this changes

- **The path condition is now localized.** In Theorem 374.3 the symmetrization only has to touch the scales where the twist is nonzero. Scales below $\eta^*(s)\approx0.35$ turn spacings are invisible, with the sharp constant. So the obstacle lives only where v374 saw it: twist between about $s^{1+1/a}$ and the comb radius.
- **A first class of genuinely asymmetric points is protected:** points whose twist vanishes on one long enough annulus. This goes beyond v370–v372, which needed exact symmetry.
- **Open:** points asymmetric at every scale, such as v374's default profile. There the symmetrization must rotate turns at the comb radius, and v374 found that this does create and destroy squares of side $s$. The question is whether those events must cancel in pairs. The surgery rule (v378) counts them exactly, as graceful near-squares swept by the rotating comb.

[Code](../code/symmetric_window_v385.py) · [Tests](../tests/test_symmetric_window_v385.py) · [Data](../results/v385/window_check.json)
