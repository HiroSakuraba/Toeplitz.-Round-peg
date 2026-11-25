# v392. Big odd squares are pinned uniformly

v391 certified closed regions one label at a time, for pass offsets up to 20. That left the obvious question open: does *every* big enough odd label close, however large its offsets get? This note answers yes, with explicit thresholds, using one computer check for infinitely many labels at a time.

**The idea.** Measure everything relative to the label's own resonance. Write:
- $m=\min(|n_1|,|n_3|)$ and $M=\max(|n_1|,|n_3|)$ for the offsets;
- $\Theta^*=(2a+1)\pi mM/|j|$ for the resonance, and $s=\Theta_A/\Theta^*$ for the scaled lift;
- $\varepsilon=1/m$ and $\mu=1/M$.

Then each error term in v391's estimates becomes a function of $(s,\varepsilon,\mu)$ alone, and it only grows as $\varepsilon$ or $\mu$ grows. So a check at one value of $\varepsilon$ covers every label with offsets at least that big. It works like checking a speed limit once at the steepest point of a hill.

**Results.** These are computer-assisted, in interval arithmetic with outward rounding. They cover the twist bounds $(\Phi,L)=(0.4,2)$, $(0.8,4)$ and $(1.2,6)$: any twist with $|\varphi|\le\Phi$ and $|r\varphi'|\le L$, not just the sine twist used for testing.

1. **All levels at once (Theorem 392.1).** Every odd label with $\min(|n_1|,|n_3|)\ge c\,|j|$ has a certified region, for any level $j$, any larger offset and any signs.
   - $c=11$, 13, 16 for $\Phi=0.4$, 0.8, 1.2.
   - Restricting to levels $|j|\ge3$: $c=9$, 10, 11.
   - Restricting to $|j|\ge9$: $c=8$, 8, 9.
2. **Sharp thresholds for small levels (Theorem 392.2).** For $|j|=1$, every odd label with $\min(|n_1|,|n_3|)\ge10$, 12, 14 has a certified region. These are exactly the thresholds v391 saw in its table up to 20, now proved with no upper limit on the offsets. For the next levels:

   | $|j|$ | $\Phi=0.4$ | $\Phi=0.8$ | $\Phi=1.2$ |
   |---|---|---|---|
   | 1 | 10 | 12 | 14 |
   | 3 | 24 | 26 | 28 |
   | 5 | 37 | 40 | 42 |
   | 7 | 51 | 53 | 56 |

3. **Infinitely many squares (Corollary 392.3).** By v390, a certified region keeps its label's degree for every twist within the bounds. So **every two-arm point with such a twist carries a square at every covered label, if the symmetric spiral's degree on that region is nonzero.** There are infinitely many covered labels, and their squares have side about $|A|/m$, so they accumulate at the wild point. The degree condition is the uniqueness assumption U0 from v390, which is supported numerically but not proved.
4. **The method works for any twist size (Proposition 392.4).** For any $\Phi<\pi$ and any $L$, some constant $c(\Phi,L)$ works. As $\varepsilon\to0$ every enclosure tends to its limit: the twist slack tends to $\Phi$, the ratio $F$ tends to 1, and the side bound stays finite. The walls $s_{\rm lo}<\pi|j|/(\pi|j|+\Phi)$ and $s_{\rm hi}>\pi|j|/(\pi|j|-\Phi)$ are excluded with strict margins, and all enclosures are continuous in $\varepsilon$.
5. **Checks.**
   - v391's independent per-label certifier, run at v392's face lifts, agrees on 90 of 90 random covered labels, with offsets up to 256.
   - Plain-spiral squares at 15 big labels (offsets up to 100, levels up to 5, all sign patterns) were followed through twists ×1, ×2, ×3. All 45 runs reach the end with their local degree kept. All 37 covered runs stay strictly inside their certified region.

## The pieces

**Scaled formulas.** Relative to $\Theta_A$:
- the larger offset sits at $\varepsilon\beta_0$ and the smaller at $\mu\beta_0$, with $\beta_0=|j|/((2a+1)s)$;
- the level sits at $j|j|\varepsilon\mu/((2a+1)s)$;
- each corner's deviation is at most $e\,\varepsilon\mu\beta_0/\pi$.

The larger offset is controlled by the *smaller* one, because $\pi M/\Theta_A=|j|/((2a+1)ms)$. That is why the error terms depend on $\min(|n_1|,|n_3|)$.

v386's identity, divided by $\pi|j|$, becomes
$$1\pm\frac{\Phi+{\rm osc}+d}{\pi|j|}=\frac{1}{s}\,\frac{\zeta^{2a+1}}{\xi^{2a+2}}\Big(1\pm\frac{e\varepsilon}{\pi}\Big)\Big(1\pm\frac{e\mu}{\pi}\Big),$$
with $\zeta$ and $\xi$ the intermediate lifts over $\Theta_A$, as in v391. In the variables $\tilde\varepsilon=|j|/m$, $\tilde\mu=|j|/M$ and $\kappa=1/|j|$, the level disappears except through $\kappa$. Every term grows with $\kappa$, which is how one check covers all levels $|j|\ge j_0$.

**The identity step needs no monotonicity.** For each of the four sign patterns, the intermediate lifts are enclosed with $\varepsilon$ and $\mu$ entered as intervals (for example $\varepsilon\in[0,1/N]$). Interval arithmetic then covers every value at once.

Treating the sign patterns separately matters. When $n_1$ and $n_3$ have the same sign, all offsets lie on one side of $\Theta_A$, which shrinks the enclosure about threefold compared with a two-sided bound. A first version without this needed offsets of 36–49.

**The side bound needs two monotonicity facts.** The side map $G$ (v391) is bounded by $\sqrt{(c_M+\kappa_q)^2+(c_m+\kappa_q)^2}$. Here $c_n=a x_n/(1-x_n)^{a+1}$, the radial steps are $x_M=\varepsilon\beta_0+\varepsilon\mu e_{\max}\beta_0/\pi$ (likewise $x_m$), and $\kappa_q$ is the curvature term.
- **Lemma 392.5 (in the lift, the offsets and the side).** $G$ increases with $\varepsilon$, $\mu$ and $q$, and decreases with $s$. So a side bound at the lower wall covers every larger lift.
- **Lemma 392.6 (in the scale).** $H(Q;\varepsilon)=G(\varepsilon Q;\varepsilon)/\varepsilon$ increases with $\varepsilon$, since $c_M/\varepsilon$ and $\kappa_q(\varepsilon Q)/\varepsilon$ both do.

So a march $G(q)<q$ over $[\hat Q/N,\,0.3]$ at $\varepsilon=1/N$ gives $H(Q;\varepsilon')<Q$ for $Q\in[\hat Q,\,0.3N]$, for every $\varepsilon'\le1/N$. The remaining sides $q\in[0.3N\varepsilon',\,0.3]$ follow from one scalar inequality,
$$\sqrt2\,C_1\varepsilon/0.3+\sqrt2\,\kappa(0.3)/0.3<1,\qquad C_1=c_M/\varepsilon\ \text{ at }\varepsilon=1/N,$$
using $c_M(\varepsilon')\le\varepsilon' C_1$ and that $\kappa(q)/q$ increases. Together: every covered label has side at most $\varepsilon\hat Q\,|A|$ at every lift beyond the lower wall, with $\varepsilon$ its own $1/m$, and none has side $0.3|A|$.

**The inner side face.** The side is at least the radial step to the corner with the larger offset, minus curvature. In scaled form $q\ge\varepsilon f$, where
$$f=\frac{a(\beta_0-\mu e_B\beta_0/\pi)}{(1+\varepsilon\beta_0+\varepsilon\mu e_B\beta_0/\pi)^{a+1}}-\frac{\kappa(\varepsilon\hat Q)}{\varepsilon}.$$
This decreases with $\varepsilon$, $\mu$ and $s$, so it is evaluated at the worst corner of the family and at the upper wall. The face sits at $\varepsilon f/2$.

**Walls.** $s_{\rm lo}=\frac{\pi|j|}{\pi|j|+\Phi}(1-\delta_{\rm lo})$ and $s_{\rm hi}=\frac{\pi|j|}{\pi|j|-\Phi}(1+\delta_{\rm hi})$, with each width picked independently from a grid. The upper face also certifies the whole tail beyond it, because the right side's upper enclosure decreases in $s$ while the twist slack is uniform.

**The certificate types.**
- *Uniform*: $\varepsilon=\mu\in[0,1/N]$ covers every label of the level with both offsets at least $N$.
- *Family*: $\varepsilon=1/m$ exactly with $\mu\in[0,1/M_0]$ covers a fixed smaller offset $m$ with every $M\ge M_0$.
- *Single*: one offset pair, falling back on v391's per-label certificate if needed.
- *All levels*: $\tilde\varepsilon=\tilde\mu\in[0,1/c]$ and $\kappa\in[0,1/j_0]$, with walls at the worst level $j_0$.

The level-1 thresholds combine these. For example, at $\Phi=0.4$:
- the uniform check covers $m\ge11$;
- $m=10$ needs the family for $M\ge11$ plus the single label $(10,10)$;
- $(9,9)$ fails.

Each threshold is certified by its own checks, so nothing rests on the bisection used to find it.

## Checks

| | $\Phi=0.4$ | $\Phi=0.8$ | $\Phi=1.2$ |
|---|---|---|---|
| all levels, $\min$ offset $\ge c|j|$: $c$ | 11 | 13 | 16 |
| walls $s_{\rm lo}$, $s_{\rm hi}$ | 0.444, 1.490 | 0.319, 1.744 | 0.289, 2.103 |
| side bound $\hat Q$ ($q\le\hat Q|j|/m$), inner floor $f$ | 1.05, 0.103 | 1.57, 0.030 | 1.67, 0.016 |
| v391 agrees at the faces (random covered labels) | 30 / 30 | 30 / 30 | 30 / 30 |
| big-label runs reaching the end, degree kept | 15 / 15 | 15 / 15 | 15 / 15 |
| covered runs staying strictly inside their region | 15 / 15 | 11 / 11 | 11 / 11 |

The twist bound for each column is that of the run's twist, $m\cdot0.4\sin(5\log1/r)$, with $m=1,2,3$. Labels below a column's constant aren't covered there. They still stayed inside numerically, but that is not certified.

The squares move noticeably. At $(-16,23,-1)$ under twist ×3, the scaled lift goes from 1.004 to 1.585, well short of the certified wall at 2.103.

## What is and isn't proved

- **Proved (computer-assisted):** certified regions for every odd label above the thresholds, for the three twist bounds tabulated.
- **Proved (continuity argument):** for every $\Phi<\pi$ and $L$, some threshold works.
- **Conditional:** existence of squares at those labels needs a nonzero degree of the symmetric spiral on each region, which is U0.
- **Not proved:**
  - a direct, uniform proof of nonzero degree. The leading-order system as $\varepsilon\to0$ is explicit, so a perturbation argument looks feasible; that's the natural next step;
  - odd labels below the thresholds;
  - even labels;
  - anything outside the two-arm power-spiral model.

**Where this leaves the program.** Within the model, a bounded twist can never remove the odd squares at big labels, provided their degree starts out nonzero. Proving that degree condition uniformly would make "every bounded-twist two-arm point has infinitely many small squares" an unconditional theorem.

[Code](../code/uniform_closure_v392.py) · [Tests](../tests/test_uniform_closure_v392.py) · [Data](../results/v392/uniform_check.json)
