# v391. Odd squares can't collapse, and big ones are pinned

v390 reduced survival to one question: can a label's square be pushed out of its region, either into the wild point or away from it? This note answers it for squares with an **odd number of corners on arm two** ("odd squares"), using the exact identity from v386 and v390's anchored labels.

**The idea.** For an odd square the identity's left side contains $\pi j$ with $j$ odd. A twist of size $\Phi<\pi$ can shift it, but it cannot cancel it. The right side is roughly $(2a+1)\,\pi^2 n_1n_3/\Theta_A$. So the lift $\Theta_A$ of corner $A$ (which grows as the square moves inward) is held near the resonance $(2a+1)\pi|n_1n_3|/|j|$, whatever the twist does.

**Results.** Below, *certified* means checked in interval arithmetic with outward rounding, any failed step counting as "not certified". *Sampled* means floating-point checks on a grid, which prove nothing between or beyond the samples.

1. **No collapse (Theorem 391.1).** Fix an odd label and twist bounds $|\varphi|\le\Phi<\pi$, $|r\varphi'|\le L$. Then there is a lift $T_\ell$ above which the label has no squares of side at most $0.3|A|$, for any such twist. So no odd square can shrink into the wild point.
   - For each fixed label this follows from the limits as $\Theta_A\to\infty$, together with a tail certificate that turns one check at $T$ into exclusion on $[T,\infty)$.
   - $T_\ell$ is **certified** for all 1 600 sign-compatible labels with $|j|=1$ and offsets up to 20, at each of $\Phi=0.4,0.8,1.2$ ($L=5\Phi$).
   - The bound depends on the label. It does not stop a *sequence of different* labels from approaching the point.
2. **Big odd squares are pinned (Theorem 391.2).** A label region $\Omega_\ell=\{T_{\rm lo}<\Theta_A<T_{\rm hi},\ q_{\rm in}<|w|<0.3\}$ is **certified** when no square of the label can sit on any part of its boundary, for any admissible twist. The boundary has four parts: both lift faces, the outer side face, and the inner side face. By v390 the label's degree is then preserved along any path of such twists.
   - Every sign-compatible label with $\min(|n_1|,|n_3|)\ge N$ has a certified region, with $N=10$, 12, 14 for $\Phi=0.4$, 0.8, 1.2 (offsets up to 20).
   - Certified regions also exist at eight larger labels, out to offsets of 200 and including $|j|=3$.
   - So every bounded-twist two-arm point with these bounds carries a square at each certified label, **if** the symmetric spiral's degree on that region is nonzero. That holds under the uniqueness assumption U0, which is supported numerically but not proved.
   - Extending this to *all* big labels, which would give infinitely many squares accumulating at the point, needs a uniform argument that isn't written yet.
3. **Small odd squares survive numerically.**
   - Setup: all 156 odd squares from v382 with side at most $0.2|Z|$, followed along the twists $m\cdot0.4\sin(5\log1/r)$ for $m=1,2,3$.
   - All 156 reach the full twist at each $m$, with final residuals below $10^{-12}$, and all keep a nonzero local degree.
   - All stay below their certified collapse bound; the largest lift reached is 78% of it.

## The pieces

**Coordinates.** As in v390: $A$ sits on its arm at lift $\Theta_A$, and $w$ is the side vector in $A$'s frame divided by $|A|$. The label is $(\mathrm{arm}_A,o_B,o_C,o_D)$. Write $n_1=o_B$, $n_3=o_D$ and $j=o_C-o_B-o_D$. A label is *odd* when an odd number of its corners lie on arm two.

**The identity, with the right ranges.** Write $C'=\Theta_B+\Theta_D-\Theta_A$. Then
$$\Theta_C-C'=\pi j-\Sigma^*+\Delta\theta=(2a+1)\,\zeta^{2a+1}\xi^{-2a-2}\,\delta_B\delta_D,$$
with $\delta_B=\Theta_B-\Theta_A$ and $\delta_D=\Theta_D-\Theta_A$.
- $\xi$ comes from the second difference of $\Theta^{-2a}$ over $A$, $B$, $D$, $C'$, so it lies in their hull.
- $\zeta$ comes from a first difference, so it lies between $\Theta_C$ and $C'$.

(v386 stated a cruder range. The tight ranges are what make lower walls possible.)

**The error terms.** Here $q$ bounds $t/|A|$.
- **Lift deviations.** Each corner's lift is $\Theta_A+\pi o_X$ up to $e=\arcsin(\sqrt2q)+\max(\Phi,-L\log(1-\sqrt2q))$. That is its angle seen from $A$, plus the twist difference: at most $\Phi$ if one of the two corners is on arm one, and at most $L|\log(|X|/|A|)|$ if both are on arm two.
- **Twist sum.** For an odd label, $|\Sigma^*|\le\Phi+{\rm osc}$, with ${\rm osc}=L\log\frac{1+u}{1-u}$ and $u=(q/\sqrt2)/(1-q/\sqrt2)$. With one arm-two corner, $\Sigma^*$ is a single twist value; with three, it is a value plus a difference.
- **Angular flag.** In anchored terms $\Delta\theta=\arg\big(1-\frac{bd}{(1+b)(1+d)}\big)$, with $b=w$ and $d=iw$. So $|\Delta\theta|\le\arcsin\frac{q^2}{(1-q)^2}$.
- **Sign.** Both sides have definite signs, so squares of the label need $\operatorname{sign}j=\operatorname{sign}(n_1n_3)$. Only these sign-compatible labels are tabulated.

**A bound on the side.** Because $AB\perp AD$ with $|AB|=|AD|=t$, the radial components satisfy $t^2=\langle AB,\hat u\rangle^2+\langle AD,\hat u\rangle^2$ *exactly*. Each component is $|B|-|A|$ (controlled by the lift step) minus a curvature term at most $(1+q)(1-\sqrt{1-q^2})|A|$. This gives $q_t\le G(q_t,\Theta_A)$ for an increasing map $G$.
- **Certified version.** Take $\hat q$ slightly above $G$'s least fixed point. Check in interval arithmetic that $G(\hat q)\le\hat q$, and that $G(\text{cell})<\min(\text{cell})$ on adaptive cells covering $[\hat q,0.3]$. Then every square of the label with side at most $0.3|A|$ has $q_t\le\hat q$, and none sits exactly at $0.3|A|$.
- **Monotonicity in the lift.** $G$ decreases in $\Theta_A$, because its radial term $a\,\delta\,\Theta^a/(\Theta-\delta)^{a+1}$ has logarithmic derivative $a/\Theta-(a+1)/(\Theta-\delta)<0$. So a certificate at a lift covers every larger lift.

## The certificates

- **A lift face** $\Theta_A=T$: with the certified $\hat q$, interval enclosures of $|\pi j-\Sigma^*+\Delta\theta|$ and of the right side must be disjoint. Then no square of the label, of any orientation and any side up to $0.3|A|$, has $A$ at lift $T$.
- **The tail** $[T,\infty)$:
  - $\hat q(T)$ bounds the side for every $\Theta_A\ge T$, so the error terms are uniform.
  - The left side is at least $N_{\rm lo}=\pi|j|-\Phi-{\rm osc}-d$.
  - The right side is at most $(2a+1)D_{\rm hi}(\Theta+P)^{2a+1}/(\Theta+Q)^{2a+2}$, where $P\ge Q$ are the largest $\zeta$-offset and the smallest $\xi$-offset. This decreases in $\Theta$ for $\Theta+Q>0$, since its logarithmic derivative is $\frac{2a+1}{\Theta+P}-\frac{2a+2}{\Theta+Q}<0$.
  
  So one interval check at $T$ covers $[T,\infty)$. As $T\to\infty$, $\hat q\to0$ and the right side tends to 0 while $N_{\rm lo}\to\pi|j|-\Phi>0$. So a tail certificate exists for every fixed label (Theorem 391.1).
- **The outer side face** $|w|=0.3$: the certified side bound at $T_{\rm lo}$ covers the whole lift interval, by monotonicity.
- **The inner side face** $|w|=q_{\rm in}$: $t\ge|\langle AB,\hat u\rangle|\ge\big||B|-|A|\big|-|B|(1-\cos)$, with $\big||B|-|A|\big|\ge a(\pi|n_1|-e)\max(\Theta_A,\Theta_B)^{-a-1}$. The resulting lower bound on $q$ decreases in the lift once $\Theta_A>a(\pi|n_1|+e)$, so it is evaluated at $T_{\rm hi}$. The face is put at half of it.
- **Twists.** All bounds use only $|\lambda\varphi|\le\Phi$ and $|r(\lambda\varphi)'|\le L$, so they hold for every twist on the path $\lambda\in[0,1]$.

The sampled band (a grid of 800–6 000 lifts) only *locates* candidate walls. They are the excluded samples just outside the allowed run, and are reported separately from the allowed range. If a candidate wall has too little margin for the interval check, it is slid outward while the sampled status stays "excluded", and then certified.

## Why small labels don't close

The identity sees $\Theta_A$ only through intermediate lifts spread over about $\pi(|n_1|+|n_3|)$. Next to $\Theta_A\approx2\pi|n_1n_3|$, that is a relative uncertainty of order $1/\min(|n_1|,|n_3|)$, raised to the power $2a+2$. For offsets below about 10 this washes out the lower exclusion before the side bound gives out.

The tail needs no such precision, so collapse is ruled out for every label. Drifting outward is ruled out only where a region is certified. For the small labels the numerics show no drift, but that is evidence, not proof.

## Checks

| | $m=1$ | $m=2$ | $m=3$ |
|---|---|---|---|
| twist size $\Phi$, slope bound $L$ | 0.4, 2 | 0.8, 4 | 1.2, 6 |
| odd squares followed | 156 | 156 | 156 |
| reach the full twist (final residual) | 156 ($<10^{-12}$) | 156 | 156 |
| nonzero local degree, kept to the end | 156 | 156 | 156 |
| certified collapse bound found | 156 | 156 | 156 |
| stay below it (largest lift / bound) | 156 (0.78) | 156 | 156 |
| label has a certified region / square stays inside it | 4 / 4 | 0 / 0 | 0 / 0 |
| side / side bound, worst (squares where the bound applies) | 0.89 (101) | 0.88 (102) | 0.87 (102) |
| largest change of $|w|$ along the path | 23% | 54% | 90% |

At $m=1$, the four labels with certified regions have offsets $(n_1,n_3)=(\pm8,\mp8)$, $(7,7)$ and $(8,8)$. v382's sample stops at offsets of 8.

Closure table (sign-compatible labels, $|j|=1$, all sign patterns, $L=5\Phi$, offsets 1–20, 1 600 labels for each $\Phi$):

| $\Phi$ | certified collapse bound | sampled band with two walls | offset pairs with a certified region | every pair with $\min\ge N$ certified |
|---|---|---|---|---|
| 0.4 | 1600 / 1600 | 853 | 171 / 400 | $N=10$ |
| 0.8 | 1600 / 1600 | 649 | 119 / 400 | $N=12$ |
| 1.2 | 1600 / 1600 | 468 | 75 / 400 | $N=14$ |

Certified regions and tails also exist at all three twist sizes, in every sign pattern, for $(25,25,1)$, $(30,40,1)$, $(50,50,1)$, $(40,120,1)$, $(100,100,1)$, $(25,200,1)$, $(60,60,3)$ and $(90,150,3)$.

An example: $(12,-13,-1)$ has a certified region $\Theta_A\in(678,1299)$, $|w|\in(0.0068,0.3)$ at $\Phi=0.4$, and $(406,1946)$ at $\Phi=1.2$, around the resonance 980. Its side is certified to stay below $0.048|A|$ and $0.089|A|$ respectively.

## What is and isn't proved

- **Proved:** for each fixed odd label, no collapse under bounded twists (Theorem 391.1), with certified bounds for the 1 600 labels per twist size listed above.
- **Certified for the listed labels:** closed regions, with every boundary part checked.
- **Conditional:** existence at a certified label needs a nonzero degree of the symmetric spiral on its region. U0 would give that, and is supported numerically (v390) but not proved. A direct proof of nonzero degree would remove the assumption.
- **Not proved:**
  - that every big label has a certified region, which needs a uniform argument in the offsets;
  - outward drift of small odd squares;
  - anything about even squares, where $j$ is even, $\pi j$ can be 0, and the identity forces nothing.
- **Inputs:** v386's identity and v390's anchored labels.

**Where this leaves the program.** Within the two-arm power-spiral model, a bounded twist cannot make an odd square collapse into the wild point, and at the certified big labels it cannot make one disappear at all (given nonzero initial degree). The next steps are a uniform large-label argument and a direct proof of nonzero degree. Together those would give unconditional existence of infinitely many squares at every bounded-twist two-arm point.

## Review corrections

An external review of this PR made four points, all fixed here.
- **Wall indexing (F2).** The band's endpoints were the first and last *allowed* samples, not walls. Walls are now the excluded lifts just outside, reported separately.
- **Sampling (F3).** Every boundary part is now certified in interval arithmetic: both lift faces, both side faces, and a tail. The tail is reduced to one check by proved monotonicity. Sampled results are labelled diagnostics.
- **Tests (F8).** The tests now require 156/156, with nonzero degrees at both ends. They regenerate two continuations in CI rather than only reading stored results.
- **Coordinates and U0 (F1, F4).** Continuation and degrees use v390's anchored, chart-free map (the first version's polar-angle chart jumped by $2\pi$ across the negative real axis). U0 is stated as an assumption.

[Code](../code/pinned_squares_v391.py) · [Tests](../tests/test_pinned_squares_v391.py) · [Data](../results/v391/pinned_check.json)
