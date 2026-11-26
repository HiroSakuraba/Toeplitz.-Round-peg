# v396. Any winding rate

Everything since v386 assumed arms winding like $r^{-2}$: the power spiral $r=\Theta^{-a}$ with $a=\tfrac12$. v395 allowed slow angular deviations from that profile, but a different exponent isn't slow: $r^{-1/a'}-r^{-2}$ has unbounded log-slope. This note removes the restriction. Every certificate is now run with the exponent $a$ as an interval, so one run covers every spiral with $a$ in that interval at once. Five intervals cover $a\in[\tfrac14,4]$, i.e. winding rates from $r^{-4}$ down to $r^{-1/4}$.

(A constant factor in the winding, $\theta=Cr^{-1/a}$, is a dilation of the picture by $C^a$. Inscribed squares don't care about scale, so the exponent is the only real parameter of a power-law profile.)

**Results.** These are computer-assisted, in interval arithmetic with outward rounding, for level-one labels ($|j|=1$) with comparable offsets ($M\le2m$), as before.

1. **Theorem 396.1 (nonzero degree, any exponent).** For every $a\in[\tfrac14,4]$, every such label with $m\ge30$, in each of the eight sign-and-corner patterns, has exactly one zero of its label map in the certified region on the plain spiral, and that zero is nondegenerate. The degree is $-1$ for $j=+1$ and $+1$ for $j=-1$, the same signs as at $a=\tfrac12$.
2. **Theorem 396.2 (existence, any exponent).** Take a two-arm curve
   $$\theta_1(r)=r^{-1/a}+e_1(r),\qquad\theta_2(r)=r^{-1/a}-\pi+e_1(r)+\varphi(r)$$
   with $a$ in one of the six exponent pieces, $|\varphi|\le\Phi$, $|r\varphi'|\le L$ and $|re_1'|\le L_1$. For each of the four tabulated $(\Phi,L,L_1)$ and each piece, a certified threshold $c$ is given such that the curve inscribes a square of every label with $m\ge c$, $M\le2m$, $|j|=1$. So **infinitely many squares accumulate at the point**, whatever the winding exponent in $[\tfrac14,4]$. All 24 entries satisfy the degree conditions: the region's lower lift is at least $\tfrac14$ and the region contains $B^*$.
3. **Checks.** Real plain-spiral squares were tested at five labels (offsets 30 to 70, all four sign patterns) for each of five exponents, $a=0.3,0.55,1,1.7,3$.
   - All 25 lie inside the certified box for their parameters.
   - Their local degrees all match the certified sign.

## What changed in the degree step

v393's zoomed system had closed-form remainders that used $a=\tfrac12$: the radius factor $|1+w|^{-1/a}=|1+w|^{-2}$ is a rational function. For general $a$ it is $(1+\rho)^{-p}$ with $\rho=2\operatorname{Re}w+|w|^2$ and $p=1/(2a)$, and the remainder has to be enclosed rather than written down.

**Taylor with integral remainder.** Exactly,
$$(1+\rho)^{-p}=1-p\rho+\tfrac{p(p+1)}2\rho^2K(\rho),\qquad K(\rho)=2\int_0^1(1-t)(1+t\rho)^{-p-2}\,dt .$$
$K$ is an average of $(1+t\rho)^{-p-2}$ with positive weights, so it lies in the hull of that function over $t\in[0,1]$. The same holds for its derivatives, and for $G=(K-1)/\rho$ and $H=(G+\frac{p+2}3)/\rho$, the next two remainders. Each is enclosed by one power of the interval $[\min(1,1+\rho),\max(1,1+\rho)]$, which is cheap and works for an interval $p$.

**Why the limit system doesn't depend on $a$.** The leading terms of the three scaled equations come from the quadratic part of the expansion, with coefficient $p(p+1)$, times the square of the zoom scale $c=a|j|/(2a+1)$. These combine as
$$p(p+1)\,a^2=\frac{(2a+1)}{4},$$
which cancels the $(2a+1)$ in the scale exactly. So the limit system is v393's $(sU-1,\ sV-1,\ \sigma_1\sigma_3sUV-\operatorname{sign}j)$ for every exponent, with the single simple zero $(1,1,1)$. What depends on $a$ is only the size of the corrections.

**Keeping the interval in $a$ from blowing up.** Written naively, the third equation subtracts large $a$-dependent terms, and with $a$ an interval those don't cancel (each copy of $a$ ranges independently). Two exact rearrangements fix this:
- the zoom scale is factored out of every term (variables $\hat x=-\sigma_1t_1U$, $\hat y=\sigma_3t_3V$), so $a$ appears only through $p$ and the small parameter $\varepsilon c$;
- the cubic terms of $R_C-R_B-R_D$ are combined by hand into $S=2(\hat x^3-\hat y^3)+(4p+6)\hat x\hat y(\hat x-\hat y)$, with prefactor exactly $\tfrac14$.

After that, the width of each equation shrinks with the width of the $a$-interval, and Krawczyk's test succeeds on pieces of the parameter box $(\varepsilon,r,a)$.

**Speed.** The search spends nearly all its time evaluating the map on boxes it is trying to rule out, where only values are needed. The forward-mode class now skips derivative work in that case, and the hull powers are cached. One box evaluation went from about 65 ms to about 3 ms.

## The four certified steps

1. **Fine step** (Krawczyk plus localization around $(1,1,1)$), on $\varepsilon=1/m\in[0,\tfrac1{30}]$, $r=m/M\in[\tfrac12,1]$, and each $a$-piece. All 40 combinations (five pieces, eight patterns) are certified, using about 2,400 parameter pieces in total. The degree signs are $-1$ throughout $j=+1$ and $+1$ throughout $j=-1$. The $a$-interval has to be cut fine ($m\ge30$ rather than v393's $m\ge11$) for Krawczyk's contraction to survive the width in $a$.
2. **Plain-spiral tail.** For each $a$-piece, no plain zero has scaled lift $s$ above a certified $s_{\rm tail}$ (v394's monotone tail, rerun with the interval $a$). The tail starts at $s_{\rm tail}=1.3$ for every piece.
3. **Coarse localization.** Every plain zero with $s\in[\tfrac14,s_{\rm tail}]$ lies in $B^*=(1\pm0.15)^3$, for every $a$ in the piece. The search region in $(U,V)$ comes from the uniform side bound $Q_0$, which grows with $a$: 1.81, 2.71, 3.56 and 4.26 on the first four pieces, then 4.30 and 4.01 on $[2,3]$ and $[3,4]$. All 48 combinations (six pieces, eight patterns) are certified, using about 2.1 million boxes and roughly 9 hours of CPU.
   - The last fine piece $[2,4]$ is split in two here. Over all of $[2,4]$, the interval dependency on $a$ inflates $Q_0$ to 6.4, and the search region becomes too big to search.
4. **Regions** (v392/v395). v392's estimates were already written for general $a$ and accept an interval. The table below gives, per $a$-piece, the threshold $c$ above which the twisted label's region is certified and contains $B^*$.

Two snags in step 3, both about the remainder formula, which needs $1+\rho>0$ on the whole box:
- **The first search box.** For the larger regions (from $a\approx0.4$ up), the very first box is so big that its estimate of $1+\rho$ reaches zero, though the true value stays far from it. Such a box is simply split, like any box that can't yet be excluded.
- **Far corners at large $a$.** In $\rho=\varepsilon c(2u+\varepsilon c|w|^2)$, the parameter $\varepsilon c$ appears twice, and the search splits only $(s,U,V)$, not the parameters. So at far corners the estimate of $1+\rho$ can reach zero however small the box gets: at one corner the true value is about 0.2. The same quantity written as $|1+w|^2=(1+\varepsilon cu)^2+(\varepsilon cv)^2$ is never negative and much tighter. It is used only when the first form fails, so every earlier evaluation is unchanged.

## The table

Thresholds $c$ (smallest offset $m$), with gap $\delta=\pi-\Phi$ to contact. Each $a$-piece column holds for every exponent in that piece at once.

| $\Phi$ | $L$ | $L_1$ | $[\frac14,0.4]$ | $[0.4,0.7]$ | $[0.7,1.2]$ | $[1.2,2]$ | $[2,3]$ | $[3,4]$ | $a=\frac12$ (v394/v395) |
|---|---|---|---|---|---|---|---|---|---|
| 0.4 | 2 | 0 | 30 | 30 | 30 | 30 | 30 | 30 | 11 |
| 1.2 | 6 | 0 | 30 | 30 | 30 | 30 | 30 | 30 | 16 |
| $\pi-0.3$ | 5 | 0 | 54 | 69 | 88 | 103 | 91 | 75 | 32 |
| $\pi-0.1$ | 5 | 2 | 113 | 159 | 204 | 241 | 240 | 223 | 103 |

- The first two rows sit at the floor $m\ge30$ set by the degree step, not by the regions.
- Near contact, the thresholds are two to three times the single-exponent values. Most of that is the price of covering a whole interval of $a$ at once, not of the exponent itself: the $[0.4,0.7]$ column contains $a=\tfrac12$ and still needs 69 against 32. Narrower pieces would bring each column down.
- Thresholds peak for $a$ between about 1.2 and 3 and come down again on $[3,4]$.

## What this leaves of the reduction

After v395 and v396, a radially monotone wild point is covered if its shared winding follows *some* power law $r^{-1/a}$, $a\in[\tfrac14,4]$, up to a slowly varying angle, and its twist stays a fixed angle short of contact at a bounded rate. Remaining:

1. **Radial monotonicity.** Arms that swing back out in radius.
2. **Winding outside power laws**, for example $\log(1/r)$-type or $r^{-1/a}\log(1/r)$ winding. These differ from every power law by an angle of unbounded slope relative to that law, so neither v395 nor v396 covers them. Exponents outside $[\tfrac14,4]$ are a matter of running more pieces.
3. **Twist of unbounded slope**, and arms approaching contact.
4. **Non-winding and rough-everywhere points.**

**In short.** The winding exponent was a modelling choice, not a real restriction: the degree argument and the region certificates go through for every exponent from $\tfrac14$ to $4$, with the same limit system and the same degree signs.

[Code](../code/any_exponent_v396.py) · [Tests](../tests/test_any_exponent_v396.py) · [Data](../results/v396/exponent_check.json)
