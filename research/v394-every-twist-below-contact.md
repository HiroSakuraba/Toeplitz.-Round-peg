# v394. Every twist below contact, at every slope, and only where each square lives

The plan was to push v393 to twists that grow without limit. That case turns out not to exist. In the two-arm model, the twist is automatically trapped below π. So this note does the version that does exist: every twist short of contact, at any slope, with the bounds needed only on the annulus where each square lives. It also says how fast the arms may approach contact before the method runs out.

**Results.** These are computer-assisted, in interval arithmetic with outward rounding, for the power spiral with $a=1/2$ and level-one labels with comparable offsets ($|j|=1$, $M\le2m$), as in v393.

1. **The twist can't grow without limit (Lemma 394.1).** Each arm has exactly one point at each radius $r$. Arm one sits at polar angle $\Theta$ and arm two at $\Theta-\pi+\varphi(r)$, with $\Theta=r^{-1/a}$. The two meet exactly when $\varphi(r)\equiv\pi\pmod{2\pi}$. A Jordan curve can't meet itself and $\varphi$ is continuous, so $|\varphi|<\pi$ at one radius forces $|\varphi|<\pi$ at every radius.
   - "Unbounded twist" means the arms crossing, which can't happen.
   - What can happen is milder: **arms approaching contact** ($\sup|\varphi|=\pi$, never attained), and **fast twisting** (large $|r\varphi'|$).
2. **One degree certificate covers every twist bound (Theorem 394.2).** v393's fine certificates never involved the twist. Only its coarse localization did, through each region's shape. Here the coarse step is redone once, for the plain spiral:
   - v392's tail certificate with zero twist rules out plain-spiral squares with scaled lift $s\ge1.5$;
   - branch and bound puts every plain-spiral zero with $s\in[\tfrac14,1.5]$ into v393's box $B^*$.
   
   So every v392-style region with lower wall at $s\ge\tfrac14$, offsets at least 11 and $B^*$ inside it has degree $\pm1$, whatever twist bound it was built for. All 8 sign-and-corner patterns are certified.
3. **Thresholds up to contact (Theorem 394.3).** For every pair in a table running from $\Phi=\pi-1.5$ to $\pi-0.01$ and slopes $L$ from 0 to 50: every two-arm power spiral with $|\varphi|\le\Phi$ and $|r\varphi'|\le L$ has an inscribed square of **every** label with $m\ge c(\Phi,L)$ (and $m\ge11$), $M\le2m$, $|j|=1$. That is infinitely many squares, accumulating at the point. All 70 entries are certified. With $\delta=\pi-\Phi$ the gap left before contact, the thresholds follow $c\approx(L+3.7)/\delta$ closely once $\delta\le0.05$.
4. **Only local bounds are needed (Theorem 394.4).** Each label's certificates look at the twist only at corners of squares whose corner $A$ lies in that label's region. So $|\varphi|\le\Phi$ and $|r\varphi'|\le L$ need to hold only on the label's annulus,
$$\big[(s_{\rm hi}\Theta^*)^{-a}(1-0.3\sqrt2),\ (s_{\rm lo}\Theta^*)^{-a}(1+0.3\sqrt2)\big].$$
   A twist can be wild at one scale and tame at another, and the squares at the tame scales still exist.
5. **Examples.**
   - **Fast oscillation.** $\varphi=2.9\sin(10\log1/r)$ comes within 0.24 of contact, with slope up to 29. Every such label with $m\ge194$ is covered.
   - **Arms approaching contact.** $\varphi=\pi-0.5r^{0.2}$ is covered by local bounds at every sampled label from $m=69$ to $2\times10^6$.
   - **Continuations.** Plain squares at $(60,61,1)$ and $(-120,-130,1)$ were followed through the twists $2.0\sin(4\log1/r)$ and $2.6\sin(2\log1/r)$. They travel down to scaled lift 0.73 and 0.55, stay strictly inside their certified regions, and keep their degree.

## Why it works for every twist bound at once

v393 proved the degree in two steps:
1. a fine step: one nondegenerate zero inside a small box, for all big labels;
2. a coarse step: no other zero in v392's region.

The fine step only involves the plain spiral. The coarse step used each region's side bound $\hat Q$ and lift range $[s_{\rm lo},s_{\rm hi}]$, which depend on the twist bound. As the bound nears contact, $s_{\rm hi}\approx\pi/\delta$ grows without limit, so redoing the coarse step for every region would never end.

The fix is to localize the plain spiral's own zeros instead, once:
- **Above $s=1.5$:** v392's upper-face certificate, run with zero twist, already excludes every plain-spiral square with scaled lift $\ge1.5$. It covers the whole tail by monotonicity, uniformly for offsets at least 11.
- **Below:** v392's side bound with zero twist at $s=\tfrac14$ gives side at most $2.39\,|A|/m$. With that, branch and bound shows every plain zero with $s\in[\tfrac14,1.5]$ lies in $B^*=(1\pm0.15)^3$.
- **Excision:** a region of any twist bound contains no plain zeros outside $B^*$, provided its lower wall is at $s\ge\tfrac14$. The table enforces that by capping the lower wall's width.

## The table

$c(\Phi,L)$: every level-one, comparable-offset label with smaller offset $m\ge c$ is certified. Gap $\delta=\pi-\Phi$:

| $\delta$ \ $L$ | 0 | 2 | 5 | 10 | 20 | 30 | 50 |
|---|---|---|---|---|---|---|---|
| 1.5 | 15 | 16 | 18 | 21 | 27 | 34 | 51 |
| 1.0 | 20 | 20 | 20 | 23 | 34 | 45 | 69 |
| 0.5 | 26 | 26 | 27 | 35 | 58 | 76 | 124 |
| 0.3 | 30 | 30 | 32 | 51 | 85 | 126 | 185 |
| 0.2 | 42 | 43 | 43 | 65 | 115 | 168 | 278 |
| 0.1 | 59 | 59 | 66 | 114 | 225 | 318 | 527 |
| 0.05 | 89 | 89 | 115 | 213 | 424 | 634 | 1054 |
| 0.03 | 143 | 143 | 179 | 354 | 705 | 1055 | 1757 |
| 0.02 | 190 | 190 | 267 | 530 | 1056 | 1582 | 2635 |
| 0.01 | 374 | 374 | 530 | 1056 | 2110 | 3163 | 5269 |

- The minimum offset 11 from v393 applies on top of every entry. Entries never need fewer labels for harder bounds, since certification only gets harder as $\Phi$ or $L$ grows.
- **Why the $1/\delta$.** For an odd label the identity's left side is at least $\pi-|\varphi|-{\rm osc}-d\approx\delta-O(1/m)$. The error terms, of size $1/m$, must fit inside the gap.
- **Why the $L$.** The twist's variation across a square, ${\rm osc}\approx L\cdot\text{(side)}/|A|\approx L/(3m)$, eats into the same gap.

So $m\gtrsim(L+\text{const})/\delta$ is the natural shape, and the table matches it with constant about 3.7.

## How fast may the arms approach contact?

Say the gap left at radius $r$ is $g(r)=\pi-|\varphi(r)|$.
- **What a label needs.** Its region reaches inward to radius about $r_{\rm in}\approx0.13\,\delta^{1/2}/m$, because $s_{\rm hi}\approx\pi/\delta$. So it needs $g(r_{\rm in})\ge\delta$, with $\delta\approx(L+3.7)/m$.
- **For $g(r)=\kappa r^\beta$.** The condition reads $m^{1-3\beta/2}\gtrsim\text{const}$.
  - It holds for all large $m$ when $\beta<2/3$, for any $\kappa>0$.
  - It fails for this method when $\beta>2/3$.

This is a **heuristic**, because it extrapolates the table's scaling below $\delta=0.01$. The scans agree with it within the table's range:
- the $r^{0.2}$ example is certified at every sampled label up to $m=2\times10^6$;
- the $r^{0.5}$ example only qualifies beyond the table, so no sampled label is covered.

Making the contact-rate statement rigorous needs a certificate uniform in $\delta$, written in scaled variables much as v392 did for the offsets.

## Checks

| | result |
|---|---|
| plain-spiral tail at zero twist | certified for $s\ge1.5$ |
| plain-spiral side bound at $s=\tfrac14$ | $\hat Q_0=2.39$ |
| coarse localization, 8 patterns | certified: 16.6k–93k boxes per pattern for $j=+1$, 46k–135k for $j=-1$ |
| table entries certified, with the degree step applying ($s_{\rm lo}\ge\tfrac14$, $B^*$ inside) | 70 / 70 |
| examples followed to the full twist, strictly inside their region, degree kept | 2 / 2 |

For the scans, both contact profiles and their slopes are monotone in $r$, so their extremes on an annulus sit at its ends, and the sampled check is exact. The oscillating example has exact global bounds $|\varphi|\le2.9$ and $|r\varphi'|\le29$.

## What is and isn't proved

- **Proved (computer-assisted):** existence of infinitely many inscribed squares at every two-arm power-spiral point ($a=\tfrac12$) whose twist obeys one of the 70 tabulated bounds, at least on the annuli of infinitely many labels.
- **Proved (continuity, as in v392):** some threshold works for every bound with $\Phi<\pi$ and $L<\infty$. With Lemma 394.1, that covers every embedded two-arm power spiral whose arms stay a fixed angle apart and twist at a bounded rate.
- **Open:**
  - **arms approaching contact faster than about $r^{2/3}$;**
  - **twist with unbounded slope** (rough twisting);
  - a certificate uniform in $\delta$ that would turn the contact-rate heuristic into a theorem;
  - higher levels, lopsided labels, other exponents, and more than two arms.

[Code](../code/contact_twist_v394.py) · [Tests](../tests/test_contact_twist_v394.py) · [Data](../results/v394/contact_check.json)
