# v393. The starting degree is ±1, so the big squares exist outright

v392 showed that a bounded twist can never push a big odd label's square out of its certified region. Turning that into existence needed one more fact: on the symmetric (untwisted) spiral, the label's degree on that region must be nonzero. Until now that was the uniqueness assumption U0, supported only by Newton searches. This note proves it, uniformly, for every big level-one label whose two offsets are comparable.

**The idea.** Zoom in on the square at the scale its label sets. As the offsets grow, the label's equations, in those zoomed coordinates, converge to a simple explicit system with exactly one solution. A rigorous interval Newton test (Krawczyk) then pins down one solution for all big labels at once. A rigorous search shows the certified region holds no other solution. One nondegenerate solution means degree ±1.

**Results.** These are computer-assisted, in interval arithmetic with outward rounding, for the power spiral with $a=1/2$.

1. **Theorem 393.1 (nonzero degree).** Take a label with level $|j|=1$, smaller offset $m\ge N$ and larger offset $M\le2m$, in any compatible sign pattern, with either neighbour of $A$ carrying the smaller offset. On the plain spiral, its label map has exactly one zero in v392's certified region, and that zero is nondegenerate. So the degree there is $\pm1$.
   - Here $N=11$, 13, 16 for the regions certified at $(\Phi,L)=(0.4,2)$, $(0.8,4)$, $(1.2,6)$.
2. **Corollary 393.2 (unconditional existence).** Every two-arm point whose twist satisfies $|\varphi|\le\Phi$ and $|r\varphi'|\le L$, for one of those three pairs, has a square with every such label. The proof chains three steps:
   - v390: a label's degree can't change while no square reaches its region's boundary;
   - v392: no square reaches the boundary under any such twist;
   - Theorem 393.1: the starting degree is $\pm1$.
   
   In particular, **such a point carries infinitely many inscribed squares, accumulating at the wild point**: for example, one with label $(m,m,1)$ for every $m\ge N$. Their sides are about $|A|/(3m)$. U0 is no longer needed for these labels.
3. **Checks.**
   - 19 real plain-spiral squares (offsets 11 to 79, all eight sign-and-corner patterns) all lie inside the certified box for their parameters.
   - Their local degrees all agree with the certified sign: $-1$ for every $j=+1$ pattern and $+1$ for every $j=-1$ pattern.

## Zooming in

Write $\varepsilon=1/m$, $r=m/M\in[\tfrac12,1]$ and $\Theta^*=(2a+1)\pi mM/|j|$. The offsets are $n_1=\sigma_1k_1$ and $n_3=\sigma_3k_3$, where $k_1,k_3$ are the two offset sizes and $t_i=k_i/M$. New unknowns $(s,U,V)$ come from
$$\Theta_A=s\,\Theta^*,\qquad \operatorname{Re}w=-\varepsilon\sigma_1\tfrac{a|j|}{2a+1}t_1U,\qquad \operatorname{Im}w=\varepsilon\sigma_3\tfrac{a|j|}{2a+1}t_3V .$$
These are the leading-order guesses with correction factors $U$, $V$ and $s$. On the plain spiral, v390's residuals are $R_X=\Theta_A(|1+wm_X|^{-1/a}-1)-\arg(1+wm_X)-\pi o_X$, and A's arm drops out. The scaled map is
$$\Psi=\Big(\frac{R_B}{\sigma_1\pi k_1},\ \frac{R_D}{\sigma_3\pi k_3},\ \frac{R_C-R_B-R_D}{\pi|j|}\Big)=\big(sU-1,\ sV-1,\ \sigma_1\sigma_3\,sUV-\operatorname{sign}j\big)+O(\varepsilon/r).$$
- The first two components say each neighbour's radial step matches its pass offset.
- The third is v386's identity, the small mismatch $\pi j$ balanced by the square's curvature.
- For compatible signs ($\operatorname{sign}j=\sigma_1\sigma_3$) the limit system has the single zero $(1,1,1)$, with Jacobian determinant $-\sigma_1\sigma_3$.

$\Psi$ is the label map composed with an invertible diagonal change of variables and an invertible linear recombination of its components. Degrees therefore agree up to a fixed sign.

**Making the remainders computable.** Every remainder is written in closed form with $\varepsilon$ factored out, so interval arithmetic can evaluate it on $\varepsilon\in[0,1/N]$, including $\varepsilon=0$.
- For $a=\tfrac12$: $|1+w|^{-2}-1+2\operatorname{Re}w=\varepsilon^2\hat h(\omega,\varepsilon)$, where $w=\varepsilon\omega$.
- The angle defect is $\arg\big(1-\tfrac{iw^2}{(1+w)(1+iw)}\big)$, which is $O(\varepsilon^2)$.
- One step had to be done by hand. The quadratic parts in the third component cancel exactly to $-8xy$, and $-8xy\cdot(2a+1)s/(j^2r)=\sigma_1\sigma_3sUV$, with $\hat h=(3x^2-y^2)+\varepsilon\hat h_1$. Left to interval arithmetic, the three nearly equal terms lose that cancellation, and the enclosure is useless.

## The certificates

- **Krawczyk (existence, uniqueness, sign).**
  - For a small parameter box (a range of $\varepsilon$ and $r$), take $y$ near that box's zero, $C\approx J(y)^{-1}$, and a box $X$ around $y$.
  - If $K(X)=y-C\Psi(y)+(I-CJ(X))(X-y)$ lies in the interior of $X$, then for every parameter value $\Psi$ has exactly one zero in $X$. Every Jacobian in $J(X)$ is nonsingular, with $\det J$ of the sign of $\det C$.
  - The Jacobian enclosures come from forward-mode automatic differentiation over intervals.
  - Parameter boxes are bisected until the test passes, at box half-widths of 0.015–0.06.
- **Fine localization.** Every zero in the coarse box $B^*=(1\pm0.15)^3$ lies in that small box. This is a branch-and-bound search: each sub-box outside the target has some component of $\Psi$ bounded away from 0.
- **Coarse localization.** Every zero in v392's twist-0 region lies in $B^*$, by the same search.
  - The region means $s\in[s_{\rm lo},s_{\rm hi}]$, with v392's side bound giving $|\omega|\le|j|\hat Q$ and hence $|U|\le4\hat Q/t_1$ and $|V|\le4\hat Q/t_3$.
  - This step depends on the twist size, through v392's region. The fine step doesn't: it involves only the plain spiral, so one run with $\varepsilon\le1/11$ serves all three twist sizes.
- **Excision.** $B^*$, mapped back to $(\Theta_A,w)$, lies inside v392's region: $s\in[0.85,1.15]$ is inside the walls, and its side ratio lies between the inner face and $0.3$. So the label map's degree on the region equals its degree on $B^*$, which is $\pm1$.

## Checks

| | $j=+1$ | $j=-1$ |
|---|---|---|
| sign patterns × which corner has the smaller offset | 4 | 4 |
| fine pieces (Krawczyk + local localization), interval boxes | 691, 223 223 | 162, 45 424 |
| certified sign of $\det J$ | $-1$ in all four | $+1$ in all four |
| coarse localization boxes at $\Phi=0.4$ / 0.8 / 1.2 | 34 481 / 49 350 / 50 434 | 34 726 / 24 118 / 27 314 |
| plain squares inside their certified box, local degree as certified | 11 / 11 | 8 / 8 |

The plain squares sit at $s\approx0.96$–$1.05$ and $U,V\approx0.94$–$1.07$, closer to $(1,1,1)$ as the offsets grow. For example, $(60,60,1)$ sits at $(1.008,0.988,0.988)$.

## What is and isn't proved

- **Proved (computer-assisted):** nonzero degree, and hence unconditional existence under the three twist bounds, for $|j|=1$ labels with $m\ge11,13,16$ and $M\le2m$. That is infinitely many labels per point.
- **Not covered:**
  - **Very lopsided labels.** When $M$ is much larger than $m^2$, the small offset's radial step is dominated by curvature, and the limit system changes shape. These labels still have certified regions (v392), but their degree isn't computed here. Extending the ratio beyond 2 is straightforward up to $M\sim m$; beyond that it needs a second scaling.
  - **Higher levels** $|j|\ge3$: the same method should apply, from $m\ge11|j|$.
  - **Other exponents** $a$: the closed forms here are written for $a=\tfrac12$.
  - **Small or even labels**, and anything outside the two-arm power-spiral model with bounded twist.

**Where this leaves the program.** For the class of wild points this line of work set out to understand, two power-spiral arms twisted by any smooth bounded amount within the certified bounds, the point is not a counterexample. It inscribes infinitely many squares accumulating at it, with no unproved assumption left in the chain. The proof is computer-assisted: interval certificates in v391–v393, resting on v386's identity and v390's degree argument. The natural next questions are:
- which other twist bounds and exponents follow by rerunning the same certificates;
- whether "two arms, bounded twist" can be widened toward the general accumulating wild points of v384.

[Code](../code/uniform_degree_v393.py) · [Tests](../tests/test_uniform_degree_v393.py) · [Data](../results/v393/degree_check.json)
