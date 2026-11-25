# v391. Odd squares can't collapse, and big ones are pinned

v390 reduced survival to one question: can a label's square be pushed out of its region, either into the wild point or away from it? This note answers it for squares with an **odd number of corners on arm two** ("odd squares"), using the exact identity from v386.

**The idea.** For an odd square the identity's left side contains $\pi j$ with $j$ odd. A twist of size $\Phi<\pi$ can shift it, but it cannot cancel it. The right side is roughly $(2a+1)\,\pi^2 n_1n_3/\Theta_A$. So the lift $\Theta_A$ of corner $A$ (which grows as the square moves inward) is held near the resonance $(2a+1)\pi|n_1n_3|/|j|$, whatever the twist does.

**Results.**
1. **No collapse (Theorem 391.1).** For every odd label and every twist with $|\varphi|\le\Phi<\pi$ and $|r\varphi'|\le L$, the label's small squares (side at most $0.3|A|$) have $\Theta_A\le\Theta_{\rm hi}(\ell,\Phi,L)$, a finite bound. So no odd square can shrink into the wild point under any bounded twist.
   - This is proved for every label, because the exclusion only improves as $\Theta_A\to\infty$.
   - The bound $\Theta_{\rm hi}$ itself is computed. It is found for all 3 200 labels with offsets up to 20 at each of $\Phi=0.4,0.8,1.2$.
2. **Big odd squares are pinned (Theorem 391.2).** When both pass offsets are large enough, the identity also excludes a lift below the band. The band is then a region whose boundary no square of the label can touch, for any twist up to that size. By v390 the label's degree is preserved, so every such curve has a square with that label.
   - With $|j|=1$ and $L=5\Phi$, every label with $\min(|n_1|,|n_3|)\ge N$ closes, where $N=10$, 12, 14 for $\Phi=0.4$, 0.8, 1.2 (checked up to offsets of 20).
   - Spot checks out to offsets of 200, including $|j|=3$, all close.
   - Conditional on uniqueness on the plain spiral (U0, supported by v390's checks), **every bounded-twist two-arm point carries infinitely many squares, accumulating at the wild point.**
3. **Small odd squares survive numerically.**
   - Setup: all 156 odd squares from v382 with side at most $0.2|Z|$, followed along the twists $m\cdot0.4\sin(5\log1/r)$ for $m=1,2,3$.
   - All 156 reach the full twist at each $m$, and all keep their local degree.
   - The largest lift any square reaches is 65–78% of its collapse bound.
   - The side bound holds for every square (worst case 0.89 of the bound).

## The pieces

**The identity, with the right ranges.** Write $C'=\Theta_B+\Theta_D-\Theta_A$. Then
$$\Theta_C-C'=\pi j-\Sigma^*+\Delta\theta=(2a+1)\,\zeta^{2a+1}\xi^{-2a-2}\,\delta_B\delta_D,$$
with $\delta_B=\Theta_B-\Theta_A$ and $\delta_D=\Theta_D-\Theta_A$.
- $\xi$ comes from the second difference of $\Theta^{-2a}$ over $A$, $B$, $D$, $C'$, so it lies in their hull.
- $\zeta$ comes from a first difference, so it lies between $\Theta_C$ and $C'$.

v386 stated the cruder range $[2\Theta_{\rm lo}-\Theta_{\rm hi},\,2\Theta_{\rm hi}-\Theta_{\rm lo}]$ for both. The tight ranges are what make lower walls possible at all.

**The terms.** Each corner's lift is $\Theta_A+\pi o_X$ up to $e=W+\max(\Phi,{\rm osc})$, where:
- $o=(0,n_1,n_1+n_3+j,n_3)$ are the pass offsets;
- $W$ is the window width;
- the twist difference between two corners is at most $\Phi$ if one is on arm one, and at most ${\rm osc}=2L\operatorname{atanh}u$ if both are on arm two.

For an odd label, $|\Sigma^*|\le\Phi+{\rm osc}$: with one arm-two corner it is one twist value, and with three it is one value plus a difference. Together with $|\Delta\theta|\le2\arcsin(u^2)$, this gives $|\pi j-\Sigma^*+\Delta\theta|\ge\pi-\Phi-{\rm osc}-d>0$.

**A bound on the side that is actually a bound.** All of $W$, osc and $d$ depend on the side. Because $AB\perp AD$ with $|AB|=|AD|=t$, the radial components satisfy $t^2=\langle AB,\hat u\rangle^2+\langle AD,\hat u\rangle^2$ *exactly*. Each component is $|B|-|A|$ (controlled by the lift step) minus a curvature term at most $(1+q)(1-\sqrt{1-q^2})|A|$, with $q=t/|A|$. This gives $q\le G(q)$ for an increasing map $G$.
- The code finds $G$'s least fixed point and steps slightly above it to $\hat q$ with $G(\hat q)\le\hat q$.
- It then checks $G(q)<q$ on $[\hat q,0.3]$, cell by cell.

So every square of the label with side at most $0.3|A|$ has $q\le\hat q$, and none sits exactly at $0.3|A|$. Where any step fails the lift is reported as *unknown*, never as excluded.

My first version used a first-order side estimate with a 10% margin. Real squares beat it by up to 38%. The rigorous version holds for all of them, with the worst case at 0.89.

**Walls.** Each lift is classed as excluded by the identity, allowed, or unknown.
- **Upper wall:** every lift above $\Theta_{\rm hi}$ is excluded. It always exists: as $\Theta_A\to\infty$, $\hat q$, osc, $d\to0$ and $e\to\Phi$. So the left side stays at least $(\pi|j|-\Phi)/\big((\pi|n_1|+\Phi)(\pi|n_3|+\Phi)\big)$ relative to $\delta_B\delta_D$, while the right side's factor tends to 0.
- **Lower wall:** the lift just below the band is excluded, so the band can't be left by drifting outward either.

The region of v390 is $\Omega_\ell=\{\Theta_A\in[\Theta_{\rm lo},\Theta_{\rm hi}]\}\times S^1\times\{0<t\le0.3|A|\}$. The orientation has no boundary, and the side can't reach either end (a label with $n_1\ne0$ has $t\ge\big||B|-|A|\big|>0$). So a closed band means no square of the label ever touches $\partial\Omega_\ell$. Bands only grow with $\Phi$, so the band at the full twist covers the whole path.

## Why small labels don't close

The identity sees $\Theta_A$ only through intermediate lifts spread over about $\pi(|n_1|+|n_3|)$. Next to $\Theta_A\approx2\pi|n_1n_3|$, that is a relative uncertainty of order $1/\min(|n_1|,|n_3|)$, raised to the power $2a+2$. For offsets below about 10, this washes out the lower exclusion before the side bound gives out.

The upper wall needs no such precision, so collapse is ruled out for every label. Drifting outward is ruled out only for big labels. For the small labels the numerics show no drift: side changes stay below 37% at $m=1$, and the squares stay well inside their region. But that is evidence, not proof.

## Checks

| | $m=1$ | $m=2$ | $m=3$ |
|---|---|---|---|
| twist size $\Phi$, slope bound $L$ | 0.4, 2 | 0.8, 4 | 1.2, 6 |
| odd squares followed | 156 | 156 | 156 |
| reach the full twist | 156 | 156 | 156 |
| keep their local degree | 156 | 156 | 156 |
| stay below the collapse bound (max lift / $\Theta_{\rm hi}$) | 156 (0.78) | 156 (0.73) | 156 (0.65) |
| side / proven side bound, worst | 0.89 | 0.88 | 0.87 |
| largest side change along the path | 37% | 88% | 159% |
| in a closed band | 4 | 0 | 0 |

At $m=1$, the four squares in closed bands have labels $(8,8,1)$, $(-8,8,-1)$, $(8,-8,-1)$ and $(7,7,1)$. v382's sample stops at offsets of 8, so the closed labels sit just beyond most of it.

Closure table ($|j|=1$, all sign patterns, $L=5\Phi$, offsets 1–20):

| $\Phi$ | labels with an upper wall | offset pairs whose band closes | every pair with $\min\ge N$ closes |
|---|---|---|---|
| 0.4 | 3200 / 3200 | 165 / 400 | $N=10$ |
| 0.8 | 3200 / 3200 | 115 / 400 | $N=12$ |
| 1.2 | 3200 / 3200 | 71 / 400 | $N=14$ |

Spot checks at $(25,25,1)$, $(30,40,1)$, $(50,50,1)$, $(40,120,1)$, $(100,100,1)$, $(25,200,1)$, $(60,60,3)$ and $(90,150,3)$ close at all three twist sizes. An example band: $(12,-13,-1)$ is pinned to $\Theta_A\in[683,1290]$ at $\Phi=0.4$ and $[409,1932]$ at $\Phi=1.2$, around the resonance 980.

## A numerical fix along the way

v383's `window` takes the polar angle of the square's centre from `np.angle`, which jumps by $2\pi$ when the centre crosses the negative real axis. Inside a Newton step this shifts all four residuals by exactly $2\pi$. In the first run, two thirds of the continuations stalled for this reason alone. v391 unwraps the angle continuously from the previous step.

v390's checks are unaffected: they only evaluate at found squares, with perturbations far too small to cross the axis.

## What is and isn't proved

- **Proved, given v386's identity bounds:** no odd square collapses under a bounded twist (Theorem 391.1), and the closed bands of Theorem 391.2 are never left.
- **Conditional:** existence at each closed label needs the label's degree on the plain spiral to be $\pm1$, which is uniqueness U0 (supported numerically: 73/73 in v390). It also needs the plain square to sit in the band, which it does near the resonance.
- **Computed, not proved for all labels:** that *every* big enough label closes. The table goes to offsets of 20 and the spot checks to 200. Every error term shrinks like $1/\min(|n_1|,|n_3|)$, so an asymptotic proof looks routine, but it isn't written.
- **Open:** outward drift of small odd squares, and anything about even squares. For even squares $j$ is even, $\pi j$ can be 0, and the identity forces nothing.

**Where this leaves the program.** A bounded twist cannot destroy the odd squares near a smooth two-arm wild point. The big ones provably persist, so a bounded-twist two-arm point always has infinitely many small squares (given U0). A counterexample of that type would need the twist to be unbounded, or the curve to leave the two-arm class.

[Code](../code/pinned_squares_v391.py) · [Tests](../tests/test_pinned_squares_v391.py) · [Data](../results/v391/pinned_check.json)
