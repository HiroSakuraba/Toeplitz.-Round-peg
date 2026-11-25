# v390. Labels are conserved: why squares survive unless they collapse

The arc picture from the last discussion said that a counterexample at the end of a symmetrizing path needs squares to be swallowed by the wild point. Writing it down carefully, the parity version turns out not to work. The version that works uses the square's **label**, the passes its corners sit on, instead of the parity.

**The idea.**
- A square moving continuously along a family of curves keeps its label.
- For a region of square configurations defined by a label, the mod-2 degree of the square map cannot change unless a square crosses the region's boundary.

So a square of the symmetric curve can disappear only by leaving its label region: shrinking into the wild point, or drifting out of its band of radii. Nothing else can remove it. Not the twist's events, not folds, not partners.

**Results.**
1. **Parity bookkeeping alone is empty (Remark 390.1).** For any path from $\gamma_0$ to $\gamma_1$, the number of graceful events at cutoff $s$ is $L(\gamma_0,s)+L(\gamma_1,s)$ mod 2. It says nothing new.
   - It also cannot rule out a cascade: a chain of branches, each joining one square to a smaller square's half-turn partner, could carry an odd count down through every scale with no single branch collapsing.
   - The earlier sketch's "every arc joins two symmetric squares, so the count is even" overlooked that the global squares of $\gamma_0$, unlike the small ones near $p$, are not paired.
2. **Labels are conserved (Lemma 390.2).** Along a continuous family of curves near a two-arm point, a square's corners stay on the same passes, as long as no corner passes through the wild point. A pass is one sweep of an arm at polar angle $\theta$, with lift $\theta+\pi k-\varphi^*$, and it moves continuously with the twist.
3. **Label degree is invariant (Theorem 390.3).** Fix a label $\ell$ and a compact region $\Omega_\ell$ of configurations (corner $A$ in an annulus, any orientation, side in an interval), with corners on the passes $\ell$ prescribes. Suppose no curve $\gamma_\lambda$ of the path has a square with label $\ell$ on $\partial\Omega_\ell$. Then the degree of the four comb-phase equations on $\Omega_\ell$ is the same for every $\lambda$.
4. **No cascades (Corollary 390.4).** Half-turn partners carry different labels: the half-turn swaps the arms, so every pass parity flips. So no branch can join partners, and the chain of Remark 390.1 is impossible. Suppose the symmetric curve has exactly one nondegenerate square with label $\ell$ in $\Omega_\ell$, and the label's squares never reach $\partial\Omega_\ell$. Then **every curve on the path has a square with label $\ell$**. The only ways out are collapse and drift.
5. **Checks on the plain spiral:**
   - all 81 squares of a v380 search are nondegenerate zeros of their label's phase map (local degree $+1$ for 31, $-1$ for 50);
   - none has the same label as its half-turn partner;
   - an independent search with another seed finds 73 squares with 73 distinct labels, all on their v382 branch, so it gives one square per label.

## Remark 390.1: why parity bookkeeping cannot do it

Run the twist path $\gamma_\lambda$ and let $E(s)$ count graceful squares whose side passes $s$. The path principle (v373) gives $L(\gamma_1,s)\equiv L(\gamma_0,s)+E(s)$, and that is all the parity says.

The arc picture refines this. Each square traces an arc in (twist, configuration) space, and an arc can end in one of three places: at $\lambda=0$, at $\lambda=1$, or at the wild point. Counting arcs whose ends straddle $s$ gives $E(s)\equiv\#\{\text{ends above }s\}$. But the $\lambda=0$ ends include the large, global squares of $\gamma_0$, which the half-turn does not pair. Their number above $s$ is odd, which is exactly why $L(\gamma_0,s)=1$. So arcs joining $\lambda=0$ squares to each other can straddle every cutoff an odd number of times: a chain from a global square, through a small square $t_1$, its partner $t_1'$ (same size), a smaller square $t_2$, its partner, and so on. No single arc needs to collapse. Parity cannot exclude this. Labels can.

## Lemma 390.2 and Theorem 390.3

A square's corner on pass $k$ satisfies one comb-phase equation, $|z|^{-1/a}-\theta-\pi k+[k\text{ odd}]\lambda\varphi(|z|)=0$, with $\theta$ the polar angle unwrapped around the square's centre. Along a continuous family the four corners move continuously. A corner can reach a different pass only by passing through the wild point, where all passes meet, or through a point where two passes meet. On a two-arm point passes never meet away from $p$, because the gap $\pi\pm\varphi$ stays in $(0,2\pi)$. So, as long as the radii stay positive, the label is constant.

For a compact region $\Omega_\ell$ in the four configuration variables, the map $F_\lambda$ of the four phase residuals is continuous in $(\lambda,v)$. If $F_\lambda\ne0$ on $\partial\Omega_\ell$ for every $\lambda$, homotopy invariance gives $\deg(F_\lambda,\Omega_\ell)=\deg(F_0,\Omega_\ell)$. ∎

The degree uses no transversality and no genericity. Folds, births and deaths inside $\Omega_\ell$ are invisible to it. Only the boundary matters, and for an annulus-times-interval region the boundary means the square's corner $A$ at the inner radius (moving toward the point), at the outer radius (drifting away), or its side at an end of the interval.

## Corollary 390.4

On the plain spiral, v382's branches give one square per label, at each resonance. v387 checked that the three-corner solve has no second solution, and the checks below confirm one square per label in an independent sample. So $\deg(F_0,\Omega_\ell)=\pm1$ for a small enough region around that square. If the label's squares stay strictly inside $\Omega_\ell$ for every twist on the path, the degree stays $\pm1$, and every $\gamma_\lambda$, the twisted curve included, has a square with label $\ell$.

What remains is a confinement statement: for which labels does the label's square stay inside a fixed compact region, whatever the twist? That is v391.

## Checks

On the plain spiral, from v380's random search (6 000 seeds):

| check | result |
|---|---|
| squares that are nondegenerate zeros of their label map | 81 / 81 |
| local degree $+1$ / $-1$ | 31 / 50 |
| squares whose half-turn partner has the same label | 0 / 81 |
| independent search: squares, distinct labels, on their v382 branch | 73, 73, 73 |

The local degree's sign depends on the label's orientation conventions and is not meaningful by itself. What matters is that it is never 0.

## What this changes

- **The obstacle is now a single, local question:** can a label's square be pushed out of a compact region of radii? That means into the wild point, or away from it.
- **Parity and partners drop out.** The degree argument needs neither the sector rule nor the half-turn pairing, and it applies to crossed squares as much as graceful ones.
- **It replaces the path condition.** Turning off the twist may create squares of side $s$ anywhere (v389). That no longer matters, as long as some label's square cannot leave its region.

[Code](../code/label_degree_v390.py) · [Tests](../tests/test_label_degree_v390.py) · [Data](../results/v390/label_check.json)
