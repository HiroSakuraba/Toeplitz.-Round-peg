# v380. Comb laws for any arms

v375–v379 worked with one shape of wild point: two arms winding in as $r=\Theta^{-a}$, only ever moving inward. This note drops the shape assumption. It proves a comb law for any arc, smooth or polygonal, with nothing assumed about how the arms wind. It then tests that law, and v379's selection rule, on spirals whose radius wobbles, including wobbles strong enough to make the arms back up.

The law is simple. A short stretch of arc is too straight to hold three corners of a square, so some side of every square must jump across a quarter-turn of the arc. So the side is at least the distance from a corner to the nearest point a quarter-turn away. The tests then show exactly when this law is attained and when squares sit much deeper. The deciding feature is whether the gap between neighbouring turns has stationary points. That one criterion also explains the polygon kink squares of v379.

One correction to the last list of directions. A point of a Jordan curve always has exactly two arms, so "points where more than two arms meet" do not exist. The nearest real case is wild points accumulating at a limit point; that is left for later.

**Results.**
1. **General comb law (Theorem 380.1).** For any $C^1$ or polygonal arc, every inscribed square has at least two sides whose corners are joined along the arc by a subarc turning through at least $\pi/4$. Hence $t\ge\min_{\text{corners}}G(z)$, where $G(z)$ is the distance from $z$ to the nearest arc point a quarter-turn away along the arc.
2. **Aligned squares sit at stationary gaps (Principle 380.2, numerical).** Squares with exactly two such sides straddle two neighbouring turns. They sit where the gap between those turns is stationary, with side equal to that gap. Where the gap is monotone they do not exist.
3. **Tests on a wobbling family** $r=\Theta^{-1/2}\big(1+\frac\kappa\Theta\sin\omega\Theta\big)$, which is Jordan for $\kappa<\pi/4$ and backs up radially when $\kappa\omega>1/2$:
   - Theorem 380.1 holds on all 1 867 squares found, and on every small square of four polygon drawings;
   - without stationary gaps (plain spiral; $\kappa=0.2,\omega=2$) every square is tilted and v379's smooth law holds;
   - with stationary gaps ($\kappa,\omega=0.3,1$; $0.5,3$; $0.7,2$) aligned squares appear and come within 0–13% of the general law ($t/G$ as low as $1.000$, with $G$ evaluated on 60 squares per family). They reach 2.2–3.6 times v379's smooth radius;
   - every aligned square sits at a stationary point of its gap ($|g'|\le0.6\%$ of the turn's maximum, median $0.01\%$), with side equal to the local gap within 4%.

## The general comb law

Let $\alpha$ be an arc, either $C^1$ or polygonal. For a subarc $\sigma$, let ${\rm turn}(\sigma)$ be the range (max minus min) of the lifted tangent angle along $\sigma$. For a polygon, take the edge directions. A subarc through a wild point at the end of both arms turns without bound.

**Lemma 380.0.** If ${\rm turn}(\sigma)<\pi/2$, no two chords of $\sigma$ are perpendicular.

*Proof.* A chord is the integral of the unit tangent (or a positive combination of edge vectors), so its direction lies in the arc of tangent directions. Two directions in an arc shorter than $\pi/2$ cannot differ by $\pi/2$. ∎

**Theorem 380.1.** Let a square have its four corners on $\alpha$. Then at least two of its sides join corners $u,v$ with ${\rm turn}(\alpha[u,v])\ge\pi/4$. In particular
$$t\ \ge\ \min_{z\ \text{a corner}}G(z),\qquad G(z)=\inf\{|w-z|:\ w\in\alpha,\ {\rm turn}(\alpha[z,w])\ge\pi/4\}.$$

*Proof.* List the corners in arc order $w_1<w_2<w_3<w_4$ and call a pair **short** if its subarc turns by less than $\pi/4$.

A short pair $(w_i,w_k)$ with $k\ge i+2$ has a subarc containing three corners. Three corners of a square include two perpendicular sides, which Lemma 380.0 forbids.

If $(w_1,w_2)$ and $(w_2,w_3)$ were both short, then ${\rm turn}(\alpha[w_1,w_3])<\pi/2$, with three corners inside. That is the same contradiction, and likewise for $(w_2,w_3),(w_3,w_4)$.

So the short pairs are among $(w_1,w_2),(w_3,w_4)$. At most two of the four sides are short. ∎

For the spirals of v375 this recovers the comb law up to a constant: the nearest point a quarter-turn away is on the next turn. But nothing here uses spirals, radial monotonicity, or even smoothness. It applies to polygonal arms directly, without v376's sag correction. The constant $\pi/4$ is the best this argument gives, since two short links must add up to less than $\pi/2$.

## Testing it: wobbling and backing-up spirals

The family is half-turn symmetric: arm one is $r(\Theta)e^{i\Theta}$ with $r=\Theta^{-1/2}(1+\frac\kappa\Theta\sin\omega\Theta)$, and arm two is its negative. Neighbouring turns stay ordered, so the arc is Jordan, when $\kappa<a\pi/2=\pi/4$. The radius backs up somewhere on every turn when $\kappa\omega>a=1/2$. Inscribed squares were found by Newton's method, with each corner assigned to the turn nearest its seed: 6 000 seeds per family, with corners restricted to radii up to $0.98\,\theta_0^{-1/2}$.

| $\kappa,\omega$ | radius backs up | stationary gaps ($\Theta\in[6,120]$) | squares | sides crossing a quarter-turn: 2 / 3 / 4 | max centre / smooth radius | min $t/G$ |
|---|---|---|---|---|---|---|
| 0, – | no | 0 | 73 | 0 / 0 / 73 | 1.007 | 3.95 |
| 0.2, 2 | no | 0 | 70 | 0 / 0 / 70 | 1.010 | 3.92 |
| 0.3, 1 | no | 36 | 396 | 62 / 0 / 334 | 2.54 | 1.001 |
| 0.7, 2 | yes | 73 | 353 | 203 / 2 / 148 | 2.16 | 1.125 |
| 0.5, 3 | yes | 109 | 975 | 243 / 15 / 717 | 3.55 | 1.000 |

The theorem never fails: no square has fewer than two quarter-turn sides. The values of $t/G$ just below 1 are within the sampling error of $G$. On the four polygon drawings of v374 ($\lambda=0,1$; $s=1/20,1/30$) every small square has two or three quarter-turn sides, mostly two: these are the kink squares of v379.

## When the law is attained: stationary gaps

The two-side squares are **aligned**. Two corners sit on each of two neighbouring turns, and the square's sides run along and across them. For two nearly parallel arcs with gap $g$, such a square fits where the arcs are locally parallel at distance $t$, that is, where $g'=0$, up to a correction from the spiral's pitch.

**Principle 380.2 (numerical).** An aligned square straddling neighbouring turns sits at a stationary point of the gap $g(\psi)=r(\psi+L_1)-r(\psi+L_2)$ between them, and its side equals the gap there. If the gap is monotone along the comb, there are no aligned squares, and every square is tilted and governed by v379's selection rule.

The evidence: across the three families with stationary gaps, all 508 aligned squares have $|g'(\psi_c)|$ at most $0.6\%$ of the largest $|g'|$ over the surrounding turn (median $0.01\%$), and $t/g(\psi_c)\in[0.98,1.04]$. The two families without stationary gaps have no aligned squares at all.

This explains several earlier findings at once:
- **The plain spiral (v379).** Its gap shrinks steadily inward, so there are no aligned squares, and squares sit deep at the resonant radii.
- **Even versus odd wobble.** With $\omega$ odd, neighbouring turns (which are half a turn apart) wobble out of phase, so the gap oscillates strongly and has stationary points. With $\omega=2$ they wobble in phase and the gap stays monotone for small $\kappa$. For larger $\kappa$ the $1/\Theta$ factor still produces stationary points, as at $\kappa=0.7,\omega=2$.
- **Polygons (v379).** At a vertex the gap function has a corner, which acts as a stationary point for every slope in between. So aligned kink squares appear at every vertex pair across neighbouring turns.
- **Backing up is not what matters.** The arm backing up is neither needed nor decisive. The gentle wobble $\kappa=0.3,\omega=1$ never backs up but has aligned squares. What decides is the gap between neighbouring turns, not the motion of a single arm.

## What this changes

- **The comb law no longer depends on the arm shape.** Theorem 380.1 gives, for any two-arm wild point, a confinement scale set by the quarter-turn gap $G$. It is attained when neighbouring turns have stationary gaps, and it is loose (by the v379 mechanism) when they do not.
- **Where the path condition must reach depends on the gap profile.** The radius that the v374 path condition must reach is set by the confinement scale of the symmetrized curve. For cores with stationary gaps, including all polygons, that is the $G$-scale. For smooth cores with monotone gaps it is the deeper resonant scale of v379.
- **Still open.** Principle 380.2 is numerical. Proving it needs the same broken-family (implicit function) argument as the existence half of v379. The selection rule of v379 should carry over to any smooth phase with lattice levels, since the identity $\pi j=\iint p^{\mathsf T}\nabla^2\Phi\,q$ holds for any $C^2$ phase. Wild points accumulating at a limit point remain untreated.

[Code](../code/general_comb_v380.py) · [Tests](../tests/test_general_comb_v380.py) · [Data](../results/v380/general_comb_check.json)
