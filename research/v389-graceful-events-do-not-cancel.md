# v389. Graceful events along the symmetrizing path do not cancel

This is the third of three follow-ups. The open case of the path condition is a point twisted at every scale. Symmetrizing it, by turning its twist down from $\lambda=1$ to $\lambda=0$, makes squares of side exactly $s$ appear and disappear. The question raised in v385 and v386 was whether those events must cancel in pairs. This note counts them on smooth arms, for v374's twist $\varphi=0.4\sin(5\log(1/r))$. The answer is no.

The events do not cancel, and they are not avoidable either:
- each square's side barely moves along the whole path;
- but there are thousands of squares in every side window;
- so at a typical cutoff some twenty graceful squares cross it, and the parity of that count behaves like a coin flip.

**Results.**
1. **Bookkeeping (from v373 and v386).** Let $\gamma_\lambda$ be the twisted smooth spiral with twist $\lambda\varphi$. Since $L(\gamma_0,s)=1$ for the symmetric point,
$$L(\gamma_1,s)\equiv1+\#\{\text{graceful squares whose side passes }s\text{ as }\lambda\text{ runs from 0 to 1}\}\pmod2 .$$
2. **Partners split.**
   - At $\lambda=0$, graceful squares come in half-turn pairs with equal sides (v370).
   - As the twist turns on, the two partners move independently.
   - A pair whose legs both survive to $\lambda=1$ and move monotonically contributes to the parity exactly when $s$ lies between their two end sides.
   - So cancellation would need partners to move together. They don't: median split $0.07\%$ and $0.4\%$ of the side, at $s_0=0.01$ and $0.02$, maximum $11\%$ and $26\%$.
3. **Sides barely move.** Along the whole path, a square's side changes by a median $0.05\%$ ($s_0=0.01$) and $0.2\%$ ($s_0=0.02$), and at most $6\%$ and $18\%$. On smooth arms, rotating arm two mostly slides it along itself. This is unlike v374's polygon drawings, where kink squares at the comb radius made the events.
4. **But there are very many squares.**
   - Two independent searches with 1 500 seeds each, for squares at $\lambda=0$ with side in $[0.9,1.6]\,s_0$, overlap in **0** of 65 and 52 pairs at $s_0=0.01$, and in **5** of 86 and 84 at $s_0=0.02$.
   - Directly observed: 117 distinct pairs (71 with both legs complete) at $s_0=0.01$, and 165 (129 complete) at $s_0=0.02$.
   - A capture–recapture *heuristic* puts the population at about 1 400 pairs at $s_0=0.02$. It assumes every pair is equally likely to be found, which random Newton searches need not satisfy. At $s_0=0.01$ there is no overlap, so it gives no estimate at all. The product of the sample sizes, about 3 400, is a scale, not a lower bound.
   - With the measured moves, the heuristic population gives roughly 20–28 graceful events at a typical cutoff.
5. **The parity looks random, in the sample.** Within the combined census at $s_0=0.02$ (165 pairs, 58 graceful), $51\%$ of cutoffs see an odd number of graceful events and only $23\%$ see none. This is a property of the sample. A sample's parity says nothing certain about the full population's, and discrete continuation can miss crossings between samples.
6. **Partners always share a sector.** All 282 found pairs have both partners graceful or both crossed, as v370's pairing requires. This check exposed a tie bug in the sector rule, fixed in v378's and v386's code: float keys like $10^{18}-\Theta$ round to $10^{18}$. v378's tests and v386's results are unchanged. v388's graceful counts moved slightly and were updated there.

## Setting

The arms are $\Theta^{-a}e^{i\Theta}$ and $\Theta^{-a}e^{i(\Theta+\pi+\lambda\varphi(r))}$ with $a=1/2$, smooth everywhere. The path turns the twist down uniformly. $\gamma_0$ is the symmetric spiral, whose large parity is 1 at every small cutoff (Theorem 386.4). By the path principle (v373), $L$ changes along the path exactly when a graceful square crosses side $s$ transversally. Crossed squares never matter.

**Census method.**
- Seeds are placed on both arms at radii between $1.2\times$ the band's low end and $1.3\,r_{\rm smooth}$ of its high end, with random orientation and side in the band. Corners go to their nearest passes, and Newton's method solves.
- Each square found at $\lambda=0$ is completed by its exact half-turn partner.
- Each leg is continued in $\lambda$ over $[0,1]$ in 40 steps, with step halving, passes held fixed.
- Crossings of side $s$ are counted on the stored samples.

A leg that stops before $\lambda=1$ ended at a fold or where continuation failed. 46 of 117 pairs at $s_0=0.01$ and 36 of 165 at $s_0=0.02$ have such a leg. Only pairs with both legs complete enter the event counts.

## Why cancellation fails

Fix a pair with common side $t_0$ at $\lambda=0$ and end sides $t_1,t_2$ at $\lambda=1$, both legs monotone. Leg $i$ crosses $s$ iff $s$ lies between $t_0$ and $t_i$. So the pair contributes:
- $0$, if $s$ is outside both intervals, or inside both (the legs moved the same way past $s$);
- $1$, if $s$ lies between $t_1$ and $t_2$.

Pairing at $\lambda=0$ is therefore no help once the twist splits the partners. The parity becomes the number of partner windows $(t_1,t_2)$ containing $s$. These windows are tiny, but there are on the order of a thousand of them per window of sides, so they overlap many times. Fold-born pairs, which appear together at one side and then separate, behave the same way.

| $s_0$ | pairs found (A, B, overlap) | distinct pairs observed (complete) | population (heuristic) | graceful share | median move | expected graceful events at a cutoff (heuristic) | odd share in census |
|---|---|---|---|---|---|---|---|
| 0.01 | 65, 52, 0 | 117 (71) | none (no overlap; size-product scale 3 380) | 35 / 71 | 0.05% | ≈ 28 | 16% (census very incomplete) |
| 0.02 | 86, 84, 5 | 165 (129) | ≈ 1 400 | 58 / 129 | 0.2% | ≈ 20 | 51% |

The expected count is (population) × (graceful share) × 2 × (mean graceful move) / (band width), with the heuristic population. Pairs are matched between runs by their side rounded to $10^{-9}$. That can merge two squares but never split one, so the overlaps are upper bounds on the true overlaps. The mean move is much larger than the median, because a few legs move by 10–20%. With some twenty roughly independent windows at each cutoff, the parity is close to a fair coin.

## What this means for the path condition

- **v374's route cannot work by counting events locally,** at least for twists like v374's. The symmetrizing path passes squares of side $s$ at almost every cutoff, and nothing pairs them. So condition (b) of Theorem 374.3 fails generically. And "events cancel" is not a theorem waiting to be proved: it is false here.
- **The parity is still odd about half the time.** So for this twisted point, $L(\gamma_1,s)=1$ at roughly half of all cutoffs, and any one of them gives a square. But certifying a specific cutoff would need a complete census, and that is out of reach: thousands of squares per window, accumulating toward the point.
- **The obstacle is now sharp.** v385–v388 settled every case in which the asymmetry leaves a symmetric window, or turns on slowly beyond one. For asymmetry at every scale, protection has to come from something other than a path to a symmetric point. That could be a global invariant that sees all these crossings at once, or an argument that $L(\gamma,s)$ cannot vanish at every cutoff simultaneously.

## Caveats

The census is a sample: population estimates come from two runs, and partial legs are excluded, which can bias it. The population and event-rate figures are heuristics. What the census shows directly is evidence *against* pairwise cancellation: partners split, and many distinct pairs exist. That is not a proof that the complete event count fails to cancel. Only the twist $0.4\sin(5\log(1/r))$ was tested.

## Review corrections

An external review (F7) noted that the capture–recapture figures were presented too strongly. The output fields are renamed accordingly: `lincoln_petersen_heuristic`, `population_scale_heuristic`, `zero_overlap_scale_heuristic` and `estimated_events_per_cutoff_heuristic`, with the values unchanged. The misnamed `population_lower_bound` is dropped, and the directly observed pair counts with completeness are reported. The conclusions are restated as evidence against pairwise cancellation.

[Code](../code/event_census_v389.py) · [Tests](../tests/test_event_census_v389.py) · [Data](../results/v389/event_census.json)
