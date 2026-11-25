# v390. Labels are conserved: why squares survive unless they collapse

The arc picture from the last discussion said that a counterexample at the end of a symmetrizing path needs squares to be swallowed by the wild point. Writing it down carefully, the parity version turns out not to work. The version that works uses the square's **label**, the passes its corners sit on, instead of the parity.

**The idea.**
- A square moving continuously along a family of curves keeps its label.
- For a region of square configurations defined by a label, the degree of the square map cannot change unless a square crosses the region's boundary.

The label has to be defined without a polar-angle chart. Absolute pass indices change by 2 whenever an angle crosses its branch cut, although nothing geometric happens. So corner $A$ is placed on its arm by its own lift, and the other corners are given by their pass offsets from $A$ (Lemma 390.2).

So a square of the symmetric curve can disappear only by leaving its label region: shrinking into the wild point, or drifting out of its band of radii. Nothing else can remove it. Not the twist's events, not folds, not partners.

**Results.**
1. **Parity bookkeeping alone is empty (Remark 390.1).** For any path from $\gamma_0$ to $\gamma_1$, the number of graceful events at cutoff $s$ is $L(\gamma_0,s)+L(\gamma_1,s)$ mod 2. It says nothing new.
   - It also cannot rule out a cascade: a chain of branches, each joining one square to a smaller square's half-turn partner, could carry an odd count down through every scale with no single branch collapsing.
   - The earlier sketch's "every arc joins two symmetric squares, so the count is even" overlooked that the global squares of $\gamma_0$, unlike the small ones near $p$, are not paired.
2. **Labels are conserved (Lemma 390.2).** The label is $(\text{arm of }A,\ o_B,o_C,o_D)$: $A$'s arm and the other corners' pass offsets from $A$. Along a continuous family of curves near a two-arm point, it stays fixed as long as no corner passes through the wild point and the side stays below $|A|/\sqrt2$.
3. **Label degree is invariant (Theorem 390.3).** Fix a label $\ell$ and the region $\Omega_\ell=\{\Theta_{\rm lo}<\Theta_A<\Theta_{\rm hi},\ q_{\rm in}<|w|<q_{\rm max}\}\subset\mathbb R^3$, where $w$ is the side vector in $A$'s frame divided by $|A|$, and $q_{\rm max}<1/\sqrt2$. Suppose no curve $\gamma_\lambda$ of the path has a square with label $\ell$ on $\partial\Omega_\ell$. Then the Brouwer degree of the three-equation label map on $\Omega_\ell$ is the same for every $\lambda$.
4. **No cascades (Corollary 390.4).** Half-turn partners carry different labels: the half-turn keeps $\Theta_A$, $w$ and the offsets, and swaps $A$'s arm. So no branch can join partners, and the chain of Remark 390.1 is impossible. Suppose the label map of the symmetric curve has nonzero total degree on $\Omega_\ell$, and the label's squares never reach $\partial\Omega_\ell$. Then **every curve on the path has a square with label $\ell$**. The only ways out are collapse and drift.
   - Nonzero total degree follows if the symmetric curve has exactly one zero of the label in $\Omega_\ell$ and it is nondegenerate. This is the uniqueness hypothesis U0.
   - U0 is *assumed, not proved*. The searches below find no second root, but Newton runs from random starts cannot exclude one, and an unseen root could carry the opposite sign.
5. **Checks on the plain spiral:**
   - every one of the 81 squares of a v380 search, in all four choices of corner $A$ (324 configurations), is a nondegenerate zero of its label map;
   - every half-turn partner is a zero of its own, different label;
   - an independent search with another seed finds 73 squares whose 292 anchored labels are all different, so no label is carried by two squares;
   - moving a square across the negative real axis changes the anchored residuals by $2\cdot10^{-10}$, and the old absolute ones by exactly $2\pi$.

## Remark 390.1: why parity bookkeeping cannot do it

Run the twist path $\gamma_\lambda$ and let $E(s)$ count graceful squares whose side passes $s$. The path principle (v373) gives $L(\gamma_1,s)\equiv L(\gamma_0,s)+E(s)$, and that is all the parity says.

The arc picture refines this. Each square traces an arc in (twist, configuration) space, and an arc can end in one of three places: at $\lambda=0$, at $\lambda=1$, or at the wild point. Counting arcs whose ends straddle $s$ gives $E(s)\equiv\#\{\text{ends above }s\}$. But the $\lambda=0$ ends include the large, global squares of $\gamma_0$, which the half-turn does not pair. Their number above $s$ is odd, which is exactly why $L(\gamma_0,s)=1$. So arcs joining $\lambda=0$ squares to each other can straddle every cutoff an odd number of times: a chain from a global square, through a small square $t_1$, its partner $t_1'$ (same size), a smaller square $t_2$, its partner, and so on. No single arc needs to collapse. Parity cannot exclude this. Labels can.

## Lemma 390.2 and Theorem 390.3

**Coordinates.** Corner $A$ sits on arm $\mathrm{arm}_A\in\{0,1\}$ at its own lift $\Theta_A$:
$$A=\Theta_A^{-a}\exp\!\big(i(\Theta_A-\pi\,\mathrm{arm}_A+\mathrm{arm}_A\,\lambda\varphi(\Theta_A^{-a}))\big).$$
The other corners are $X=A(1+w\,m_X)$, with $m_B=1$, $m_C=1+i$, $m_D=i$. Here $w=(t/|A|)e^{i\alpha}$ is the side vector in $A$'s frame. Corner $X$'s residual is
$$R_X=|X|^{-1/a}-\Theta_A-\arg(X/A)-\pi o_X+\lambda\big(\mathrm{arm}_X\varphi(|X|)-\mathrm{arm}_A\varphi(|A|)\big),\qquad \mathrm{arm}_X\equiv\mathrm{arm}_A+o_X\pmod 2.$$
- This is the difference of $X$'s and $A$'s comb-phase equations. The polar angles enter only through $\arg(X/A)$.
- That argument is the principal one. It is continuous because $|w\,m_X|\le\sqrt2|w|<1$ keeps $X/A$ off the negative real axis.

So the label map $F_\lambda:(\Theta_A,\operatorname{Re}w,\operatorname{Im}w)\mapsto(R_B,R_C,R_D)$ is continuous on the closure of the bounded open set $\Omega_\ell\subset\mathbb R^3$, jointly in $\lambda$. There is no chart and no wrapping.

**Lemma 390.2.** Along a continuous family of curves, a square's $\Theta_A(\lambda)$ and $w(\lambda)$ move continuously. $A$ cannot change arm without passing through the wild point, because away from $p$ the arms are disjoint. Each offset $o_X$ is an integer. It is determined by $|X|^{-1/a}-\Theta_A-\arg(X/A)$ up to the twist terms, which stay well below $\pi$, and that quantity moves continuously while $|w|<1/\sqrt2$. So the label is constant. ∎

**Theorem 390.3.** If $F_\lambda\ne0$ on $\partial\Omega_\ell$ for every $\lambda$, homotopy invariance of the Brouwer degree gives $\deg(F_\lambda,\Omega_\ell)=\deg(F_0,\Omega_\ell)$. ∎

The degree uses no transversality and no genericity. Folds, births and deaths inside $\Omega_\ell$ are invisible to it. Only the boundary matters, and it has three parts:
- the lift faces $\Theta_A=\Theta_{\rm lo}$ (drifting outward) and $\Theta_A=\Theta_{\rm hi}$ (moving toward the point);
- the outer side face $|w|=q_{\rm max}$;
- the inner side face $|w|=q_{\rm in}$.

The orientation is part of $w$, so it adds no boundary.

**Relabeling.** Each geometric square gives four configurations, one for each choice of corner $A$. They lie in four different label regions (the offsets shift cyclically), and each is counted once in its own region.

## Corollary 390.4

Suppose $\deg(F_0,\Omega_\ell)\ne0$ and no square of the label reaches $\partial\Omega_\ell$ for any twist on the path. Then the degree stays nonzero, and every $\gamma_\lambda$, the twisted curve included, has a square with label $\ell$.

The simplest route to $\deg(F_0,\Omega_\ell)=\pm1$ is U0: one zero of the label in the region, and a nondegenerate one. On the plain spiral, v382's branches give one square per resonance label, and the checks below find no label carried by two squares. But that is evidence, not proof, and a direct proof of nonzero total degree would make U0 unnecessary.

What remains is a confinement statement: for which labels does the label's square stay inside a fixed compact region, whatever the twist? That is v391.

## Checks

On the plain spiral, from v380's random search (6 000 seeds), each square taken with all four choices of corner $A$:

| check | result |
|---|---|
| configurations that are zeros of their label map (max residual) | 324 / 324 ($1.0\cdot10^{-11}$) |
| nondegenerate, local degree $+1$ / $-1$ | 324 / 324: 162 / 162 |
| half-turn partner is a zero of its own label, and that label differs | 324 / 324 |
| configurations inside the cut-free domain $|w|<1/\sqrt2$ | 322 / 324 |
| independent search: squares, anchored labels, labels carried by two squares, on their v382 branch | 73, 292, 0, 73 |
| residual jump when a square crosses the negative real axis: anchored / absolute | $1.9\cdot10^{-10}$ / $2\pi$ |

The two configurations outside the domain are large squares, with side up to $0.83|A|$. The local degree's sign depends on the orientation conventions and is not meaningful by itself. What matters is that it is never 0.

## What this changes

- **The obstacle is now a single, local question:** can a label's square be pushed out of a compact region of radii? That means into the wild point, or away from it.
- **Parity and partners drop out.** The degree argument needs neither the sector rule nor the half-turn pairing, and it applies to crossed squares as much as graceful ones.
- **It replaces the path condition.** Turning off the twist may create squares of side $s$ anywhere (v389). That no longer matters, as long as some label's square cannot leave its region.

## Review corrections

An external review of this PR found that the first version's labels were absolute pass indices read off a polar-angle chart. Those are not globally conserved: crossing the chart's branch cut changes every index by 2. So the four-equation map was not shown to be continuous on an annular region. This version replaces them with the anchored coordinates above, which need no chart.

It also states U0 as an assumption. The checks now cover all four anchorings, partner zeros, and an explicit chart-crossing test.

[Code](../code/label_degree_v390.py) · [Tests](../tests/test_label_degree_v390.py) · [Data](../results/v390/label_check.json)
