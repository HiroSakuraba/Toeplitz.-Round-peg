# v383. Small squares on the plain spiral touch at least three turns

[v379](v379-squares-on-a-smooth-spiral.md) found, only numerically, that the plain spiral has no aligned squares (Proposition 379.3). These are squares with two corners on one turn and two on another. This note proves it, and more. On the plain two-arm spiral, no small square has all four corners on one or two passes of the arms. So every small square is of the tilted kind that [v382](v382-tilted-squares-exist.md) builds at the resonances.

The tool is an elementary fact, the **British flag theorem**. For any rectangle $P_1P_2P_3P_4$ and any point $O$,
$$|OP_1|^2+|OP_3|^2=|OP_2|^2+|OP_4|^2 .$$
Take $O$ at the wild point. Each corner's squared distance from it is fixed by which pass the corner sits on. So the theorem becomes an exact bookkeeping rule on the passes.

Here is the whole idea in one line. Squared radius jumps by a definite amount from one pass to the next. The flag sums can only balance if the corners are spread over the passes in a balanced way. When two adjacent corners sit on one pass and the other two on another (the aligned case), the sums do balance to first order. One order further, though, they balance only if the two passes have the same radial slope there, which means the gap between them is stationary. On the plain spiral the gap shrinks steadily, so that never happens.

**Results.**
1. **Slope identity (Theorem 383.1).** Take two passes $\rho=R_{\rm in}(\theta)$ and $\rho=R_{\rm out}(\theta)$ that are polar graphs over the square's angular window. Suppose $P_1,P_2$ lie on the inner pass and $P_4,P_3$ on the outer one, with $P_1P_2$ and $P_4P_3$ opposite sides. Then
$$R_{\rm in}'(\eta_1)\,(\theta_2-\theta_1)(\rho_1+\rho_2)=R_{\rm out}'(\eta_2)\,(\theta_3-\theta_4)(\rho_3+\rho_4),$$
   for some $\eta_1$ between $\theta_1,\theta_2$ and some $\eta_2$ between $\theta_4,\theta_3$.
2. **Aligned squares need matching slopes (Corollary 383.2).** The equal chords $P_1P_2=P_4P_3$ pin the slope ratio $Q=R_{\rm in}'(\eta_1)/R_{\rm out}'(\eta_2)$ to an explicit band around 1. The band has second-order width, in the window and the spiral's pitch. This holds for any two polar graphs. It is the radial form of v381's localization: the gap $R_{\rm out}-R_{\rm in}$ must be nearly stationary under the square.
3. **No square on one or two passes (Theorem 383.3).** On the plain half-turn spiral $\rho=(\theta+\pi k)^{-a}$, consider squares of side $t\le0.678\,|Z|$, where $Z$ is the centre, so the corners span at most 1 radian of angle. Every such square whose outermost pass is beyond $\Theta=L_0(a)$ has its corners on at least three different passes. For $a=1/2$, $L_0=6.78$ (radius 0.384). The cases are:
   - all four corners on one pass: impossible on a convex arc that turns by less than $\pi$;
   - a 3–1 split, or a 2–2 split with each pair diagonal: the flag sums cannot balance;
   - a 2–2 split with each pair adjacent (aligned): $Q\le\kappa^{a+1}<1$ on the spiral, while the chords force $Q$ above that.
4. **One three-pass pattern also excluded (Proposition 383.4).** Take two adjacent corners on one pass and the other two on two different passes. This is impossible once the corners' radii differ by less than a factor $(\pi-1)^{a/(2a+1)}$. For $a=1/2$ that factor is 1.21, which holds when $t\le0.134\,|Z|$. What remains are the patterns 1-1-1-1 and 2-1-1 with the doubled pass on a diagonal. These are exactly the labels of v382's tilted squares.
5. **Checks:**
   - Plain spiral, random search: 77 squares, all of pattern 1-1-1-1 or 2-1-1 diagonal. All have window at most 0.71 and outermost pass beyond 9.2, so they are inside the theorem's range.
   - Newton's method with the passes prescribed, 600 seeds for each excluded pattern: no solution at all on the plain spiral. On a wobbling spiral the same seeds find 143 aligned squares.
   - Wobbling spirals: all 514 aligned squares lie inside the slope band of Corollary 383.2, with $|1-Q|\le0.73\%$. Wobbles also allow the 2-1-1 adjacent pattern (14 squares), which the plain spiral forbids.
   - All 807 labelled squares of v382 have pattern 1-1-1-1 or 2-1-1 diagonal.

## The flag identity on passes

Near a two-arm spiral point, the arms $z=\pm\Theta^{-a}e^{i\Theta}$ together form a family of polar graphs, one per pass:
$$\rho=R_k(\theta)=(\theta+\pi k)^{-a}.$$
Arm one gives even $k$ and arm two odd $k$. A point at polar angle $\theta$ on pass $k$ has lift $\Theta=\theta+\pi k$ and squared radius $F(\Theta)=\Theta^{-2a}$. Larger $k$ is further in. $F$ is decreasing and convex.

A square of side $t$ with centre $Z$ fits in the disk of radius $t/\sqrt2$ about $Z$. So its four polar angles lie in one arc of width
$$W\le2\arcsin\!\big(t/(\sqrt2|Z|)\big),$$
which is at most 1 once $t\le\sqrt2\sin\tfrac12|Z|=0.678|Z|$. Measure each corner's angle $\theta_X$ in that arc. Its pass $k_X$ is then well defined, and $\Theta_X=\theta_X+\pi k_X$. The flag theorem reads
$$F(\Theta_A)+F(\Theta_C)=F(\Theta_B)+F(\Theta_D)\qquad(\ast)$$
for corners $A,B,C,D$ in cyclic order.

## Unbalanced splits

**Diagonal 2–2.** Say $A,C$ are on pass $k+n$ and $B,D$ on pass $k$, with $n\ge1$. Then every inner lift is at least $\theta_{\min}+\pi(k+n)$, and every outer lift is at most $\theta_{\max}+\pi k$. Since $W<\pi$, every outer lift is smaller than every inner lift. $F$ is decreasing, so the right side of $(\ast)$ is strictly larger. ∎

**3–1, lone corner outside.** Say $A,B,C$ are on pass $k+n$ and $D$ on pass $k$. Let $\Theta_1=\theta_{\min}+\pi(k+n)$. Then $(\ast)$ needs $F(\Theta_A)+F(\Theta_C)-F(\Theta_B)=F(\Theta_D)$. The left side is at most $2F(\Theta_1)-F(\Theta_1+W)$, and $\Theta_D\le\Theta_1+W-\pi$. Convexity gives
$$F(\Theta_1+W-\pi)-F(\Theta_1)\ge(\pi-W)|F'(\Theta_1)|>W|F'(\Theta_1)|\ge F(\Theta_1)-F(\Theta_1+W)$$
whenever $W<\pi/2$. So $F(\Theta_D)$ is too large. ∎

**3–1, lone corner inside.** Now $D$ is on pass $k+n$ and the triple on pass $k$, with floor $\Theta_0=\theta_{\min}+\pi k$. Then $(\ast)$ needs $F(\Theta_D)+F(\Theta_B)=F(\Theta_A)+F(\Theta_C)\ge2F(\Theta_0+W)$, while the left side is at most $F(\Theta_0+\pi)+F(\Theta_0)$. Convexity bounds $F(\Theta_0)-F(\Theta_0+W)\le W|F'(\Theta_0)|$ and $F(\Theta_0+W)-F(\Theta_0+\pi)\ge(\pi-W)|F'(\Theta_0+\pi)|$. So the case is impossible once
$$(\pi-W)\Big(\frac{\Theta_0}{\Theta_0+\pi}\Big)^{2a+1}>W .$$
The left side increases with $\Theta_0$. For $W=1$ this holds beyond $\Theta_0=6.78$ when $a=1/2$, and beyond $10.87$ when $a=1$. This case sets $L_0$. ∎

**One pass.** On a pass, the curvature has the sign of $\Theta^{-2a}\big(1+(a^2-a)/\Theta^2\big)>0$ for $\Theta>1/2$. Over the window the tangent turns by $W$ plus the change in the angle between radius and tangent, which is less than $\pi/2$. So the total turning is less than $1+\pi/2<\pi$. Points on a strictly convex arc that turns by less than $\pi$ are in convex position in arc order. So that order is the square's cyclic order, and three consecutive sides of the square are chords of the arc. Each chord points along a tangent of the arc, so their directions lie in an interval shorter than $\pi$. But the first and third sides are antiparallel. ∎

**2-1-1 adjacent (Proposition 383.4).** Say $A,B$ are on one pass, and $C,D$ on two different passes. By $(\ast)$, $F(\Theta_C)-F(\Theta_D)=F(\Theta_B)-F(\Theta_A)$. The left side jumps across at least $\pi-W$ in lift. The right side moves within one window, at most $W$ in lift. With $|F'(\Theta)|=2a\rho^{(2a+1)/a}$, they cannot be equal if $(\pi-W)\rho_{\min}^{(2a+1)/a}>W\rho_{\max}^{(2a+1)/a}$, where $\rho_{\min},\rho_{\max}$ are the smallest and largest corner radii. Since $\rho_{\max}/\rho_{\min}\le(|Z|+t/\sqrt2)/(|Z|-t/\sqrt2)$, this holds for $t\le0.134|Z|$ when $a=1/2$ and $W\le1$. ∎

## The aligned case

Label the square so that $P_1,P_2$ are on the inner pass $R_{\rm in}$ and $P_4,P_3$ on the outer pass $R_{\rm out}$, with $P_1P_2$ and $P_4P_3$ opposite sides. Write $\rho_X=|P_X|$, $\varphi_{\rm in}=\theta_2-\theta_1$ and $\varphi_{\rm out}=\theta_3-\theta_4$.

**Theorem 383.1 (slope identity).** The flag theorem gives $\rho_2^2-\rho_1^2=\rho_3^2-\rho_4^2$. By the mean value theorem, $\rho_2-\rho_1=R_{\rm in}'(\eta_1)\varphi_{\rm in}$ and $\rho_3-\rho_4=R_{\rm out}'(\eta_2)\varphi_{\rm out}$. Substituting gives the identity in the summary. ∎

On the spiral both slopes are negative, so $\varphi_{\rm in}$ and $\varphi_{\rm out}$ have the same sign. The identity says
$$Q:=\frac{R_{\rm in}'(\eta_1)}{R_{\rm out}'(\eta_2)}=x\cdot\frac{\rho_3+\rho_4}{\rho_1+\rho_2},\qquad x=\frac{|\varphi_{\rm out}|}{|\varphi_{\rm in}|}.$$

**The chords pin $Q$ near 1.** The two equal sides give
$$4\rho_1\rho_2\sin^2\tfrac{\varphi_{\rm in}}2+(\rho_2-\rho_1)^2=t^2=4\rho_3\rho_4\sin^2\tfrac{\varphi_{\rm out}}2+(\rho_3-\rho_4)^2 .$$
Put $p=|\rho_3-\rho_4|/t$. This is the outer pass's pitch, a first-order small quantity. Then
$$\frac{\sin(|\varphi_{\rm out}|/2)}{\sin(|\varphi_{\rm in}|/2)}\ge s:=\sqrt{1-p^2}\sqrt{\frac{\rho_1\rho_2}{\rho_3\rho_4}} .$$

Two facts about angles up to $W$ turn this into a bound on $x$, with $g=\frac{W/2}{\sin(W/2)}$ and $T=\frac{\tan(W/2)}{W/2}$:
- Always $x\ge s/g$, because $xu\ge\sin xu\ge s\sin u\ge su/g$.
- If $s<1$, either $x\ge1$, or the mean value theorem on $\sin$ gives $1-x\le(1-s)\tan u/u\le T(1-s)$. Here $u=|\varphi_{\rm in}|/2$.

Hence
$$Q\ \ge\ \max\big(1-T(1-s),\ s/g\big)\cdot\frac{\rho_3+\rho_4}{\rho_1+\rho_2}.$$
The upper bound comes from the same argument with the roles swapped. This is Corollary 383.2. For concentric circles, $p=0$, $s=\rho_{\rm in}/\rho_{\rm out}$ and $Q$ is exactly 1 to first order. The correction is $T-1\approx W^2/12$ times a first-order quantity.

**Theorem 383.3, aligned part.** On the plain spiral $R_k'(\theta)=-a\Theta^{-a-1}$. Let $L$ be the window floor on the outer pass, $\theta_{\min}+\pi k$, and suppose the passes are $n\ge1$ apart. The inner lifts are then at least $L+\pi n$, and the outer lifts at most $L+W$. So
$$Q=\Big(\frac{\Theta_{\rm out}(\eta_2)}{\Theta_{\rm in}(\eta_1)}\Big)^{a+1}\le\kappa^{a+1},\qquad\kappa=\frac{L+W}{L+\pi n}<1 .$$

On the other side, each factor above can be bounded in terms of $\kappa$ and $L$:
- $(\rho_3+\rho_4)/(\rho_1+\rho_2)\ge\kappa^{-a}$;
- $s\ge\sigma\kappa^a$, with $\sigma=\sqrt{1-p^2}\,(1+W/L)^{-2a}$;
- $p\le a\,g\,(1+W/L)^a/L$.

The condition "chord bound beats spiral bound" becomes
$$\phi(\kappa):=T-(T-1)\kappa^{-a}-\kappa^{a+1}\ >\ T(1-\sigma)\qquad\text{or}\qquad\sigma/g>\kappa^{a+1}.$$
Here $\phi$ is concave with $\phi(1)=0$, and $1-\sigma\le p^2+2aW/L$. Fix a cut $\kappa_1$:
- For $\kappa\le\kappa_1$, the second inequality holds once $\sigma/g>\kappa_1^{a+1}$.
- For $\kappa_1\le\kappa<1$, concavity gives $\phi(\kappa)\ge\phi(\kappa_1)(1-\kappa)/(1-\kappa_1)$, and $1-\kappa\ge(\pi-W)/(L+\pi)$ for every $n\ge1$. So the first inequality holds once
$$\frac{\pi-W}{L+\pi}\cdot\frac{\phi(\kappa_1)}{1-\kappa_1}>T\Big(p^2+\frac{2aW}{L}\Big).$$

Both conditions are monotone in $L$ (multiply the second by $L$), so each holds beyond a single computed threshold, uniformly in $n$. With $W=1$ and $\kappa_1$ optimized:

| $a$ | aligned case from $L=$ | 3–1 case from $L=$ | $L_0$ | radius $L_0^{-a}$ |
|---|---|---|---|---|
| 1/4 | 1.08 | 4.75 | 4.75 | 0.68 |
| 1/2 | 2.45 | 6.78 | 6.78 | 0.38 |
| 1 | 6.17 | 10.87 | 10.87 | 0.092 |
| 2 | 18.8 | 19.1 | 19.1 | 0.0027 |

The margin is far from tight. At $a=1/2$, $L=10$, $n=1$, the chords force $Q\ge0.891$ while the spiral allows at most $0.766$. ∎

**Remark (where v379's rule comes from).** Run $(\ast)$ on a tilted square with passes $0,n_1,n_1+n_3+j,n_3$ and expand $F$ to second order in the pass offsets. This gives $\pi jF'+\tfrac12\pi^2F''\cdot2n_1n_3\approx0$, that is, $j\approx(2a+1)\pi n_1n_3/\Theta$. With $K=a r^{1/a}=a/\Theta$, this is exactly v379's leading-order selection rule $j=(2+1/a)\pi K n_1n_3$. The window offsets enter only at the next order. So the flag identity is the selection rule's first-order skeleton. The aligned case is its degenerate corner, $n_1=j=0$, where one has to go one order further, as above.

## Checks

The test families are half-turn symmetric spirals $r=\Theta^{-1/2}(1+\frac\kappa\Theta\sin\omega\Theta)$ from v380. The plain spiral is $\kappa=0$. Squares come from v380's random Newton search, 6 000 seeds per family. Each square's corners are assigned window angles and passes as above.

| $\kappa,\omega$ | squares | 1-1-1-1 | 2-1-1 diag | 2-2 adjacent | 2-1-1 adjacent | aligned in band | max $\lvert1-Q\rvert$ |
|---|---|---|---|---|---|---|---|
| 0, – (plain) | 77 | 66 | 11 | 0 | 0 | – | – |
| 0.2, 2 | 79 | 59 | 20 | 0 | 0 | – | – |
| 0.3, 1 | 404 | 267 | 71 | 65 | 1 | 65/65 | 0.29% |
| 0.5, 3 | 1013 | 569 | 194 | 239 | 11 | 239/239 | 0.73% |
| 0.7, 2 | 344 | 114 | 18 | 210 | 2 | 210/210 | 0.66% |

- The flag identity holds on every square to $10^{-15}$ relative.
- No square in any family has pattern 4, 3–1 or 2–2 diagonal. For the wobbles this is not covered by the theorem, since their $F$ is not monotone along a pass, but it is what the flag rule predicts.
- On the plain spiral, every square has window $W\le0.71$ and outer floor $\ge9.2$, inside the theorem's range.
- **Prescribed passes.** Newton's method was run with each corner held on a prescribed pass, 600 seeds per excluded pattern, at $\Theta\in[8,60]$ with pass gaps 1–3 in both directions. On the plain spiral it finds nothing for any pattern. On the wobble $\kappa,\omega=0.5,3$ it finds 143 aligned squares and nothing in the other four patterns. So the search can find these squares when they exist.
- **v382.** Every one of v382's 807 labelled squares has pattern 1-1-1-1 (692) or 2-1-1 diagonal (115).

## What this changes

- **Proposition 379.3 is now a theorem.** On the plain spiral, for $L\ge L_0(a)$ and squares spanning at most a radian, there are no aligned squares, and none of the other one- or two-pass patterns either. Corollary 379.4 (the smooth comb law, $|j|\ge1$) no longer needs the numerical assumption in that range.
- **The classification near a plain spiral point is complete in kind.** A small square there must have pattern 1-1-1-1 or 2-1-1 diagonal, which is v382's labelling with $n_1n_3\ne0$. v382 shows these exist at every resonance. Still open is whether every such square lies on v382's branches, which would need the three-corner solution to be unique.
- **Aligned squares on any spiral sit at matched slopes.** Corollary 383.2 holds for any pair of polar graphs. It is the polar companion of v381's frame-based localization, and it is what made v380's Principle 380.2 work with the radial gap.
- **Still open:** the v374 path condition, wild points that accumulate, and a sharp radius formula beyond v379's first order.

[Code](../code/two_turn_squares_v383.py) · [Tests](../tests/test_two_turn_squares_v383.py) · [Data](../results/v383/two_turn_check.json)
