# v382. Tilted squares exist at every resonance

[v379](v379-squares-on-a-smooth-spiral.md) showed that on a smooth spiral, small squares are tilted and obey a selection rule. With the corners' turn offsets labelled $(n_1,n_3)$, the alternating phase sum $j$ satisfies $j\approx(2+\frac1a)\pi K(r)\,n_1n_3$, so squares can only sit at special "resonant" radii. Whether a square actually exists at each resonance was left open. This note proves that it does, with the same one-unknown idea as [v381](v381-squares-between-two-curves.md).

Put one corner $A$ on an arm. Its two neighbours can then be placed exactly on the turns $n_1$ and $n_3$ away, which leaves the opposite corner free. That corner lands on the spiral exactly when one real number $h(A)$ is an integer, and that integer is $j$. As $A$ moves inward, $h$ shrinks steadily to $0$. So $h$ passes through every integer between its starting value and $0$, and each crossing is an inscribed square.

**Results.**
1. **Reduction (Theorem 382.1).** Along any continuous branch of three-corner near-squares (corner $A$ on an arm, $B$ and $D$ on the turns $n_1,n_3$ away), inscribed squares are exactly the points where $h=(\Phi(C)-\Phi(A))/\pi-n_1-n_3$ is a nonzero integer $j$. That $j$ is v379's alternating sum.
2. **Existence (Theorem 382.2).** If $h$ runs from below $1$ to above $J$ along a branch, there are squares with $j=1,\dots,J$. To leading order $h=(2+\frac1a)\pi K(r)n_1n_3\to0$ as $A$ moves in. So every resonance $(n_1,n_3,j)$ whose radius lies inside the arc carries a square.
3. **Exact checks on the half-turn spiral ($a=1/2$):**
   - all 252 branches with $|n_1|,|n_3|\le8$: $h$ is monotone on every one, and every integer in its range is a square, 807 labelled squares in all, with arc residuals at most $4\times10^{-12}$;
   - these are 448 distinct squares, since a square gets a different label from each corner;
   - all 78 of v379's randomly found squares in this range lie on these branches.
4. **The radius rule is first order only (correction to v379).** The predicted radius depends on which corner is labelled $A$. Labels with $n_1n_3<0$ predict to a median 1.4% for $t/r<0.1$; labels with $n_1n_3>0$ are off by about 10–30%. An earlier claim of "median 0.5%" (from 7 squares) and one mislabelled example in v379 are corrected.

## The reduction

On the half-turn spiral, the arc is $\{\Phi\in\pi\mathbb Z\}$ with $\Phi(z)=|z|^{-1/a}-\arg z$, the argument lifted from $A$. Put $A$ on arm one at lift $\Theta_A$, and set $B=A+te^{i\theta}$, $D=A+ite^{i\theta}$. The **three-corner conditions**
$$\Phi(B)-\Phi(A)=\pi n_1,\qquad\Phi(D)-\Phi(A)=\pi n_3$$
are two equations in $(\theta,t)$. Near the straight-comb solution, with radial components $-\pi K r n_1$ and $-\pi K r n_3$ (the phase grows inward), the Jacobian is a near-similarity. This is the same structure as Lemma 378.2. So the implicit function theorem gives a smooth branch $(\theta(\Theta_A),t(\Theta_A))$ while $t/r$ and $K$ are small.

**Theorem 382.1.** On such a branch, with $C=B+D-A$,
$$h(\Theta_A)=\frac{\Phi(C)-\Phi(A)}{\pi}-n_1-n_3$$
is continuous, and $ABCD$ is an inscribed square iff $h\in\mathbb Z\setminus\{0\}$. Then $j=h$ is the alternating sum $k_A-k_B+k_C-k_D$ of v379.

*Proof.* $A,B,D$ lie on the arc by construction, so $C$ lies on it iff $\Phi(C)\in\pi\mathbb Z$, that is, iff $h\in\mathbb Z$. Since $k_B=k_A+n_1$ and $k_D=k_A+n_3$, we get $h=k_C-k_A-n_1-n_3=k_A-k_B+k_C-k_D$. The value $0$ is excluded because v379 rules out $j=0$ for tilted squares; on these branches $h$ never reaches it. ∎

**Theorem 382.2.** If $h$ takes values below $1$ and above $J\ge1$ on a branch (or the mirror statement for $n_1n_3<0$), there are inscribed squares on it with $j=1,\dots,J$. By Theorem 379.2,
$$h=(2+\tfrac1a)\pi K(r)\,n_1n_3\,\big(1+O(t/r+K)\big),$$
which tends to $0$ as $A\to p$ and grows outward. So every $(n_1,n_3,j)$ with $j\,n_1n_3>0$ and $|j|<(2+\frac1a)\pi K(r_{\max})|n_1n_3|(1-o(1))$ is realized. In particular, with $j=1$ and $|n_1|,|n_3|\to\infty$, resonant squares exist at every scale down to the wild point, with $t/r\to0$.

*Proof.* Intermediate value theorem on the continuous $h$. ∎

This is the existence half of the smooth comb law. Squares sit only at resonances (v379), and every resonance has one.

## Checks

**Branches.** For each $(n_1,n_3)$ with $1\le|n_1|,|n_3|\le8$ and $|n_1|+|n_3|\ge3$ (252 pairs), the branch is followed from deep inside ($\Theta_A=900$) out to $\Theta_A=6$. It is restricted to the regime of the theorem: side at most $0.3|A|$, every corner within a quarter turn of $A$. Beyond that, the phase lift around $A$, and so the turn count, stops being meaningful. Each crossing of an integer is refined by bisection and checked on the arc.

- $h$ is monotone on all 252 branches. So each label gives exactly one square per branch.
- 807 labelled squares, all with arc residual at most $4\times10^{-12}$, all with $j\,n_1n_3>0$. No integer in any branch's range is missed.
- Collapsing labels gives 448 distinct squares: 265 seen under one label in this range, 90 under two, 10 under three, 83 under all four.
- Of v379's 200 randomly found squares, 78 have $|n_1|,|n_3|\le8$ and $t/r\le0.3$. All 78 lie on these branches; one needs the size cutoff relaxed from $0.30$ to $0.32$, since v379 measures $t/r$ from the centre and the scan from $A$. So the random search found nothing outside the branches.

**Radius rule.** The leading-order rule evaluates $K$ at one radius, but a square spans a range of radii. With the same signs ($n_1n_3>0$) the corners spread radially by up to $t\sqrt2$; with opposite signs they hardly spread at all. That shows in the errors:

| labels | $t/r<0.1$ | $0.1\le t/r<0.2$ | $0.2\le t/r<0.3$ |
|---|---|---|---|
| $n_1n_3<0$: median (max) error of predicted radius | 1.4% (4.4%) | 5.1% (16%) | 6.8% (21%) |
| $n_1n_3>0$: median (max) | 9.6% (13.5%) | 20% (31%) | 32% (40%) |

Every square has labels of both kinds, so the opposite-sign label is the one to use. Evaluating the factors at the side midpoints and the centre roughly halves the error but leaves it first order. A sharp radius formula would need the next term.

**Corrections to v379.**
- Its summary claimed the predicted radius was accurate to a median 0.5% for $t/r<0.15$. That came from 7 squares, two of them off by 18%.
- Its worked example labelled a square $(2,3,1)$. That square's label is $(-2,3,-1)$, with the same predicted radius; the $(2,3,1)$ square is at $0.136$.

Both are fixed in the v379 note. None of v379's theorems or other numbers change.

## Where the comb picture stands

On a two-arm spiral point the location of small squares is now fully described, for both kinds:
- **Aligned squares** exist exactly at stationary gaps between turns (v381).
- **Tilted squares** exist exactly at the resonances of the selection rule (v379 and v382).
- **Outer bound:** both obey the general comb law (v380).

What remains open is unchanged: v379's exclusion of aligned squares on the plain spiral is numerical, and the path condition of v374 is still open.

[Code](../code/tilted_squares_v382.py) · [Tests](../tests/test_tilted_squares_v382.py) · [Data](../results/v382/tilted_check.json)
