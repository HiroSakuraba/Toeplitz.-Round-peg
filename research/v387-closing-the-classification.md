# v387. Closing the gaps in the square classification

Two gaps were left in the map of small squares near a smooth two-arm point, and both v386 and the path-condition work lean on them:
- v383's exclusion of one- and two-pass squares was proved for the plain spiral only;
- nobody had checked that every allowed square lies on one of v382's families.

This note closes the first for twisted arms, with explicit hypotheses, and gives strong evidence for the second.

**Results.**
1. **Twisted arms (Theorem 387.1).** Let arm two be arm one turned by $\pi+\varphi(r)$, with $|\varphi|\le\Phi$ and $|r\varphi'|\le L$. Consider small squares whose corners span at most $W$ in angle, so $t\le\sqrt2\sin(W/2)\,|Z|$. On twisted arms, v383's flag arguments go through with two substitutions:
   - the pass spacing $\pi$ becomes $\pi-\Phi$;
   - the window $W$ becomes $W_1=W+2L\,\operatorname{artanh}(\sin\frac W2)$, which adds the twist's oscillation across the square.

   The aligned case also picks up the twist's slope: on arm two, $dR/d\theta=-a\Theta^{-a-1}/(1+\mu)$ with $|\mu|\le aL/\Theta$. So it needs the margin
$$(2a+1)(\pi-\Phi-W)\ >\ T(W)\,a\,(\pi+\Phi+W)+aL .$$
   For v374's twist $\varphi=0.4\sin(5\log(1/r))$ ($\Phi=0.4$, $L=2$), with $a=1/2$ and $W=0.3$ (so $t\le0.21|Z|$), every condition holds once the lift is at least $6.41$, that is radius $\le0.39$. So:
   - no small square there has its corners on one or two passes;
   - the pattern 2-1-1 adjacent is excluded for $t\le0.126|Z|$;
   - what is left is 1-1-1-1 and 2-1-1 diagonal, as on the plain spiral.
2. **Steep twists break it.** The aligned margin is $2.95-0.5L$ at $W=0.3$, so the argument stops at $L\approx5.9$. This is the same gradient effect v386 found for graceful squares.
3. **v382's families are complete, numerically (Proposition 387.2).** Every square found by an independent search on the plain spiral has a corner whose two neighbours are off its pass. At that corner's lift, v382's three-corner solve, started both from the square itself and from the generic guess, returns exactly that square. So each found square is the branch square of its label. This held for 81 of 81 squares.
4. **Checks on twisted arms:**
   - Newton's method with prescribed one- or two-pass labels, 300 seeds for each of five patterns: no solution;
   - a random search finds 42 squares, all of pattern 1-1-1-1 or 2-1-1 diagonal;
   - arm two is convex at every lift checked from 3 to $10^8$;
   - at three times the twist ($\Phi=1.2$, $L=6$), where the criterion no longer applies, the prescribed-label search still finds no aligned squares. So the criterion looks conservative.

## Lifts on twisted arms

A corner on pass $k$ has lift $\Theta=\theta+\pi k-\varphi^*(\rho)$, with $\varphi^*=\varphi$ on odd passes (arm two) and $0$ on even ones, and squared radius $F(\Theta)=\Theta^{-2a}$ on both arms. For a square whose angles lie in an arc of width $W$:
- **Within a pass,** corner lifts lie in a window of width at most $W_1$. The angles add $W$. The twist adds at most $L\log(R/r)\le2L\operatorname{artanh}(u)$, where $u=t/(\sqrt2|Z|)\le\sin(W/2)$ bounds the radius ratio.
- **Between neighbouring passes,** the floor of one window is at least $\pi-\Phi$ above the floor of the next, and at least $\pi-W-\Phi$ above its top. Only one of two neighbouring passes is on arm two.

## Theorem 387.1

**Unbalanced splits.** Write $G=\pi-\Phi$. v383's three flag arguments apply verbatim, with every lift window widened to $W_1$ and every floor-to-floor gap shrunk to $G$:
- **2–2 diagonal:** impossible when $W_1<G$.
- **3–1, lone corner outside:** impossible when $G>2W_1$ (the convexity argument needs top-to-floor gap above window width; this is a slightly stronger, simpler form).
- **3–1, lone corner inside:** impossible when $(G-W_1)\big(\Theta_0/(\Theta_0+G)\big)^{2a+1}>W_1$, where $\Theta_0$ is the triple's floor. The left side grows with $\Theta_0$.
- **2-1-1 adjacent:** impossible when the corner radii differ by less than $((G-W_1)/W_1)^{a/(2a+1)}$.

**One pass.** Arm two is the polar graph $\theta=\Theta+\pi+\varphi(\Theta^{-a})$, with $d\theta/d\Theta=1+O(aL/\Theta)>0$. Its curvature has the sign of $r^2+2r_\theta^2-rr_{\theta\theta}$, which is checked positive for the profile at every lift from 3 to $10^8$. The tangent turns by less than $W_1+\pi/2<\pi$ across the window, and v383's convexity argument applies.

**Aligned (2–2 adjacent).** v383's slope identity holds for any polar graphs. The chord bound for $Q=R_{\rm in}'(\eta_1)/R_{\rm out}'(\eta_2)$ is unchanged, except that the outer pitch bound picks up a factor $(1-aL/\Theta)^{-1}$. The spiral side becomes
$$Q\le\Big(\frac{\Theta_{\rm out}^{\max}}{\Theta_{\rm in}^{\min}}\Big)^{a+1}\frac{1+\mu_{\rm out}}{1-\mu_{\rm in}} .$$
The lift windows carry the $\Phi$ shift of whichever pass is on arm two. To leading order this gives the margin in the summary for odd pass gaps, where one pass is on arm two. Even gaps put both passes on the same arm, so $\Phi$ cancels down to the oscillation, and the margin is much larger.

For the profile at $W=0.3$:

| condition | requirement | holds from lift |
|---|---|---|
| diagonal, 3–1 outside | $W_1=0.90<G=2.74$, and $2W_1<G$ | all |
| 3–1 inside | monotone in the floor | 6.41 |
| aligned, every gap $n\le200$ | chord bound beats spiral bound | 2.44 (grid) |
| 2-1-1 adjacent | radius ratio below 1.195 | for $t\le0.126\,|Z|$ |

The aligned threshold comes from a check on a grid of 1 500 floors up to $10^9$, with the leading-order margin positive for every gap and parity. The monotonicity argument that makes v383's threshold exact carries over term by term but is not written out here. Overall, Theorem 387.1 holds for lifts $\ge6.41$ (radius $\le0.39$). ∎

## Completeness of v382's families

v382 showed that every resonance carries a square. The converse is what makes the classification complete: every small square is one of those. v382's branch for a label $(n_1,n_3)$ is built by solving, at each lift $\Theta_A$, for the side and angle that put $B$ and $D$ on their passes. So the converse needs that solve to have a unique solution in the small-square regime. Then any square with that label and that $\Theta_A$ is the branch square.

**Proposition 387.2 (numerical).** For each of the 81 squares from a v380 random search on the plain spiral (6 000 seeds), the following holds:
- choose a corner $A$ with both neighbours off its pass;
- map it to arm one by the half-turn if needed, and read off $(n_1,n_3)$;
- solve v382's three-corner system at its lift, once from the square and once from v382's generic straight-comb guess.

Both solves return the square's own side and angle, to $10^{-7}$. So in this sample, every square is its label's branch square, and the solve has no second solution near it.

A proof would show that the two-equation map $(\theta,t)\mapsto(\Phi(B)-\Phi(A),\Phi(D)-\Phi(A))$ is injective on the regime. Its Jacobian is a rotation-scaling by $1/(K(r)r)$ plus $O(t/r)$, so a Kantorovich-type argument should do it. That is left open.

## What this changes

- **v386's twisted bounds now rest on a proved classification** for twists of moderate gradient. That covers v374's default profile in the radius range above.
- **The gradient limit appears twice.** The aligned exclusion fails at about $L\approx5.9$, and v386's graceful confinement degrades as $(Ls)^{a/(1+a)}$. Both say the same thing: a twist that varies quickly across scales is where new squares can come from.
- **Open:** a written proof of uniqueness (Proposition 387.2), and of the monotone aligned threshold.

[Code](../code/twisted_classification_v387.py) · [Tests](../tests/test_twisted_classification_v387.py) · [Data](../results/v387/classification_check.json)
