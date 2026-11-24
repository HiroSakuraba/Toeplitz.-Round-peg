# v377. Invisibility collapses near spiral points

[v374](v374-asymmetric-wild-points.md) reduced an asymmetric two-arm point to a path condition: symmetrize the point out to radius about $8\rho(s)$ without passing a square of side $s$. [v375](v375-the-comb-law.md) and [v376](v376-the-comb-law-for-polygons.md) pinned $\rho(s)\approx(s/(ag))^{a/(1+a)}$. This note takes stock of the other end of the scale: how small a change near the wild point already counts at cutoff $s$.

The answer is: almost none can be ignored. Near a spiral point there are squares of side exactly $s$ with three corners on the spiral and the fourth within about one turn spacing of the point. So the radius below which deformations are invisible at cutoff $s$ is at most about $d(s)=a\,g\,s^{1+1/a}$. That is far below $s$, and further still below the comb radius. Two consequences follow. Finite drawings cannot certify the wild curve's parity at a fixed cutoff unless they reach depth of order $d(s)^{-1/a}$. And the path condition cannot be localized: side-$s$ events can occur anywhere from radius $d(s)$ out to $\rho(s)$.

**Results.**
1. **Near-squares at the wild point (Proposition 377.1, proof sketch).** For every small $s$ there is a square of side $s$ with three corners on the arc and the fourth within $\kappa\,d(s)$ of $p$, where $\kappa\approx(2\pi-g)/(\sqrt2\,g\,(1-2^{-1/(2a)}))$. That is $\sqrt2$ for the half-turn spiral with $a=1/2$.
2. **Invisibility radius (Corollary 377.2).** If a Jordan curve containing the arc has no square of side $s$, then $\varepsilon_{\rm inv}(s)\le\kappa\,d(s)(1+o(1))$. This is much smaller than the general cap $s/\sqrt2$ from v374.
3. **The scale ladder.**
$$\varepsilon_{\rm inv}(s)\lesssim s^{1+1/a}\ \ll\ s\ \ll\ \rho(s)\sim s^{a/(1+a)},$$
and the symmetrizing path of Theorem 374.3 must reach about $8\rho(s)$.
4. **Exact and numerical checks:**
   - Newton's method finds the near-squares of item 1 on both arcs at $s=0.05,0.02,0.01$, with free corner between $0.32\,d(s)$ and $1.16\,d(s)$ from $p$ and arc residuals below $2\times10^{-12}$;
   - on 24 asymmetric drawings, the exact large graceful parity at a fixed cutoff changes between depths 5 and 7 in 9 cases, and is odd at depth 7 in exactly half. At these depths the parity describes the drawing, not the wild curve.

## Setting

Use the arcs of v375: $\gamma_1(\Theta)=\Theta^{-a}e^{i\Theta}$ and $\gamma_2(\Theta)=\Theta^{-a}e^{i(\Theta+\pi+\varphi)}$ for $\Theta\ge\theta_0$, with $|\varphi|\le\bar\varphi<\pi$ and $g=\pi-\bar\varphi$. Write
$$d(r)=a\,g\,r^{1+1/a},$$
the (leading-order) turn spacing at radius $r$. The comb phase $\Phi(z)=|z|^{-1/a}-\arg z$ (lifted) satisfies: $z\in\alpha$ iff $\Phi(z)\equiv0$ or $\Phi(z)\equiv-\pi-\varphi(|z|)\pmod{2\pi}$. Its gradient is dominated by the radial term,
$$\nabla\Phi=-\frac1{K(r)\,r}\,e_r-\frac1r\,e_\psi,\qquad K(r)=a\,r^{1/a},$$
so $\Phi$ changes by $2\pi$ over a radial step of about $2\pi K r$. The zero set has gaps of at most $2\pi-g$ in $\Phi$.

## Near-squares at the wild point

**Proposition 377.1.** For all sufficiently small $s$, there is a square of side $s$ whose corners $B,C,D$ lie on $\alpha$ and whose fourth corner $A$ satisfies $|A|\le\kappa\,d(s)(1+o(1))$, with
$$\kappa=\frac{2\pi-g}{\sqrt2\,g\,\big(1-2^{-1/(2a)}\big)}.$$

*Proof sketch.* Fix an angle $\theta$ and put $B=A+se^{i\theta}$, $D=A+ise^{i\theta}$, $C=A+(1+i)se^{i\theta}$. Use rotating coordinates $A=(x+iy)e^{i\theta}$. For $|A|\le\eta\ll s$:

*Step 1 (two corners).* Up to relative errors $O(\eta/s+K(s))$, $\Phi(B)$ depends on $x$ alone and $\Phi(D)$ on $y$ alone, each with slope $-1/(K(s)s)$. So $A\mapsto(\Phi(B),\Phi(D))$ is a near-similarity of the $\eta$-disk onto a disk of radius about $\eta/(K(s)s)$. Once that radius exceeds the gap $2\pi-g$, some $A$ with $|A|=O(K(s)s)=O(d(s))$ puts both $B$ and $D$ on $\alpha$. Following one such solution as $\theta$ varies gives a branch $A(\theta)$.

*Step 2 (drift).* On the branch, $\arg B$ increases with $\theta$, so $|B|$ must decrease to keep $\Phi(B)$ fixed: $\Delta x=-K(s)s\,\Delta\theta$, and likewise $\Delta y=-K(s)s\,\Delta\theta$. So $|A|$ drifts at rate $\sqrt2K(s)s=\sqrt2\,d(s)/g$ per radian, and a branch starting near $p$ stays within $\eta$ for a $\theta$-interval of length about $\sqrt2\,g\,\eta/d(s)$.

*Step 3 (third corner).* Along the branch $|C|-s\sqrt2\approx(x+y)/\sqrt2$ decreases at rate $\sqrt2K(s)s$, so $|C|^{-1/a}$ increases at rate $K(s)/K(s\sqrt2)=2^{-1/(2a)}$, while $\arg C$ increases at rate 1. Hence
$$\frac{d\Phi(C)}{d\theta}=2^{-1/(2a)}-1+o(1)\neq0.$$
Over the branch, $\Phi(C)$ moves through $(1-2^{-1/(2a)})\sqrt2\,g\,\eta/d(s)$. Once that exceeds the gap $2\pi-g$, $C$ crosses $\alpha$. This happens for $\eta=\kappa\,d(s)(1+o(1))$. ∎

The sketch's error terms are $O(\eta/s+K(s))=O(s^{1/a})$. Writing them out is routine but not done here. The point of the sketch is the scale $d(s)$. The constant is a worst case: the numerics find the free corner closer than $\kappa\,d(s)$.

**Corollary 377.2.** Let $\gamma$ be a Jordan curve with $\gamma\cap\bar B(p,R)=\alpha$ and no inscribed square of side $s$. Then $\varepsilon_{\rm inv}(s)\le\kappa\,d(s)(1+o(1))$.

*Proof.* For $\varepsilon>|A|$, reroute $\gamma$ inside $B(p,\varepsilon)$ through $A$. The result is a closed curve in $\gamma\cup B(p,\varepsilon)$ with an inscribed square of side $s$. So $\varepsilon\ge\varepsilon_{\rm inv}(s)$. ∎

## What this means

**For the path condition.** Theorem 374.3 needs a symmetrizing path supported out to about $8\rho(s)$ that passes no square of side $s$. Blanket invisibility reaches only $\kappa\,d(s)$. For $a=1/2$ the two radii scale as $s^{1/3}$ and $s^3$. Between them lies the whole comb, and side-$s$ events can occur anywhere in it, including at radius about $s$ (squares straddling the point). So a proof of the path condition has to control events globally: either count them and show they pair up, or find a different invariant. It cannot argue that the path is too small to matter.

**For computations.** By Lemma 374.2 and Corollary 377.2, a depth-$n$ drawing has the wild curve's large degree at cutoff $s$ once its inner truncation lies inside $B(p,\varepsilon_{\rm inv}(s))$, that is, at radius $\lesssim d(s)$. For $s=1/20$ and $a=1/2$ this means $r\lesssim2\times10^{-4}$, or $\Theta$ beyond $10^7$. The exact drawings used so far reach $\Theta\approx25$. The symmetric results of v370–v372 are unaffected: there the parity is forced by pairing at every depth. But for asymmetric drawings, a fixed-cutoff parity is a property of the drawing only.

**The depth scan confirms this.** Asymmetric drawings with four choices of $\varphi$ and six cutoffs $s=1/16,\dots,1/36$, at depths 5 and 7. Each entry shows the parity at depth 5, then at depth 7:

| $\varphi$ | 1/16 | 1/20 | 1/24 | 1/28 | 1/32 | 1/36 |
|---|---|---|---|---|---|---|
| $0.4\sin(5\log\frac1r)$ | 1 1 | 1 1 | 0 0 | 1 1 | 1 1 | 1 1 |
| $1.0\sin(5\log\frac1r)$ | 0 0 | 0 1 | 0 0 | 1 0 | 1 1 | 0 0 |
| $1.0\sin(\log\frac1r+1)$ | 0 0 | 1 0 | 1 1 | 1 0 | 1 1 | 0 1 |
| $0.7$ (constant) | 1 0 | 1 0 | 0 1 | 0 1 | 0 0 | 0 0 |

The parity changes with depth in 9 of 24 cases, and at depth 7 it is odd in 12 of 24. Two further choices, $2\sin(\log\frac1r+1)$ and the constant $2$, give drawings whose outer closure crosses the arm, so the simplicity check rejects them. Crossed-sector squares, absent at depth 5 for small cutoffs, appear in quantity at depth 7. That suggests they reach toward the inner end of the drawing.

## Where this leaves asymmetric points

The program now has three proved scales near a two-arm spiral point: the invisibility radius $\lesssim s^{1+1/a}$, the general cap $s/\sqrt2$, and the comb radius $\sim s^{a/(1+a)}$. Symmetric points escape the problem because pairing works at every scale at once. For asymmetric points, the open step is a global count of the side-$s$ squares that a symmetrizing path creates and destroys across the comb. Proposition 377.1 suggests some structure for that count. Along a branch, near-squares are organized by a single phase that moves monotonically with the angle, so side-$s$ squares near the point come in regular families rather than at random.

[Code](../code/invisibility_v377.py) · [Tests](../tests/test_invisibility_v377.py) · [Data](../results/v377/invisibility_check.json)
