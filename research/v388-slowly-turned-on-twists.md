# v388. Symmetric windows on smooth arms, with a slowly turned-on twist

The second of three follow-ups. [v386](v386-smooth-cores.md) found that on smooth arms the twist's size only frees crossed squares. Graceful squares, the only ones the parity sees, escape only through the twist's gradient $L=\sup|r\varphi'|$. This note turns that into a protection theorem, by redoing v385's symmetric-window argument with the smaller reach that smooth arms allow. It also tests the graceful bound for the first time: v386 had found no graceful twisted squares to check it on.

**Results.**
1. **Graceful confinement radius.** On twisted smooth arms with $|\varphi|\le\Phi$ and $|r\varphi'|\le L$, every graceful square of side $t$ has a corner within
$$r_g(t)=\max\Big(r_{\rm smooth}(t)\,(1-\Phi/\pi)^{-a/(2a+1)},\ \big(L\,t/c_L\big)^{a/(1+a)}\Big),\qquad c_L=\frac{(2a+1)\pi^2}{2\sqrt2},$$
   to leading order (Proposition 386.5, with the classification of v387).
2. **Symmetric window, smooth arms (Theorem 388.1).** Suppose:
   - the twist vanishes on $[\rho,\varepsilon]$;
   - outside, it has size at most $\Phi$ and gradient at most $L$;
   - for some cutoff, $8\,r_g(s)\le\varepsilon$ and $\rho<\min(s/2,\eta^*(s))$.

   Then the point inscribes a square of side greater than $s$.
3. **What the gradient costs (Corollary 388.2).** Take the largest allowed cutoff.
   - If $L$ is below a crossover $L^*\asymp(\varepsilon/8)^{1/(2a)}$ (about $1.05\,\varepsilon$ for $a=1/2$), it is $s\asymp\varepsilon^{(2a+1)/(2a)}$, and the window may start as deep as $\rho\asymp\varepsilon^{(1+1/a)(2a+1)/(2a)}$. That is $\varepsilon^6$ for $a=1/2$, against $\varepsilon^9$ for polygonal arms (v385).
   - For larger $L$ it is $s\asymp\varepsilon^{(1+a)/a}/L$, giving $\rho\asymp\varepsilon^9/L^3$. Same exponent as polygons, but a better constant.
4. **Checks.** Newton searches find graceful squares on smooth twisted arms at one, two and three times v374's twist (13, 9 and 15 of them).
   - Every one lies inside $r_g$: largest ratio $0.81$, $0.61$, $0.53$.
   - In fact every one lies inside even the **plain** spiral's $r_{\rm smooth}$ (largest ratio $0.91$, $0.79$, $1.000$).
   - The crossed squares of the same searches reach $1.40$–$1.53$ times $r_{\rm smooth}$.

   So graceful squares are confined like the symmetric spiral's, and the twist's size does free only crossed ones.

## Theorem 388.1

**Setting.** The arms are $C^2$ away from $p$, as in Theorem 386.4, with twist $\varphi$. The twist is 0 on $[\rho,\varepsilon]$, and $|\varphi|\le\Phi<\pi/2$, $|r\varphi'|\le L$ on $[\varepsilon,R]$. These hypotheses include the stretch where the twist rises from 0, so a slow turn-on means a small $L$ there.

*Proof.*
- If $\gamma$ has a square of side $s$, we are done.
- Otherwise, run the gap path of v385, supported inside $B(p,\rho')$ with $\rho<\rho'<\min(s/2,\eta^*(s))$. By Corollary 385.2 it passes no square of side $s$; its end curve $\gamma_1$ is symmetric on $B(p,\varepsilon)$.
- Theorem 386.4 needs every graceful square of side below $s$ on $\gamma_1$, and on its smooth-cap truncations, to have a corner in $B(p,\varepsilon/8)$. Squares inside the symmetric ball obey Theorem 386.2. Squares reaching into the twisted part obey Proposition 386.5, whose pattern constraints v387 supplies in the range considered. Both give a corner within $r_g(s)\le\varepsilon/8$.
- Theorem 386.4, which is Theorem 374.3 on smooth arms, finishes the proof. ∎

The leading-order form of $r_g$ is the same level of precision as v386's Proposition 386.5. The explicit inequalities behind it are written out there.

## Corollary 388.2

$r_g$ increases with $t$, so the best cutoff is the largest $s$ with $8r_g(s)\le\varepsilon$, and $\rho_{\max}(\varepsilon)=\min(s/2,\eta^*(s))$. For $a=1/2$ and $\Phi=0.4$:

| $\varepsilon$ | $L$ | best $s$ | $\rho_{\max}$ (smooth) | $\rho_{\max}$ (polygons, v385) | $\log\rho_{\max}/\log\varepsilon$ |
|---|---|---|---|---|---|
| 0.1 | 0–0.1 | $1.3\times10^{-4}$ | $1.2\times10^{-12}$ | $1.6\times10^{-17}$ | 11.9 |
| 0.1 | 2 | $6.8\times10^{-6}$ | $1.8\times10^{-16}$ | $1.6\times10^{-17}$ | 15.8 |
| 0.01 | 0–0.01 | $1.3\times10^{-6}$ | $1.2\times10^{-18}$ | $1.6\times10^{-26}$ | 8.96 |
| 0.01 | 1 | $1.4\times10^{-8}$ | $1.4\times10^{-24}$ | $1.6\times10^{-26}$ | 11.9 |
| $10^{-5}$ | 0 | $1.3\times10^{-12}$ | $1.2\times10^{-36}$ | $1.6\times10^{-53}$ | 7.18 |
| $10^{-5}$ | 0.01 | $1.4\times10^{-15}$ | $1.4\times10^{-45}$ | $1.6\times10^{-53}$ | 8.97 |
| $10^{-5}$ | 2 | $6.8\times10^{-18}$ | $1.8\times10^{-52}$ | $1.6\times10^{-53}$ | 10.4 |

With a vanishing gradient the exponent tends to 6 (it is already 7.18 at $\varepsilon=10^{-5}$). Any fixed $L>0$ eventually brings back the polygonal exponent 9, with a gain of order $L^{-3}$. The crossover is at $L^*=\frac{c_L}{c_2\sqrt{1-\Phi/\pi}}\,(\varepsilon/8)^{1/(2a)}$, which is $L^*\approx1.05\,\varepsilon$ for $a=1/2$ and $\Phi=0.4$. For example, at $\varepsilon=0.01$ the rows for $L=0$ and $L=0.01$ coincide, and $L=0.1$ is past the crossover.

## Graceful squares on twisted arms

v386's search was aimed at $j\approx0$, where crossed squares live. Here a random search places seeds on both arms, with random orientation and side between 0.7 and 8 turn spacings at radius 0.08–0.3. Corners are assigned to their nearest pass and solved by Newton's method. The twist is $\lambda\cdot0.4\sin(5\log(1/r))$ for $\lambda=1,2,3$, so $\Phi=0.4\lambda$ and $L=2\lambda$.

| $\lambda$ | squares | graceful | max $r_{\min}/r_g$ (graceful) | max $r_{\min}/r_{\rm smooth}$, graceful | same, crossed |
|---|---|---|---|---|---|
| 1 | 32 | 13 | 0.81 | 0.91 | 1.40 |
| 2 | 24 | 9 | 0.61 | 0.79 | 1.51 |
| 3 | 42 | 15 | 0.53 | 1.000 | 1.53 |

Every graceful square sits within the plain spiral's smooth radius. The largest, at $\lambda=3$, is at $0.99997$, essentially on it. The crossed squares go well beyond. So in this sample even the gradient term of $r_g$ is never needed, and $r_g$ looks conservative. The sample is small: 37 graceful squares, from 9 000 seeds.

## What this changes

- **A new class of protected points:** smooth two-arm points that are symmetric on one window and whose twist outside it turns on slowly. For a slow enough twist, the window may start at $\varepsilon^6$ (for $a=1/2$), not $\varepsilon^9$.
- **The data suggest more.** Graceful squares on twisted smooth arms behaved exactly like the symmetric spiral's, even at three times v374's twist. If that holds in general, the gradient term is an artifact of the proof, and the path condition's reach on smooth arms is $r_{\rm smooth}(s)$ whatever the twist.
- **Still open:** points twisted at every scale, with no symmetric window. There the path itself creates and destroys graceful squares of side $s$, which is the subject of the third follow-up.

[Code](../code/slow_twist_v388.py) · [Tests](../tests/test_slow_twist_v388.py) · [Data](../results/v388/slow_twist_check.json)
