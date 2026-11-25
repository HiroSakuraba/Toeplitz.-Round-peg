# v379. Squares on a smooth spiral: the selection rule

[v375](v375-the-comb-law.md) proved that every small square near a spiral point has a corner within the comb radius $\rho(t)\approx(t/(ag))^{a/(1+a)}$. The drawings suggested that squares really do reach out that far; that was the "lower half" left open. This note settles it, with a surprise. **On a smooth spiral the lower half is false.** Small squares sit much deeper, at radius about $t^{2a/(2a+1)}$ ($\sqrt t$ for $a=1/2$, against $t^{1/3}$). The squares that reach the comb radius in the drawings are all held in place by polygon vertices.

The reason is a selection rule. Every inscribed square obeys an exact identity: the alternating sum of the comb phase over its four corners equals an integral of the phase's Hessian over the square. On a spiral, that Hessian has a definite shape. Squares aligned with the turns would need the sum to vanish, which curvature forbids. Tilted squares need it to be a nonzero multiple of $\pi$, which happens only at special radii fixed by three integers. The smallest such square at radius $r$ has side about $r\sqrt{2\pi K(r)/(2+1/a)}$, far larger than the turn spacing.

**Results.**
1. **Exact identity (Lemma 379.1).** For an inscribed square $A,B=A+p,C=A+p+q,D=A+q$ on the half-turn spiral, with $\Phi(Z)=\pi k_Z$,
$$\pi j:=\Phi(A)-\Phi(B)+\Phi(C)-\Phi(D)=\int_0^1\!\!\int_0^1p^{\mathsf T}\,\nabla^2\Phi(A+up+vq)\,q\,du\,dv .$$
2. **Selection rule (Theorem 379.2).** With $n_1=k_B-k_A$ and $n_3=k_D-k_A$,
$$j=(2+\tfrac1a)\,\pi K(r)\,n_1n_3\,\big(1+O(t/r+K)\big).$$
So $j$ has the sign of $n_1n_3$, and $j=0$ only for squares aligned with the turns.
3. **No aligned squares (Proposition 379.3, numerical; proved in [v383](v383-small-squares-touch-three-turns.md)).** For aligned configurations the alternating sum keeps a fixed sign and stays strictly between consecutive integers. Checked for turn gaps 1–3 over $\Theta\in[8,90]$ with $a=1/2$.
4. **Smooth comb law (Corollary 379.4).** A square at radius $r$ has side at least $r\sqrt{2\pi K(r)/(2+1/a)}\,(1-o(1))$. Equivalently, squares of side $t$ sit within
$$r_{\rm smooth}(t)=\Big(\frac{(2+1/a)\,t^2}{2\pi a}\Big)^{a/(2a+1)}(1+o(1)),$$
which is $(4t^2/\pi)^{1/4}$ for $a=1/2$, and far inside $\rho(t)$.
5. **Exact and numerical checks:**
   - 200 inscribed squares on the smooth spiral, found by Newton's method: none aligned, all with $jn_1n_3>0$, none more than 0.7% beyond $r_{\rm smooth}$, all within $0.64\,\rho(t)$. The predicted radius is accurate only to first order in $t/r$. It depends on which corner is labelled $A$: labels with $n_1n_3<0$ do well, labels with $n_1n_3>0$ are off by about $1.3\,t/r$. (Corrected in [v382](v382-tilted-squares-exist.md); an earlier version quoted a median 0.5% from 7 squares.);
   - all 186 small squares of the three symmetric v374 drawings lie beyond $r_{\rm smooth}$, and every one has a corner within 10% of an edge's length from a polygon vertex (about half within 2%).

## The identity

On the half-turn spiral $\alpha=\{z:\Phi(z)=|z|^{-1/a}-\arg z\in\pi\mathbb Z\}$ (both arms of $r=\Theta^{-a}$), every corner of an inscribed square has $\Phi=\pi k$. Lift $\arg$ continuously across the square, so the $k$ are well defined. Since $C=A+p+q$, the alternating sum is the mixed second difference $\Delta_p\Delta_q\Phi(A)$. That equals the double integral of $p^{\mathsf T}\nabla^2\Phi\,q$ exactly (Lemma 379.1).

In the polar frame $(e_\psi,e_r)$, with $K(\rho)=a\rho^{1/a}$, $\nabla\Theta=-e_r/(K\rho)$ and $\nabla^2\rho=e_\psi e_\psi^{\mathsf T}/\rho$:
$$\nabla\Phi=-\frac{1}{K\rho}e_r-\frac1\rho e_\psi,\qquad
\nabla^2\Phi=\frac{1}{K\rho^2}\Big(-e_\psi e_\psi^{\mathsf T}+\big(1+\tfrac1a\big)e_re_r^{\mathsf T}\Big)+\frac1{\rho^2}\big(e_re_\psi^{\mathsf T}+e_\psi e_r^{\mathsf T}\big).$$
The first part is the curvature of the circles together with the change of turn spacing with radius. The second part comes from the angle.

## The selection rule

**Theorem 379.2.** Let a square of side $t$ have its corners on $\alpha$ at radii in $[r,r+\sqrt2t]$, with $t/r$ and $K(r)$ small. Then
$$\pi n_1=-\frac{t}{K r}(\sin\beta+K\cos\beta)(1+O(t/r)),\qquad\pi n_3=-\frac{t}{Kr}(\cos\beta-K\sin\beta)(1+O(t/r)),$$
$$\pi j=\frac{t^2}{r^2}\Big(\frac{(2+1/a)\sin2\beta}{2K}+\cos2\beta\Big)\big(1+O(t/r)\big),$$
where $\beta$ is the angle of $p$ from $e_\psi$ at $A$. Hence $j=(2+1/a)\pi K(r)\,n_1n_3\,(1+O(t/r+K))$.

*Proof outline.* The first two lines are $\pi n=\int_0^1\nabla\Phi(A+up)\cdot p\,du$ with the gradient above. Its relative variation across the square is $O(t/r)$, because $\nabla\Phi$ changes on the scale $r$. The third line is the identity with the Hessian above, using $p=t(\cos\beta,\sin\beta)$ and $q=t(-\sin\beta,\cos\beta)$; $\beta$ rotates by $O(t/r)$ across the square. Substituting $\sin\beta\approx-\pi Krn_1/t$ and $\cos\beta\approx-\pi Krn_3/t$ gives $\sin2\beta\approx2\pi^2K^2r^2n_1n_3/t^2$. The $\cos2\beta$ term is smaller by a factor $K$. ∎

The rule is a resonance condition. On a comb of straight, evenly spaced lines, squares with corners on lines $n_1$ and $n_3$ apart form continuous families, sliding along the lines. Curvature breaks each family, and a square survives only where the curvature term $(2+1/a)\pi K(r)n_1n_3$ equals an integer $j$. Along each arm that picks out one radius per $(n_1,n_3,j)$.

**Proposition 379.3 (no aligned squares; numerical).** If $j=0$ the rule forces $n_1n_3=0$, so one pair of sides runs along the turns. For such configurations put $A$ and $B$ on the same turn and $D=A+i(B-A)$ exactly on the turn $n_3$ levels away, leaving $C$ free. Then the alternating sum $m=j$ must be an integer. Numerically, $m$ keeps a fixed sign in each of the four orientations. It tends to $0$ as $r\to0$ but never reaches it, and $|m|<0.18$ throughout. This holds for $n_3=1,2,3$ and 500 values of $\Theta_A\in[8,90]$, about 5 700 configurations. So no aligned square closes up. The leading term $t^2\cos2\beta/r^2$ is outweighed by third-order terms of size $t^3/(Kr^3)$. Their sign is what the numerics record; a symbolic third-order computation would turn this into a proof. *(v383 proves it by a different route, the British flag identity, for outer lift $\Theta\ge6.78$ and squares spanning at most a radian: no square on the plain spiral has its corners on only one or two passes.)*

**Corollary 379.4 (smooth comb law).** Granting Proposition 379.3, every square has $|j|\ge1$. Then $n_1n_3\ge(1-o(1))/((2+1/a)\pi K)$, and
$$t^2=\pi^2K^2r^2(n_1^2+n_3^2)(1+O(\cdot))\ \ge\ 2\pi^2K^2r^2n_1n_3\ \ge\ \frac{2\pi K(r)\,r^2}{2+1/a}(1-o(1)).$$
Inverting gives $r_{\rm smooth}(t)$ above. Since $\frac{2a}{2a+1}>\frac a{1+a}$, this is far inside $\rho(t)$ as $t\to0$.

**The lower half.** For the smooth spiral the correct statement is now the existence of squares at the resonant radii. The numerics show them where predicted: for example the label $(n_1,n_3,j)=(-2,3,-1)$ predicts centre radius $0.1629$, and Newton finds a square with that label at centre $0.1621$ and side $0.0232$. (An earlier version gave this label as $(2,3,1)$; the square with label $(2,3,1)$ sits at $0.136$, see v382.) Squares near the minimal side exist whenever $n_1\approx n_3$ can be chosen. Existence at each resonance is proved, along explicit branches, in [v382](v382-tilted-squares-exist.md).

## Why the drawings looked different

A polygon has vertices, where the direction jumps. At a vertex, a corner of a square can take any direction in between. That is exactly the freedom an aligned square lacks on a smooth curve. So polygons carry aligned "kink squares" with side near the turn spacing, which reach out to the comb radius $\rho(t)$. In the three symmetric v374 drawings, every one of the 186 small squares has a corner within $0.102$ of a vertex in edge parameter, and 94 within $0.02$. Every one also sits beyond $r_{\rm smooth}$. The drawings stop at radius $0.25$, while the smooth squares of side $\le1/20$ live inside $r_{\rm smooth}(1/20)=0.24$.

So the v375 bound is sharp for polygons, as the drawings showed, but not for smooth spirals.

## What this changes

- **v375 stands as an upper bound.** It is attained up to a constant by polygonal spirals and is loose for smooth ones.
- **The data of v370–v378 are kink data.** The parity theorems for symmetric points do not depend on where squares sit, so they are unaffected. But the asymmetric drawings' small squares, their events (v374) and the surgery trials (v378) all involve kink squares. The surgery rule itself is general.
- **Regularity matters near a wild point.** A smooth spiral core confines squares of side below $s$ to radius about $s^{2a/(2a+1)}$. A polygonal core with vertices reaches $s^{a/(1+a)}$. Theorem 370.2 assumes the curve is polygonal away from $p$. Extending it to piecewise-smooth curves would let the smooth law shrink the radius the v374 path condition must reach, from about $8\rho(s)$ to about $8\,r_{\rm smooth}(s)$. That radius is still much larger than $s$, but it is smaller.

[Code](../code/smooth_squares_v379.py) · [Tests](../tests/test_smooth_comb_v379.py) · [Data](../results/v379/smooth_comb_check.json)
