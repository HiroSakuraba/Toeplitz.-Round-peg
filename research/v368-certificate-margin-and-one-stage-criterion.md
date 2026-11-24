# v368. The margin of a sweep certificate, and a one-stage criterion

[v366](v366-sweeps-read-from-the-curve.md) reported a *clearance* for each certificate but left open how large a redraw it tolerates. This note defines an exact margin, proves that every registered Jordan redraw below the margin keeps the certificate, and computes the margin with rigorous rational brackets. Combined with [v359](v359-fixed-arcs-through-refinement.md), it yields a criterion that needs only one stage of an approximation, and the criterion is sharp. A curve satisfies it exactly when the curve itself has an odd three-stretch certificate whose box boundary carries no square. The circle never satisfies it.

**Results.**
1. **Tolerance (Lemma 368.1).** Every registered Jordan redraw within uniform distance $\delta<\mu$ keeps the certificate's degree.
2. **The margin is 1-Lipschitz (Lemma 368.2)** under registered redraws.
3. **One-stage criterion (Theorem 368.3).** Suppose a polygon within $\delta$ of a Jordan curve $\gamma$ has an odd certificate with $\mu>\delta$. Then $\gamma$ has a square with three corners in the three stretches, and its side is bounded below by the stretch distances on $\gamma$. Conversely, an odd zero-free certificate box on $\gamma$ gives such polygons at every scale below half its margin.
4. **Circle (Corollary 368.4).** On every approximation of a circle, every odd certificate has $\mu\le\delta$.
5. **Exact numbers.** The rigorous brackets are $[0.0096337,0.0096341]$ for the roof square, $[0.0107747,0.0107749]$ for the decagon, and $[0.027782,0.027784]$ for the v367 union box. On the bumped octagon the margin is about $0.062$ times the bump, which is far below the distance to the circle.

## The margin

Use the notation of v366: stretches $A,B,E$, orientation $\sigma$, box $\mathcal B=A\times B\times E$, and
$$
F_1=\gamma(b)+\sigma J(\gamma(b)-\gamma(a))-\gamma(e),\qquad D=\gamma(a)+\sigma J(\gamma(b)-\gamma(a)),\qquad F=(F_1,\rho(D)).
$$
Define
$$
\mu=\min_{\partial\mathcal B}\ \max\!\Big(\frac{|F_1|}{4},\ \frac{|\rho(D)|}{2+\sqrt2}\Big).
$$
It is positive exactly when $\partial\mathcal B$ carries no zero of $F$.

**Lemma 368.1 (tolerance).** Let $\gamma'$ be a Jordan curve with $|\gamma'(t)-\gamma(t)|\le\delta<\mu$ for every $t$. Then $F'$ has no zero on $\partial\mathcal B$, and its degree equals that of $F$.

*Proof.* Each corner moves by at most $\delta$, so $|F_1'-F_1|\le\delta+2\delta+\delta=4\delta$. Since $D=(I-\sigma J)\gamma(a)+\sigma J\gamma(b)$ and $\lVert I-\sigma J\rVert=\sqrt2$, the fourth corner moves by at most $(1+\sqrt2)\delta$.

*Step 1.* Move the corners along $\gamma_s=(1-s)\gamma+s\gamma'$, and keep the signed distance $\rho$ of the original curve. The curves $\gamma_s$ need not be embedded, but $\rho$ is fixed, so this does not matter. At a point of $\partial\mathcal B$, either $|F_1|>4\delta$, so the first two coordinates never vanish, or $|\rho(D)|>(2+\sqrt2)\delta$. In the second case $|\rho(D_s)|>\delta>0$, since $\rho$ is 1-Lipschitz.

*Step 2.* With the corners on $\gamma'$, replace $\rho$ by $(1-s)\rho+s\rho'$. Where the first two coordinates can vanish, $|\rho(D')|>\delta$. A point at distance more than $\delta$ from $\gamma$ is never crossed by $\gamma_s$, so its winding number is the same for $\gamma$ and $\gamma'$. Thus $\rho$ and $\rho'$ have the same strict sign there, and the convex combination does not vanish.

The boundary stays zero-free throughout, so the degree is unchanged. ∎

**Lemma 368.2.** Under the same redraw, $|\mu'-\mu|\le\delta$.

*Proof.* We have $|F_1'|\ge|F_1|-4\delta$. Also $\operatorname{dist}(D',\gamma')\ge\operatorname{dist}(D,\gamma)-(1+\sqrt2)\delta-\delta$. After dividing by 4 and by $2+\sqrt2$, each term changes by at most $\delta$. The same holds for their maximum, and then for the minimum over $\partial\mathcal B$. The roles of the two curves are symmetric. ∎

## Computing it exactly

On a polygon, each face of $\mathcal B$ fixes one corner at a stretch end and lets the other two run over edge pieces. On each product of pieces, $F_1$ and $D$ are affine in the two local edge parameters. On a rectangle of those parameters the checker bounds:
- $|F_1|_2$ from below by its value at the centre, in the Euclidean or the max norm, minus an $\ell^1$ bound on its variation;
- $\operatorname{dist}(D,\gamma)$ from below by the exact squared distance at the centre, through an integer square root, minus an $\ell^1$ bound on how far $D$ moves.

It divides by $2+141422/100000>2+\sqrt2$ and refines the worst rectangle first. The minimum over the leaves is a rigorous lower bound. The value at the centre of the worst leaf, and at every outer locus end (where $F_1=0$), gives an upper bound.

| certificate | stretches | $\mu$ bracket after 8000 refinements |
|---|---|---|
| v364 roof square | bins 1, 7, 11 | $[0.0096337,\ 0.0096341]$ |
| v247 decagon | bins 7, 3, 10 | $[0.0107747,\ 0.0107749]$ |
| v367 union box on the decagon | $\{6,7\},\{1,2,3\},\{10,11,12\}$ | $[0.027782,\ 0.027784]$ |

The minimum is not at a locus end. On the decagon it sits where both terms are moderate: $|F_1|\approx0.043$ and $D$ about $0.037$ from the curve. The locus ends alone would give $0.095/(2+\sqrt2)\approx0.028$, so the end clearance overstates the tolerance. The tests move every decagon vertex at random so that the registered distance stays below the lower bound. The certificate keeps degree one in every trial.

**v367 in these terms.** On the final greedy polygon, the union box has margin about $0.0252$. The bump that destroyed it moves the curve by $0.4216$ in registered distance, about 17 times the margin. The bumps that only moved parity changed the union box's faces not at all. The lemma is sufficient, not necessary: the start and final polygons are $0.19$ apart, and the union box kept its degree anyway.

## A criterion from one stage

**Theorem 368.3.** Let $\gamma$ be a registered Jordan curve and $\Gamma$ a registered polygon with $\sup_t|\Gamma(t)-\gamma(t)|=\delta$. Suppose a certificate on $\Gamma$ has stretches with pairwise disjoint closures, odd degree and margin $\mu_\Gamma>\delta$. Then:
- $\gamma$ contains a square with $P,Q,R$ in the three stretches;
- its side is at least the largest of $d_\gamma(A,B)$, $d_\gamma(B,E)$ and $d_\gamma(A,E)/\sqrt2$, which is positive;
- every Jordan curve within $\mu_\Gamma-\delta$ of $\gamma$ also contains such a square.

Conversely, suppose $\gamma$ has an odd certificate box whose boundary carries no square. Then its margin $\mu_\gamma$ is positive. Every registered polygon within $\delta<\mu_\gamma/2$ of $\gamma$ has the same certificate, with odd degree and margin greater than $\delta$.

*Proof.* Lemma 368.1 applies to any Jordan redraw, $\gamma$ included. The side bound uses that $P,Q$ and $Q,R$ are adjacent and $P,R$ are opposite, together with injectivity of $\gamma$.

For the converse: on the compact set $\partial\mathcal B$ the two terms are continuous, nonnegative and never vanish together, so $\mu_\gamma>0$. Lemma 368.2 then gives $\mu_\Gamma\ge\mu_\gamma-\delta>\delta$, and Lemma 368.1 gives equal degrees. ∎

So, earlier I said the margin must beat *twice* the error at fine stages. It only needs to beat the error once. But the theorem also says this route proves exactly the curves that already carry an odd, zero-free certificate box. A curve whose inscribed square is a regular zero of $F$ qualifies, with a small box around that square. Curves whose squares are all degenerate in this sense need something else.

**Corollary 368.4 (circle).** Let $\Gamma$ be any registered polygon within $\delta$ of a circle. Every odd certificate on $\Gamma$ has $\mu_\Gamma\le\delta$.

*Proof.* Otherwise Theorem 368.3 would give an odd box with zero-free boundary on the circle. On the circle, the squares with $P$ at parameter $a$ form a single continuous family over all $a$. Disjoint stretches exclude collapsed squares. A family with a square inside the box leaves it as $a$ leaves $A$, so it crosses $\partial\mathcal B$. That is [v361](v361-blind-three-window-selection.md)'s argument. ∎

On the symmetric octagon with one vertex pushed out by a factor $1+b$, the bumped vertex is $b$ from the unit circle. So every registration has $\delta\ge b$. The best certificate margin is about $0.062\,b$ for $b=1/64,1/256,1/1024$, consistent with the corollary by a wide margin.

## What this changes

Noncollapse along refinements is now a property of the target alone. This method succeeds on $\gamma$ exactly when $\gamma$ has an odd three-stretch box whose boundary carries no square. It fails exactly when every square lies in a family that crosses the boundary of every such box, as on the circle. That is precisely the situation [v362](v362-relabeling-path-noncollapse.md) is built for: there a family that comes back to itself under a change of starting corner forces separation. The natural next lemma is a dichotomy. On every Jordan curve, either some three-stretch box is odd with zero-free boundary, or some square family is invariant under a nontrivial relabeling. Either alternative would then give a square. Neither alternative is proved for general curves, and the dichotomy itself is open.

[Code](../code/certificate_margin_v368.py) · [Tests](../tests/test_certificate_margin_v368.py)
