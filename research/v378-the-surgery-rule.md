# v378. The surgery rule and the near-square lattice

[v377](v377-invisibility-collapses-near-spirals.md) found that near a spiral wild point there are squares of side $s$ with three corners on the spiral and the fourth within about one turn spacing of the point. It also observed that these near-squares come in regular families. This note makes both observations exact, and turns them into a working tool.

The tool is a rule for how the large graceful parity changes under a small local change of the curve. Change the curve inside a ball smaller than $s/2$. The parity then changes by the number of graceful near-squares whose free corner the change sweeps over. Nothing else matters. The rule is proved in general and checked exactly on 42 bumps of the v374 drawings, including the cases where a swept near-square should not count and where two should cancel.

Near a spiral point, the near-squares turn out to form a lattice: the image of $2\pi\mathbb Z^2$ under an explicit, nearly affine map. Their number within a given radius follows a simple area law, and exactly half are graceful, by a one-line rule. Together with the surgery rule, this gives an exact description of what symmetrizing the innermost part of an asymmetric point does to the parity. It also shows why no cancellation is forced there.

**Results.**
1. **Surgery rule (Theorem 378.1).** If Jordan curves $\gamma,\gamma'$ agree outside a ball of radius below $s/2$, then
$$L(\gamma',s)+L(\gamma,s)\equiv\#\{\text{graceful near-squares }n:\ A_n\in\operatorname{int}\gamma\,\triangle\,\operatorname{int}\gamma'\}\pmod 2.$$
Lemma 374.2 is the special case of an empty count.
2. **Lattice of near-squares (Lemma 378.2, Theorem 378.3).** At a two-arm spiral point, near-squares of side $s$ with free corner within $M\,K(s)\,s$ of $p$ correspond one-to-one to lattice points of a map $P$ with $DP\approx\begin{pmatrix}-1&1\\-(1-c)&c\end{pmatrix}$, where $c=2^{-1/(2a)-1}$. Their number is $\frac{2|1-2c|}\pi M^2(1+o(1))$, which is $M^2/\pi$ for $a=1/2$. This does not depend on the asymmetry $\varphi$.
3. **Sector rule (Lemma 378.4).** A near-square with free corner at the wild point is graceful exactly when the free corner's two neighbours lie on different arms. So half are graceful.
4. **Symmetrizing the core (Corollary 378.5).** Turning arm two inside $B(p,r_1)$, $r_1<s/2$, changes the parity by the number of graceful lattice near-squares lying in the thin spiral bands that arm two sweeps. Those bands are thinner than the lattice spacing, so their parity is not forced.
5. **Exact checks:**
   - 42 bumps on three drawings: 15 sweep one graceful near-square and flip the parity; 7 sweep one crossed near-square and do not; 7 sweep two and do not; 1 sweeps three and flips; 12 random bumps sweep none and do not flip. All 42 agree with the rule;
   - lattice enumeration on the continuous arcs: 322–326 near-squares within $M=32$ against $325.9$ predicted, for $s=0.02,0.01$, symmetric and asymmetric. The sector rule holds for every one, and the nearest ones agree with v377's independent solver.

## The surgery rule

Fix $s$. A **near-square** for a ball $B$ of radius $r_1<s/2$ is a square of side $s$ with three corners $B',C,D$ on $\gamma\setminus B$ and the fourth, the free corner $A$, in $B$. Suppose $\gamma\cap\bar B$ is a single arc $\beta$. Then all points of $\beta$ occupy one interval of the curve's cyclic order. So whether the square is aligned, reversed or crossed, with $A$ placed on $\beta$, depends only on $B',C,D$. Call the near-square **graceful** if it is aligned or reversed.

**Theorem 378.1 (surgery rule).** Let $\gamma,\gamma'$ be Jordan curves that agree outside a ball $B$ of radius $r_1<s/2$, each meeting $\bar B$ in a single arc with the same endpoints. Suppose neither has a square of side exactly $s$, and the near-squares of $\gamma$ for $B$ are nondegenerate (the generic case). Then
$$L(\gamma',s)-L(\gamma,s)\equiv\#\{n\text{ graceful}:A_n\in\operatorname{int}\gamma\,\triangle\,\operatorname{int}\gamma'\}\pmod2.$$

*Proof.* Join the two inner arcs by a generic homotopy rel endpoints inside $\bar B$. The intermediate curves may be non-simple, which is allowed for $L$ (v371, v373).

Two corners of a side-$s$ square are at least $s>2r_1$ apart, so at most one corner can lie in $B$. A side-$s$ square on an intermediate curve with no corner in $B$ would be a square of $\gamma$, which is excluded. So every such square is a near-square with its free corner on the moving arc. The near-squares do not move, since their other corners lie on the fixed part.

Near-square $n$ therefore becomes an inscribed square exactly when the moving arc passes over $A_n$. For a generic homotopy each passage is transverse. Nondegeneracy then makes the inscribed square cross the cutoff transversally. By v373 each such crossing changes $L$ if the square is graceful, and leaves it alone if it is crossed.

The number of passages over $A_n$ has the parity of the change in the winding number of the curve around $A_n$. For Jordan curves that is the inside/outside indicator. ∎

Two remarks. Lemma 374.2 is the case where no graceful near-square has its free corner in the swept region. And the rule is local in both directions. Any deformation of a curve can be cut into moves supported in balls of radius below $s/2$, so the total change of $L$ is a sum of such counts.

## The near-square lattice at a spiral point

Use the arcs of v375 and write $K(s)=a\,s^{1/a}$, so that $K(s)\,s=d(s)/g$. Put the free corner at $A=(x+iy)e^{i\theta}$, with $B'=A+se^{i\theta}$, $D=A+ise^{i\theta}$, $C=A+(1+i)se^{i\theta}$. Scale $X=x/(K(s)s)$ and $Y=y/(K(s)s)$.

For a corner $Q$ and arm $j$, let $h_Q^{(j)}$ be the comb phase $|Q|^{-1/a}-\arg Q$ computed in the rotating frame, plus $\pi+\varphi(|Q|)$ if $j=2$. Then $Q$ lies on arm $j$ iff $h_Q^{(j)}-\theta\in2\pi\mathbb Z$. The angle $\theta$ enters all three conditions in the same way, with coefficient $-1$, and the radii do not depend on it. This is the regular-family structure of v377.

**Lemma 378.2 (lattice map).** For each arm triple $(j_B,j_D,j_C)$, near-squares are in bijection with solutions of
$$P(X,Y)=\big(h_{B'}-h_D,\ h_{B'}-h_C\big)\in2\pi\mathbb Z^2,\qquad\theta\equiv h_{B'}\pmod{2\pi}.$$
On $X^2+Y^2\le M^2$,
$$DP=\begin{pmatrix}-1&1\\-(1-c)&c\end{pmatrix}+O\big((M+L+1)K(s)\big),\qquad c=2^{-1/(2a)-1},\qquad\det=1-2c\neq0.$$

*Proof.* Eliminate $\theta$ using the $B'$ condition. In the rotating frame $|B'|\approx s+x$, $|D|\approx s+y$ and $|C|\approx\sqrt2s+(x+y)/\sqrt2$. Differentiating $r^{-1/a}$ gives slopes $-1$, $-1$ and $-c$ (times $1/(K(s)s)$) for $h_{B'},h_D,h_C$ in the directions shown. The frame angles, the $\varphi$ terms and the corrections to the radii contribute $O(MK(s))$, $O(LK(s))$ and $O(K(s))$ to the scaled derivatives. ∎

**Theorem 378.3 (count).** For $M\ge1$ with $(M+L+1)K(s)$ small, the number of near-squares with $|A|\le M\,K(s)\,s$ and $\theta\in[0,2\pi)$ is
$$N(M)=\frac{2|1-2c|}{\pi}M^2\,(1+o(1))+O(M).$$

*Proof.* By the lemma, $P$ is injective on the disk and its image has area $|1-2c|\pi M^2(1+o(1))$, with boundary length $O(M)$. Each of the 8 arm triples contributes the lattice points of $2\pi\mathbb Z^2$ (shifted) in that image: area$/(4\pi^2)+O(M)$. Add the 8 contributions. ∎

The leading term does not involve $\varphi$. Asymmetry moves the lattice but does not thin it. In terms of the turn spacing, $|A|\le\eta$ means $M=g\eta/d(s)$, so on average there are about $(2|1-2c|/\pi)(g\eta/d(s))^2$ near-squares within $\eta$ of $p$. This makes Proposition 377.1 quantitative.

**Lemma 378.4 (sector rule).** If $|A|\ll s$ and $A$ is placed at the wild point in the curve order, the near-square is graceful iff $B'$ and $D$ lie on different arms.

*Proof.* Along the curve, arm one runs inward with $\Theta$ increasing, then comes the wild point, then arm two runs outward with $\Theta$ decreasing. The opposite corner $C$ is at radius about $\sqrt2s$, beyond $B'$ and $D$, so it has the smallest $\Theta$ on whichever arm it lies on.

If $B'$ and $D$ are on the same arm, $C$ is adjacent to $A$ in the cyclic order. That is a crossed order. If they are on different arms, the order starting from $A$ is $A,D,C,B'$ or $A,B',C,D$. Both are square orders. ∎

## Symmetrizing the core

**Corollary 378.5.** Let $\gamma$ have a two-arm spiral point at $p$, with no square of side $s$, and let $\gamma'$ turn arm two inside $B(p,r_1)$, $r_1<s/2$, to the half-turn position. Then $L(\gamma',s)+L(\gamma,s)$ is the parity of the number of lattice near-squares with $B'$ and $D$ on different arms whose free corner lies in the region swept by arm two.

That region is a family of thin spiral bands, one per turn. A band at radius $\rho$ has width about $|\varphi(\rho)|K(\rho)\rho$, which is $|\varphi|/g$ times the turn spacing $d(\rho)$. For $\rho<s/2$ that is far below the lattice spacing $K(s)s$. The total swept area is about $|\bar\varphi|r_1^2/2$, so the expected number of near-squares in it is of order $(r_1/K(s)s)^2$, which is large. Whether each lattice point falls in a band is decided at scale $d(\rho)\ll K(s)s$. Nothing in the local geometry forces the count to be even.

**What a counterexample must do.** If $\gamma$ had no inscribed square at all, then $L(\gamma,s)=0$ for every $s$. That is because $L$ can change only at a square's side length, and it vanishes for $s$ above the diameter. Meanwhile, the curve symmetrized out to about $8\rho(s)$ is protected (v374). So every symmetrizing path must cross an odd number of graceful side-$s$ squares, at every scale $s$ at once. The surgery rule splits that number into contributions from local moves. The innermost contribution is the lattice-band count of Corollary 378.5. No contradiction follows, but this is now a concrete arithmetic condition rather than an unstructured one.

## Exact checks

**Surgery rule on drawings.** Each bump replaces a short piece of one edge by two segments through an apex, inside a ball of radius below $s/2$. The parity before and after is computed exactly from the small-square counts and Lemma 370.1. The near-squares are found in closed form: three linear conditions plus $|c|=1$ give at most two solutions per edge triple.

| drawing | bump | trials | swept (graceful) | predicted flip | actual flip |
|---|---|---|---|---|---|
| asymmetric, $s=1/20$ | targeted at one near-square | 15 | 1 (1) | 1 | 1 |
| asymmetric, $s=1/20$ | random | 12 | 0 | 0 | 0 |
| asymmetric, $s=1/20$ and $1/25$ | targeted at a crossed near-square | 7 | 1 (0) | 0 | 0 |
| symmetric and asymmetric, $s=1/20,1/25$ | over two near-squares | 7 | 2 (2) | 0 | 0 |
| asymmetric, $s=1/25$ | over a cluster | 1 | 3 (3) | 1 | 1 |

On the drawings, crossed near-squares near the curve occur only in the asymmetric drawings: 2 of 20 candidates at $s=1/20$, 5 of 35 at $1/25$, none of 16 in the symmetric drawing.

**Lattice on the continuous arcs.** Lattice points of $P$ are enumerated from the affine guess and refined by Newton's method, for all 8 arm triples.

| $M$ | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|
| predicted $M^2/\pi$ | 1.3 | 5.1 | 20.4 | 81.5 | 325.9 |
| found ($\varphi=0$, $s=0.02$) | 2 | 6 | 22 | 82 | 322 |
| found ($\varphi=0$, $s=0.01$) | 2 | 4 | 20 | 82 | 326 |
| found (asymmetric, $s=0.02$) | 2 | 4 | 20 | 82 | 326 |
| found (asymmetric, $s=0.01$) | 2 | 4 | 19 | 83 | 326 |

At $M=32$ the graceful share is between 0.497 and 0.506, and the sector rule of Lemma 378.4 holds for every near-square. The nearest near-squares, at $0.38\,d(s)$ ($\varphi=0$) and at $0.317$ and $0.569\,d(s)$ (asymmetric) for $s=0.01$, are the ones v377 found by a different method.

## Where this leaves asymmetric points

The surgery rule turns the path condition into bookkeeping: count the graceful near-squares swept by each local move. At the core of a spiral point the near-squares form an explicit lattice, and the moves sweep explicit thin bands. The open question is now arithmetic: can the parities of these lattice-in-band counts, summed over all the moves of a symmetrizing path, be shown to be even at some scale? Or, the other way round, does the requirement that they be odd at every scale constrain the curve enough to rule it out?

[Code](../code/surgery_v378.py) · [Tests](../tests/test_surgery_v378.py) · [Data](../results/v378/surgery_and_lattice.json)
