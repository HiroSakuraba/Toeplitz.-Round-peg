# v373. The path sweep: a square-free curve forces squares of every small size along any path to it

[v370](v370-half-turn-wild-points-cannot-absorb-parity.md) and [v372](v372-mirror-symmetric-wild-points.md) gave two classes of *protected* curves: curves whose large graceful square count is odd at every small cutoff. [v371](v371-stability-of-the-half-turn-mechanism.md) showed that this parity can change only when some curve on a path has a square of side exactly the cutoff. This note turns that into a statement about any hypothetical counterexample, and checks it on an exact path.

**Results.**
1. **Sweep theorem (Theorem 373.1).** Deform a protected curve continuously into a curve with no inscribed square. Then every size below the protection scale occurs as the side of an inscribed square on some intermediate curve. When the crossings of a size are finitely many and transverse, their number is odd. The intermediate curves need not be simple.
2. **Collapse must pass through every size.** A square-free curve can be reached only through curves carrying squares of every small size. The parity can never simply vanish into tiny squares on the way down.
3. **Exact check on a twist path.** Along the twist path of the v370 double spiral, the large graceful parity at cutoff $1/20$ changes from odd to even. Exact bisection locates the change in a twist interval of width about $3\cdot10^{-6}$. Across it, one reversed square's side goes from $0.05000008$ to $0.04999994$, while its corners move by about $1.6\cdot10^{-6}$.

## Protected curves

For a registered closed curve $C$ and $s>0$ with no square of side exactly $s$, let $L(C,s)$ be the mod-two degree of the large aligned plus large reversed squares, those of side above $s$, as in v371. It is defined whether or not $C$ is simple, because zeros on the boundary of the row space are collapsed and never large.

Call $C$ **protected on $(0,s_0)$** if $L(C,s)=1$ for every $s\in(0,s_0)$ at which $C$ has no square of side exactly $s$.

Corollary 371.3 and the proof of Theorem 372.3 give protected curves:
- every curve of Theorem 370.2 (half-turn wild point), with $s_0$ its explicit scale;
- every curve of Theorem 372.2 (mirror-cone wild point), since there the large aligned degree is odd and the large reversed degree is even;
- every curve of Theorem 372.3 (finitely many such points);
- every Jordan polygon, by Lemma 370.1, with $s_0$ below its small-square scale $\tau/\sqrt2$ (Lemma 369.6).

## The sweep theorem

**Theorem 373.1.** Let $\gamma_\lambda$, $\lambda\in[0,1]$, be a uniformly continuous family of registered closed curves, with $\gamma_0$ protected on $(0,s_0)$. Let
$$
\Sigma=\{s>0:\ \text{some }\gamma_\lambda\text{ has a graceful inscribed square of side exactly }s\}.
$$
Then for every $s\in(0,s_0)\setminus\Sigma$, the curve $\gamma_1$ has a graceful inscribed square of side greater than $s$.

*Proof.* If $s\notin\Sigma$, no curve on the path has a square of side exactly $s$. In particular $L(\gamma_0,s)=1$. By Corollary 371.2, $L(\gamma_1,s)=1$, so $\gamma_1$ has a zero in the large region, which is a square of side greater than $s$ with distinct corners. ∎

**Corollary 373.2 (sweep).** If a Jordan curve $\gamma_1$ has no inscribed square, then $(0,s_0)\subseteq\Sigma$ for every path to $\gamma_1$ from a curve protected on $(0,s_0)$.

Moreover, fix $s\in(0,s_0)$ such that the curves on the path with a side-$s$ square occur at finitely many parameters, each with a single square crossing the cutoff transversally. Then the number of such crossings is odd.

*Proof.* The first part is Theorem 373.1. For the second, $L$ starts at 1 and must end at 0, and it changes exactly at transverse crossings of the cutoff wall. ∎

So a square-free curve cannot be reached from a polygon, or from a protected wild curve, by letting the squares quietly collapse. Every size below $s_0$ must be crossed, an odd number of times when crossings are clean. The large region never meets the collapsed configurations. So these crossings are genuine squares of positive size $s$, not tiny squares disappearing into a point.

This is the classical deformation idea, made precise for the large region and for wild starting curves. It does not by itself exclude a counterexample. It says what a counterexample must look like from every protected starting point.

## An exact twist path

The path turns arm two of the v370 double spiral (radii $\theta^{-1/2}$, depth three half-turns, 24 vertices per turn) by $\pi+\phi$, with rational parameter $u=\tan(\phi/2)$. Every drawing on it is an exact Jordan polygon. By Lemma 370.1 the large graceful parity is $1+(\text{small aligned})+(\text{small reversed})$.

| $u$ | small graceful squares (side $<1/20$) | $L(\cdot,1/20)$ |
|---|---|---|
| $0$ | 10 | 1 |
| $1/200$ | 6 | 1 |
| $1/20$ | 13 | 0 |

Between $u=1/200$ and $u=1/20$ the parity changes, so by Theorem 373.1 some drawing in between has a square of side exactly $1/20$. Fourteen exact bisection steps narrow the change to
$$
u\in\left[\tfrac{38029}{3276800},\ \tfrac{19019}{1638400}\right],
$$
an interval of width about $3\cdot10^{-6}$.

At the upper end, the largest small square has side $0.04999994$. At the lower end, the same square has corners within $1.6\cdot10^{-6}$ and side $0.05000008$, just above the cutoff. It lies in the reversed sector at both ends. So the parity change is exactly one square crossing the cutoff, as the theorem requires. It sits on the outer turns, as v371's data showed, not at the wild point.

## Where this points

Every result in this strand now reduces to one question: along some path from a protected curve, can the side-$s$ crossings be controlled? The next note leaves symmetry behind. It studies wild points whose arms approach the point steadily but wind or swing without any symmetry, and asks what replaces the pairing there.

[Code](../code/path_sweep_v373.py) · [Tests](../tests/test_path_sweep_v373.py)
