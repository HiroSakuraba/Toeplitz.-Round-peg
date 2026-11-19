# v253: Boundary-corrected square parity and an exact GLOBAL maximum-square bottleneck

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact-rational theorem and regression-checked computational certificate for the genuinely moving convex nine-gon from v250. Builds on the **complete** large-square orbit census of v252 and establishes both (i) a boundary-flux-corrected **integer** count and (ii) the **global** maximum-square side as an algebraic size envelope over the *entire* interval. The result is explicitly restricted to one certified family and cutoff, and **does not** solve universal Big Peg or Toeplitz.

## The actual mathematical advance

v251 found two different square configurations that exchange relative size. At that time it was not legitimate to conclude that the **globally largest square** changes identity: an unknown third square might always be larger.

v252 removed that problem for the interval

\[
W=\left[\frac12,\frac{501}{1000}\right]
\]

and the strict squared-side cutoff

\[
\eta=\frac{94257}{10000}=9.4257.
\]

It exhaustively classified **all 126 possible** large-square edge configurations: 87 impossible by global pair-edge geometry, 36 impossible by exact segment-membership signs, one below the cutoff, and exactly two surviving square branches. The pre-existing v223 complete geometric-orbit census independently verifies the endpoint counts, including a square touching a polygon vertex, while v251 proves at least one of the two squares stays strictly above \(\eta\) at every real time.

The two surviving configurations have edge assignments

\[
A=(0,3,5,7),\qquad B=(1,3,5,8).
\]

Write their squared sides as exact rational functions \(S_A(\lambda)\) and \(S_B(\lambda)\), obtained from their four-by-four Cramer matrices.

**The global maximum is therefore EXACTLY**

\[
\boxed{s_{\max}(P(\lambda))^2=
\max\big(S_A(\lambda),S_B(\lambda)\big)
\qquad\forall\lambda\in W.}
\]

This is a theorem about the true **global maximum over every inscribed square**, not merely the maximum over two selected witnesses. The exhaustive v252 exclusion and v251 full-time existence theorem are essential premises.

## 1. Unique algebraic minimizer of the global square maximum

The side functions are rational and continuous on the entire window because their Cramer determinants are strictly nonzero. Write a given side function as

\[
S(\lambda)=\frac{N_u^2(\lambda)+N_v^2(\lambda)}
{D(\lambda)^2}=\frac{H(\lambda)}{D(\lambda)^2}.
\]

Its derivative has the exact sign of the **polynomial**

\[
D(\lambda)\big(H'(\lambda)D(\lambda)-2H(\lambda)D'(\lambda)\big),
\]

because the denominator in the resulting quotient is \(D^4>0\).

The v253 verifier uses exact rational Bernstein coefficient positivity on the *complete closed interval* \(W\) to prove

\[
\boxed{S_A'(\lambda)<0,\qquad S_B'(\lambda)>0
\quad\forall\lambda\in W.}
\]

The independently verified exact side-comparison polynomial has one and only one zero in \(W\), with \(S_A>S_B\) before it and \(S_B>S_A\) after it. Call this algebraic real time \(\lambda_*\).

It follows that

\[
M(\lambda):=\max(S_A,S_B)=
\begin{cases}
S_A(\lambda),&\lambda\le\lambda_*,\\
S_B(\lambda),&\lambda\ge\lambda_*,
\end{cases}
\]

which is **strictly decreasing before** and **strictly increasing after** \(\lambda_*\).

**Thus \(\lambda_*\) is the unique GLOBAL MINIMIZER of the maximum attainable squared side over the entire homotopy window.**

Its exact defining equation is the polynomial identity

\[
(N_{A,u}^2+N_{A,v}^2)D_B^2
-
(N_{B,u}^2+N_{B,v}^2)D_A^2=0,
\]

taken at the unique root of that polynomial in \(W\). No numerical root approximation is substituted into a proof.

### Very tight certified enclosures

The root is rigorously isolated between the **exact rationals**

\[
\boxed{\frac{50054929969}{10^{11}}<
\lambda_*<
\frac{50054929970}{10^{11}}}.
\]

Let \(L,U\) denote these two fractions. By strict monotonicity and the equality \(S_A(\lambda_*)=S_B(\lambda_*)\),

\[
\max(S_A(U),S_B(L))
<
M(\lambda_*)
<
\min(S_A(L),S_B(U)).
\]

Every term outside the unknown root is a **rational number**. The exact verifier checks the stronger rational enclosure

\[
\boxed{
9.4257027585<
\min_{\lambda\in W}s_{\max}(P(\lambda))^2
<
9.4257027586.
}
\]

This is a certified interval of width \(10^{-10}\) for the **true global bottleneck**, and not just a bound for a convenient square.

Since the whole moving polygon lies in a disk of radius three, its globally largest interior disk has \(r_{\rm in}\le3\). Therefore the entire window also satisfies the new restricted-family Big Peg bound

\[
\boxed{J(P(\lambda))>
\frac{\sqrt{9.4257027585}}{3}
>1.02337.}
\]

Neither result asserts an optimal constant for general polygons.

## 2. The corrected square-count theorem

Now return to the **strict size filter** \(S>\eta\). The complete v252 census and the nonvanishing, transverse cutoff derivatives establish exactly two signed events:

- \(\lambda_B\): square branch \(B\) crosses **up** through the cutoff, contributing \(+1\) to the raw count;
- \(\lambda_A\): square branch \(A\) crosses **down** through the cutoff, contributing \(-1\).

Their exact roots are known to be unique and are isolated in disjoint rational brackets

\[
\lambda_B\in\Big(\frac{50053}{100000},\frac{50054}{100000}\Big),\qquad
\lambda_A\in\Big(\frac{500551}{10^6},\frac{50056}{100000}\Big),
\]

with \(\lambda_B<\lambda_*<\lambda_A\).

Define \(N_>(t)\) to be the number of **geometric** inscribed square orbits whose squared sides are **strictly greater** than \(\eta\).

Set the cumulative signed size-boundary flux

\[
\Phi(t)=
\underbrace{{\bf1}_{\,t>\lambda_B}}_{\text{entry counted only after equality}}
-
\underbrace{{\bf1}_{\,t\ge\lambda_A}}_{\text{exit counted already at equality}}.
\]

The different endpoint conventions are not arbitrary: for a strict \(>\) filter, a square is still **excluded** at its upward-crossing equality and becomes excluded **at** its downward-crossing equality.

**The boundary-corrected integer index** is

\[
\boxed{\mathcal I(t)=N_>(t)-\Phi(t)=1
\quad\text{for every real }t\in W.}
\]

Thus its reduction modulo two is also constantly one:

\[
\boxed{(N_>(t)-\Phi(t))\bmod2=1.}
\]

This equality is demonstrated using all possible above-cutoff square orbits, not by assuming an individual square remains the selected witness. The positive square's size at the crossover \(\lambda_*\) is well above the cutoff, so this size-equality event **does not** contribute any flux.

### Exact event ledger

| Parameter region, including equality conventions | Above-cutoff orbits \(N_>\) | Signed flux \(\Phi\) | Corrected integer index \(\mathcal I\) |
|---|---:|---:|---:|
| \(1/2\le t\le\lambda_B\) | 1 | 0 | 1 |
| \(\lambda_B<t<\lambda_A\) | 2 | +1 | 1 |
| \(\lambda_A\le t\le501/1000\) | 1 | 0 | 1 |

For every **rational** query time, the Python program decides these counts using exact signs of the two cutoff-polynomials. **It never rounds the two algebraic event times to rationals.** At the algebraic equality times, the strict cutoff convention and unique monotonicity establish the corresponding count symbolically.

The integer conservation formula is the standard finite-event accounting statement: raw counts change by precisely the signed number of transversely crossing objects. The nontrivial work in this instance is establishing a **complete orbit census**, proving exactly two event roots, and handling the endpoint label duplication. This is **not** a novel unfiltered topological invariant and **not** a degree that survives arbitrary singular orbit births or deaths.

## 3. Independent validation and adversarial controls

The exact v253 verifier *replays* the v252 complete 126-chart census, the independent v223 exact geometric-orbit endpoint checks, the v251 full-time high-size chart atlas, and the v250 all-time Jordan family verification before declaring the new global maximum theorem.

The 18 new regression tests include:
- Every-time squared-side monotonicity certified by derivative-sign Bernstein coefficients.
- Both sides of a \(10^{-11}\)-wide exact isolating bracket for the algebraic bottleneck time.
- A \(10^{-10}\)-wide exact enclosure for the global minimum of the **true maximum** square side squared.
- Boundary event order, plus integer and mod-two corrected count tests at representative rational times and dozens more rational times inside the algebraic event brackets.
- Independent v223 **unfiltered** exact geometric-square census at rational times, checking the global maximum identity from a separate algorithm.
- Strict-cutoff equality conventions for entry/exit, and rejection of unsupported polygon families, times, floating inputs, and arbitrary cutoffs.

Reproduce with Python 3.12+ and no nonstandard library:

~~~sh
python code/big_peg_relative_parity_v253.py
python -m unittest -v tests.test_big_peg_relative_parity_v253
~~~

The dedicated CI also rechecks the full v252 cutoff census and v251 high-size witness existence theorem.

## 4. Limitations and the genuinely unresolved problem

A corrected integer count **cannot disappear mysteriously** when every event is classified: this is finite 1-parameter bookkeeping. It does not establish that an above-cutoff square exists for an *arbitrary* Jordan curve or that the corrected index has useful absolute meaning in every family.

The obstruction remains precisely what v252 exposed: even the **parity** of above-threshold squares changes at the size boundary. To obtain a positive universal Big Peg bound, one needs to show that all large-square branches **cannot simultaneously escape** through size thresholds arbitrarily close to zero, uniformly in polygon complexity. This is **not proved** here.

The result does not handle:
- A complete orbit census for arbitrary affine polygon deformations, including rank-deficient incidence matrices.
- Birth/death of unfiltered square orbits; these require their own signed/singular event terms.
- Boundary curves with infinitely many arbitrarily small corners or fractal geometry.
- A universal positive Big Peg constant or the unrestricted Toeplitz conjecture.

## References

- Benjamin Matschke, [*A Survey on the Square Peg Problem*](https://www.ams.org/notices/201404/rnoti-p346.pdf), *Notices of the AMS* 61 (2014). Existing parity and configuration-space approaches provide context; the present corrected count is **not claimed as new topological parity theory**.
- [v251 exact size-aware two-square atlas](v251-size-aware-multibranch-atlas.md).
- [v252 complete size-filtered 126-chart ledger](v252-complete-large-square-event-ledger.md).

**v254 priority:** Investigate genuine **unfiltered square-orbit birth/death** transitions in a controlled moving-polygon family, including singular four-edge incidence charts. Any proposed corrected invariant must explicitly record *both* interior orbit events and the cutoff-boundary flux; the present v253 program has only verified the latter.
