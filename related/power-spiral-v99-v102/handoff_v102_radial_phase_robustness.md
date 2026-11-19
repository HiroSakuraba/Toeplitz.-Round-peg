# Project Handoff v102: A Quantitative Robustness Theorem for Power-Like Radial Phase Profiles

**Author:** Benjamin John Schulz

## Executive result

v100/v101 proved, with independently reproduced face certificates, that every exact power spiral

\[
\gamma_p(r)=r e^{i/r^p},\qquad p>0,
\]

contains a non-graceful square of radial temporal order \(0132\).

v102 extracts the first theorem that is **not restricted to an exact power law**.

The decisive object is the *normalized distribution of angular turning across logarithmic radius*.

Let

\[
\Gamma_\Theta(r)=r e^{i\Theta(r)},\qquad 0<r\le1,
\]

with \(\Theta\) decreasing.  On a radial window \([r_0,\rho r_0]\), define

\[
\boxed{
P_{\Theta;r_0,\rho}(s)
=
\frac{\Theta(r_0)-\Theta(r_0\rho^s)}
{\Theta(r_0)-\Theta(\rho r_0)},
\qquad 0\le s\le1.
}
\tag{0.1}
\]

For an exact power law \(\Theta(r)=r^{-p}\),

\[
\boxed{
P_p(s)=G_a(s)
=
\frac{1-e^{-as}}{1-e^{-a}},
\qquad
 a=p\log\rho.
}
\tag{0.2}
\]

The v101 face margins imply:

> **Theorem A — quantitative power-profile robustness.**  Fix \(p>0\).  Suppose a strictly radial phase law \(\Theta\) realizes continuously all windows in the appropriate v99/v101 Miranda rectangle for exponent \(p\), with the same endpoint radius ratio and the same total angular drop.  If on every such window
> \[
> \boxed{
> \sup_{0\le s\le1}
> \left|P_{\Theta;r_0,\rho}(s)-G_{p\log\rho}(s)\right|
> \le10^{-6},
> }
> \tag{0.3}
> \]
> then \(\Gamma_\Theta\) contains an inscribed Euclidean square.  The forced square is in radial temporal order \(0132\).

Thus the v100 result has a genuine open neighborhood in normalized phase-profile space.

The constant \(10^{-6}\) is deliberately conservative.  It is obtained only from the *stopping margins* in v101, not from the much larger apparent true face margins.

---

# 1. Radial windows and forced square geometry

Take two proposed opposite square vertices on a strictly radial spiral,

\[
z_0=\Gamma_\Theta(r_0),
\qquad
z_2=\Gamma_\Theta(r_3),
\qquad
r_3=\rho r_0,
\quad \rho>1.
\]

Rotate and scale so that

\[
z_0=1.
\]

The v100 boxes prescribe a lifted endpoint angular drop

\[
\boxed{
\Delta=3\pi-\varepsilon.
}
\tag{1.1}
\]

Thus

\[
z_2=-\rho e^{i\varepsilon}.
\]

The other two Euclidean square vertices are forced by the diagonal transform,

\[
Q_+=qz_0+pz_2,
\qquad
Q_-=pz_0+qz_2,
\]

where

\[
p=\frac{1+i}{2},
\qquad
q=\frac{1-i}{2}.
\]

Their radii determine normalized logarithmic positions

\[
\boxed{
s_\pm=\frac{\log |Q_\pm|}{\log\rho}\in(0,1).
}
\tag{1.2}
\]

All of this is pure Euclidean geometry.  It does not depend on the radial phase law.

The fixed target arguments likewise determine lifted angular drops \(\delta_\pm\).  Therefore for an arbitrary radial phase profile the two square residuals are

\[
\boxed{
\widetilde E_\pm
=
\delta_\pm-\Delta P_{\Theta;r_0,\rho}(s_\pm).
}
\tag{1.3}
\]

For the power reference,

\[
E_\pm
=
\delta_\pm-\Delta G_a(s_\pm).
\]

Consequently

\[
\boxed{
|\widetilde E_\pm-E_\pm|
\le \Delta\eta
<3\pi\eta,
}
\tag{1.4}
\]

whenever

\[
\|P_\Theta-G_a\|_\infty\le\eta.
\]

Equation (1.4) is the entire robustness mechanism.

---

# 2. What it means to realize the same Miranda window

A general radial phase law does not automatically have the exact total angular drop required by a point in the v100 boxes.

For a given radius ratio \(\rho\), define

\[
D_\rho(r)
=
\Theta(r)-\Theta(\rho r).
\]

A parameter point with target \(\Delta\) is **realized** if there is

\[
0<r_0\le\frac1\rho
\]

such that

\[
\boxed{D_\rho(r_0)=\Delta.}
\tag{2.1}
\]

Theorem A assumes these solutions can be chosen continuously across the relevant parameter rectangle.

This is not an artificial condition.  There is a simple sufficient criterion.

## Lemma 2.1 — turn-completeness criterion

Assume \(\Theta\in C^1(0,1]\) and put

\[
\boxed{V(r)=-r\Theta'(r).}
\tag{2.2}
\]

Suppose:

1. \(V(r)>0\);
2. \(V\) is strictly decreasing as \(r\) increases;
3. \(V(r)\to\infty\) as \(r\downarrow0\);
4. at the outermost admissible radius,
   \[
   D_\rho(1/\rho)<\Delta
   \]
   for every target \((\rho,\Delta)\) in the chosen box.

Then (2.1) has a unique solution, varying continuously with \((\rho,\Delta)\).

### Proof

Differentiate with respect to \(\log r\):

\[
\frac{d}{d\log r}D_\rho(r)
=
-r\Theta'(r)+\rho r\Theta'(\rho r)
=
-V(r)+V(\rho r)<0.
\]

Thus \(D_\rho\) is strictly decreasing with \(r\).

Also

\[
D_\rho(r)
=
\int_r^{\rho r}V(u)\,\frac{du}{u}
\ge
(\log\rho)V(\rho r)
\to\infty
\]

as \(r\downarrow0\).  The endpoint inequality and the intermediate value theorem give a unique root.  Strict monotonicity gives continuous dependence. \(\square\)

Exact power laws satisfy this criterion.

---

# 3. The three reference regimes

The proof uses the same parameter regions already established in v99/v101.

## Regime A: \(p\ge2\)

Put

\[
q_p=\frac1{1+p},
\]

and use

\[
1\le c\le\frac{19}{10},
\qquad
\frac9{20}\le d\le\frac{18}{25},
\]

with

\[
\rho=e^{cq_p},
\qquad
\varepsilon=dq_p.
\]

The Poincare--Miranda coordinates are

\[
(E_-,E_+).
\]

## Regime B: \(1/200\le p\le2\)

Use

\[
\frac34\le c\le\frac32,
\qquad
\frac25\le d\le\frac{18}{25},
\]

again with

\[
\rho=e^{cq_p},
\qquad
\varepsilon=dq_p.
\]

The coordinates are

\[
\boxed{
E_-,
\qquad
G=E_++\frac{57}{400}E_-.
}
\tag{3.1}
\]

## Regime C: \(0<p\le1/200\)

Use the v99 rectangle

\[
\frac{119}{50}\le\rho\le\frac{62}{25},
\qquad
\frac{27}{50}\le\varepsilon\le\frac{61}{100}.
\]

Its coordinates are

\[
\boxed{
H_1=E_+-E_-,
\qquad
H_2=\frac{E_++E_-}{2}.
}
\tag{3.2}
\]

The three regimes overlap at their endpoints and cover every \(p>0\).

---

# 4. Conservative certified margins

v101's independent Arb calculation certifies all eight v100 faces on the exact rational boxes.

For the robustness calculation we deliberately round those margins *downward* to much simpler rational numbers.

## Regime A

\[
\begin{array}{c|c|c}
\text{face}&\text{v101 Arb margin}&\text{margin used in v102}\\ \hline
c=1 & 9.451231\times10^{-2} & 9\times10^{-2}\\
c=19/10 & 5.255473\times10^{-3} & 5\times10^{-3}\\
d=9/20 & 6.767832\times10^{-5} & 6\times10^{-5}\\
d=18/25 & 3.620809\times10^{-5} & 3\times10^{-5}
\end{array}
\]

## Regime B

\[
\begin{array}{c|c|c}
\text{face}&\text{v101 Arb margin}&\text{margin used in v102}\\ \hline
c=3/4 & 1.277418\times10^{-2} & 1.2\times10^{-2}\\
c=3/2 & 3.930228\times10^{-3} & 3\times10^{-3}\\
d=2/5 & 6.283801\times10^{-5} & 6\times10^{-5}\\
d=18/25 & 2.183953\times10^{-5} & 2\times10^{-5}
\end{array}
\]

## Regime C

The exact v99 rational proof, after allowing the full drift from \(p=0\) to \(p=1/200\), retained margins greater than

\[
0.0285,\quad0.0195,\quad0.000564,\quad0.001564.
\]

v102 rounds these down again to

\[
0.028,\quad0.019,\quad0.0005,\quad0.0015.
\]

---

# 5. Exact robustness budget

Take

\[
\boxed{\eta=10^{-6}.}
\tag{5.1}
\]

Use the rigorous upper bound

\[
\pi<\frac{355}{113}.
\]

Since every box has

\[
\Delta<3\pi,
\]

(1.4) gives

\[
\boxed{
|\delta E_\pm|
<
3\frac{355}{113}\,10^{-6}
=
9.424778761\ldots\times10^{-6}.
}
\tag{5.2}
\]

### Regime A

Each coordinate is one \(E\)-residual, so the maximum possible displacement is (5.2).

The smallest conservative face margin is

\[
3\times10^{-5},
\]

more than three times the perturbation budget.

### Regime B

On the two \(c\)-faces, the coordinate remains \(E_-\), so (5.2) applies.

For

\[
G=E_++\frac{57}{400}E_-,
\]

\[
|\delta G|
\le
3\pi\left(1+\frac{57}{400}\right)\eta.
\]

Hence

\[
\boxed{
|\delta G|
<1.076780974\times10^{-5}.
}
\tag{5.3}
\]

The smallest conservative \(G\)-face margin is

\[
2\times10^{-5}.
\]

So both signs survive.

### Regime C

Because

\[
H_1=E_+-E_-,
\]

\[
|\delta H_1|<1.884955753\times10^{-5}.
\]

Because

\[
H_2=(E_++E_-)/2,
\]

\[
|\delta H_2|<9.424778762\times10^{-6}.
\]

The corresponding conservative face margins are at least

\[
0.019
\]

and

\[
0.0005.
\]

Again the signs survive with large slack.

Therefore every Poincare--Miranda degree used in v99/v101 is unchanged by (0.3).  This proves Theorem A.

---

# 6. Angular-turning measure formulation

The normalized phase profile has a useful differential interpretation.

Let

\[
V(r)=-r\Theta'(r)>0.
\]

On a realized window, put

\[
r(s)=r_0\rho^s.
\]

Then

\[
\Theta(r_0)-\Theta(r_0\rho^s)
=
(\log\rho)
\int_0^s V(r_0\rho^u)\,du.
\]

Thus

\[
\boxed{
P_\Theta'(s)
=
\frac{V(r_0\rho^s)}
{\int_0^1 V(r_0\rho^u)\,du}.
}
\tag{6.1}
\]

For a power law,

\[
V_p(r)=p r^{-p},
\]

and therefore

\[
\boxed{
G_a'(s)
=
\frac{a e^{-as}}{1-e^{-a}}.
}
\tag{6.2}
\]

Equations (6.1)--(6.2) identify the structural quantity:

> the probability distribution of angular turning over logarithmic radius.

This gives an immediate corollary.

> **Corollary B — \(L^1\) turn-density robustness.**  Under the window-realization hypotheses of Theorem A, suppose on every Miranda window
> \[
> \boxed{
> \int_0^1
> \left|
> \frac{V(r_0\rho^s)}{\int_0^1V(r_0\rho^u)\,du}
> -
> \frac{a e^{-as}}{1-e^{-a}}
> \right|ds
> \le10^{-6}.
> }
> \tag{6.3}
> \]
> Then \(\Gamma_\Theta\) contains a non-graceful \(0132\) square.

Indeed, the difference of the two cumulative distribution functions is bounded by the \(L^1\) distance in (6.3).

This formulation is more useful than pointwise closeness of \(\Theta\) itself: additive phase constants disappear, total winding speed is normalized away, and only the *shape* of angular-turn distribution matters.

---

# 7. A simple multiplicative speed cone

Suppose on each realized window there is some positive constant \(K\) such that

\[
\boxed{
V(r_0\rho^s)
=K e^{-as}(1+e(s)),
\qquad
|e(s)|\le\epsilon.
}
\tag{7.1}
\]

After normalization, the cumulative phase profile differs from \(G_a\) by at most

\[
\boxed{
\frac{2\epsilon}{1-\epsilon}.
}
\tag{7.2}
\]

Therefore the explicit choice

\[
\boxed{\epsilon=4\times10^{-7}}
\]

implies

\[
\frac{2\epsilon}{1-\epsilon}
<8.000004\times10^{-7}
<10^{-6}.
\]

So Theorem A applies.

This relative cone is extremely narrow.  It should be read as a first rigorous openness statement, not as the expected optimal robustness radius.

---

# 8. A quantitative obstruction for a square-free radial spiral

Take the contrapositive.

A strictly radial square-free spiral satisfying the turn-completeness assumptions **cannot** remain inside the v102 power-profile tube throughout the relevant Miranda family.

For every reference exponent \(p>0\), at least one required boundary window must satisfy

\[
\boxed{
\sup_s
|P_\Theta(s)-G_{p\log\rho}(s)|
>10^{-6}.
}
\tag{8.1}
\]

Equivalently, under differentiability, the logarithmic angular-turn distribution cannot stay within \(10^{-6}\) in cumulative-distribution distance of the corresponding truncated exponential on all required windows.

This is a genuine structural restriction on any radial-monotone local model of a hypothetical Toeplitz counterexample.

It says that the exact exponent \(p\) is not the essential object.  What the square degree sees is the distribution of angular turning across the radial scale.

---

# 9. Why the numerical tolerance is so small

The \(10^{-6}\) constant is **certificate-limited**, not evidently geometry-limited.

The Arb verifier stops subdividing a cell as soon as its enclosure lies on the correct side of zero.  Therefore its reported smallest leaf margin can be arbitrarily smaller than the true minimum of the face function.

Direct high-precision sampling of the exact power faces gives apparent true minima of order

\[
10^{-1},
\]

not \(10^{-5}\).

For example, numerical scans put the eight signed face minima approximately at

\[
0.451,
0.227,
0.162,
0.116,
0.113,
0.337,
0.100,
0.097.
\]

These numbers are **not certificates** and are not used in Theorem A.

They suggest that a margin-targeted Arb run could plausibly enlarge the rigorous profile radius from \(10^{-6}\) to roughly \(10^{-3}\) or even a few \(10^{-3}\).

That is the highest-value computational follow-up.

---

# 10. Relation to the Haar-area program

v94 identified the dyadic midpoint-triangle tree simultaneously as:

1. the generator of winding data;
2. the increment tree for the symplectic-area primitive.

v102 identifies a different local statistic relevant to the non-graceful spiral branch:

\[
\boxed{
\text{normalized angular-turn distribution over logarithmic radius}.
}
\]

The next bridge should therefore not try to prove merely

\[
\text{area divergence}\Rightarrow\text{power law}.
\]

That is false in spirit and unnecessarily rigid.

The more plausible target is a compactness/dichotomy statement of the form

\[
\boxed{
\begin{gathered}
\text{persistent localized Haar-area defect}\\
\Downarrow\\
\text{after choosing a radial center and logarithmic scale, either}\\
\text{(i) the turn distribution enters a v102 power-like tube, or}\\
\text{(ii) the arc becomes graph-like/oscillatory in another controlled way.}
\end{gathered}
}
\]

Branch (i) is now excluded from any square-free counterexample by Theorem A.

The task is to understand branch (ii).

---

# 11. Claim boundary

## Established in v102

- a quantitative stability lemma for the radial two-phase residuals;
- an explicit all-exponent profile tolerance
  \[
  \boxed{\eta=10^{-6}};
  \]
- persistence of the v99/v101 Poincare--Miranda square degree for any continuously realized radial window family inside that profile tube;
- the equivalent \(L^1\) logarithmic angular-speed formulation;
- a simple multiplicative speed cone sufficient for the theorem;
- a quantitative obstruction that any square-free radial spiral must leave the power-profile tube on at least one required boundary window.

## Inputs inherited from earlier versions

- v99's exact small-exponent face margins;
- v101's independent Arb face certificates and their certified margins;
- the exact radial two-target geometry from v96.

## Not established

- that an arbitrary radial spiral is turn-complete in the sense of Section 2;
- that an arbitrary radial spiral lies in the \(10^{-6}\) power-profile tube;
- a useful large robustness radius such as \(10^{-2}\);
- a theorem deriving v102 turn-distribution control from Haar-area primitive failure;
- the universal Toeplitz square-peg conjecture.

---

# 12. Recommended next work

1. **Margin-targeted Arb certification.**  Instead of stopping when a leaf has the correct sign, require a prescribed face margin, initially \(0.05\).  If all eight faces certify at that level, the profile robustness radius jumps to approximately
   \[
   \frac{0.05}{3\pi(1+57/400)}\approx4.6\times10^{-3}.
   \]
2. **Perturbation experiments.**  Test radial laws whose logarithmic angular speed is an exponential times smooth oscillatory factors, and measure the actual failure threshold of the \(0132\) degree.
3. **Haar-to-turn compactness.**  On nested intervals carrying persistent signed Haar area, search for centers for which the induced turning measures have tight weak limits.  Classify those limits rather than assuming a power law.
4. **Keep Arb primary.**  Use the v101 Decimal layer as an independent regression checker, not as the publication-critical transcendental implementation.
