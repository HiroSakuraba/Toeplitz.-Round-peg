# Project Handoff v100: Every Power Spiral Contains a Non-Graceful Square

**Author:** Benjamin John Schulz

## Executive result

For every real exponent

\[
\boxed{p>0},
\]

the embedded strictly radial spiral

\[
\boxed{\gamma_p(t)=t\,e^{i/t^p},\qquad 0<t\le1}
\]

contains an inscribed Euclidean square.

Because \(|\gamma_p(t)|=t\) is strictly increasing, no gracefully ordered square can occur.  The square supplied here has temporal order

\[
\boxed{0132}
\]

(up to the reflected chirality).

This is a **computer-assisted theorem for the entire one-parameter power-spiral family**.  It is not a proof of the Toeplitz conjecture for arbitrary Jordan curves.

The proof combines:

1. the exact small-exponent theorem from v99;
2. one Poincare--Miranda box for \(1/200\le p\le2\);
3. one Poincare--Miranda box for \(p\ge2\).

No root-finding output is used in the theorem.

---

# 1. Radial two-phase reduction

Write

\[
q=\frac1{1+p},
\qquad
p=\frac{1-q}{q}.
\]

Normalize the nearest proposed square vertex to

\[
z_0=1.
\]

Use

\[
\rho=e^{cq},
\qquad
\varepsilon=dq,
\qquad
\Delta=3\pi-\varepsilon,
\]

and

\[
z_2=-\rho e^{i\varepsilon}
\]

for the farthest opposite vertex.

The other two Euclidean square vertices are forced by the diagonal transform.  Denote them by \(Q_+\) and \(Q_-\), ordered so that

\[
R_+=|Q_+|>|Q_-|=R_-.
\]

Their squared radii are

\[
\boxed{
R_\pm^2
=
\rho\bigl(\cosh(cq)\pm\sin(dq)\bigr).
}
\tag{1.1}
\]

For each forced radius, strict radiality gives exactly one point of the spiral.  The square condition is therefore only two phase-matching equations.

---

# 2. Analytic compactification

Let

\[
x=cq,
\qquad
y=dq.
\]

Define the entire functions

\[
S(y)=\frac{\sin y}{y},
\qquad
C(x)=\frac{\cosh x-1}{x^2},
\]

with their removable values at zero, and

\[
L(z)=\frac{\log(1+z)}{z}.
\]

Set

\[
\boxed{
w_\pm=c^2q\,C(x)\pm d\,S(y).
}
\]

Then

\[
\cosh x\pm\sin y=1+q w_\pm,
\]

so the normalized radial logarithms are

\[
\boxed{
s_\pm
=\frac12\left(1+\frac{w_\pm L(qw_\pm)}{c}\right).
}
\tag{2.1}
\]

This formula is analytic at \(q=0\).

Put

\[
a=c(1-q)=p\log\rho.
\]

Define

\[
E(u)=\frac{1-e^{-u}}u,
\qquad
\boxed{
g_a(s)=s\frac{E(as)}{E(a)}}.
\tag{2.2}
\]

Then

\[
\boxed{
g_a(s_\pm)
=\frac{1-R_\pm^{-p}}{1-\rho^{-p}}.
}
\tag{2.3}
\]

This representation remains regular when \(p\downarrow0\) or \(p\to\infty\).

---

# 3. Forced target angles

On both v100 boxes we have

\[
0<dq<\frac\pi4
\]

and \(d<c\).  Hence the target imaginary parts never vanish and one can use fixed argument branches.

Define

\[
A_+
=
\arctan
\frac{
\rho(\cos y-\sin y)-1
}{
1+\rho(\cos y+\sin y)
},
\]

\[
A_-
=
\arctan
\frac{
\rho(\cos y+\sin y)-1
}{
1+\rho(\cos y-\sin y)
}.
\]

The required lifted angular drops are

\[
\boxed{\delta_+=\frac{5\pi}{2}+A_+},
\]

\[
\boxed{\delta_-=\frac{3\pi}{2}-A_-}.
\]

Therefore define the two lifted phase residuals

\[
\boxed{
E_+=\delta_+-\Delta\,g_a(s_+),
}
\tag{3.1}
\]

\[
\boxed{
E_-=\delta_--\Delta\,g_a(s_-).
}
\tag{3.2}
\]

An exact square on the selected lifted sheet is equivalent to

\[
\boxed{E_+=E_-=0.}
\tag{3.3}
\]

---

# 4. Radial-extreme condition is automatic

For both v100 boxes,

\[
0<d<c.
\]

For every \(q>0\),

\[
\sin(dq)<dq<cq<\sinh(cq).
\]

But the radial-extreme condition for the forced square targets is exactly

\[
\sin(dq)<\sinh(cq).
\]

Thus throughout both boxes,

\[
\boxed{1<R_-<R_+<\rho.}
\tag{4.1}
\]

Hence a simultaneous phase zero produces four distinct spiral vertices in temporal order

\[
\boxed{0,1,3,2.}
\]

---

# 5. Box A: every \(p\ge2\)

For

\[
p\ge2
\]

we have

\[
0<q\le\frac13.
\]

Use the fixed rectangle

\[
\boxed{
1\le c\le\frac{19}{10},
\qquad
\frac9{20}\le d\le\frac{18}{25}.
}
\tag{5.1}
\]

Use \((E_-,E_+)\) as the Poincare--Miranda coordinates.

The validated face enclosures certify uniformly in

\[
0\le q\le\frac13
\]

that

\[
\boxed{c=1:\quad E_->0,}
\tag{5.2}
\]

\[
\boxed{c=\frac{19}{10}:\quad E_-<0,}
\tag{5.3}
\]

\[
\boxed{d=\frac9{20}:\quad E_+>0,}
\tag{5.4}
\]

\[
\boxed{d=\frac{18}{25}:\quad E_+<0.}
\tag{5.5}
\]

The adaptive enclosure verifier regenerates these face certificates with zero unresolved cells.  Its certified leaf margins are:

\[
\begin{array}{c|c}
\text{face}&\text{certified margin}\\ \hline
c=1&E_->5.22248\times10^{-2}\\
c=19/10&-E_->1.80874\times10^{-2}\\
d=9/20&E_+>6.42061\times10^{-6}\\
d=18/25&-E_+>1.75117\times10^{-5}
\end{array}
\]

Poincare--Miranda therefore gives

\[
E_+=E_-=0
\]

for every fixed \(q\in(0,1/3]\), hence every

\[
\boxed{p\ge2.}
\]

---

# 6. Box B: \(1/200\le p\le2\)

This exponent interval is

\[
\frac13\le q\le\frac{200}{201}.
\]

Use

\[
\boxed{
\frac34\le c\le\frac32,
\qquad
\frac25\le d\le\frac{18}{25}.
}
\tag{6.1}
\]

Keep \(E_-\) as the first coordinate and define

\[
\boxed{
G=E_++\frac{57}{400}E_-.
}
\tag{6.2}
\]

The linear change

\[
(E_-,E_+)\mapsto(E_-,G)
\]

is invertible, so

\[
E_-=G=0
\iff
E_-=E_+=0.
\]

The validated face enclosures certify uniformly in the entire \(q\)-interval:

\[
\boxed{c=\frac34:\quad E_->0,}
\tag{6.3}
\]

\[
\boxed{c=\frac32:\quad E_-<0,}
\tag{6.4}
\]

\[
\boxed{d=\frac25:\quad G>0,}
\tag{6.5}
\]

\[
\boxed{d=\frac{18}{25}:\quad G<0.}
\tag{6.6}
\]

Certified leaf margins:

\[
\begin{array}{c|c}
\text{face}&\text{certified margin}\\ \hline
c=3/4&E_->1.06004\times10^{-3}\\
c=3/2&-E_->8.40199\times10^{-4}\\
d=2/5&G>1.35811\times10^{-6}\\
d=18/25&-G>1.48986\times10^{-6}
\end{array}
\]

Thus Poincare--Miranda gives a square for every

\[
\boxed{\frac1{200}\le p\le2.}
\]

---

# 7. Small exponents from v99

v99 proved, by an exact rational Poincare--Miranda argument around the logarithmic-spiral endpoint, that

\[
\boxed{0<p\le\frac1{200}}
\]

also contains a non-graceful \(0132\) square.

Combining Sections 5--7 covers every positive exponent.

---

# 8. Main theorem

> **Theorem 8.1 — Power-spiral square theorem.**  For every real
> \[
> p>0,
> \]
> the spiral
> \[
> \gamma_p(t)=t e^{i/t^p},\qquad0<t\le1,
> \]
> contains an inscribed Euclidean square.  Because the radius is strictly increasing, the square cannot be graceful.  The construction above supplies the radial temporal order
> \[
> \boxed{0132}.
> \]

This theorem is much stronger than the original v94 model claim, which only had a numerical square for \(p=2\).

---

# 9. The vertices remain on \(0<t\le1\)

Let \(T=t_0^{-p}\) be the first angular coordinate.  The far endpoint satisfies

\[
t_3^p=\frac{\rho^p}{T}.
\]

Using

\[
T=\frac{\Delta}{1-\rho^{-p}},
\]

we get

\[
t_3\le1
\iff
\rho^p-1\le\Delta.
\]

On Boxes A and B,

\[
c(1-q)\le1.9,
\]

so

\[
\rho^p=e^{c(1-q)}<e^2<9.
\]

Hence

\[
\rho^p-1<8.
\]

Meanwhile

\[
\Delta=3\pi-dq>9-0.72=8.28.
\]

Therefore

\[
\boxed{t_3<1.}
\]

The other three radii are smaller.

---

# 10. Validation method

`verify_v100_power_spiral_all_exponents.py` regenerates every face certificate from scratch.

The verifier does **not** use Newton roots or sampled extrema.

Its analytic evaluator uses:

- outward-rounded decimal arithmetic;
- explicit alternating Taylor bounds for \(\sin\), \(\cos\), \(\arctan\), and \(E(x)=(1-e^{-x})/x\);
- an explicit positive Taylor enclosure for \(e^x\);
- the identity
  \[
  \log(1+z)=2\operatorname{atanh}\frac{z}{2+z}
  \]
  to avoid a black-box logarithm;
- explicit tail estimates for each series;
- adaptive bisection only on cells whose current enclosure does not yet determine the required sign.

The final adaptive certificate has **zero unresolved cells** on all eight faces.

For endpoint safety, the machine boxes slightly over-cover the exact rational \(q\)-intervals:

\[
0\le q\le0.333333333333333334,
\]

and

\[
0.333333333333333333
\le q\le
0.995024875621890548.
\]

These strictly contain \(1/3\) and \(200/201\).

### Audit note

This is a custom outward-rounded enclosure verifier rather than a SageMath interval calculation.  All transcendental enclosures are explicit in the source.  Before publication, an independent SageMath/Arb reproduction of the eight face signs would still be prudent.

---

# 11. What this says about the broader Toeplitz program

The original spiral lead asked whether geometry invisible to graceful-square detection could force a square in another order.

For the entire power family, the answer is now yes:

\[
\boxed{
\text{strict radial monotonicity excludes graceful squares, but every }t e^{i/t^p}\text{ has a }0132\text{ square.}
}
\]

This also shows that the symplectic-area threshold at \(p=2\) is **not** a square-existence threshold.  Squares persist in both the area-tame range \(0<p<2\) and the area-active range \(p\ge2\).

The next useful general theorem should therefore isolate which features of the power-law phase function make the two-phase degree unavoidable, rather than using area divergence alone as the local hypothesis.

A natural target is a perturbative radial-phase theorem: prove that a strictly radial graph whose logarithmic angular profile lies in an appropriate cone/neighborhood of a power law inherits the nonzero \(0132\) phase degree.

---

# Claim boundary

## Proved/computer-assisted in v100

- every \(p>0\) power spiral \(t e^{i/t^p}\) contains a square;
- the square can be taken in temporal order \(0132\);
- one uniform Miranda box covers every \(p\ge2\);
- a second uniform Miranda box covers \(1/200\le p\le2\);
- v99 covers \(0<p\le1/200\);
- all eight face certificates regenerate with zero unresolved adaptive cells.

## Still not proved

- the analogous theorem for every monotone radial spiral;
- that divergent \(\int r^2|d\theta|\) forces a square for arbitrary radial spirals;
- a bridge from failure of the Haar-area primitive to the power-law/radial phase geometry;
- the universal Toeplitz square-peg conjecture.
