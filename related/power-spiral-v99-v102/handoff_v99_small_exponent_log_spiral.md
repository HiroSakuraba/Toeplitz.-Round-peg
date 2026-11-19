# Project Handoff v99: The Logarithmic-Spiral Endpoint and an Explicit Small-Exponent Square Theorem

**Author:** Benjamin John Schulz

## Executive state

v97 proved a non-graceful square on

\[
\gamma_2(t)=t e^{i/t^2}.
\]

v98 extended this to an explicit interval around \(p=2\) and to every sufficiently large exponent.

v99 closes the opposite endpoint.

> **Theorem A — Logarithmic-spiral limit square.**  
> There exists a logarithmic spiral
> \[
> \boxed{
> L_k(r)=r e^{-ik\log r}
> }
> \]
> containing an inscribed Euclidean square in temporal radial order
> \[
> \boxed{0132}.
> \]

More precisely, the v98 power-spiral branch has a nondegenerate \(p\to0^+\) limit inside the explicit box

\[
\boxed{
\frac{119}{50}\le \rho\le\frac{62}{25},
}
\]

\[
\boxed{
\frac{27}{50}\le\varepsilon\le\frac{61}{100}.
}
\]

> **Theorem B — Explicit small-exponent power spirals.**  
> For every
> \[
> \boxed{
> 0<p\le\frac1{200},
> }
> \]
> the power spiral
> \[
> \boxed{
> \gamma_p(t)=t e^{i/t^p},
> \qquad 0<t\le1,
> }
> \]
> contains an inscribed non-graceful square of temporal order \(0132\).

The same fixed Miranda box works for the entire interval.

Numerically the limiting logarithmic square is located at

\[
\boxed{
\rho_0
\approx
2.4290082042304468262,
}
\]

\[
\boxed{
\varepsilon_0
\approx
0.5742561839941456240,
}
\]

with logarithmic pitch

\[
\boxed{
k_0
=
\frac{3\pi-\varepsilon_0}{\log\rho_0}
\approx
9.9726096175632322563.
}
\]

The numerical root is only a locator.  Existence is proved by Poincare--Miranda with rational inequalities.

This means the selected \((-1,-1)\) radial square branch is now rigorous at all three regimes:

1. \(p\to0^+\);
2. a neighborhood of \(p=2\);
3. \(p\to\infty\).

The remaining power-family problem is compact in exponent.

---

# Part I. The \(p\to0\) normalized curve

For

\[
\gamma_p(t)=t e^{i/t^p},
\]

choose normalized opposite endpoint geometry as in v98:

\[
z_0=1,
\qquad
z_2=\rho e^{-i\Delta},
\]

where

\[
\Delta=3\pi-\varepsilon.
\]

The first endpoint angular coordinate is

\[
T_p
=
\frac{\Delta}{1-\rho^{-p}}.
\]

For a normalized intermediate radius \(R\),

\[
T_p(R^{-p}-1)
=
-\Delta
\frac{1-R^{-p}}{1-\rho^{-p}}.
\]

As \(p\downarrow0\),

\[
\frac{1-R^{-p}}{1-\rho^{-p}}
\longrightarrow
\frac{\log R}{\log\rho}.
\]

Hence the normalized power spiral converges to

\[
\boxed{
R
\longmapsto
R
\exp\left(
-i\frac{\Delta}{\log\rho}\log R
\right).
}
\]

This is a logarithmic spiral with pitch

\[
\boxed{
k=\frac{\Delta}{\log\rho}.
}
\]

---

# Part II. Forced target geometry

Define

\[
A=\frac{1+\rho^2}{2},
\qquad
S=\rho\sin\varepsilon,
\]

\[
B=A^2-S^2.
\]

The two forced square-target radii satisfy

\[
R_+^2=A+S,
\qquad
R_-^2=A-S.
\]

On the box used below the radial-extreme condition holds, so

\[
1<R_-<R_+<\rho.
\]

Use the same target argument combinations as v97:

\[
\phi_+-\phi_-=-\pi-\beta,
\]

where

\[
\boxed{
\beta
=
\arctan
\frac{\rho^2-1}{2\rho\cos\varepsilon},
}
\]

and

\[
\phi_++\phi_-=\varphi,
\]

where

\[
\boxed{
\varphi
=
\arctan
\frac{\rho^2\sin2\varepsilon}
{1+\rho^2\cos2\varepsilon}.
}
\]

All argument branches are fixed throughout the box because

\[
0<\varepsilon<\frac\pi4
\]

and the same quadrant inequalities used in v97 hold with larger margins.

---

# Part III. Limiting phase functions

Set

\[
L=\log\rho.
\]

Because

\[
\log\frac{R_+}{R_-}
=
\operatorname{atanh}\frac SA,
\]

the \(p=0\) phase coordinates are

\[
\boxed{
H_1^{(0)}
=
\pi+\beta
-
\frac{\Delta}{L}
\operatorname{atanh}\frac SA.
}
\tag{3.1}
\]

Also

\[
R_+R_-=\sqrt B,
\]

so

\[
\boxed{
H_2^{(0)}
=
2\pi
-
\frac{\Delta}{4L}\log B
-
\frac{\varphi}{2}.
}
\tag{3.2}
\]

The equation

\[
H_1^{(0)}=H_2^{(0)}=0
\]

is exactly the selected logarithmic-spiral square sheet

\[
f_+=f_-=-2\pi.
\]

---

# Part IV. Explicit Miranda rectangle

Use

\[
\boxed{
\rho_-
=
\frac{119}{50},
\qquad
\rho_+
=
\frac{62}{25},
}
\]

\[
\boxed{
\varepsilon_-
=
\frac{27}{50},
\qquad
\varepsilon_+
=
\frac{61}{100}.
}
\]

Thus

\[
2.38\le\rho\le2.48,
\]

\[
0.54\le\varepsilon\le0.61.
\]

---

# Part V. Monotonicity of \(H_1^{(0)}\)

## 5.1 Increasing in \(\rho\)

Let

\[
u=
\frac{2\rho\sin\varepsilon}{1+\rho^2}.
\]

For

\[
\rho>1,
\]

\(u\) decreases strictly with \(\rho\).

Hence

\[
\operatorname{atanh}u
\]

decreases while

\[
\log\rho
\]

increases.

Therefore

\[
\frac{\operatorname{atanh}u}{\log\rho}
\]

decreases strictly.

Meanwhile

\[
\frac{\rho^2-1}{2\rho\cos\varepsilon}
=
\frac{\rho-\rho^{-1}}{2\cos\varepsilon}
\]

increases strictly.

Thus

\[
\boxed{
\partial_\rho H_1^{(0)}>0.
}
\tag{5.1}
\]

No numerical estimate is needed.

## 5.2 Decreasing in \(\varepsilon\)

Put

\[
F=
\frac{\operatorname{atanh}u}{\log\rho}.
\]

Then

\[
\partial_\varepsilon H_1^{(0)}
=
\beta_\varepsilon
+
F
-
\Delta F_\varepsilon.
\]

On the full box, elementary rational Taylor bounds give

\[
\boxed{
\Delta F_\varepsilon
-
F
-
\beta_\varepsilon
>
4.66.
}
\tag{5.2}
\]

Therefore

\[
\boxed{
\partial_\varepsilon H_1^{(0)}<0.
}
\tag{5.3}
\]

The verifier evaluates the lower margin using exact rational arithmetic only.

---

# Part VI. Monotonicity of \(H_2^{(0)}\)

Define

\[
G=
\frac{\log B}{4\log\rho}.
\]

Then

\[
H_2^{(0)}
=
2\pi-\Delta G-\frac\varphi2.
\]

## 6.1 Decreasing in \(\rho\)

We have

\[
\varphi_\rho>0.
\]

Also

\[
G_\rho>0.
\]

Indeed,

\[
G_\rho
=
\frac{
(\log\rho)B_\rho/B
-
(\log B)/\rho
}{
4(\log\rho)^2
}.
\]

On the whole rectangle the exact rational lower estimate is

\[
\boxed{
(\log\rho)\frac{B_\rho}{B}
-
\frac{\log B}{\rho}
>
0.0979.
}
\tag{6.1}
\]

Therefore

\[
\boxed{
\partial_\rho H_2^{(0)}<0.
}
\tag{6.2}
\]

## 6.2 Increasing in \(\varepsilon\)

Because

\[
B_\varepsilon
=
-2\rho^2\sin\varepsilon\cos\varepsilon<0,
\]

we have

\[
G_\varepsilon<0.
\]

On our box

\[
\cos2\varepsilon>0,
\]

and the exact target-angle formula gives

\[
\frac{\varphi_\varepsilon}{2}<1.
\]

The verifier obtains

\[
\boxed{
G
+
\Delta|G_\varepsilon|
-
1
>
0.698.
}
\tag{6.3}
\]

Hence

\[
\boxed{
\partial_\varepsilon H_2^{(0)}>0.
}
\tag{6.4}
\]

---

# Part VII. Exact corner signs

Only two corners need to be checked.

Use:

\[
\frac{333}{106}<\pi<\frac{355}{113},
\]

the alternating Taylor bounds for \(\sin,\cos,\arctan\), and

\[
\log x
=
2\operatorname{atanh}
\frac{x-1}{x+1}.
\]

Every estimate below is verified by exact rational arithmetic.

At

\[
\left(
\frac{119}{50},
\frac{27}{50}
\right)
\]

we obtain

\[
\boxed{
H_1^{(0)}>0.046,
}
\]

\[
\boxed{
H_2^{(0)}>0.018.
}
\]

At

\[
\left(
\frac{62}{25},
\frac{61}{100}
\right)
\]

we obtain

\[
\boxed{
H_1^{(0)}<-0.037,
}
\]

\[
\boxed{
H_2^{(0)}<-0.019.
}
\]

The monotonicity results turn these into the four full side signs.

Poincare--Miranda therefore gives an exact zero in the rectangle.

This proves Theorem A.

---

# Part VIII. Persistence for positive \(p\)

The v98 exponential-interpolation formula gives

\[
g_a(s)
=
\frac{1-e^{-as}}{1-e^{-a}},
\qquad
a=p\log\rho.
\]

Its universal derivative bound is

\[
0<\partial_ag_a(s)\le\frac14.
\]

Therefore throughout the v99 box

\[
\boxed{
|\partial_pH_i|
\le
\frac{\Delta\log\rho}{4}
<
\frac{7881}{2260}
<3.488,
}
\tag{8.1}
\]

for

\[
i=1,2.
\]

If

\[
0<p\le\frac1{200},
\]

the total displacement from the \(p=0\) side values is less than

\[
\boxed{
\frac{7881}{452000}
\approx0.017436.
}
\]

The four conservative \(p=0\) margins retain signs:

\[
0.046-\frac{7881}{452000}>0,
\]

\[
0.037-\frac{7881}{452000}>0,
\]

\[
0.018-\frac{7881}{452000}
=
\boxed{
\frac{51}{90400}>0,
}
\]

\[
0.019-\frac{7881}{452000}>0.
\]

Thus the same Miranda box works simultaneously for every

\[
0\le p\le\frac1{200}.
\]

For \(p>0\) this proves Theorem B.

---

# Part IX. Numerical limiting square

Newton iteration on the limiting equations locates

\[
\rho_0
=
2.4290082042304468262203585655\ldots,
\]

\[
\varepsilon_0
=
0.5742561839941456239982256775\ldots.
\]

The corresponding logarithmic pitch is

\[
k_0
=
\frac{3\pi-\varepsilon_0}{\log\rho_0}
=
9.9726096175632322562823305191\ldots.
\]

The two intermediate normalized radii are approximately

\[
R_-
=
1.45964989752102648\ldots,
\]

\[
R_+
=
2.18391919101547275\ldots.
\]

So the radial order is

\[
1<R_-<R_+<\rho_0,
\]

hence temporal square order

\[
\boxed{0132}.
\]

Again, these decimals are not theorem inputs.

---

# Part X. State of the power family

The selected square branch is now rigorous in:

### Small exponent

\[
\boxed{
0<p\le1/200.
}
\]

### Critical interval

v98 proved

\[
\boxed{
97/50\le p\le103/50.
}
\]

### Large exponent

v98 proved existence for all sufficiently large \(p\).

Numerically the branch remains nondegenerate for every tested

\[
p>0.
\]

Thus the conjecture

\[
\boxed{
\forall p>0,\quad
t e^{i/t^p}
\text{ contains a non-graceful square}
}
\]

has now been reduced from two singular endpoint problems to a compact intermediate continuation problem.

---

# Part XI. Literature connection

Cole Hugelmeyer's 2023 relation-avoiding-path program independently singles out exponential/logarithmic spiral normal forms.  In one complex dimension the paper proves its Spiral Conjecture and proposes stronger monotonic versions as a route toward inscribed-square questions.

The v99 logarithmic endpoint is not an application of Hugelmeyer's theorem, but it is compatible with the same structural theme: logarithmic spirals are natural normal forms at the end of relation-avoidance/renormalization problems.

---

# Claim boundary

## Proved in v99

- the exact logarithmic-spiral \(p=0\) phase equations;
- fixed monotonicity signs on the explicit v99 box;
- exact rational side signs on that box;
- existence of a logarithmic spiral containing a \(0132\)-order square;
- persistence of the same Miranda degree for every
  \[
  0<p\le1/200.
  \]

## Numerical only

- the decimal location of the limiting root and its pitch;
- apparent continuation through the remaining intermediate exponents.

## Not proved

- the power-spiral square theorem for every \(p>0\);
- a finite cover of the remaining compact exponent interval;
- the area-active radial square conjecture for arbitrary radial spirals;
- the Haar-area-to-radial-phase bridge;
- universal square-peg existence.


---

# Post-v99 audit: global phase monotonicity is false

After completing the small-exponent proof, the four v97/v99 sign derivatives were tested on the broader strict radial-extreme domain.

The tempting universal statement

\[
H_{1,\rho}>0,\qquad
H_{1,\varepsilon}<0,\qquad
H_{2,\rho}<0,\qquad
H_{2,\varepsilon}>0
\]

is false globally.

A numerical counterexample occurs near

\[
p=0.0269650,\qquad
\rho=3.4689804,\qquad
\varepsilon=0.0554029,
\]

where the configuration is still strict radial-extreme, but

\[
\boxed{
H_{2,\varepsilon}<0.
}
\]

The other three derivative signs retain their expected directions at that sample.

This is a numerical structural warning, not a certified theorem input.

## Consequence

Do **not** try to prove the full power-family theorem by establishing the v97 four-sign monotonicity on the whole radial-extreme domain.

The better next target is a **branch corridor theorem**:

1. numerically trace the selected \((-1,-1)\) square branch on a compact \(p\)-interval;
2. construct a tubular region around that branch;
3. prove the four Miranda side signs only inside that corridor;
4. cover the compact exponent interval by finitely many overlapping certified corridors.

The endpoint theorems v98--v99 mean such a finite corridor argument would be enough to prove the power-spiral square conjecture for every \(p>0\).
