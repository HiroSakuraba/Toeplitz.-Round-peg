# v236 — A provably square-free logarithmic spiral, and a certified exceptional pitch

**Research program:** Benjamin John Schulz  
**19 November 2026**  
**Status:** Elementary analytical theorem, supported by exact rational inequality regressions; separate historical v99 Poincaré–Miranda existence certificate rerun locally. Not a Toeplitz counterexample, not a classification of all logarithmic-spiral pitches, and not a result about a closed Jordan curve.

## 1. Research trigger and corrected premise

The user-supplied independent brainstorm, *Toeplitz attack ideas.html* (19 November 2026), proposes testing whether an open logarithmic spiral has no inscribed squares, *except perhaps at isolated openness values*. The **exception clause is essential**: the project's recovered v99 manuscript and exact verifier already prove that **some** logarithmic spiral contains an inscribed square, with numerical pitch `k ≈ 9.972609617563232`, so openness `a=1/k ≈ 0.1002747`.

We now prove the complementary statement that **all sufficiently open logarithmic spirals are square-free**, with a simple fully explicit threshold **`a >= 4`**. Thus the proposed local model is mathematically substantive, but the headline 'log spirals are square-free' would be false. The classification of `0<a<4` remains open here, including whether the square-bearing pitches are isolated.

## 2. Main theorem: no square if `a>=4`

For a real `a>0`, define the **compactified open logarithmic spiral**

\[
S_a=\{e^{(-a+i)t}:t\ge0\}\cup\{0\}.
\]

This is an injective arc with an accumulation endpoint, **not** a Jordan curve. An inscribed square means four distinct points of `S_a` forming a nondegenerate Euclidean square, in any order; the square need not lie inside any domain.

> **Theorem 1 (large-openness square exclusion).** For every real `a>=4`, the compactified open logarithmic spiral `S_a` contains **no nondegenerate inscribed square**.

### Proof, with all constants explicit

Set `epsilon=1/4`. Suppose by contradiction that the four corners of a square of side `s>0` lie on `S_a`. Let `R` be the **largest** of their distances from the origin, so `R>0`. Rotate the plane so that the farthest corner lies on the positive real axis. Parameterize by radial coordinate `r=exp(-at)`; within the annulus `epsilon R <= r <= R`, the relative polar angle is

\[
\vartheta(r)=\frac1a\log\frac Rr\in[0,\delta],\quad
\delta=\frac{\log 4}{a}.
\]

The rational Taylor inequality

\[
e^{3/2}>1+\frac32+\frac98+\frac9{16}=\frac{67}{16}>4
\]

implies `log 4 < 3/2`. Hence for `a>=4`,

\[
0\le\delta<\frac{3}{8}.
\tag{1}
\]

Call a square corner *inner* if `r<epsilon R`, and *outer* otherwise. The farthest corner is outer. Since the distance between any two square corners is at least `s` and at most `sqrt(2)s`, the distance from the farthest to any inner corner yields

\[
s\ge\frac{(1-\varepsilon)R}{\sqrt2}=\frac{3R}{4\sqrt2}.
\tag{2}
\]

**Case A: at least two inner corners.** Their mutual distance is less than `2 epsilon R=R/2`, whereas the square-side lower bound (2) is strictly greater than `R/2` because `2sqrt(2)<3`. Contradiction.

**Case B: exactly one inner corner.** The other three corners constitute a right-isosceles triangle with two legs of length `s`. Its minimum width across any pair of parallel lines, in *any* direction, equals the altitude to its hypotenuse, `s/sqrt(2)`. In particular its width perpendicular to the outermost radial direction is at least

\[
\frac{s}{\sqrt2}\ge\frac{3R}{8}.
\tag{3}
\]

But every outer corner has `r<=R` and angle in `[0,delta]`. Its transverse coordinate is therefore in `[0,R delta]`; all three outer corners fit in a strip of width **strictly less than** `3R/8` by (1). Contradiction.

**Case C: no inner corners.** All four corners lie on the same radial subarc `epsilon R <=r<=R`. In coordinates aligned with the outermost corner, differentiating the spiral with respect to radius gives

\[
\frac{dz}{dr}=e^{i\vartheta(r)}\left(1-\frac{i}{a}\right).
\]

The tangent makes an angle of magnitude at most `delta+atan(1/a)` with the positive radial axis. Since `atan(x)<x` for `x>0`,

\[
\delta+\arctan(1/a)<\frac38+\frac14=\frac58<\frac\pi4,
\]

using `pi>3`. Thus `dx/dr>0` and `|dy/dr|<dx/dr` throughout this *entire* radial interval. Integration shows that any two points on this arc obey the **strict global graph inequality** `|Delta y|<|Delta x|`. But a square has two perpendicular adjacent edges; for any orthogonal pair of nonzero vectors, **at least one** has `|Delta y|>=|Delta x|`. This contradicts placing all four corners on the arc.

The three exhaustive cases rule out every nondegenerate square. The proof includes the limiting point `0` as a possible corner in the inner category. QED.

### Corollaries and interpretation

1. Every subarc of `S_a` is also square-free if `a>=4`.
2. Logarithmic spirals are **not uniformly square-free over all openness parameters**: v99 supplies at least one special square-bearing pitch.
3. If an open log spiral `S_a` **does** contain a square, self-similarity produces a one-parameter family of smaller inscribed squares: shifting all four spiral parameters by the same `h>=0` multiplies every corner by `e^{(-a+i)h}`. The square locus therefore cannot be treated as isolated in the four raw source parameters. This is a warning against naively applying the ordinary square-branch implicit-function theorem to infer an open interval of square-bearing `a` values.
4. **Do not infer anything about a Jordan closure** from the square-free property of this single arc. A closing connector or second arm can support squares that do not lie on the first arm. The unrestricted Toeplitz conjecture remains open.

## 3. Independent positive control: v99 special pitch

The preserved `related/power-spiral-v99-v102/handoff_v99_small_exponent_log_spiral.md` formulates two real phase equations on the rational rectangle

\[
\frac{119}{50}\le\rho\le\frac{62}{25},\qquad
\frac{27}{50}\le\epsilon\le\frac{61}{100}.
\]

It proves the required boundary signs using exact rational inequalities and Poincaré–Miranda. At its zero, the logarithmic-spiral pitch is `k=(3pi-epsilon)/log rho` and the openness is `a=1/k=log rho/(3pi-epsilon)`. The numerical location `k≈9.972609617563232` is **not** a proof input. Elementary bounds show it lies in `1/20<a<1/4`: because `rho>2`, `log rho>1/2` and `Delta=3pi-epsilon<10`; also `rho<3`, `log rho<2` and `Delta>8`. The v99 checker was rerun separately in this development session; it printed **PASS**, with exact endpoint inequalities. A separate restored CI workflow already verified v99 and the later v101 Arb proof. None of the v99 assertions implies a whole interval of openness values containing squares.

This is a counterexample to **all** log spirals being square-free; it is entirely consistent with the artifact's *possible isolated exceptional parameters* caveat. Between this certified special pitch and the `a>=4` exclusion theorem, no classification is claimed.

## 4. Reproduction

Run

```bash
python code/log_spiral_exclusion_v236.py
python -m unittest discover -s tests -p 'test_log_spiral_exclusion_v236.py' -v
```

The implementation uses Python standard-library `fractions.Fraction` for every theorem-critical numerical inequality (`e^{3/2}>4`, the two-inner separation, the three-outer strip contradiction and the graph tangent bound). The 15 new tests include independent numerical geometry *sanity* tests, clearly nonessential to the proof. The mathematical width and Lipschitz lemmas above, not the unit tests, establish the theorem. Re-running the old v99 existence checker is a distinct complementary check and does not establish all intermediate pitches.

## 5. Next research tasks

The strongest remaining log-spiral question is not the proved large-`a` exclusion. It is **classifying which intermediate pitches actually admit square configurations**. This could be approached by exact interval isolation of the two normalized phase constraints, with robust handling of special self-similar sheets, and with a decisive distinction between *absence on the open arc* and *inscriptions in a closed double-log-spiral Jordan curve*.

For general Toeplitz strategy, the local-negative model shows why a proof cannot assume each sufficiently irregular *open arc* must produce its own square. It does not say that all proofs must use global topology, and cannot rule out robust local arguments on other classes of arcs.

**Scientific status:** Unrestricted Toeplitz square-peg remains open.
