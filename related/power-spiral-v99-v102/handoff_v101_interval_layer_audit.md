# Project Handoff v101: Audit and Replacement of the v100 Interval Layer

**Author:** Benjamin John Schulz

## Scope

This revision does not touch the mathematics of v100.  The reduction of
Section 1, the compactification of Section 2, the forced target angles of
Section 3, the radial-extreme argument of Section 4, the two Poincare-Miranda
boxes of Sections 5 and 6, and the endpoint containment argument of Section 9
are unchanged.  Theorem 8.1 is unchanged.

What changes is the certified-computation layer underneath the eight face
certificates.  v100 evaluated every transcendental with a hand-written
outward-rounded Decimal enclosure.  That layer has been audited, three defects
were found, and it has been replaced.  The eight face signs have then been
reproduced a second time in Arb ball arithmetic on a disjoint code path.

Section 10 of the v100 handoff asked for exactly this before treating the
computer-assisted component as publication-final.  It is now done, with one
correction that the v100 certificate needed.

---

# 1. Audit findings on the v100 layer

## 1.1 Unsound enclosure: `Efun_scalar` at negative argument

`Efun_scalar(x)` encloses

\[
E(x)=\frac{1-e^{-x}}{x}=\sum_{n\ge0}\frac{(-1)^n x^n}{(n+1)!}
\]

by an alternating-series bracket: the lower bound is the largest odd-index
partial sum, the upper bound the smallest even-index partial sum.  That
bracket is valid only for \(x\ge0\).

For \(x<0\) every term is positive, the partial sums increase, and the
smallest even-index partial sum is \(v_0=1\).  Since \(E(x)>1\) for \(x<0\),
the returned upper bound lies below the true value.

Instrumentation of the eight v100 faces shows the branch is reached.  The
most negative argument encountered is

\[
x=-0.140898577957246973939300440805816705\ldots
\]

at which `Efun_scalar` returns the interval \([1.073878167822,\,1.000000000000]\)
while the true value is \(1.073878167822\).  The returned interval does not
contain the value, and its endpoints are inverted.

The negative argument arises because \(a s_\pm\) is evaluated with \(s_\pm\)
an outward-rounded enclosure whose lower endpoint can fall below zero on a
coarse cell, not because \(a s_\pm\) is mathematically negative.

**Impact.**  Of the 9644 certified leaves across the eight v100 faces, exactly
one passed its sign test through this branch: face `B_c_low`, cell

\[
d\in[0.56,\,0.72],
\qquad
q\in\left[\tfrac13,\,0.6641791044776119405\right],
\qquad
\text{depth }2 .
\]

On that cell the corrected layers give

| layer | enclosure of \(E_-\) |
|---|---|
| v100 (defective) | \([1.0135815599,\,4.9547286393]\) |
| v101 Decimal | \([0.9879463673,\,4.9582334812]\) |
| v101 Arb | \([0.5226679697,\,5.2233090802]\) |

The v100 enclosure is narrower than the sound one, as expected from an
upper bound that was too small.  The required sign \(E_->0\) holds in every
case with margin at least \(0.52\), so the conclusion on that cell is not in
doubt.  The defect is real and was load-bearing for one leaf, but it did not
change any face sign.

## 1.2 Accuracy failure: `atan_pos_scalar` near argument one

`atan_pos_scalar` reduces \(x>1\) by \(\arctan x=\pi/2-\arctan(1/x)\) and then
uses the alternating Gregory series with 80 terms.  Across the eight faces the
reduced argument reaches

\[
|x|=0.9999895610136561,
\]

where the first omitted term is

\[
\frac{|x|^{159}}{159}=6.28\times10^{-3}.
\]

The enclosure is sound but three orders of magnitude wider than the certified
margins.  This is the reason the two \(d\)-faces of Box B needed several
thousand adaptive cells.

## 1.3 Latent domain gaps

The guard in `evalE` tests \(z_-\) only:

```
if zm.lo <= D('-0.995'): raise
```

\(z_+\) is unguarded.  If \(z_+\) fell below \(-1\), the negative branch of
`L1_scalar` would form `tail = power/(N+1)/(1-a)` with \(1-a\le0\) and return
an interval with inverted endpoints rather than raising.  Instrumentation
shows \(z_\pm\) never goes below \(-0.620527\) on the eight faces, so this is
unexercised rather than active.

Separately, the stated domains in the source are wrong.  The comments claim
\(x\in[0,2]\) for `Efun_scalar` and \(|z|<0.7\) for `L1_scalar`.  The measured
arguments are

| evaluator | comment | measured range over the eight faces |
|---|---|---|
| `dec_exp_bounds` | \(0\le x\le1.6\) | \([0,\,1.492538]\) |
| `sin_scalar`, `cos_scalar`, `sinc_scalar` | \(\lvert x\rvert\le0.8\) | \([0,\,0.716418]\) |
| `coshC_scalar` | none | \([0,\,1.492538]\) |
| `L1_scalar` | \(\lvert z\rvert<0.7\) | \([-0.620527,\,2.050882]\) |
| `Efun_scalar` | \(x\in[0,2]\) | \([-0.140899,\,2.224946]\) |
| `atan_scalar` | none | \([-0.408983,\,5.547701]\) |

For \(2<x<3\) the `Efun_scalar` bracket happens to remain valid, because the
only non-monotone step is at \(n=0\) and \(v_0=1\) is still an upper bound
when \(x>0\).  This was luck, not design.

## 1.4 Rounding discipline

Every scalar series in v100 accumulates its partial sums in the ambient
Decimal context, which is round-half-even at 60 digits, and then widens by a
single `next_minus` or `next_plus`.  With 40 to 180 terms the accumulated
rounding error can exceed one unit in the last place, so the widening does
not provably cover it.  At \(10^{-58}\) against margins of \(10^{-6}\) this is
a gap in the proof rather than in the numbers, but no v100 evaluator's
enclosure was actually established.

---

# 2. The replacement layer: `iv101.py`

Design rules, each aimed at one finding above.

1. Every operation touching an enclosure endpoint runs under an explicit
   directed rounding context.  No step relies on one ulp of widening covering
   half-even error.
2. Every series is a partial sum accumulated in interval arithmetic plus an
   explicitly bounded remainder interval, with the remainder bounded
   geometrically by magnitude.  No evaluator uses an alternating-bracket
   argument, so none can fail the way `Efun_scalar` did.
3. Every evaluator asserts its validated argument domain and raises on exit.
   The asserted domains are the measured ranges of Section 1.3 with margin.
4. \(\pi\) is certified in-module by Machin's formula
   \(\pi=16\arctan\frac15-4\arctan\frac1{239}\)
   using the same arctan evaluator, instead of the rational literals
   \(333/106\) and \(355/113\).
5. `atan` applies four half-angle reductions
   \[
   \arctan x=2\arctan\frac{x}{1+\sqrt{1+x^2}}
   \]
   before the series, sending any \(x\ge0\) into \([0,\tan(\pi/32)]\subset[0,0.0985]\).

Working precision is 40 decimal digits.

## 2.1 Evaluator validation

Each evaluator was checked against mpmath at 100 decimal digits on 404 points
covering its asserted domain, including both endpoints, zero, and 400 random
interior points.

| evaluator | max enclosure width | containment |
|---|---|---|
| `exp` | \(1.412\times10^{-33}\) | PASS |
| `sin` | \(4.572\times10^{-25}\) | PASS |
| `cos` | \(1.314\times10^{-23}\) | PASS |
| `sinc` | \(5.715\times10^{-25}\) | PASS |
| `coshC` | \(6.672\times10^{-25}\) | PASS |
| `L1` | \(2.932\times10^{-34}\) | PASS |
| `Efun` | \(2.646\times10^{-34}\) | PASS |
| `atan` | \(1.2\times10^{-38}\) | PASS |
| \(\pi\) | \(7.5\times10^{-38}\) | PASS |

The arctan worst case improves from \(6.28\times10^{-3}\) to \(1.2\times10^{-38}\).

## 2.2 Face certificates with the new layer

Boxes, residuals, coordinate choice, \(\alpha=57/400\), and maximum depth 18
are all as in v100.

| face | nodes | leaves | maxdepth | unresolved | certified margin |
|---|---|---|---|---|---|
| `A_c_low` | 9 | 5 | 4 | 0 | \(5.235007\times10^{-2}\) |
| `A_c_high` | 47 | 24 | 8 | 0 | \(1.821355\times10^{-2}\) |
| `A_d_low` | 1959 | 980 | 11 | 0 | \(2.151460\times10^{-4}\) |
| `A_d_high` | 2745 | 1373 | 13 | 0 | \(1.121684\times10^{-5}\) |
| `B_c_low` | 25 | 13 | 8 | 0 | \(1.185224\times10^{-3}\) |
| `B_c_high` | 47 | 24 | 8 | 0 | \(9.611056\times10^{-4}\) |
| `B_d_low` | 5439 | 2720 | 14 | 0 | \(5.525805\times10^{-6}\) |
| `B_d_high` | 8865 | 4433 | 14 | 0 | \(7.574730\times10^{-6}\) |

Zero unresolved cells on all eight faces.  The worst margin improves from
\(1.35811\times10^{-6}\) in v100 to \(5.525805\times10^{-6}\).

---

# 3. Independent Arb reproduction

`verify_v101_arb.py` shares no code with the Decimal path.  Differences are
deliberate:

- cell endpoints are exact `Fraction` values, subdivided exactly, so the
  boxes are the true rational intervals \(q\in[0,1/3]\) and
  \(q\in[1/3,200/201]\) with no decimal over-cover;
- all enclosures are Arb balls at 160-bit working precision;
- every transcendental is Arb's own implementation: `exp`, `sin`, `cos`,
  `sinh`, `log1p`, `expm1`, `atan`, `pi`.  No series in this project is used;
- removable singularities are handled by closed forms that do not cancel, in
  particular
  \[
  C(x)=\frac{\cosh x-1}{x^2}=\frac12\left(\frac{\sinh(x/2)}{x/2}\right)^2 ,
  \qquad
  E(x)=\frac{-\operatorname{expm1}(-x)}{x}.
  \]

Two adjustments were needed and both are worth carrying back into the Decimal
driver.

First, Arb ball arithmetic overestimates \(c^2\) badly on a full-width cell:
for \(c\in[3/4,3/2]\) the ball square gives \([0.28,2.25]\) against the true
\([0.5625,2.25]\).  On the root cell of `B_d_high` that pushed \(z_-\) below
\(-1\).  An endpoint-based square removes the blowup at its source.

Second, a cell whose enclosure is too coarse to evaluate is now treated as a
subdivision request rather than a hard failure.  This is sound because the
two children cover the parent, and it turns a spurious abort into ordinary
branch-and-bound behaviour.  A cell still unevaluable at maximum depth is
reported as unresolved.

Result, 1.3 seconds total:

| face | nodes | leaves | maxdepth | unresolved | certified margin |
|---|---|---|---|---|---|
| `A_c_low` | 11 | 6 | 4 | 0 | \(9.451231\times10^{-2}\) |
| `A_c_high` | 53 | 27 | 8 | 0 | \(5.255473\times10^{-3}\) |
| `A_d_low` | 2027 | 1014 | 12 | 0 | \(6.767832\times10^{-5}\) |
| `A_d_high` | 2811 | 1406 | 13 | 0 | \(3.620809\times10^{-5}\) |
| `B_c_low` | 35 | 18 | 9 | 0 | \(1.277418\times10^{-2}\) |
| `B_c_high` | 57 | 29 | 8 | 0 | \(3.930228\times10^{-3}\) |
| `B_d_low` | 5513 | 2757 | 14 | 0 | \(6.283801\times10^{-5}\) |
| `B_d_high` | 9133 | 4567 | 14 | 0 | \(2.183953\times10^{-5}\) |

Margins are not comparable face by face between the two tables, since the
adaptive partitions differ.  What is comparable is the sign of every face and
the absence of unresolved cells, and those agree.

---

# 4. Revised claim boundary

## Established

- The eight face signs of v100 Boxes A and B hold, certified twice on
  independent arithmetic stacks with zero unresolved cells.
- Each of the eight transcendental evaluators in `iv101` is validated for
  containment on its asserted domain against a 100-digit reference.
- The Arb reproduction runs on the exact rational boxes, so the decimal
  over-cover of v100 Section 10 is no longer part of the argument.
- The one v100 leaf that depended on the defective `Efun_scalar` branch is
  re-certified with margin at least \(0.52\).

## Not established

- Machine-checked correctness of the interval layer itself.  Both stacks are
  tested, not verified.  A Lean or Coq formalisation of the enclosure
  arithmetic remains open, as does a translation of the face certificates
  into a proof assistant.
- The remainder bounds in `iv101` are proved on paper in the source comments
  and checked numerically, but not formally.
- Everything listed as unproved in the v100 claim boundary stays unproved.
  The audit changes the standing of the computation, not of the theorem's
  reach.

## Recommended before publication

- Port the endpoint-square and subdivide-on-unevaluable changes back into the
  Decimal driver, so both stacks use the same branch-and-bound policy.
- Add the \(z_+\) domain guard, which is currently unexercised but only by
  accident of the box choice.
- Consider dropping the Decimal layer entirely in favour of the Arb path.  It
  is 200 times faster, runs on exact rational cells, and does not require the
  project to own its own transcendental implementations.  The Decimal layer's
  only advantage is that it has no compiled dependency.

---

# Appendix A: reproduction

```
python test_iv101.py                  # evaluator containment against mpmath
python verify_v101_faces.py           # eight faces, Decimal layer, about 6 min
python verify_v101_arb.py             # eight faces, Arb layer, about 1.3 s
python instrument_ranges.py           # measured argument ranges of Section 1.3
python probe_bugs.py                  # reproduces findings 1.1 and 1.2
python audit_efun_impact.py           # locates the one tainted v100 leaf
```

`verify_v101_arb.py` requires `python-flint`.  The other scripts require only
the standard library, plus `mpmath` for the containment tests.

# Appendix B: files

| file | role |
|---|---|
| `iv101.py` | replacement interval layer |
| `verify_v101_faces.py` | eight face certificates on `iv101` |
| `verify_v101_arb.py` | independent Arb reproduction |
| `test_iv101.py` | evaluator containment tests |
| `instrument_ranges.py` | argument-range instrumentation of the v100 layer |
| `probe_bugs.py` | minimal reproductions of the two defects |
| `audit_efun_impact.py` | impact analysis of finding 1.1 |
| `_core_v100.py` | v100 evaluator core, split off for instrumentation |
| `v101_verifier_output.txt` | captured output of the three verification runs |
