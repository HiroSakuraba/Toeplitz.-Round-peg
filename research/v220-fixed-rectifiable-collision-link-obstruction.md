# v220 — A fixed rectifiable Jordan curve with infinitely many collapsing squares: obstruction to a naive relative collision link

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** An elementary, explicit geometric countermodel to a proposed *collision-link regularity assumption*, with exact-rational checks of every finite stage. **Not a counterexample to the square-peg conjecture and not a proof of it.**

## Executive result

The v219 program sought a relative mod-2 configuration-space invariant after resolving cyclic symmetry. The next indispensable hypothesis was that the scale-normalized square-equation section has a meaningful, preferably nonvanishing, restriction to the total-collision link. For a fixed regular C¹ curve, v219 justified such a restriction by tangent estimates. For general continuous curves it was unresolved.

**Theorem (fixed-curve obstruction).** There exists a **single, rectifiable Jordan curve** \(\Gamma_\infty\subset\mathbb R^2\), of exact length \(37/6\), containing a sequence of pairwise distinct, gracefully ordered inscribed squares of side \(4^{-(n+1)}\to0\). At the accumulation parameter, two sequences of ordered four-point configurations tend to the **same naive parameter-space blow-up point**, with identical normalized successive gap ratios \((1/3,1/3,1/3)\), while the square residual divided by their common total parameter span has **different limits**: one zero, one nonzero.

Consequently **no continuous extension of this particular scale-normalized residual to the *bare parameter-gap collision chart* exists** for \(\Gamma_\infty\). Even a shape-enriched compactification that retains every limiting physical square shape cannot have a *universally zero-free* collision restriction: this curve exhibits actual small squares converging to its collision stratum. A relative obstruction must either accommodate such boundary zeroes, add still more data, or force the existence of an independent macroscopic witness. It cannot simply import the regular-\(C^1\) link argument.

This is stronger than v219's **moving-polygon** shrinking-notch example: all squares here lie on the **same fixed curve**, which also has finite length. It does not challenge the known rectangular-peg results for rectifiable curves: the curve visibly contains many squares.

## 1. Explicit infinite-notch curve

Begin with the counterclockwise boundary of \([0,2]\times[0,1]\), parametrized from \((0,0)\) to the right along the bottom. Put

\[
 a_n=2^{-n},\qquad s_n=4^{-(n+1)}=2^{-2n-2},\quad n\ge1.
\]

Replace the bottom interval from \((a_n,0)\) to \((a_n+s_n,0)\) by the three-edge outward rectangular detour

\[
P_{n,0}=(a_n,0),\quad
P_{n,1}=(a_n,-s_n),\quad
P_{n,2}=(a_n+s_n,-s_n),\quad
P_{n,3}=(a_n+s_n,0).
\]

All coordinates are rational. For every \(n\ge1\), the intervals are separated:

\[
a_n-3s_n>a_{n+1}+s_{n+1},
\]

because with \(u=2^{-n}\le1/2\), their difference is
\(u(8-13u)/16>0\). In particular the notches do not touch, and the bottom path is strictly advancing in \(x\) between each vertical detour. As \(a_n,s_n\to0\), their endpoints and detours accumulate only at \((0,0)\). The outward notches stay below \(y=0\) and do not touch the remaining rectangle boundary. The path is an injective, continuous, closed parametrized curve, so \(\Gamma_\infty\) is a **Jordan curve**.

The \(N\)-notch simple polygon \(\Gamma_N\) is obtained by using notches \(N,N-1,\dots,1\) in increasing \(x\), then finishing the rectangle. Their omitted detours have depth at most \(s_{N+1}\), giving the geometric Hausdorff bound

\[
 d_H(\Gamma_N,\Gamma_\infty)\le s_{N+1}=4^{-(N+2)}\to0.
\]

Because the total added lengths are summable, the infinite curve admits a continuous arclength parametrization. Its finite-stage and limiting lengths are exactly

\[
L_N=6+2\sum_{n=1}^N s_n
=6+\frac16(1-4^{-N}),\quad
\boxed{L_\infty=6+\frac16=\frac{37}{6}}.
\]

Its enclosed area is

\[
A_N=2+\sum_{n=1}^N s_n^2
=2+\frac1{240}(1-16^{-N}),\quad
\boxed{A_\infty=\frac{481}{240}}.
\]

These formulae also prove rectifiability. The exact checker independently verifies polygon simplicity, area and perimeter for finite \(N\), in addition to the algebraic argument for the limit.

## 2. Exact squares at every scale on the *fixed* limit curve

The four corners \((P_{n,0},P_{n,1},P_{n,2},P_{n,3})\) all belong to \(\Gamma_\infty\) and occur in that cyclic order. With \(J(x,y)=(-y,x)\), their first edge is \(e=(0,-s_n)\) and the next is \(Je=(s_n,0)\). Thus they are a square of orientation \(+1\), with exact squared side length \(s_n^2\) and squared diameter \(2s_n^2\).

For the square residual

\[
F_+(p_0,p_1,p_2,p_3)=
\bigl(p_2-p_1-J(p_1-p_0),\ p_3-p_0-J(p_1-p_0)\bigr)\in\mathbb R^4,
\]

we have \(F_+(P_{n,0},\ldots,P_{n,3})=0\) **exactly for every \(n\)**. Their vertices all converge to \((0,0)\), but remain distinct for every finite \(n\). This shows that **a fixed rectifiable Jordan curve need not have a positive minimum side length among its inscribed squares**. It does *not* mean it lacks other, macroscopic squares.

## 3. Two incompatible link limits at the same parameter-gap direction

Parametrize \(\Gamma_\infty\) by normalized arclength \(t\in[0,1)\), starting at \(\gamma(0)=(0,0)\), and set \(L=L_\infty=37/6\). The portion of curve to the left of notch \(n\) contains the earlier-in-traversal notches \(k>n\), whose **additional** vertical length is

\[
2\sum_{k>n}s_k=\frac23s_n.
\]

The four notch corners correspond to parameters starting at

\[
t_n^{\rm sq}=\frac{a_n+2s_n/3}{L},
\]

with all three successive gaps equal to \(s_n/L\). The total parameter span is

\[
\varepsilon_n=t_{n,3}-t_{n,0}=\frac{3s_n}{L}\to0.
\]

For the **flat control**, take four points on the unmodified bottom segment immediately preceding the notch:

\[
Q_{n,j}=(a_n-(3-j)s_n,0),\qquad j=0,1,2,3.
\]

The separation inequality from Section 1 shows that this complete interval is flat and contains no smaller notch. These points have the **same three arclength gaps** \(s_n/L\) and the same total span \(\varepsilon_n\), starting at

\[
t_n^{\rm flat}=\frac{a_n-3s_n+2s_n/3}{L}.
\]

Both first parameters tend to zero, and their normalized interior gap ratios agree exactly:

\[
\left(\frac{g_0}{\varepsilon_n},\frac{g_1}{\varepsilon_n},\frac{g_2}{\varepsilon_n}\right)
=(1/3,1/3,1/3).
\]

Thus the two configuration sequences approach **the same point** \((t=0;a=b=c=1/3)\) on the elementary (and the unnested top-level) parameter collision-link chart.

But the normalized square residuals differ. On the notch-square sequence,

\[
\frac{F_+(P_{n,0},\ldots,P_{n,3})}{\varepsilon_n}=(0,0,0,0).
\]

On the flat control, with \(e=(s_n,0)\), direct substitution yields \(F_+(Q_n)=(s_n,-s_n,3s_n,-s_n)\). Therefore

\[
\boxed{\frac{F_+(Q_n)}{\varepsilon_n}=
\left(\frac L3,-\frac L3,L,-\frac L3\right)
=\left(\frac{37}{18},-\frac{37}{18},\frac{37}{6},-\frac{37}{18}\right)\ne0.}
\]

These are incompatible limits of *the same explicitly defined function* \(F_+/\varepsilon\) at the *same parameter blow-up point*. Hence that function has **no continuous extension there**. In particular, no universally positive collision-link lower bound is possible for that normalization and this class of Jordan curves.

**Scope of the conclusion.** The argument targets a compactification built from collision positions and relative **parameter** gaps alone. It does **not** prove that every conceivable ambient-shape-enriched compactification lacks a continuous residual extension: an ambient shape blow-up can distinguish the flat control from the notch square. But it must then permit an **actual zero of the normalized square section on the boundary**, since the square-shape sequence has residual zero at every scale. That is the new relative-boundary obstruction.

## 4. Compactness: what is sufficient, and what isn't

A basic compactness fact remains valid. If continuous embeddings \(\gamma_m:S^1\to\mathbb R^2\) converge **uniformly** to an embedding \(\gamma\), and each \(\gamma_m\) contains a square with side \(\ge\delta>0\) (one fixed positive \(\delta\)), then after taking a subsequence the four parameter values converge. Their images form a nondegenerate square on \(\gamma\), because square equations are closed and the side-length lower bound survives the limit.

The infinite-notch example shows that **merely choosing any inscribed square** on each approximant supplies no such \(\delta\), even if the *limit curve itself is rectifiable*. The selected square can lie in the smallest notch of \(\Gamma_m\), with side \(s_m\to0\). It does not imply that all squares on \(\Gamma_m\) shrink, nor does it contradict the fact that the limiting curve has other squares.

Therefore a relative parity or Floer approach must supply a **macroscopic protected witness**, or use a homological/spectral mechanism that cannot lose *all* nontrivial information through the collision stratum. Merely recording the normalized parameter gaps is inadequate.

## 5. Consequences for the research program

| Proposed step | v220 conclusion |
|---|---|
| Assume every *fixed rectifiable* Jordan curve has a uniform square-size floor | **False**, by the infinite-notch curve. |
| Extend v219's C¹ collision-link section using just rectifiability and gap ratios | **False** for the specific normalized residual, by the two-sequence calculation. |
| Declare a source-parameter cyclohedral collision boundary universally zero-free | **Unsupported**, and false for any shape-retaining residual that preserves these small actual-square zeros. |
| Retain extrinsic normalized shape information | **Potentially useful**, but boundary zeros still occur and must be treated. |
| Infer a square on a C⁰ limit from approximating squares whose sides have a common positive floor | **Valid**, by compactness. |
| Infer Toeplitz merely from *existence* of squares on every approximating polygon | **Not justified**; selected witnesses may collapse. |
| Conclude that the general Toeplitz conjecture is false | **No.** The constructed limit contains infinitely many explicit squares. |

**Recommended v221:** Treat the collision-zero set as a **boundary stratum that can carry incoming/outgoing mod-2 incidence**, rather than attempting to force it away. Build a finite-dimensional collar model around the collision boundary, determine which incidence counts are intrinsic and which depend on scale cutoffs, and compare with known cyclohedral configuration-space approaches. Demand a relative cycle/continuation law that can tolerate the explicit infinite-notch case and the v210 rotating-scaling packet. A degree or parity computed without such a law has no universal limiting content.

## 6. Reproduction and sources

```bash
python -m unittest discover -s tests -v
python code/accumulating_notches_v220.py
```

The implementation uses only rational arithmetic and a direct geometric simplicity check for each finite polygon. The **infinite**-curve Jordan and rectifiability arguments are proved above, not inferred from finitely many tests. Tests also check normalized arclength-gap identities and the two different residual limits at arbitrarily large indices.

- v219 starting point: [equivariant quotient and C¹ collision link](v219-equivariant-quotient-collision-link.md).
- Siniša Vrećica and Rade Živaljević, *Fulton–MacPherson compactification, cyclohedra, and the polygonal pegs problem*, [arXiv:0810.1439](https://arxiv.org/abs/0810.1439). The cyclohedral framework is a precedent, not a theorem asserting our desired C⁰ collision regularity.
- Jason Cantarella, Elizabeth Denne and John McCleary, *Transversality for Configuration Spaces and the “Square-Peg” Theorem*, [arXiv:1402.6174](https://arxiv.org/abs/1402.6174), and *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, [arXiv:2103.07506](https://arxiv.org/abs/2103.07506). Their transversality statements include boundary-disjointness hypotheses not supplied by a general C⁰ curve.

**Originality caution:** Summable geometric notches and their small squares are elementary constructions; no claim is made that the construction is the first of its kind in the literature. The contribution to this research program is the **exact incompatibility of two collision-link residual limits on one fixed rectifiable curve**, and the resulting falsification of an overly strong proposed extension assumption.