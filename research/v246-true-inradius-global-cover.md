# v246: Global inradius bounds and a true improvement in the Big Peg ratio

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Completed reproducible **exact-rational global covering certificates** for the inradii of both the original Pochuev ten-gon (v244) and our modified v245 ten-gon. Together with exact complete maximum-square enumerations from v243, they prove that **v245 has a strictly smaller *actual* Big Peg ratio**, not just a smaller convenient-disk upper bound.

## Theorem

Let \(P_0\) be the original ten-vertex Pochuev polygon from v244 and \(P_1\) the modified rational ten-gon from v245. Let \(r_{\mathrm{in}}(P)\) be the largest radius of a disk contained in the bounded polygonal region, and \(s_{\max}(P)\) the largest side length of any square inscribed in its boundary. Define

\[
J(P)=\frac{s_{\max}(P)}{r_{\mathrm{in}}(P)}.
\]

**Certified result:**

\[
1.40898288313781 < J(P_0) < 1.40898289277737,
\]

\[
1.40869527908545 < J(P_1) < 1.40869528263098.
\]

Each asserted decimal interval is intentionally rounded *outward* for display. All comparisons in the verification program use exact rational squares; the decimal approximations are not proof inputs.

Consequently,

\[
\boxed{J(P_1)<J(P_0)}
\]

with a proven separation greater than 0.0002876005 (in the real, not squared, ratio). The new polygon truly improves upon Pochuev's original ratio. This resolves the comparison limitation explicitly flagged in v245.

## Exact inradius enclosures

The rational lower centres were found through numerical optimization but are then **verified**, not assumed. For each centre we compute the minimum exact squared distance to every closed polygon segment and check the centre belongs to the interior using exact ray crossing.

| Polygon | Rational interior centre \((x,y)\) | Exact squared lower radius | Exact upper radius |
|---|---|---|---|
| Pochuev original \(P_0\) | \((45637/10^8,\,10485/10^8)\) | \(246132835013175286089001/245641050000000000000000\) | \(100100053/10^8\) |
| v245 polygon \(P_1\) | \((302622629/10^{12},\,2143102/10^{13})\) | \(68382681809394993823195766186089/68233625000000000000000000000000\) | \(100109166/10^8\) |

The lower and upper radii differ by less than \(7\times 10^{-9}\) for the original and \(3\times 10^{-9}\) for v245.

The upper bounds are *global*, not local optimality conclusions based on three presumed active edges.

## How an exact box cover proves a GLOBAL upper bound

Let \(B\) be any closed rational rectangle containing the polygon. The implementation chooses the smallest axis-aligned bounding rectangle of all polygon vertices, which contains the entire bounded polygonal region.

For a subrectangle \(Q=[x_l,x_h]\times[y_l,y_h]\), let its rational midpoint be \(m\), and put

\[
H=\frac{x_h-x_l+y_h-y_l}{2}.
\]

For every \(x\in Q\), Euclidean distance \(\|x-m\|_2 \le H\), because Euclidean norm is bounded by the \(L^1\) norm.

Write \(d(m)\) for the exact Euclidean distance of \(m\) from the polygon's **entire closed boundary**, computed by taking the minimum over distances to *all* closed segments. Each squared distance and each comparison is rational.

The algorithm declares one of two terminal certifications:

**A. Whole box is close to the boundary.** If \(H\le R\) and the exact rational inequality
\[
d(m)^2 \le (R-H)^2
\]
holds, then the 1-Lipschitz property of distance to any nonempty closed set proves \(d(x)\le d(m)+H\le R\) for *every point* \(x\in Q\). No disk centered in this box can have radius greater than \(R\).

**B. Whole box is exterior.** If \(m\) is outside the open polygonal interior (checked exactly) and
\[
d(m)^2>H^2,
\]
then no point in the box can belong to that interior. If there were an interior point, the line segment joining it to the exterior midpoint would cross the polygon boundary at distance at most \(H\), contradicting the strict inequality.

If neither rule applies, split \(Q\) into four congruent rational boxes. The search starts with the global bounding rectangle and processes every child. Any unclassified box at depth 48, or resource-limit exhaustion, **throws an exception and refuses to certify**. There is no probabilistic completion step, unverified geometry assumption or floating-point branch choice.

### Recorded complete certificate statistics

| Polygon | Processed rational boxes | Maximum subdivision depth |
|---|---:|---:|
| Pochuev original | 885 | 29 |
| Modified v245 | 917 | 30 |

These are small, deterministic, executable global covers. Both box-completion identities and the exact predicates are regression-tested. The entire computation uses Python's standard library: fractions.Fraction and exact integer arithmetic.

## Exact separation using maximum inscribed squares

The independently certified **maximum squared sides** are retained from v244/v245:

\[
S_0^2=\frac{56381275521625791352241234309}
{28343589374868261752462760000},\quad
S_1^2=\frac{2263154862866037715637393248853633089}
{1137974336077532830058267809922880000}.
\]

Let \(L_0^2,L_1^2\) be the two rational lower radius-squared values in the table and \(U_0,U_1\) their rational upper radius bounds. The complete global covers establish

\[
\frac{S_i^2}{U_i^2}\le J(P_i)^2\le \frac{S_i^2}{L_i^2}.
\]

The key theorem is the **strict exact rational inequality**

\[
\boxed{\frac{S_1^2}{L_1^2}<\frac{S_0^2}{U_0^2}}.
\]

The program verifies it by Fraction comparison, not rounded decimals. Therefore every possible actual inradius consistent with the independent certificates satisfies \(J(P_1)<J(P_0)\).

The alternative hypotheses that the v245 result merely improved a suboptimal disk choice, or that a larger undiscovered square invalidated the ratio, are both eliminated: inradii are globally bounded and square maxima were exhausted by v243's edge-chart method.

## Reproduction and controls

Run from the repository root:

    python code/big_peg_inradius_v246.py
    python -m unittest -v tests.test_big_peg_inradius_v246

The v246 workflow also reruns v244 and v245's full exact maximum-square verification commands. The dedicated tests check the rational contact centres, exact squared boundary distances, both complete global covers, the exact J comparison, fail-closed depth and node limits, incorrect radii, and a square polygon benchmark.

**What is NOT proved:** Exact closed-form symbolic algebraic values of either true inradius, or optimality of v245 among all ten-gons, or existence/nonexistence of a positive universal Big Peg constant. The unrestricted Toeplitz conjecture remains open. These certificates give **narrow rigorous intervals** for both true inradii, strong enough to establish the actual strict ratio ordering.
