# v251: A size-aware atlas for multiple, competing square branches

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact rational *existential* size-cutoff certificate for every real parameter of a changing Jordan polygon, plus a unique verified crossover between **two distinct** inscribed square configurations. Does **not** prove a universal Big Peg constant, generic homotopy parity or Toeplitz.

## The structural obstacle after v250

The v250 exact construction showed that a square can survive a genuine non-collinear polygon corner, changing its four-edge chart without shrinking. However, its verifier still followed one continuous geometric square through the handoff.

That is stronger than necessary for the Big Peg strategy. To prove that some large square exists throughout a deformation, **the square serving as witness may change discontinuously**. It is the existence of at least one square above a size threshold that matters, not whether we can continuously follow one named square.

We therefore distinguish:
- **Continuous root atlas:** each edge chart provides one square branch over a parameter interval, and the two selected square roots agree at each seam. This was v249/v250.
- **Existential size atlas:** each edge chart provides *some* inscribed square with squared side strictly above a common rational cutoff throughout its time interval. Closed intervals cover all parameter times. At chart boundaries, the chosen squares may be **different**. This new v251 criterion is enough to show that the **maximum** square side exceeds the cutoff at every time, without tracking the maximum itself.

## General finite-atlas theorem

Let \(P(\lambda)\), \(\lambda\in[a,b]\), be a family of polygonal Jordan curves with affine rational vertex coordinates. Assume independent geometry checks establish the Jordan property for all times. For an ordered quadruple \(E=(e_0,e_1,e_2,e_3)\) of polygon edges, solve the square's four supporting-line equations

\[
M_E(\lambda) z=b_E(\lambda),\qquad
z=(p_x,p_y,u,v).
\]

The entries are rational polynomials in \(\lambda\). When \(D_E(\lambda)=\det M_E(\lambda)\) has constant nonzero sign, the unique square candidate is the continuous rational function \(z_E=N_E/D_E\), where \(N_E\) denotes the Cramer determinant vector.

An **exact chart certificate** proves for every real \(\lambda\) in a specified **closed rational interval** \(I\):
1. \(D_E(\lambda)\ne0\).
2. All four square corners lie on their **closed** chosen polygon edge segments, verified by exact parameter-sign polynomial inequalities.
3. Its squared side \(u^2+v^2\) is strictly greater than a given rational cutoff \(\eta>0\), verified by the rational polynomial inequality
\[
N_{E,u}^2+N_{E,v}^2-\eta D_E^2>0.
\]

These statements are each decided by the exact rational Bernstein coefficient/subdivision verifier developed in v249. Indeterminate signs fail closed.

**Existential size-atlas theorem:** If finitely many such chart intervals have union containing the full deformation interval \([a,b]\), then **every** \(P(\lambda)\) has an inscribed square with squared side **strictly greater than \(\eta\)**. It does **not** matter whether the chart roots coincide at transition times. The proof simply chooses a valid chart whose closed time interval contains the given real parameter value. The proof is **pointwise existential**; no unproven continuous choice is required.

The verifier checks every provided chart and the entire closed interval cover, rejecting uncovered open intervals, missing endpoints, invalid rational input, singular charts and false side cutoffs.

## A quantified competition between TWO different large squares

Reuse v250's rigorously verified, actually changing, **strictly convex non-collinear nine-gon family**, with a moving kink and positive nonadjacent-edge separation. v250 established the Jordan property and that every polygon remains inside a disk of radius three for all \(\lambda\in[0,1]\).

Focus on the short interval

\[
W=\left[\frac12,\frac{501}{1000}\right].
\]

The following two *different* inscribed square candidates both exist for **every** parameter in \(W\). Their ordered supporting-edge assignments are

\[
A=(0,3,5,7),\qquad B=(1,3,5,8).
\]

Exact Cramer/Bernstein checks prove the two squares lie on the polygon and have squared side greater than \(9\) throughout \(W\).

The key discovery is that their **relative size order changes**.

At \(\lambda=1/2\), the exact comparison \(s_A^2>s_B^2\) holds; at \(\lambda=501/1000\), the exact opposite strict inequality \(s_A^2<s_B^2\) holds. Both inequalities are checked by exact rational evaluation, not rounded calculations.

For illustration only:

| parameter | square \(A\) side² | square \(B\) side² |
|---|---:|---:|
| \(1/2\) | about 9.42595291 | about 9.42556944 |
| \(501/1000\) | about 9.42549764 | about 9.42581282 |

A continuous rational polynomial equal in sign to \(s_A^2-s_B^2\) on the interval is

\[
G=(N_{A,u}^2+N_{A,v}^2)D_B^2
 -(N_{B,u}^2+N_{B,v}^2)D_A^2.
\]

Because \(D_A^2D_B^2>0\), the sign of \(G\) is exactly the relative side-length ordering. We certify \(G(1/2)>0\), \(G(501/1000)<0\), and **\(G'<0\)** on the entire interval with rational Bernstein sign checks. Therefore there is **exactly one** time when these two distinct configurations have equal square size. It is isolated in the tight rational interval

\[
\frac{500549}{10^6}<\lambda_{\rm crossover}
<\frac{10011}{20000}.
\]

This is a **competition between two certified square candidates**, not a proven exchange of the **globally largest** square. There may be other square orbits at these times.

## A high-size guarantee that requires changing the selected branch

Define the strict squared-side threshold

\[
\eta=\frac{94257}{10000}=9.4257.
\]

Neither candidate \(A\) nor \(B\), individually, remains above \(\eta\) for the whole interval \(W\): by the right endpoint, \(s_A^2<\eta\), while at the left endpoint \(s_B^2<\eta\).

Nevertheless, the rational time

\[
\lambda_*=\frac{10011}{20000}=0.50055
\]

provides a successful **two-chart existence atlas**:

| Closed real time interval | Certified inscribed square | Exact strict side bound |
|---|---|---|
| \([1/2,10011/20000]\) | Square \(A\), edges \((0,3,5,7)\) | side² > 94257/10000 |
| \([10011/20000,501/1000]\) | Square \(B\), edges \((1,3,5,8)\) | side² > 94257/10000 |

At \(\lambda_*\), **both** squares are valid, both exceed the cutoff and their first square corners are **different**, so the selected witness can jump legitimately. Exact Bernstein tests establish the size bounds over every real instant, not just at endpoints or samples.

It follows that

\[
\boxed{\max_{\text{boundary-inscribed squares }Q\text{ on }P(\lambda)}
\operatorname{side}(Q)^2>\frac{94257}{10000}}
\qquad \forall\,\lambda\in W.
\]

This is a **certified lower bound on the actual maximum**, proved by two alternative candidate branches, neither of which alone would suffice to prove that particular high cutoff on all of \(W\).

## Whole-family demonstration without continuous witness selection

For the complete homotopy \(\lambda\in[0,1]\), a four-chart existential atlas uses:

1. \([0,1/2]\), edges \((1,3,5,7)\);
2. \([1/2,\lambda_*]\), edges \(A=(0,3,5,7)\);
3. \([\lambda_*,501/1000]\), edges \(B=(1,3,5,8)\);
4. \([501/1000,1]\), edges \(A=(0,3,5,7)\).

Each square has **side² > 9**, with no uncovered time, and the selected witness changes to a different square at two atlas seams. The polygon remains within radius three, so its true inradius is at most three. Consequently,

\[
\boxed{J(P(\lambda))=\frac{s_{\max}(P(\lambda))}{r_{\rm in}(P(\lambda))}>1}
\quad\text{for every }\lambda\in[0,1].
\]

This restricted-family result was also achievable by tracking one continuous square in v250. The **new conceptual gain** is that our certificate **does not require such tracking**, and can prove a stronger high-cutoff inequality on a subinterval by *combining competing squares*.

## What the theorem does not establish

- **Not the global maximum's crossover:** The two certified squares are not proven to exhaust every square orbit. The equal-size event is between these two particular candidates.
- **Not a genuine birth/death of square orbits:** Both candidates exist throughout the observation window. Their **relative sizes** cross, not their existence domains.
- **Not a homotopy-invariant parity or Floer class:** The existential finite cover is a sufficient certificate, not a conserved invariant for every polygon deformation.
- **Not a uniform Big Peg theorem:** The size threshold and chosen charts depend on the specific nine-gon family. There is no bound independent of edge count, configuration singularities or feature size.
- **Not a Toeplitz proof:** Arbitrary Jordan curves remain outside the established family.

Generic odd-parity arguments for smooth and generic polygon curves (e.g., Benjamin Matschke, *A Survey on the Square Peg Problem*, Notices AMS 2014, DOI 10.1090/noti1100) motivate a future size-filtered **parity** study, but generic parity alone does not prevent every square shrinking toward a rough curve as polygonal approximations are refined. We make no new assertion of a universally conserved cutoff parity in this PR.

## Reproduce

Pure Python standard library; all formal comparisons are exact Fractions:

~~~sh
python code/big_peg_size_atlas_v251.py
python -m unittest -v tests.test_big_peg_size_atlas_v251
~~~

The 17 new tests check full closed-time coverage, both strict high-cutoff charts, a discontinuous selected-witness handoff, reversed side ordering, exact isolated crossover, unique-event derivative sign, bad time gaps, unverified endpoints, wrong rational type and invalid individual square charts.

**Next:** Develop a complete **square-orbit event ledger** for a moving polygon that enumerates all branches and cutoff crossings, including rank-deficient incidence charts, rather than choosing a sufficient subset of square witnesses. This could support future size-aware parity studies, but eliminating all geometry-dependent scales still requires a new universal idea.
