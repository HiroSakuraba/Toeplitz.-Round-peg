# v265: Generic exact gluing of two positive squares at an irrational polygon vertex

**Research program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** A new reusable **sufficient certificate** for a single, nonsingular, irrational **vertex-induced two-square annihilation**, with an exact relative-degree-zero kink atlas and an optional completion of the v264 whole-polygon event ledger. It reproduces the v254/v263 complete 3→2→1 example without hardcoding the polygon's event polynomial into its verifier. **The unrestricted Toeplitz conjecture remains open.**

## 1. The precise improvement over v264

[v264](v264-reusable-convex-square-event-atlas.md) produced an exhaustive event-free chamber census for a rational-affine motion of a strictly convex polygon. Its principal missing case was **the orbit count at an irrational event time**. Samples on either side could not determine what happens *at* the singularity; the v254 moving pentagon, with three orbits before and one afterward, has **two** orbits exactly at the collision instant.

The new [v265 implementation](../code/square_irrational_vertex_gluing_v265.py) supplies a **generic event resolver**, under specifically verified hypotheses, rather than a polygon-specific collision lemma.

Consider rational endpoint convex polygons \(P_0,P_1\), moving affinely as

\[
P_i(t)=(1-t)(P_0)_i+t(P_1)_i,\qquad t\in[0,1].
\]

Let \(A=(a_0,a_1,a_2,a_3)\) and \(B=(b_0,b_1,b_2,b_3)\) be increasing supporting-edge assignments which differ in exactly one corner \(k\), with

\[
a_k=j-1,\qquad b_k=j,\qquad a_i=b_i \ (i\ne k).
\]

The \(k\)-th square corner crosses the **actual polygon vertex** \(P_j(t)\) from edge \(j-1\) to edge \(j\). The other three corners stay on the same three polygon edges.

The supported case has four distinct supporting edges on each side, a fixed nonzero Cramer determinant for both charts, a positive square side bound, and a single irrational contact time. A determinant-singular event, tangency, wraparound vertex, simultaneous multi-event or nonconvex/nongraceful configuration is deliberately not certified by this first solver.

## 2. Deriving the event polynomial, rather than supplying it

The two exact affine square-defect systems yield Cramer data

\[
D_A(t)=\det M_A(t),\quad u_{A,i}(t)=N_{A,i}(t)/D_A(t);
\qquad D_B(t),N_{B,i}(t).
\]

The *critical polygon-vertex contacts* are

\[
u_{A,k}(t)=1,\qquad u_{B,k}(t)=0.
\]

The verifier first certifies the two determinants have constant opposite signs \(s_A,s_B\in\{\pm1\}\) at every real time. It then constructs the signed contact polynomials

\[
H_A=s_A(D_A-N_{A,k}),\qquad
H_B=s_B N_{B,k}.
\]

It computes their **exact monic greatest common polynomial divisor**

\[
\boxed{Q(t)=\gcd(H_A(t),H_B(t))}
\]

using rational polynomial arithmetic. There is **no hardcoded v254 quadratic in this step**.

The built-in exact Sturm procedure proves there is precisely **one irrational root** \(\alpha\in(0,1)\) of \(Q\), supplied with a strict rational isolating bracket. Further exact division establishes

\[
H_A=QG_A,\qquad H_B=QG_B,
\]

where both quotient factors \(G_A,G_B\) have **strict constant nonzero sign** throughout the complete interval. Exact endpoint signs additionally establish

\[
H_A(0)>0,\ H_B(0)>0,\qquad
H_A(1)<0,\ H_B(1)<0.
\]

Since the root of \(Q\) is unique, the branch occupancy is fully determined: **two genuine roots before the event**, both reach the vertex at \(\alpha\), and **neither is on its assigned polygon edge afterward**.

This establishes a *death*, not a birth or an unproved tangency. Reversing the time direction fails the death-specific checker.

## 3. Five exact algebraic identities prove a single geometric square at the singular instant

It is insufficient to prove that the two roots have their critical corners at a common polygon vertex. They could describe **different** squares sharing that vertex.

For every one of the three other corner coordinates, v265 computes

\[
R_i(t)=N_{A,i}(t)D_B(t)-N_{B,i}(t)D_A(t)
\qquad(i\ne k)
\]

and requires that **the complete polynomial \(Q\) divides \(R_i\) exactly**. These are three exact rational polynomial congruences. The two critical coordinates are checked separately by the exact divisibilities

\[
Q\mid(D_A-N_{A,k}),\qquad Q\mid N_{B,k}.
\]

Together, these are the **five exact divisibility obligations** proving that, at the irrational root \(\alpha\),

- the other three parameterized corners agree exactly, including their actual Euclidean locations;
- the fourth lies at exactly \(P_j(\alpha)\) on both adjacent edges;
- both Cramer determinants are nonzero, so each chart has a unique extended root;
- both square-side cutoff polynomials remain strictly positive, so the shared square has **positive side length**.

This proves the merging square is **one geometric orbit** at the singular time, not two differently labeled roots and not a collapsed zero-size square.

## 4. A continuous 4D kink chart and an all-time external boundary certificate

Write the parameter at the crossing corner as

\[
h=
\begin{cases}
u_{A,k}-1,&h\le0 \quad\text{(previous edge)},\\
u_{B,k},&h\ge0 \quad\text{(next edge)}.
\end{cases}
\]

The boundary point is exactly \(P_j(t)\) at \(h=0\) from either side, so the **actual square-defect map** is continuous, although its derivative generally jumps there.

Use one glued four-dimensional rectangular region:

\[
U=\{\,0<u_i<1\ (i\ne k),\quad-1<h<1\,\}.
\]

For each of the two affine halfcells, the verifier proves by **full-time exact rational Bernstein inequalities** that its *extended* Cramer root remains **strictly inside the corresponding outer box**, even after the root leaves the physically admissible halfcell. This matters because the outer-boundary proof must apply to *both affine extensions* throughout the entire interval, not just the currently actual squares.

For each determinant, the code gives an exact positive Bernstein lower bound \(d\) on \(|D|\). Because every matrix coefficient is affine in \(t\), its absolute values are bounded by their rational endpoint maxima \(B\). The adjugate bound gives

\[
\|M(t)^{-1}\|_\infty\le\frac{24B^3}{d}.
\]

The sixteen strict root-to-face inequalities (eight per chart) give a rational source-coordinate margin \(\rho>0\). It follows that every outer-boundary configuration satisfies a **uniform strictly positive square-defect bound**

\[
\boxed{
\|F_t(u)\|_\infty
\ge \min_{\text{two halfcells}}
\left(\frac{\rho\,d}{24B^3}\right)>0.
}
\]

The continuous glued map is therefore free of zeros on the **external** box boundary for *every real time*, including \(\alpha\).

The local orientations are \(s_A=-s_B\). Before the event the two positive square roots have opposite local Brouwer indices and their sum is **zero**. Degree homotopy invariance and outer-boundary nonvanishing then imply

\[
\boxed{\deg(F_t,U,0)=0\qquad\forall t\in[0,1],}
\]

including the singular time where the individual halfcell degrees are undefined because their common internal interface contains a zero.

This is an exact mathematical gluing certificate, **not** a numerical count at an approximating rational time. Degree **zero** does not force an inscribed square by itself.

The strictly positive residual floor also supplies a quantitative budget for sufficiently small continuous perturbations of the source-arc images to remain boundary zero-free. The perturbations must refer to a **single consistent glued family**; no claim is made that the selected pair's zero-degree guarantees a surviving square after perturbation.

## 5. Completing a whole-polygon irrational event when the global atlas permits it

The function `complete_v264_vertex_event` takes the exact same moving polygon, the two adjacent-edge assignments and the contact corner index.

It first runs **the complete generic v264 polygon atlas**, rather than trusting a preselected branch inventory. It requires all of the following:

1. The polygon stays strictly convex throughout the interval, and the v264 checker supplies a strictly positive time-uniform **geometric** square-size floor.
2. **No unresolved chart** occurs in the global weakly increasing support-edge census.
3. There is exactly **one** candidate event cluster, and its complete list of factors consists of precisely the two proposed adjacent-edge contact factors, with no other determinant, edge membership or cutoff event.
4. The event contact root lies in that global event cluster, and the two v264 independent rational chamber counts differ by **exactly two**.
5. All *other* positive square branches are individually certified **strict interior stable charts**, in bijection with the geometric square orbits after the pair death. The checker now **rejects permanently vertex-locked third-party branches**: they could coincide with the singular square without generating another time-wall factor, so treating them as distinct is unsafe.
6. The oriented Jacobian index sums of the two regular chambers agree, as the disappearing pair contributes zero signed index.

With these **extra global checks**, all nonparticipating squares persist as distinct strict-interior roots, while the dying pair contributes two geometric orbits before the event, **one** at the event and none after. If the pre-event count is \(N\), then

\[
\boxed{N_{\rm before}=N,\quad N_{\rm at}=N-1,\quad N_{\rm after}=N-2.}
\]

This is a **complete whole-polygon three-slice result** on this restricted single-event class. It is not inferred from merely knowing the chamber counts, and the checker fails rather than guess when any premise is absent.

## 6. Reproduction on the genuine v254 moving pentagon

The verifier receives the v254 pentagon's **two rational endpoint polygons**, the two chart edge assignments

\[
A=(0,1,3,4),\qquad B=(0,2,3,4),
\]

and the crossing corner \(k=1\). It does **not** receive the collision polynomial.

In normalized time \(\tau\in[0,1]\), it independently derives the exact monic quadratic

\[
\boxed{
Q(\tau)=
\tau^2-\frac{1007331884387}{1811702}\tau
+\frac{61413740234}{905851}.
}
\]

Its unique irrational event is isolated in

\[
\boxed{\frac{255713}{2097152}<\tau<
\frac{2045705}{16777216}.}
\]

This event polynomial is exactly proportional to the v254 polynomial under the affine change of variables \(t=a+(b-a)\tau\). The regression suite checks the divisibility both ways, not just a floating-point approximation.

The two regular local indices are precisely

\[
\boxed{-1,\quad+1}.
\]

The code proves an exact whole-box root-to-exterior margin and an explicit outer square-defect lower bound:

\[
\rho\ \ge\ \frac{1547989225020375}{18880292598431446}>0,
\]

\[
\boxed{\|F\|_{\infty,\partial U}>
\frac{22266996110283865829996360231}
{204511802668101322500268416384000}.}
\]

The reusable v264 atlas independently recovers the five geometric categories: 25 three-collinear exclusions, five adjacent-edge exclusions, 37 strict off-edge exclusions, one stable interior square and two candidate event charts. The new resolver closes its formerly unknown singular slice, yielding

| Time relative to the unique irrational event | **Whole-polygon** positive square orbits | Paired-region Brouwer degree |
|---|---:|---:|
| Before | **3** | 0 |
| At event | **2** | 0 |
| After | **1** | 0 |

The persistent third square supplies the nonzero quotient-relative mod-2 degree of v262; the disappearing pair itself has degree zero. As in v263, there are **no hidden square events** in this fixed interval.

## 7. New tests and negative controls

The [new regression suite](../tests/test_square_irrational_vertex_gluing_v265.py) exercises exact equivalence to the v254 event polynomial, all five divisibilities, exact irrational root isolation, opposite Jacobian signs, full-time outer-boundary nonvanishing, the complete 3→2→1 global count, positive-size inheritance, and affine equivariance under **nonconstant translation** and **positive scaling**.

The fail-closed tests include: misaligned chart assignments, wrong crossing corner, reverse-time square birth instead of death, nonconvex input, missing common event, uncertified cutoff, perturbed polygon, a polynomial with a root falsely presented as strictly positive, and an adversarial artificial **third-party vertex-locked square branch** injected into the complete atlas.

The dedicated CI additionally replays the v264, v263, and v260 regression suites.

```sh
python code/square_irrational_vertex_gluing_v265.py
python -m unittest -v tests.test_square_irrational_vertex_gluing_v265
python -m unittest -v tests.test_convex_square_event_atlas_v264 tests.test_pentagon_complete_ledger_v263 tests.test_square_relative_degree_v260
```

## 8. Limits and the next logical step

This milestone solves an important *specific* verification problem: turning some irrational two-branch polygon-vertex events into complete, exact, **geometric orbit counts with certified local degree**, using a generic proof algorithm. It is applicable to other strictly convex rational-affine polygon families when they satisfy the same inequalities; a new example with a different algebraic event should be developed as an independent demonstration.

**It does not resolve**:

- event times with three or more interacting square orbits;
- determinant-rank singularities or true positive-dimensional square families;
- tangential or multi-parameter catastrophe events;
- nongraceful nonconvex square-configuration sectors;
- zero-size escape under arbitrary \(C^0\) approximations of Jordan curves.

The research path now has two high-value branches: (a) broaden the generic algebraic singularity resolver beyond simple two-chart vertex pair deaths, and (b) test whether any **mesh-independent** relative quotient degree/noncollapse argument can survive genuinely different approximating Jordan curves. The second is necessary for an unrestricted Toeplitz proof.

**The unrestricted square-peg conjecture remains open.**

## References

- [v218: cyclic configuration symmetry and orientation obstruction](v218-global-configuration-symmetry-and-collision-audit.md).
- [v223: complete geometric-square census on rational polygons](v223-complete-polygon-orbit-census-and-homotopy-audit.md).
- [v224: exact moving-polygon candidate event isolation](v224-algebraic-moving-polygon-event-covers.md).
- [v254: genuine square pair-death at a pentagon vertex](v254-genuine-square-pair-death.md).
- [v260: exact degree-zero gluing through the singular event](v260-relative-degree-pair-death.md).
- [v262: cyclic quotient-relative parity through the pair death](v262-cyclic-quotient-relative-parity.md).
- [v263: independent complete whole-pentagon all-real-time census](v263-complete-all-real-time-moving-pentagon-census.md).
- [v264: reusable exact convex-polygon event atlas](v264-reusable-convex-square-event-atlas.md).
