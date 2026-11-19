# v221 — Positive-size-filtered square incidence: exact local wall crossing and an accumulating-cutoff stress test

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Scientific status:** A conditional, elementary mod-2 cobordism identity; an exact-rational single-branch solver; and a fixed-curve tagged-witness negative control. **Not an unrestricted Toeplitz proof, not a homological construction on all Jordan curves, and not a global square-orbit parity invariant.**

## Executive result

The v220 infinite-notch Jordan curve demonstrates that genuine square zeros can accumulate at the collision boundary even for one fixed rectifiable curve. Consequently, an unmodified zero-free collision-link argument cannot be used. An alternative is to introduce a **positive squared-side cutoff** \(\kappa>0\), then record how square counts change when a witness crosses either the size cutoff or a marked source-arc boundary.

This batch obtains a correct **local incidence law**. Under explicit regularity and transversality hypotheses, the mod-2 number of filtered square solutions on the two ends of a one-parameter deformation differs by the mod-2 number of size-wall and source-wall crossings. Exact rational computation validates the identity in models with zero, one, and two cutoff crossings, an irrational cutoff time, and a source-boundary exit. We also show that the *explicitly designated* notch-corner squares on v220's **one fixed curve** do not have a cutoff-stabilizing tagged parity as \(\kappa\downarrow0\).

This is useful **bookkeeping**, not a guarantee of a nonzero end count. The missing global problem is to produce a protected nontrivial class whose boundary flux is controlled uniformly as the cutoff is removed.

## 1. A conditional mod-2 incidence theorem

Let \(A_i(t,u_i)\in\mathbb R^2\), \(i=0,1,2,3\), be continuous in \(t\in[t_-,t_+]\), differentiable as required for transversality, with \(u_i\in[0,1]\) labeling **four specified physical source arcs**. If working in a periodic covering, their physical deck lifts are fixed in advance; they are not silently summed over. Fix a cyclic square orientation \(\sigma\in\{+1,-1\}\) and the residual

\[
F_\sigma(t,u)=\bigl(A_2-A_1-\sigma J(A_1-A_0),\;A_3-A_0-\sigma J(A_1-A_0)\bigr)\in\mathbb R^4,
\]

where \(J(x,y)=(-y,x)\). Define squared side \(q(t,u)=|A_1-A_0|^2\) and the **filtered locus**

\[
\mathcal Z_\kappa=\{(t,u): F_\sigma(t,u)=0,\ u\in(0,1)^4,\ q(t,u)>\kappa\},\qquad \kappa>0.
\]

**Proposition (relative mod-2 incidence).** Assume that the closure of the filtered zero set in the compact prism \([t_-,t_+]\times[0,1]^4\) is a compact, one-dimensional manifold with boundary or a finite disjoint union thereof; the relevant zero locus is transverse to the two time slices, each external source face \(u_i=0,1\), and the size wall \(q=\kappa\); and that no simultaneous or corner boundary events occur. Then, writing all counts modulo two,

\[
\boxed{N_\kappa(t_-)+N_\kappa(t_+)
\equiv E_{\mathrm{source}}+E_{\mathrm{size}}\pmod2.}
\]

Here \(N_\kappa(t)\) counts solutions **with the indicated source labels, cyclic orientation, and physical lifts**, while \(E_{\mathrm{source}}\) and \(E_{\mathrm{size}}\) count transverse crossings of the two kinds of exposed boundary. The formula does **not** automatically count geometric square **orbits** in the full configuration space. It is a standard compact one-manifold boundary identity, not a new global topological obstruction.

**Proof.** A compact one-manifold is a finite union of circles and closed intervals. Each interval has two endpoints, each circle has none. In this filtered prism, every endpoint of the zero locus lies on exactly one of its external faces: initial time, terminal time, source boundary, or squared-side cutoff. Summing their cardinalities gives zero modulo two, which rearranges to the displayed formula. Because \(q>\kappa>0\), a solution in the filtered locus is a **nondegenerate square** with side \(>\sqrt\kappa\), so total-collision zeros are excluded for that fixed cutoff. This does not yield a cutoff-independent count as \(\kappa\to0\). \(\square\)

**Important precision.** The theorem assumes a suitable compact, regular zero locus. For an arbitrary merely continuous Jordan embedding, those hypotheses have *not* been proved. At a nongeneric isolated tangency or simultaneous boundary crossing, the simple transverse event counting is inapplicable; such events require perturbations and a more general intersection-chain formalism.

## 2. Exact affine-translation specialization: a complete one-branch algorithm

The implementation covers the substantial but finite-dimensional case

\[
A_i(t,u_i)=a_i+d_i u_i+v_i t,\qquad 0\le u_i\le1,
\]

where all vectors have **rational** coordinates and \(d_i\ne0\) are fixed segment directions. Translating an entire segment by \(v_i t\) does not change its direction.

The four scalar square equations are

\[
M u=b_0+t b_1.
\]

If the rational \(4\times4\) matrix \(M\) is invertible, the unique algebraic root branch is

\[
u(t)=u_0+t u_1,\quad u_0=M^{-1}b_0,\quad u_1=M^{-1}b_1.
\]

In particular, both the square edge \(e(t)=e_0+t e_1\) and its squared length

\[
q(t)=|e_1|^2t^2+2\langle e_0,e_1\rangle t+|e_0|^2
\]

are explicit rational polynomials. A branch is *physically admissible* precisely when \(0<u_i(t)<1\) for every \(i\). The implementation enumerates every potential source-face crossing \(u_i(t)=0\) or 1 **exactly** as a rational time. Between those times the source-interior condition is constant.

On each admissible open time interval, it counts the **simple real roots** of \(q(t)-\kappa\) from rational endpoint signs and, for a quadratic, the exact rational critical time \(t_*=-B/(2A)\). The actual crossing time can be irrational; no floating-point root finding is used. A quadratic with two simple roots contributes two size-wall events. A double tangency is rejected rather than misreported as a transverse crossing. Simultaneous source/cutoff crossings, endpoint events, and a singular \(M\) are also explicitly rejected.

The branch's two end indicators \(N_\kappa(t_\pm)\in\{0,1\}\) and the wall counts satisfy the proposition in every certified case. The checker proves *only this one-branch result*; a multi-cell PL homotopy may have additional simultaneous square branches, births, and knot-interface events.

## 3. Exact test models

### A. One positive-size cutoff crossing

Choose four short transverse segments centered at the moving vertices of the oriented square

\[
p_0=(1/2,0),\quad e(t)=(t/8,t/8),\quad
p_1=p_0+e,\ p_2=p_0+e+Je,\ p_3=p_0+Je.
\]

Their rational directions are \((1,-2),(1,3),(1,7/2),(1,-1/3)\), and each segment has half-width \(1/1024\) in its chosen direction coordinate. They therefore carry **one transverse labeled square** at parameters \(u_i=1/2\) for every \(t\). Let \(t_-=1/4\), \(t_+=3/4\) and \(\kappa=1/128\). Then

\[
q(t)=t^2/32,\qquad q(t)=\kappa\iff t=1/2
\]

on the specified positive time interval. The exact ledger is

\[
(N_\kappa(t_-),N_\kappa(t_+))=(0,1),\quad
(E_{\rm source},E_{\rm size})=(0,1).
\]

Thus the size filter alone changes the parity; **a parity that ignores this wall is not homotopy invariant**.

### B. Two cutoff crossings, still even total

Take \(e(t)=(1/4-t/2,0)\), \(t\in[1/4,3/4]\), and \(\kappa=1/256\). The two threshold times are \(t=3/8,5/8\), separated by a zero-sized square at \(t=1/2\). The filtered branch has a witness at both endpoint slices, no source boundary crossing, and **two** positive-size cutoff crossings. Consequently the parity change is zero, consistent with the formula. The zero-size middle square is excluded by \(\kappa>0\).

### C. An irrational cutoff time with exact counting

Take \(e(t)=(1/8,t/8)\), \(t\in[0,1]\), and \(\kappa=3/128\). Then

\[
q(t)-\kappa=(t^2-1/2)/64,
\]

with **one positive crossing** at \(t=1/\sqrt2\). The time is irrational, but exact rational sign and vertex tests return one crossing and endpoint counts \((0,1)\). Nothing is inferred from a floating-point sample grid.

### D. Source-window escape without a cutoff crossing

Keep an initial square of edge \(e=(1/8,1/8)\) and translate **source segment 1** vertically at velocity \((0,1)\), holding the other three segments fixed. Their half-width remains \(1/1024\), so the unique solution parameters move affinely. With \(\kappa=1/1024\) and \(t\in[0,1/5]\), the branch crosses its **source-parameter face** at exactly

\[
\boxed{t=55/19968}. 
\]

At this instant \(u_1=0\), with the other three parameters strictly interior and \(q(t)>\kappa\). There are no size-cutoff crossings in the active source chamber. The ledger gives

\[
(N_\kappa(0),N_\kappa(1/5))=(1,0),\quad
(E_{\rm source},E_{\rm size})=(1,0).
\]

Reversing time exchanges exit with entry and preserves the modulo-two boundary identity. Rational rotations and translations preserve the filtered ledger.

### E. Nontransverse controls

A positive-size **tangential** contact with the cutoff, a root exactly at the initial time, an identically-on-cutoff branch, or a singular residual matrix is flagged as **non-generic**. They are not silently counted as transverse parity events. An auxiliary quadratic routine also checks pairs of irrational cutoff roots by exact rational sign logic.

## 4. V220 fixed-curve stress test: *tagged* cutoff parity does not stabilize

Recall v220's one fixed rectifiable infinite-notch Jordan curve. Its **designated, explicit notch-corner squares** have squared side lengths

\[
\kappa_n=s_n^2=16^{-(n+1)},\qquad n\ge1.
\]

At any cutoff \(\kappa>0\), only finitely many of these **named** squares satisfy \(s_n^2\ge\kappa\). Define the *tagged* count (which is **not** the total square-orbit count on the curve)

\[
T(\kappa)=\#\{n\ge1:s_n^2\ge\kappa\}.
\]

For \(n\ge1\), set \(\kappa_n=16^{-(n+1)}\) and \(\eta_n=(\kappa_{n-1}+\kappa_n)/2\), using the formal threshold \(\kappa_0=1/16\) above the largest named notch. Then

\[
\boxed{T(\kappa_n)=n,\qquad T(\eta_n)=n-1.}
\]

In particular, **the parity of these tagged witnesses flips infinitely often as the positive cutoff is lowered toward zero**. The sequence of marked cutoff events accumulates at the collision boundary, even though the curve has finite length. This is consistent with and strengthens the *stress test* for a proposed relative collision incidence law. But it is crucial not to overinterpret it:

- It counts **one deliberately specified square per notch**, not all inscribed squares or all C4 orbit classes. There may be additional squares (including large ones) on the same curve.
- The exact notch corners are not claimed to be isolated or transverse square solutions. The tagged parity is **not** presented as a Brouwer degree, mod-2 intersection number, or actual global invariant.
- Therefore the alternating tagged parity **does not disprove** a possible global, symmetry-aware parity with compensating contributions from other square families. It shows only that any proof purporting to protect individual square witnesses without handling cutoff flux fails this explicit test.

The code checks \(T(\kappa_n)\) and \(T(\eta_n)\) at dozens of scales using exact fractions, and CI cross-checks the notch corner identities against the independent v220 module.

## 5. What is new for this program, and what is still missing

| Mathematical/procedural statement | Status |
|---|---|
| Parity of endpoint solutions plus all transverse boundary events is zero in a compact one-manifold | **Standard topological identity, derived conditionally** |
| Exact source-wall and squared-side-cutoff ledger for a full-rank affine-translation four-segment chart | **Implemented with rational arithmetic** |
| Exact counting of irrational quadratic cutoff times by rational signs | **Implemented; no floating-point oracle** |
| One and two size crossings, one source exit, and tangency refusals | **Exact worked models and regressions** |
| Tagged notch-corner parity stabilizes as cutoff tends to zero | **False** in the v220 fixed-curve model |
| The tagged parity equals the total geometric square-orbit parity | **Not claimed and generally unsupported** |
| A compact regular cobordism of all square-orbits for arbitrary C0 Jordan curves exists | **Not established** |
| The cutoff flux can be controlled or paired to preserve an odd global orbit count as \(\kappa\to0\) | **Open** |
| This establishes the unrestricted square-peg conjecture | **No** |

## 6. Next serious milestone (v222)

Instead of adding more local root enumerators, formulate a **source-compatible, symmetry-aware relative chain complex on a nonzero-size region**. Its generators should be **geometric square orbits**, not four cyclic labelings counted separately. Its boundary must include source-arc escape, cutoff escape, and (for PL homotopies) internal knot/bifurcation interfaces, with a certified cancellation law. The key *global* research test is whether a nontrivial class supplies a size scale bounded away from zero under arbitrary C0 approximations. The v220 accumulating notches, v218 sign cancellation, v213 wrong-deck examples, and v210 spiral obstruction remain mandatory negative controls.

A possible alternative is to recover positive scale from a symplectic/spectral quantity, but its physical-size lower bound at the **square angle** is not supplied by the cited periodic inscription theorems. A mere positive-measure family of rectangle angles cannot be substituted for a guarantee at exactly 90 degrees.

## 7. Reproduction and references

```
python -m unittest discover -s tests -v
python code/cutoff_incidence_v221.py
```

The v221 source and tests use only the Python standard library (`fractions.Fraction`, no Sage or Arb). The full repository workflow also runs v212–v220 and checks archived ZIP hashes. Local sampling-based tests supplement, rather than replace, the exact quadratic sign proof.

**Prior art, not claimed as a new global theorem:**
- Cantarella, Denne & McCleary, *Transversality for Configuration Spaces and the "Square-Peg" Theorem*, arXiv:1402.6174, https://arxiv.org/abs/1402.6174. Their boundary-disjointness/transversality hypotheses are not automatic for rough Jordan curves.
- Cantarella, Denne & McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, arXiv:2103.07506, https://arxiv.org/abs/2103.07506. Generic smooth odd-square results do not give unrestricted C0 passage to a noncollapsed limit.
- Matschke, *On the Square Peg Problem and Its Relatives*, arXiv:1001.0186, https://arxiv.org/abs/1001.0186. Prior mod-2 intersection and size-controlled criteria.
- Internal prerequisites: [v220 infinite notches](v220-fixed-rectifiable-collision-link-obstruction.md), [v218 cyclic signs](v218-global-configuration-symmetry-and-collision-audit.md), [v217 subdivision chain](v217-multiscale-relative-degree-ledger.md), and [v213 physical-lift audit](v213-exact-PL-localization-and-lift-audit.md).
