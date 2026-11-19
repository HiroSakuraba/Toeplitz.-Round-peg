# v224 — Exact algebraic event covers for linearly moving Jordan polygons

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** Conditional, finite-dimensional polynomial/event-cover theorem; exact rational Sturm implementation; independently checked polygon example and adversarial controls. **Not a global square-orbit cobordism, Floer invariant, or proof of the unrestricted Toeplitz conjecture.**

## Executive finding

v223 could enumerate all geometric square orbits of a *fixed* rational Jordan polygon, including repeated edges, vertex hits, and positive-dimensional families, but its dynamic incidence ledger applied only to global positive scaling and translation. The next unresolved operation was a polygon whose **individual vertices move affinely**.

For every ordered selection of four moving polygon edges and a fixed square orientation, v224 obtains explicit **rational univariate polynomials** whose real zeros contain every possible change in that chart's nonsingular, interior, positive-size square status. Their roots can be **exhaustively isolated, including irrational and repeated roots, with exact rational Sturm sequences**. A separate simple general-position certificate proves that the polygon stays embedded throughout the deformation, rather than merely validating its endpoints.

**The crucial distinction:** This is an *event candidate cover*, **not** a classification or count of real geometric-orbit events. For charts where square parameters or determinants vanish identically, the implementation reports **unresolved**. Even an isolated polynomial root can be spurious when the corresponding square falls outside the chosen edge intervals. A global square-orbit cobordism must still merge equivalent chart presentations and analyze singular/persistent strata.

## 1. Exact moving-edge square equations

Let the polygon's rational vertices move by

\[
V_j(t)=(1-t)V_j^0+tV_j^1,\qquad 0\le t\le1.
\]

An ordered four-edge chart chooses edges \(i_0,i_1,i_2,i_3\), **allowing repeated edges**, and for each \(r=0,\ldots,3\) represents its square point as

\[
p_r(t,u_r)=V_{i_r}(t)+u_r\bigl(V_{i_r+1}(t)-V_{i_r}(t)\bigr),
\quad 0\le u_r\le1,
\]

with polygon edge indices cyclic. The chosen oriented-square residual is

\[
F_\sigma(t,u)=\bigl(p_2-p_1-\sigma J(p_1-p_0),\quad
p_3-p_0-\sigma J(p_1-p_0)\bigr)\in\mathbb R^4,
\qquad \sigma=\pm1,\quad J(x,y)=(-y,x).
\]

The four scalar equations are **linear in the edge parameters** \(u_0,\ldots,u_3\) and **affine in deformation time in each coefficient**:

\[
\boxed{A(t)u=b(t),\quad A(t)=A_0+tA_1,\quad b(t)=b_0+tb_1.}
\]

Set

\[
D(t)=\det A(t),\qquad N_r(t)=\det A_r(t),
\]

where \(A_r\) replaces column \(r\) of \(A\) with \(b\). The determinants have rational coefficients and degree **at most four**. Where \(D(t)\ne0\), Cramer's rule gives the **unique real branch** for that chart:

\[
\boxed{u_r(t)=N_r(t)/D(t).}
\]

Thus any transition from a strictly interior to an exterior edge parameter can occur only where

\[
N_r(t)=0\quad\text{or}\quad D(t)-N_r(t)=0,
\]

or where \(D(t)=0\). The sign of these expressions must be interpreted relative to the sign of \(D\). No numerical approximation or generic-rank guess is needed for the algebraic identities.

## 2. Exact positive-size cutoff polynomial

Write \(a_r(t)=V_{i_r}(t)\) and \(d_r(t)=V_{i_r+1}(t)-V_{i_r}(t)\). The first square edge is

\[
e(t)=p_1(t,u_1(t))-p_0(t,u_0(t))=E(t)/D(t),
\]

with the polynomial vector

\[
\boxed{E(t)=D(t)(a_1(t)-a_0(t))+N_1(t)d_1(t)-N_0(t)d_0(t).}
\]

Each coordinate of \(E\) has degree at most **five**. For any rational squared-side cutoff \(\kappa\ge0\), define

\[
\boxed{H_\kappa(t)=E_x(t)^2+E_y(t)^2-\kappa D(t)^2.}
\]

This polynomial has degree at most **ten**, and, wherever \(D\ne0\),

\[
H_\kappa(t)=D(t)^2\bigl(|e(t)|^2-\kappa\bigr).
\]

Therefore an interior square can cross the size cutoff only at a real zero of \(H_\kappa\). The squared denominator is positive away from singular times, so cutoff eligibility is decided by the **exact sign of \(H_\kappa\)**.

### Theorem A: complete **candidate** cover on a nonpersistent chart

Suppose none of the ten polynomials

\[
D;\quad N_r\ (r=0,1,2,3);\quad D-N_r\ (r=0,1,2,3);\quad H_\kappa
\]

is identically zero. Between consecutive real roots of this finite family, the signs of all ten polynomials are constant. Consequently the chart either contains **one** strictly interior square of squared side greater than \(\kappa\) throughout that open interval, or contains **none** throughout. Any change in that chart status is therefore covered by the real roots of the displayed polynomials.

**Proof.** On a root-free interval, \(D\ne0\), so the unique algebraic branch \(u=N/D\) varies continuously. The signs of \(N_r/D\) and \((D-N_r)/D\) establish whether \(0<u_r<1\), while the sign of \(H_\kappa\) establishes \(|e|^2>\kappa\). All signs are fixed throughout the interval. \(\square\)

**Limit:** If one of these polynomials is identically zero, its root set is not finite. This may represent permanent singularity, a root persisting at an edge endpoint, or a cutoff equality; one must analyze that chart separately. Returning an empty root list would be incorrect. v224 returns an explicit unresolved reason instead. Likewise a candidate real root need not be a *physical* event: its root may correspond to edge parameters outside \([0,1]\), or may be a duplicate labeling of the same square at a polygon vertex. The theorem does **not** assert that counting candidate roots yields mod-2 orbit incidence.

## 3. Exact real-root isolation by rational Sturm sequences

For any nonzero rational polynomial \(p\), calculate its squarefree part \(\widehat p=p/\gcd(p,p')\), then construct the rational Sturm chain by signed Euclidean remainders. Its sign-variation difference gives the exact number of **distinct real roots** in each rational open interval. In the implementation, values at root endpoints are handled using the correct one-sided convention; exact rational roots detected at dyadic subdivision points are returned as singleton brackets.

Recursively bisect a rational interval whenever it contains at least one undiscovered real root. Every isolated irrational root is returned in a rational **open isolating bracket** of requested width, with its exact defining polynomial retained; neither approximate floating-point time nor extracted radicals are needed. Multiple roots of \(p\) are included **once** via the squarefree part. The algorithm has an explicit subdivision budget and raises an inconclusive error if exhausted.

The reproducible regression suite covers rational, irrational, and even-multiplicity roots, endpoint roots, and polynomials with no real roots. An additional **development-only independent check** matched the counts from SymPy 1.14's exact rational `Poly.intervals` on **80 seeded nonconstant integer polynomials** of degrees one through seven. **SymPy is not a dependency** of the committed checker or CI.

## 4. Certifying the polygon stays Jordan

Checking that the initial and final polygons are individually simple does **not** prove the intermediate polygons are simple. The module includes a conservative, separate certificate:

- Validate a rational starting snapshot as a simple polygon.
- For every triple of distinct polygon vertex indices, form the signed orientation determinant

  \[
  O_{ijk}(t)=\det(V_j(t)-V_i(t),V_k(t)-V_i(t)).
  \]

  It is a rational polynomial of degree at most two.
- Prove exactly that **none** of these polynomials has a root anywhere on the closed requested time interval, using Sturm isolation.

**Lemma B (sufficient simple-isotopy certificate).** Under these conditions the polygon remains embedded over the whole interval. The first failure of simplicity would require a vertex entering a nonadjacent edge, adjacent backtracking or overlapping edges, or a vertex collision. At first contact, some three distinct polygon vertices must become collinear (vertex collision implies collinearity with any third vertex). Every such event contradicts the triple-orientation certificate. \(\square\)

The certificate is intentionally **not necessary**. A harmless collinear triple for an entire homotopy causes a conservative refusal. For a general moving polygon this condition must be checked or replaced by a more capable exact simplicity proof before regarding the candidate ledger as Jordan-curve evidence.

## 5. Exact vertex-alignment event on **one strictly convex polygon**

Use the v223 convex quadrilateral

\[
P(\tau)=\bigl((0,0),(2,0),(3-\tau,3-3\tau/2),(0,2-3\tau/2)\bigr),
\quad 0\le\tau<4/3.
\]

Move linearly from \(P(0)\) to \(P(9/8)\) with normalized new parameter \(t\in[0,1]\). Then \(\tau=(9/8)t\). The v223 corner alignment at \(\tau=1\) occurs at the **exact rational time**

\[
\boxed{t=8/9.}
\]

On the repeated-edge chart \((i_0,i_1,i_2,i_3)=(0,0,1,2)\) with \(\sigma=+1\), Cramer's rule gives

\[
u(8/9)=(1/2,1,2/3,1/2).
\]

Thus the **second square vertex reaches polygon vertex \(P_1=(2,0)\)**. Its actual physical square corners are

\[
\boxed{(1,0),(2,0),(2,1),(1,1),\qquad\text{squared side }1.}
\]

The event polynomial \(D-N_1\) vanishes *exactly* at \(8/9\); the full chart uses rational coordinates and has nonzero determinant there. The independent v223 exhaustive fixed-polygon census gives **one** geometric square at \(\tau=63/64\), **two** at \(\tau=1\), and **three** at \(\tau=33/32\). These are distinct certified snapshots, not an assertion that this single event explains every change in the geometric-orbit count. The triple-orientation certificate verifies that the polygon stays a Jordan quadrilateral throughout the normalized interval \([0,1]\).

**Why repeated edges matter:** A scan restricting every square point to four **different** polygon edges produces no candidate at \(t=8/9\) in this example. That omission is material; a vertex square can use two points on one edge and one point exactly at its endpoint. The full repeated-edge search covers the event, but in this example also reports *persistently-on-edge chart identities*, leaving the full global event chain unresolved. At rational candidate times, a separate function verifies the actual four square points and checks that all four parameters lie on their indicated edges and the positive-size condition holds. It never promotes a merely algebraic root to a geometric event without this check.

## 6. Independently isolated irrational cutoff time

Start with the rational quadrilateral

\[
Q=((0,0),(2,0),(3,3),(0,2)),
\]

whose v223 census gives one isolated geometric square of squared side \(q=72/25\). Apply the genuine embedded whole-polygon homotopy \(Q(t)=(1+t)Q\), for \(t\in[0,1]\), and impose the squared-side cutoff \(\kappa=2q\). The original square persists and scales by \(1+t\), so

\[
q(t)=q(1+t)^2,\qquad q(t)=2q\iff
\boxed{t=\sqrt2-1.}
\]

The associated four-edge chart's *degree-at-most-ten cutoff polynomial* has exactly **one** real root in \([0,1]\), which the Sturm routine encloses in a certified rational bracket. The endpoint filtered counts from the independent v223 solver are \(0\) and \(1\). This checks the dynamic polynomial against a distinct complete fixed-polygon algorithm, without using floating point to infer the crossing.

## 7. Global scan boundary and what remains unresolved

The all-chart scan covers ordered edge quartets (with optional repeated edges), both orientations, and all nonsingular nonpersistent charts, and explicitly reports candidate event root intervals. For the positive convex model on \(t\in[3/4,1]\) at \(\kappa=1/100\), it visits **408** charts after the justified exclusion of three-or-more-points-on-one-straight-edge products; **128** charts contain at least one identically zero event polynomial. These are *not* 128 confirmed physical square events. The raw list also contains multiple polynomial roots for the same physical square and many candidate roots outside actual edge domains.

As a result, the implementation **does not** assert that the overall orbit count is constant on all apparently root-free chambers, nor does it assign a global mod-2 incidence index from this list. To turn candidate isolation into a full polygonal cobordism, future work must:

1. Resolve persistent edge-face and all-time singular charts with exact semialgebraic feasibility and geometric orbit deduplication.
2. Cluster algebraic roots of different polynomials into **one** physical event, including shared polygon vertices and multiple D4 labels, without comparing rounded root times.
3. Analyze genuine birth/death events where a square branch goes through \(D(t)=0\), not simply classify those as isolated crossings.
4. Distinguish and account for positive-dimensional square families such as rectangles; v223 properly refuses to assign finite parity on these.
5. Establish a useful *global*, nonzero, size-protected invariant under \(C^0\) passage to arbitrary Jordan curves. v220 proves that merely tracking selected shrinking squares is insufficient.

The last item is the main conjectural bridge; **the unrestricted Toeplitz square-peg conjecture remains open**. None of the candidate-root computations supplies a global parity theorem or a Floer continuation morphism.

## 8. Reproduction and antecedents

```bash
python -m unittest discover -s tests -v
python code/moving_polygon_events_v224.py
```

The implementation depends only on the Python standard library and the previously merged `polygon_orbit_census_v223.py`. The event-isolation theorem uses standard Cramer's rule and Sturm's theorem; neither is claimed as a new algebraic discovery. Exact-rational regression tests check all 4×4 rational Cramer numerators against the *independent* v223 linear solver at 80 seeded chart/time configurations, the known vertex square, irrational cutoff, full repeated-edge candidate audit, and the polygon-simplicity gate.

**Internal:** [v223 complete fixed-polygon square census](v223-complete-polygon-orbit-census-and-homotopy-audit.md), [v222 multibranch PL incidence](v222-multibranch-cutoff-and-orbit-audit.md), [v220 accumulating small squares](v220-fixed-rectifiable-collision-link-obstruction.md), [v218 cyclic-label orientation cancellation](v218-global-configuration-symmetry-and-collision-audit.md).  
**External context, not used as a new proof:** Cantarella–Denne–McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, arXiv:2103.07506; Matschke, *On the Square Peg Problem and Its Relatives*, arXiv:1001.0186.
