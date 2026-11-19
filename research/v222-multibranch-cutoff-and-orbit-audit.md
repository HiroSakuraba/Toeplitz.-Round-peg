# v222 — Multiple square branches, PL knot incidence, and the distinction between labeled and geometric counts

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** A conditional exact-rational theorem and an executable multichart model. The general Toeplitz square-peg conjecture remains open. This is **not** a relative Floer complex, an all-Jordan mod-2 invariant, or a proof that arbitrary rough-curve limits retain a positive square size.

## 1. Question and principal outcome

[v221](v221-positive-size-cutoff-incidence.md) implemented an exact mod-2 boundary ledger for **one** transverse square branch on four moving line segments. A polygonal arc has many line segments. The same square crosses chart boundaries as its vertex passes an internal corner, and two branches may meet and disappear at a corner. Counting each local chart separately without gluing would mistake these internal events for **true square exits**. Separately, the same geometric square has four cyclic labeled presentations for a fixed orientation, or eight if both orientations/reversals are included. Counting the eight modulo two would always give zero even when the geometric square-orbit count equals one.

v222 resolves **both finite-dimensional bookkeeping issues** within an explicitly restricted rational PL class:

1. It enumerates all regular chart branches on **four independently marked polygonal arcs translating rigidly with constant rational velocities**, with a prescribed physical deck tuple and one square orientation.
2. It exactly pairs both sides of every admissible internal PL knot in **original arc-parameter space**. A through-crossing and a two-to-zero PL fold contribute **two half-branch ends** and hence **zero exposed mod-2 boundary**.
3. It counts positive-size-filter crossings, including irrational quadratic crossing times, and actual marked source-arc exits. The resulting endpoint parity is accounted for entirely by **source** and **size** crossings.
4. It defines a **geometric D4 orbit key** for a verified square. Its eight labelings form one geometry class; the helper is not an exhaustive square-orbit solver for a Jordan curve.
5. It separately verifies the standard smooth local fold model `F_t(x,y,z,w)=(x²−t,y,z,w)`. This abstract model establishes even pair creation/destruction at a generic fold, but is **not derived from a translating-PL square-residual chart**.

The work is an exact *incidence chain auditor*, not the global nonzero mod-2 class needed for Toeplitz.

## 2. A conditional multichart incidence theorem

For `i=0,1,2,3`, let `A_i` be a compact, continuously parametrized oriented PL arc with rational vertices, each segment traversed linearly. Fix rational velocities `v_i∈Q²`; let

\[
 A_i(t,u)=A_i(0,u)+t v_i+(k_i,0),\qquad k=(0,k_1,k_2,k_3)\in\mathbb Z^4.
\]

Thus each **whole arc** translates rigidly; it does not change its internal segment directions over time. Define the square residual, with `J(x,y)=(-y,x)`,

\[
F_\sigma(t,u)=\big(A_2-A_1-\sigma J(A_1-A_0),\ A_3-A_0-\sigma J(A_1-A_0)\big)\in\mathbb R^4,
\quad \sigma\in\{\pm1\},
\]

and squared side `q(t,u)=|A_1-A_0|²`. For `a<b` and **strictly positive** `κ`, count only roots with `u∈(0,1)^4` and `q>κ` at times `a,b`.

Each choice of one of the finitely many constituent segments from each arc defines a product chart with local coordinates `w∈[0,1]^4`. On that chart, because only whole-arc translation depends on time,

\[
F_\sigma(t,w)=M w-r_0-t r_1
\]

with rational, *time-independent* `4×4` matrix `M`. If `M` is invertible, its unique algebraic branch is

\[
w(t)=M^{-1}r_0+tM^{-1}r_1,
\qquad q(t)=At^2+Bt+C.
\]

The intersection of its parameter-box constraints and `[a,b]` is a rational **closed time interval** (or empty). Its filtered part is obtained by excluding the roots of `q(t)−κ`.

**Genericity/admissibility conditions.** The following are checked rather than assumed silently by the implementation:

- Every **feasible** product chart has invertible `M`. A singular chart may be skipped only after **exact rational feasibility** in the full five-dimensional time/segment-parameter prism shows that it has **no square at all**. Singular feasible charts raise `NonGenericIncidence`.
- No entire branch is pinned to a source or internal knot face; no admissible isolated one-time chart-corner roots; no solution on a time/source/knot/cutoff corner or at a time slice's source-chart boundary.
- Every positive-size internal-knot endpoint of a regular chart branch is matched **exactly once** by the adjacent chart, at the same time and **same four global source-arc parameters** (not just the same ambient square vertex positions). The two half-branches may lie on opposite time sides (transit) or the same side (PL pair annihilation).
- Positive-size cutoff contacts on active branch intervals are **simple**. Tangential cutoff contacts and simultaneous size/source events are refused as nongeneric.

**Theorem (finite generic multichart incidence).** Under those assumptions, let `N_a,N_b` count admissible **ordered labeled** square solutions at times `a,b`, and let `E_source,E_size` be the respective **exposed** source-wall and positive-size-cutoff crossing counts over the time interval. Then

\[
\boxed{N_a+N_b\equiv E_{\mathrm{source}}+E_{\mathrm{size}}\pmod2.}
\]

Every paired internal PL interface contributes two incident half-branches and **does not appear** on the right. The square-size guard `κ>0` excludes zero-side configurations **for this fixed scale**, with no assertion of a uniform cutoff as curves vary.

**Proof.** Each regular chart has one affine-time solution branch. The four box constraints and the bounded time interval cut it into either an empty set or a closed interval. Removing points where `q≤κ` divides an admissible branch into intervals. Every resulting one-dimensional component has even cardinality of endpoints. These endpoints are on the initial/terminal time slices, the original source-arc parameter endpoints, the strictly positive cutoff `q=κ`, or an internal source chart boundary. By the exact adjacent-chart gluing hypothesis, internal endpoints pair at identical global configurations and cancel modulo two. All other endpoints are precisely those counted on the right and at the initial/terminal time slices. Summing over the finitely many branch intervals gives the displayed identity. No extra contribution is made by a PL pair birth or death at a **single internal knot**: it has exactly two incident half-branches. `□`

This is a **standard one-dimensional mod-2 boundary argument**, not a newly proved universal square-peg theorem. The code provides an exact decision/ledger on its rational input class; it does not certify all smooth or rough homotopies.

### Why the parameter-identity check matters

A source arc could revisit the same ambient point in different parameter locations. Two roots with the same **planar** ordered square but different source parameters are distinct roots in the configuration space. A gluing routine that matches only planar coordinates could silently identify them and break the incidence theorem. v222 forms global segment parameters `j+w_j`, checks corresponding planar positions, and matches only genuinely adjacent charts. This also retains the supplied integer physical lifts, without summing over periodic translates.

## 3. Exact two-root polygonal fold example

Take four short rational graph arcs centered near the oriented square with vertices

\[
 p_0=(1/2,0),\quad p_1=(5/8,1/8),\quad
 p_2=(1/2,1/4),\quad p_3=(3/8,1/8).
\]

Use horizontal half-width `h=1/32`. The linear slopes of arcs `0,1,3` are respectively `−2,3,−1/3`. Replace arc `2` by a two-segment tent with x-coordinates

\[
\ell=15/32,\ c=1/2,\ r=17/32,
\]

and initial heights `target(ℓ)−d`, `target(c)+d`, `target(r)−d`, where

\[
\operatorname{target}(x)=1/4-2(x-1/2),\qquad d=1/100.
\]

Translate only the tent **downward by `2dt`**, for `0≤t≤1`. The exact solver finds **two distinct transverse squares** at `t=0`, one in each tent segment, then their algebraic branches both reach the same **internal arc-2 knot** at `t=1/2`, where they form the planted square of squared side `1/32`. Neither branch has a positive-length continuation for `t>1/2` inside the marked arcs.

With `κ=1/1024`, both initial squares pass the filter and the result is

| Quantity | Count |
|---|---:|
| `N_0` | 2 |
| `N_1` | 0 |
| Pair of coincident knot **half-branch** endpoints | 2 (one matched internal event) |
| Exposed source-wall crossings | 0 |
| Positive cutoff crossings | 0 |

Both sides of the incidence law are zero mod two. Reversing time gives a two-square **PL birth**. This does **not** assert a smooth saddle-node bifurcation; the family is PL and has a corner at the meeting configuration.

For the different cutoff `κ=1/40`, exactly **one** initial square is admitted. One genuine cutoff crossing accounts for its loss while the two underlying branches still meet at an internal knot (the knot event contributes no exposed boundary). This tests mixed cutoff and knot bookkeeping. For `κ=1/32`, cutoff and knot meet simultaneously, so the checker refuses the non-generic event rather than assigning a false count.

An independent example inserts a collinear polygonal subdivision knot into an arc and translates that arc so one square crosses the knot at `t=1/2`. The square persists (one on each time slice), and the paired internal knot event contributes zero. Artificial static subdivision **exactly along a continuously pinned square** is deliberately rejected as a nontransverse persistent face, not double-counted.

## 4. Ordinary source/cutoff crossings and exact comparisons

When every arc consists of one segment, the new all-chart routine is compared with the **independent v221 single-branch solver**, including 36 seeded rational moving-square cases.

- **One cutoff crossing:** side `e(t)=(t/8,t/8)`, squared cutoff `κ=1/128`, time `[1/4,3/4]`; counts `(0,1)` and one size event.
- **Two cutoff crossings:** side `e(t)=(1/4-t/2,0)`, `κ=1/256`, time `[1/4,3/4]`; counts `(1,1)` and two size events.
- **Irrational event:** side `e(t)=(1/8,t/8)`, `κ=3/128`, cutoff at `t=1/√2`; count is certified by rational polynomial signs.
- **Source exit:** translating source arc `1` upward with `κ=1/1024` gives one original-arc endpoint crossing at exactly `55/19968`, counts `(1,0)` and no cutoff crossings.

All these are *four marked arcs*, not a construction extracting compatible arcs from every Jordan curve.

## 5. Cyclic symmetry, the eight labelings, and the orbit distinction

For an exact **nondegenerate** square `p=(p0,p1,p2,p3)`, the four cyclic rotations preserve the oriented square equations. Reversing cyclic order changes the orientation `σ`; another four rotations produce the opposite orientation. Therefore there are eight `D4` labelings of one physical square. The exact helper `dihedral_labels` builds those eight distinct ordered quadruples and `canonical_square_orbit` maps all eight to the same sorted set of four exact planar vertices. It refuses collapsed or nonsquare tuples.

**Crucial caution:** This is only an **orbit-key normalization** on supplied exact square witnesses, *not* a proven way to divide a general count of labeled configurations by eight. The incidence ledger's `N_a,N_b` count one **chosen labeling** on four independently designated source arcs with fixed `σ` and lifts. They are **not** geometric-square orbit counts across all six parameter-order sectors, nor do the four arcs necessarily lie on one curve. A global quotient chain would require compatibility of parameter sectors, their boundaries, geometric orbit identification and their action on a suitable homotopy.

The independent [v218](v218-global-configuration-symmetry-and-collision-audit.md) convex-quadrilateral diagnostic (one geometric square with eight labeled descriptions) is included as an optional **full-repository** cross-check. A wrong physical lift tuple is also checked: the correct `(0,−1,0,+1)` admits an exact nonzero square, while the unshifted tuple `(0,0,0,0)` does not. Reflection transforms the tent model's positive orientation `+1` into `−1`, as required.

### Smooth fold model is separate

The familiar local model

\[
G_t(x,y,z,w)=(x^2-t,y,z,w)
\]

has no root for `t<0`, one **singular** root for `t=0`, and two regular roots for `t>0`, whose Jacobian signs are opposite. The total 1-dimensional zero locus in `(t,x,y,z,w)` is the smooth parabola `t=x²`, with **no boundary** at the fold. Thus creation/destruction of a generic pair does not change mod-2 parity. This is an elementary, **abstract** normal-form verification; the current moving-PL code does not certify arbitrary smooth fold events in real square-residual families.

## 6. Validation, limits, and next move

Reproduce in a complete checkout:

```bash
python -m unittest discover -s tests -v
python code/multibranch_incidence_v222.py
```

The 25 new v222 tests run in exact `Fraction` arithmetic: **23 passed locally and two optional repository integration tests skipped** in the isolated work directory. GitHub CI runs the two integration tests (against v218 and v220), all prior v212–v221 tests, every research audit script and the original archive checksums. This passing status must be verified separately from opening the PR.

**Not proved or implemented:**

1. A complete one-parameter **generic smooth** or arbitrarily varying PL Jordan-curve orbit cobordism with creation/destruction, multiple source sectors, and collision strata. v222's arcs translate rigidly and are independent.
2. A global mod-2 cycle of **geometric** square orbits with correctly handled `D4` quotient and all six cyclic-order sectors.
3. A topological reason for an **odd** global orbit count or a nonzero relative class on an arbitrary Jordan curve.
4. Any uniform **positive size** bound or protected non-collapse class along `C0` approximation. The [v220](v220-fixed-rectifiable-collision-link-obstruction.md) fixed infinite-notch curve continues to invalidate naive assumptions about a minimum size among all square witnesses.
5. Any theorem transferring the area-balanced **four-periodic-strand** Floer results to source arcs and the physical deck lifts of **one** Jordan curve. The [v213](v213-exact-PL-localization-and-lift-audit.md) wrong-lift control remains essential.

**Recommended v223:** Rather than enlarging the restricted root solver yet again, define a legitimate **geometric orbit-level quotient incidence chain at fixed positive cutoff** for a homotopy of *one* embedded Jordan polygon. An explicit homotopy should preserve the curve's embeddedness, and its full orbit list must include chart endpoints and repeat-edge configurations. Check whether the counts can be glued across label sectors without an unaccounted collision/cutoff boundary. If not, isolate a precise counterexample and stop short of a global class claim. Any subsequent rough-curve limit still needs an independent positive-scale theorem.

### Sources

- Cantarella, Denne and McCleary, *Transversality for Configuration Spaces and the “Square-Peg” Theorem*, [arXiv:1402.6174](https://arxiv.org/abs/1402.6174): generic smooth transversality and odd-square conclusions under relevant boundary-disjointness hypotheses.
- Cantarella, Denne and McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, [arXiv:2103.07506](https://arxiv.org/abs/2103.07506): compactified configuration spaces and the boundary-disjointness qualification.
- Benjamin Matschke, *On the Square Peg Problem and Its Relatives*, [arXiv:1001.0186](https://arxiv.org/abs/1001.0186): general square-peg context and generic-curve results.
- Internal: [v221 filtered incidence](v221-positive-size-cutoff-incidence.md), [v220 collision notches](v220-fixed-rectifiable-collision-link-obstruction.md), [v218 cyclic symmetry](v218-global-configuration-symmetry-and-collision-audit.md), and [v216 exact boundary degree](v216-boundary-degree-global-arc-selection.md).

All code and tests are verification of the **stated finite model**; they are not machine-checked proofs in a proof assistant.
