# v226 — Algebraic event clustering and complete geometric-square census between events

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Status:** Exact-rational finite-model theorem and implementation, independent exact polygon-census crosschecks, and negative controls. **Not** an unrestricted Jordan-curve square-peg proof. **Not** a complete global orbit cobordism at singular times.

## Executive results

v224 produced exact *necessary* polynomial event candidates during an affine-vertex motion of one rational Jordan polygon. v225 removed 128 persistent or singular-chart flags in a specific convex-quadrilateral example, leaving 224 raw candidate entries. v226 performs the next logically necessary steps:

1. It **identifies exactly when two independently obtained polynomial roots are the same real algebraic time**, even if their rational isolating intervals overlap and the time is irrational. This uses polynomial GCDs and exact Sturm root counts, not floating-point proximity.
2. It builds an **ordered, disjoint algebraic event-time ledger** with original edge charts and equation types retained as provenance.
3. It uses v223's independent **complete geometric-square census** at rational points in every event-free open chamber. Each square is counted as one geometric orbit, not four cyclic labelings or eight dihedral presentations; repeated edges and polygon vertex hits are included.
4. It **refuses unsupported deductions**: a chart polynomial's root is not necessarily a physical square event; a singular positive-dimensional square family cannot be assigned a finite parity; and the exact square count *at an irrational event time* is not guessed from a nearby rational polygon.

Two whole-polygon controls make the distinction quantitative:

| Control | Raw chart-root entries | Distinct exact algebraic times | Geometric orbit counts on consecutive open chambers | What changes |
|---|---:|---:|---|---|
| Convex vertex-alignment motion | **224** | **1**, exactly `t=8/9` | **1, 3** | Net +2 orbits; at the event itself there are **2** orbits |
| Uniform scaling of a fixed quadrilateral | **88** | **6** | **0, 0, 0, 0, 1, 1, 1** | **One** cutoff crossing, at `t=√2−1`; five algebraic candidates do not change orbit count |

The first case has no change in the odd/even count *on generic chambers*. The second illustrates that **size cutoff walls can change parity** even when geometric squares persist continuously. Neither statement creates a protected odd global square count as the size cutoff approaches zero.

## 1. Exact algebraic root equality

Let \(p,q\in\mathbb Q[t]\setminus\{0\}\) be nonconstant polynomials. Suppose rational closed intervals \(I=[a,b]\), \(J=[c,d]\) each certify exactly **one** real root of their respective defining polynomial, with no endpoint root unless the interval has collapsed to the certified exact rational root. Write \(K=I\cap J\).

**Lemma (root identity).** The two isolated roots are the **same algebraic number** if and only if

\[
K\ne\varnothing\quad\text{and}\quad \gcd(p,q)\text{ has a real root in }K.
\]

**Proof.** If the roots agree at \(\alpha\), then \(p(\alpha)=q(\alpha)=0\), so the minimal polynomial of \(\alpha\) divides both rational polynomials and hence their GCD. Conversely, a GCD root in \(K\) is a root of both \(p\) and \(q\) in their certified one-root isolating intervals. Therefore it equals each of the unique isolated roots. \(\square\)

Sturm's theorem counts roots of the squarefree GCD in the open interior of a rational interval. Explicit rational evaluation tests both endpoints. This works for repeated roots of the original polynomials, irrational roots, and overlapping brackets that accidentally contain **different** real roots. The code groups equal roots transitively, computes their shared GCD factor (possibly reducible), then refines different roots until their returned isolating intervals are disjoint. A linear shared factor yields an **exact rational event time** even when the original dyadic isolator did not happen to hit that rational.

**Important:** An isolated irrational event is represented by a defining rational polynomial and a certified rational interval, *not* by a decimal approximation. The implementation does not claim it has extracted the irreducible minimal polynomial of every non-linear algebraic root. Its `minimal_factor` field means only a common defining factor.

## 2. Conditional chamber-count theorem

Take a rational Jordan polygon \(P(t)\) whose corresponding vertices move affinely in common time \(t\in[a,b]\subseteq[0,1]\). Fix \(\kappa>0\) as the minimum **squared** side length. Assume:

1. An exact independent certificate guarantees every intermediate polygon is an embedding, using v224's sufficient all-triple-orientations criterion (an overly conservative but valid sufficient test).
2. Every ordered four-edge chart and each square orientation has been examined; charts with three square vertices on one straight edge are geometrically impossible for a positive-side square.
3. The retained v225 strata exhaust every viable **positive-side** configuration over each open time chamber. Persistently singular feasible families and cutoff-locked charts are **reported inconclusive**, not thrown away. Canonical half-open outgoing-edge ownership is used at polygon corners. `incoming-edge-duplicate` charts are excluded only where their regular solutions are represented identically by outgoing-edge charts.
4. All nonzero Cramer determinants, source-edge endpoint polynomials, squared-side cutoff polynomials, and necessary augmented-rank obstructions have been isolated. Algebraic event roots are grouped and ordered exactly as in Section 1. A family whose exceptional-time behavior violates these hypotheses is outside this certification.
5. At a rational sample in each open chamber, the complete v223 polygon census returns a finite set, not a positive-dimensional family.

Then the **number of geometric square orbits of side squared strictly larger than \(\kappa\)** is constant on each connected open time interval between consecutive event times. The count is exactly the v223 census at *any* rational sample in that interval.

**Why this holds.** On a nonsingular four-edge chart, v224 gives the four square parameters as Cramer ratios \(u_i(t)=N_i(t)/D(t)\). The strict canonical edge-admissibility and positive-size status can only change at zeros of the finitely many chart polynomials, provided no identically-on-boundary or persistently singular positive-sized family has been overlooked. Away from this zero set, the parameters and the side length vary continuously and the chart remains admissible or inadmissible throughout the interval. At every instant of such an interval, the finite collection of all geometric configurations is obtained by taking these root branches and quotienting duplicate labels under the square's dihedral symmetries and polygon-corner ownership. Distinct nondegenerate geometric square orbits cannot merge without encountering a chart boundary or a singular equation, which would be among the isolated candidate times. Therefore the number of geometric orbits is locally constant, hence constant on each chamber. A single rational v223 census computes it without relying on labeled chart counts. \(\square\)

**What the theorem does not give:** a complete local cobordism **through** every singular algebraic event time, particularly where positive-dimensional square families or exceptional rank drops occur. For a *rational* event time, v223 can census that snapshot independently. For an **irrational** event time, v226 reports its instant count as **unknown** while preserving rigorous rational counts on the left and right. An even change \(0\to2\) or \(1\to3\) is a *net* change of two orbits, not an identification of which individual analytic branches meet at the event. Nor does an event with unchanged count necessarily have no branch reconnection.

## 3. Convex quadrilateral: 224 labels reduce to one event

Use v223's affine one-polygon deformation

\[
P(\tau)=((0,0),(2,0),(3-\tau,3-3\tau/2),(0,2-3\tau/2))
\]

and let the **normalized** homotopy parameter be \(t\) with \(\tau=9t/8\). Audit \(3/4\le t\le1\), with \(\kappa=1/100\). Its strictly positive turning determinants certify a simple, strictly convex polygon throughout this range.

v225's exhaustive chart classifier produces **408** chart records: 280 regular, 40 canonical vertex-owned, 40 redundant incoming-edge presentations, and 48 geometrically impossible incident-edge 2+2 records, with zero genuinely unresolved. The eligible records produce **224 raw polynomial-root labels**. They all identify the very same algebraic time:

\[
\boxed{t=\tfrac89\quad\Longleftrightarrow\quad\tau=1.}
\]

The shared polynomial GCD after grouping is proportional to \(9t-8\). The positive-side geometric orbit counts certified by v223 are:

- Entire open chamber \((3/4,8/9)\): **one** square orbit.
- Exactly at \(t=8/9\): **two** square orbits, including \((1,0),(2,0),(2,1),(1,1)\), an actual square with a polygon-vertex incidence.
- Entire open chamber \((8/9,1)\): **three** square orbits.

Hence the generic *net* orbit-count change is \(+2\), and its parity is unchanged. This does **not** prove an all-polygon conservation law; it proves the complete count on this certified algebraic deformation.

## 4. Similarity motion: 88 labels, six times, one real count transition

Take one embedded rational quadrilateral

\[
P=((0,0),(2,0),(3,3),(0,2))
\]

and scale the **entire polygon** by \(1+t\), \(t\in[0,1]\). V223 proves it has exactly **one** geometric square orbit, of squared side \(q_0=72/25\). Its squared side during the motion is \(q(t)=(72/25)(1+t)^2\). Select the **strict** positive-size cutoff \(\kappa=144/25=2q_0\).

The physical cutoff event is therefore exactly

\[
(1+t)^2=2\quad\Longrightarrow\quad\boxed{t=\sqrt2-1}.
\]

The v224–v225 polynomial atlas generates **88** raw candidate-root labels. v226 clusters these into **six distinct algebraic event times**, each exactly isolated. The independently certified geometric-orbit counts in the seven complementary chambers are

\[
\boxed{(0,0,0,0,1,1,1)}.
\]

Five candidate clusters do **not** change the geometric orbit count. The fourth event cluster contains the true irrational cutoff; its rational isolator is certified to contain exactly \(\sqrt2-1\). Its instant orbit count is intentionally *not* asserted by the general solver; geometrically the strict threshold excludes the orbit exactly at equality.

This is a hard control against conflating polynomial candidate roots with actual physical events. In particular, one cannot infer a parity jump from the **number of raw root labels** or the **number of clusters** alone.

## 5. Testing, implementation, and fail-closed cases

### Implementation

- `code/algebraic_orbit_events_v226.py`: exact shared-root equivalence from GCD + rational Sturm tests; disjoint root isolation; event provenance; independent complete polygon census on all rational event-free chambers; rational-time singular snapshots; summary reporting.
- `tests/test_algebraic_orbit_events_v226.py`: **24 exact tests** including close distinct irrational roots with overlapping brackets, duplicate and repeated roots, shared GCD roots, rational event grouping, the two whole-polygon examples, and explicit failures for genuinely singular rectangle square families.
- All arithmetic is from `fractions.Fraction`; no SymPy is required at runtime. The v223/v224/v225 modules remain intact and are called as independent implementations. 

### Reproduce

```sh
python -m unittest discover -s tests -v
python code/algebraic_orbit_events_v226.py
```

The standalone combined local v224+v225+v226 suite has **74 successful tests**. The GitHub workflow is extended to run all earlier research tests, including v223's independent polygon census, and verify the preserved ZIP archive SHA-256 hashes.

### Failure and scope discipline

- An **identically** cutoff-locked chart or a persistently feasible singular chart is rejected as inconclusive by the event-cover engine. The rectangle is a required regression: its continuous square families cannot be assigned an unsupported finite mod-2 parity.
- It also rejects a proposed closed time interval if a candidate event lies **exactly at an endpoint**; choose a smaller generic interval or perform separate endpoint analysis.
- A `RootBracket` is a **certified one-root isolator**. If inputs are merely overlapping approximate boxes, the equivalence lemma would not apply; the module checks that each proposed bracket is truly isolating.
- Coincident algebraic times are grouped without treating every original polynomial root as a distinct physical square. If two *different* roots initially have overlapping isolators, their defining polynomials' GCD separates them and the isolators are refined.
- A zero change between neighboring chamber counts means **no net change in the number of geometric orbits**, not necessarily absence of a branch relabeling, connection, or local interaction. No explicit branch-to-branch continuation graph is constructed here.
- Fixed positive \(\kappa\) avoids collapse at **that** cutoff. The v220 infinite-notch curve proves that an arbitrary rectifiable curve need not have a positive minimum square size across *all* its squares; no uniform positive cutoff for arbitrary Jordan limits has been obtained.

## 6. What v227 should do

The best next step is **algebraic event-germ classification**, rather than merely adding root counts: determine which geometric orbits on the two sides of one real algebraic time are continuations of the same square, which are births/deaths, and which only change chart ownership. Such a classification would need exact algebraic evaluation at irrational times, including rank drops and polygon-vertex incidence. It should be checked against *all* edge charts and both square orientations and reconciled with the orbit deduplicator, not with signed labeled counts.

Even a completed polygonal cobordism would still be a finite-model infrastructure theorem. An unrestricted Toeplitz proof also needs a nontrivial **global invariant** and a non-collapse mechanism as polygonal approximations tend to arbitrary continuous Jordan curves. This research program has **not** supplied those two ingredients.

## 7. References and provenance

Internal base: [v223 polygon orbit census](v223-complete-polygon-orbit-census-and-homotopy-audit.md), [v224 algebraic event covers](v224-algebraic-moving-polygon-event-covers.md), and [v225 singular/persistent edge charts](v225-persistent-edge-strata-rank-certificates.md).

External mathematical context (not claimed as a new all-Jordan theorem):

- Cantarella, Denne & McCleary, *Transversality for Configuration Spaces and the “Square-Peg” Theorem*, [arXiv:1402.6174](https://arxiv.org/abs/1402.6174).
- Cantarella, Denne & McCleary, *Configuration Spaces, Multijet Transversality, and the Square-Peg Problem*, [arXiv:2103.07506](https://arxiv.org/abs/2103.07506).
- Matschke, *On the Square Peg Problem and Its Relatives*, [arXiv:1001.0186](https://arxiv.org/abs/1001.0186).

**Scientific caveat:** Passing exact arithmetic and regression tests corroborates the stated **finite symbolic procedures**. It is not a Lean/Coq-certified proof and does not establish the unrestricted Toeplitz square-peg conjecture.
