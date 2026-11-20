# v262: Exact cyclic-quotient relative parity through genuine square annihilation

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Scientific status:** An exact-arithmetic, all-real-time **relative mod-2 degree** on the quotient of a selected eight-box cyclic atlas around the three tracked square configurations of the moving v254 pentagon. The quotient is protected through a genuine noncollinear-vertex square-pair collision, even though the integer signed count over all four **labeled** cyclic lifts vanishes.

**Not claimed:** an absolute Euler-class obstruction, a complete all-real-time census of every inscribed square, an all-six-cyclic-order-sector atlas, a universal positive-size lower bound, or the unrestricted Toeplitz theorem.

## 1. Why v261 was not yet an unlabeled invariant

The v261 selected atlas gave two disjoint *labeled* four-dimensional parameter regions:

- a persistent square with exact degree **+1**, on pentagon edge assignment \((0,1,2,4)\);
- a continuous piecewise-affine region containing two squares with degrees **-1** and **+1** before their collision, one shared vertex square at the event, and no roots after. Its total degree stays **0** by the v260 all-real-time outer-boundary certificate.

The disjoint union has degree \(1+0=1\). That does **not** automatically define an integer degree on unlabeled squares. In fact, v218 already proved the general obstruction: changing which corner is called the first corner reverses the four-dimensional source orientation while preserving the determinant of the residual's target representation. The cyclicly labeled integer degrees alternate \(+1,-1,+1,-1\).

v219 also proved that the absolute mod-2 Euler-class argument **alone** cannot force a square in the positive cyclic-order component. The present result addresses only a specific **relative**, zero-free-exterior problem. Its nontrivial parity comes from the geometric square-residual section with fixed external boundary data, not from the topology of an unrestricted abstract rank-four bundle.

## 2. A genuine free fourfold cyclic atlas

Write the counterclockwise boundary of the v254 pentagon using an unwrapped parameter \(q\), with edge \(e\) having local coordinate \(q=e+u\) and period \(L=5\). A positively ordered configuration satisfies

\[
0\le q_0<5,\qquad q_0<q_1<q_2<q_3<q_0+5.
\]

For the persistent square, the entire original rational box has the exact **global** parameter ranges

\[
U_P:
[0.14,0.17]\times[1.30,1.34]\times
[2.78,2.84]\times[4.27,4.31].
\]

For the colliding pair, the single **hinge-glued** box has ranges

\[
U_D:
[0.90,0.93]\times[1.999,2.001]\times
[3.43,3.46]\times[4.73,4.78].
\]

The second corner parameter of \(U_D\) crosses the actual polygon vertex at \(q_1=2\). Both regions are contained in the positive cyclic-order configuration domain for their whole parameter boxes, not only at their roots.

Let \(T\) denote cyclic relabeling:

\[
T(q_0,q_1,q_2,q_3)
=(q_1,q_2,q_3,q_0+5),
\]

followed by subtracting a common multiple of 5 so the first parameter lies in \([0,5)\). Thus \(T^4=\mathrm{id}\) on the configuration space of parameter points modulo full turns.

The exact [v262 verifier](../code/square_quotient_parity_v262.py) constructs the eight boxes

\[
U_P,TU_P,T^2U_P,T^3U_P,\quad
U_D,TU_D,T^2U_D,T^3U_D.
\]

All **eight** are pairwise disjoint. This is certified without solving a single equation, simply by checking the **strict separation of their first-coordinate intervals**:

| Region | Cyclic lift \(k\) | First parameter \(q'_0\) |
|---|---:|---|
| Persistent | 0 | \([0.14,0.17]\) |
| Pair | 0 | \([0.90,0.93]\) |
| Persistent | 1 | \([1.30,1.34]\) |
| Pair | 1 | \([1.999,2.001]\) |
| Persistent | 2 | \([2.78,2.84]\) |
| Pair | 2 | \([3.43,3.46]\) |
| Persistent | 3 | \([4.27,4.31]\) |
| Pair | 3 | \([4.73,4.78]\) |

Every translated four-box stays strictly inside the positively ordered cyclic-configuration component. Its coordinate transition is a rational affine permutation and translation. The code verifies the group law \(T^4=\mathrm{id}\), zero stabilizers, all eight full-box separation inequalities and explicit polygon-parameter agreement at exact rational time endpoints. The hinge \(q_1=2\) is separately checked on both sides and at the interface.

The free \(C_4\) quotient of this saturated atlas has exactly **two** relative chart regions, one for the persistent orbit and one for the colliding pair orbit. Each box in the base lift is a complete representative of its four cyclic images. Consequently, computing quotient mod-2 degree as a degree in one representative box **counts each geometric square once**, not four times.

## 3. Equivariance and the orientation-local-system character

Use the same four-dimensional counterclockwise square-defect map as v259–v261:

\[
F(p_0,p_1,p_2,p_3)=
\bigl(p_2-p_1-J(p_1-p_0),\,p_3-p_0-J(p_1-p_0)\bigr)
\in\mathbb R^4.
\]

Cyclic permutation of the corners changes \(F\) by a constant invertible matrix \(R\):

\[
F(Tp)=RF(p),\qquad
R=\begin{pmatrix}
-1&1&1&0\\
-1&-1&0&1\\
0&1&0&0\\
-1&0&0&0
\end{pmatrix}.
\]

The exact verifier establishes

\[
R^4=I,\quad\det R^k=+1;
\qquad
P^4=I,\quad\det P^k=(-1)^k,
\]

where \(P\) is the four-cycle permutation matrix acting on the local source coordinates. The residual identity is **linear in all eight original Cartesian corner coordinates**. The verifier tests all eight basis vectors for every \(k=0,1,2,3\), giving **32 exact universal linear identities**, plus checks on the actual moving polygon including the kink.

Hence the local signed index on cyclic sheet \(k\) is \((-1)^k\) times that on sheet 0. The character

\[
\chi(k)=(-1)^k
\]

is the orientation-local-system transition law. It does **not** provide a canonical global integer orientation of the quotient; the quotient domain remains nonorientable. Reducing the local index modulo two, however, eliminates the sign and yields a well-defined **local orbit parity** on the quotient of our chosen atlas.

**Crucial caveat:** the mod-2 degree of all the *labeled* roots is still zero, since every geometric orbit has four labels. Only the **quotient** count is the meaningful parity.

### The quotient is a bundle section, not an ordinary globally valued map

Because $F(Tq)=RF(q)$ rather than $F(Tq)=F(q)$, the quotient does **not** carry a globally defined map to a fixed $\mathbb R^4$ with this formula. Define the associated rank-four bundle

\[
E=(\mathcal U\times\mathbb R^4)/\bigl((q,v)\sim(Tq,Rv)\bigr)
\longrightarrow\mathcal U/C_4,
\]

where $\mathcal U$ is the disjoint eight-box saturation. The square residual descends to the section $s_t([q])=[q,F_t(q)]$. The local mod-2 index is the Brouwer degree in any chosen representative box, reduced modulo two. The exact source and target transition determinant calculations prove that this definition is independent of the cyclic representative. The external box boundaries are zero-free by the v260/v261 time-uniform certificates.

At a regular time, this mod-2 degree counts transverse *geometric orbits* once. At the collision, it is defined by the same relative homotopy-invariant local degree, **not** by directly counting the singular roots. On these two contractible quotient regions the degree is represented by the pair of mod-2 local indices; their sum is the atlas's relative invariant. No canonical integer trivialization is asserted globally.

## 4. Degree and orbit ledger on both sides of the irrational collision

The original v254 exact algebraic collision time \(\alpha\) lies strictly inside

\[
a=\frac{64722}{100000}<\alpha<
b=\frac{64723}{100000}.
\]

The v260 all-time zero-free outer boundary of the piecewise-affine pair box remains zero-free at the singular time when the two roots meet on its **internal** hinge. The v261 persistent box has a separate all-time zero-free boundary and integer degree \(+1\). The cyclic action maps the zero-free boundaries homeomorphically to the six other labeled sheets. Because it is a change of parameters and the target map is an invertible \(R^k\), no boundary zeros can arise after relabeling.

This proves that the relative mod-2 degree on the quotient of the **whole selected atlas** is the same for every real time, including \(\alpha\):

\[
\boxed{\deg_2(F_t,\mathcal U/C_4,0)
=(\deg F_t|_{U_P}+\deg F_t|_{U_D})\bmod2
=(1+0)\bmod2=1.}
\]

| Time | Persistent quotient region | Pair quotient region | Quotient degree |
|---|---:|---:|---:|
| \(t<\alpha\) | one root, degree 1 | two roots, indices \(-1,+1\), degree 0 | **1** |
| \(t=\alpha\) | one root, degree 1 | one **singular** vertex root, local degree 0 | **1** |
| \(t>\alpha\) | one root, degree 1 | no root, degree 0 | **1** |

An important mathematical subtlety is now explicit: **the parity of distinct geometric roots exactly at a singular event is not generally the mod-2 degree.** In these *selected* regions at \(t=\alpha\), there are **two** distinct geometric roots (the persistent one and the shared collision square), whose raw count is even; nevertheless the mod-2 relative degree is **one** because the collision square has **even local multiplicity modulo two**. Before and after the event, all roots are regular and the quotient degree equals their ordinary geometric orbit count modulo two.

At the two rational times \(a,b\), the exact verifier checks all three pre-event roots and the single post-event root directly on the actual pentagon and verifies all four distinct cyclic lifts per root: 12 labeled zeros before, four after. v254's independent complete geometric-square census at those two rational slices already established **3** and **1** whole-polygon square orbits respectively. The new verification of the \(C_4\)-saturated atlas does not elevate this to a complete whole-polygon count at every irrational time.

## 5. The cyclic-cover perturbation bound is also exact

The raw v261/v260 atlas has a rational external boundary residual lower bound

\[
\delta=\min(1/960000,1/19200000)=1/19200000.
\]

On cyclic sheet \(k\), \(F\) transforms as \(R^kF\), so

\[
\|R^kF\|_\infty
\ge\frac{\|F\|_\infty}{\|R^{-k}\|_\infty}.
\]

The program computes the **exact** \(\infty\)-norm of each integer \(4\times4\) matrix \(R^{-k}\) and takes the worst. This yields a certified strictly positive exterior bound over **all eight** translated boxes. Each component of the perturbed square residual changes by at most \(4\varepsilon\) when each Cartesian coordinate of each of the four boundary arcs is changed by at most \(\varepsilon\). The inherited budget

\[
\varepsilon=10^{-9}
\]

is strictly smaller than the resulting cyclic-cover boundary budget. Therefore the quotient mod-2 degree remains one for any common continuous perturbation of the underlying boundary parameterization that moves each of its points by at most this coordinatewise amount. The perturbed eight-sheet atlas must be constructed by cyclically relabeling that **same** perturbed curve. If each sheet is perturbed independently, the equivariant bundle section may cease to exist and no quotient-degree claim follows.

The inherited v260/v261 **positive side-length bounds** also survive each cyclic relabeling, since cyclically rotating a square does not change its side length. Hence the protected root remains a positive-size square for these perturbations.

Again, this is a quantitative theorem only near the **particular** moving pentagon and its selected source neighborhoods. It does **not** furnish a size floor for arbitrary Jordan approximations.

## 6. Limits, and what v263 must attack

The symmetry quotient and vertex gluing are rigorously resolved on this example, but several major tasks are still open:

1. **Complete the whole-pentagon event ledger for all real times.** The selected atlas contains exactly the **three tracked** square branches. The v223 exhaustive whole-polygon counts at rational slices do not preclude additional squares appearing between those slices. Exact all-real-time coverage of singular/repeated-edge charts is still needed.
2. **All cyclic-order sectors:** v218 found six labeled cyclic-order components. The present certificate only treats the positive cyclic-order sector (the appropriate one for the strictly convex example). Arbitrary nonconvex Jordan curves can admit nongraceful square configurations, and a general proof must address or legitimately discard those sectors.
3. **Collision/zero-size boundary:** the global positive-size configuration space is noncompact as the four boundary points collide. The natural compactification has unavoidable collapsed-square zeros. Neither local quotient parity nor a zero-free boundary around one existing square supplies uniform control in arbitrary \(C^0\) polygonal approximation.
4. **Do not re-prove generic polygon parity as a claimed Toeplitz breakthrough.** Generic odd-orbit parity results are part of the classical square-peg literature. A new general argument must add enough compactness or a suitable boundary invariant to handle arbitrary continuous Jordan curves.

**Bottom line:** v262 gives an exact nonzero **unlabeled, relative mod-2 degree on a saturated atlas** through a genuine square-pair annihilation. It repairs the relabeling ambiguity left by v261 *on that atlas*, and precisely separates genuine degree conservation from raw singular-root parity. It does **not** establish a new universal degree or solve Toeplitz.

## Replay

Standard-library-only exact rational arithmetic, reusing the established v260/v261 proof engines:

```sh
python code/square_quotient_parity_v262.py
python -m unittest -v tests.test_square_quotient_parity_v262
```

The [new regressions](../tests/test_square_quotient_parity_v262.py) verify exact fourfold group laws, all eight disjoint boxes, free cyclic action, universal residual identity on 32 linear basis cases, polygon coordinate compatibility across the kink, the regular root orbits, signed orientation transitions, zero sum on the labeled cover, nonzero parity on the quotient, uniform external perturbation budget and explicit fail-closed cases.

## References

- [v218: cyclic symmetry and collision-stratum obstruction](v218-global-configuration-symmetry-and-collision-audit.md).
- [v219: equivariant quotient and vanishing absolute Euler-class approach](v219-equivariant-quotient-collision-link.md).
- [v254: genuine positive-size square-pair annihilation](v254-genuine-square-pair-death.md).
- [v260: degree-zero cancellation through the genuine kink](v260-relative-degree-pair-death.md).
- [v261: oriented multi-region +1 relative degree](v261-oriented-multi-region-square-atlas.md).
- B. Matschke, *A Survey on the Square Peg Problem*, *Notices of the AMS* 61 (2014).
- D. Salamon, *Lectures on Differential Topology*, sections on degree additivity, homotopy invariance and mod-2 degree.
