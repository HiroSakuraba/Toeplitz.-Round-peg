# v218: Global ordered four-point configurations, cyclic sign cancellation, and collision boundaries

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Scientific status:** An elementary orientation/symmetry obstruction and collision-stratum analysis, with exact-rational checks and a deliberately limited polygonal model. **The unrestricted Toeplitz square-peg conjecture remains open.** These observations are related to established configuration-space approaches in the literature and are **not** claimed as a novel resolution of the global problem.

## 1. What changes after v217

[v217](v217-multiscale-relative-degree-ledger.md) proved that subdivision of four *marked* source arcs preserves the total integer Brouwer degree when no square enters or exits through the external product boundary. However, no argument supplies a nonzero total degree from *one arbitrary Jordan curve*. Attempting to sum all signed local contributions in the full four-point configuration space exposes a symmetry obstruction **before** any compactification or rough-limit question is addressed.

This note distinguishes three mechanisms that must not be conflated:

- **Cyclic relabeling** of a square changes the sign of a local Jacobian contribution in the cyclically ordered four-point domain. This makes the naive global **labeled signed** count cancel in quadruples.
- **Collision** of points on a compactified Jordan-curve configuration is the locus where squares can shrink to zero. This is not solved by canceling the cyclic signs.
- **Order sectors:** The full labeled four-point configuration space has six cyclic-order components, not only the two in which square-vertex order agrees with curve order. Any universal argument must explain why it can discard the other four, or include them.

## 2. The four-dimensional cyclic configuration chart

Let \(\gamma:S^1\to\mathbb R^2\) be a continuous *embedding*, and label four pairwise distinct parameter points in cyclic curve order by \(t_0,t_1,t_2,t_3\). In one cyclic-order component, use

\[
 (t,g_1,g_2,g_3)\in S^1\times\Delta_3^{\circ},\qquad
 g_i>0,\quad g_4=1-g_1-g_2-g_3>0,
\]

where successive parameter gaps are \(g_1,g_2,g_3,g_4\). The corresponding positions on the curve are
\(p_i=\gamma(t+\sum_{j=1}^{i}g_j)\) for \(i=1,2,3\) and \(p_0=\gamma(t)\), using periodic parameters.

Define the cyclic relabeling map

\[
T(t,g_1,g_2,g_3)=(t+g_1,g_2,g_3,1-g_1-g_2-g_3).
\]

Its local coordinate Jacobian is

\[
D T=\begin{pmatrix}
1&1&0&0\\0&0&1&0\\0&0&0&1\\0&-1&-1&-1
\end{pmatrix},\qquad \boxed{\det DT=-1}.
\]

This map generates a **free** cyclic \(C_4\) action on the space of four distinct, cyclically ordered parameters: no nontrivial cyclic relabeling fixes a labeled configuration. Since it reverses orientation, the quotient of this connected orientable open four-manifold by \(C_4\) is **nonorientable**. This statement concerns the parameter space alone; it does not assert that a chosen target bundle or Floer chain has already been constructed.

## 3. The square residual has an orientation-preserving target action

For \(\sigma\in\{+1,-1\}\), let \(J(x,y)=(-y,x)\) and define the four-component square residual

\[
F_\sigma(p_0,p_1,p_2,p_3)=
\begin{pmatrix}
p_2-p_1-\sigma J(p_1-p_0)\\
p_3-p_0-\sigma J(p_1-p_0)
\end{pmatrix}\in\mathbb R^4.
\]

A zero is a square of cyclic orientation \(\sigma\), **including the collapsed square** when all four points agree. Let \(\mathcal T(p_0,p_1,p_2,p_3)=(p_1,p_2,p_3,p_0)\). A direct expansion gives the *linear equivariance identity*

\[
F_\sigma(\mathcal T p)=B_\sigma F_\sigma(p).
\]

For \(\sigma=+1\) and \(\sigma=-1\), respectively, the exact matrices in the coordinate order \((F_x^{(1)},F_y^{(1)},F_x^{(2)},F_y^{(2)})\) are

\[
B_+=\begin{pmatrix}-1&1&1&0\\-1&-1&0&1\\0&1&0&0\\-1&0&0&0\end{pmatrix},\qquad
B_-=\begin{pmatrix}-1&-1&1&0\\1&-1&0&1\\0&-1&0&0\\1&0&0&0\end{pmatrix}.
\]

Both satisfy \(B_\sigma^4=I\) and

\[
\boxed{\det B_+=\det B_-=+1}.
\]

The implementation verifies the *entire* linear identity, not just the determinant, against 100 rational four-point inputs for each orientation.

## 4. Proposition: cancellation of a full cyclic orbit of signed local roots

**Proposition 1.** Suppose the curve parameterization is differentiable in the relevant neighborhoods and \(F_\sigma\) has an isolated regular zero at a cyclically ordered quadruple \(x\) with distinct points. Then the four roots \(x,Tx,T^2x,T^3x\) have local Jacobian signs

\[
\varepsilon,-\varepsilon,\varepsilon,-\varepsilon,
\]

and their integer sum equals zero.

**Proof.** Differentiate equivariance: \(D F_\sigma(Tx)\,D T(x)=B_\sigma\,D F_\sigma(x)\). Taking determinants and signs, with \(\det DT=-1\) and \(\det B_\sigma=+1\), gives \(\operatorname{sgn}\det D F_\sigma(Tx)=-\operatorname{sgn}\det D F_\sigma(x)\). Iterate four times. The action is free, so the four labels are distinct roots. QED.

**Consequence.** Whenever a finite global labeled signed count of these transverse zeros makes sense and respects cyclic relabeling, each orbit contributes **zero**. Thus the naive sum cannot produce the universal nonzero total sought in v217. **This does not contradict v216's nonzero degree on *fixed marked arcs*:** those marked products are not invariant under the full cyclic relabeling operation.

In a single exact convex-quadrilateral test, v216's independent boundary-degree algorithm gives \((-1,+1,-1,+1)\) for the four cyclically rotated lists of preselected vertex arcs, confirming the sign rule without using the new code's test-map linearization.

**Caution on parity.** Counting all *labeled* roots modulo two also loses information: four cyclic labelings of one square give \(4\equiv0\pmod2\). A count of **square orbits on the quotient** can instead give \(1\pmod2\) in an example, but the quotient is nonorientable and no new all-Jordan mod-2 continuation theorem has been proved here. A correct global construction would need a quotient/equivariant or local-coefficient formulation with collision control.

## 5. Six cyclic-order sectors, not one

The configuration space of *all* four distinct, **labeled** points on the parameter circle has \((4-1)!=6\) connected cyclic-order components. Recording each order from label 0, they are all permutations \((0,*,*,*)\). Exactly **two** components correspond to the labels appearing in geometric cyclic square order along the curve, either forward \((0,1,2,3)\) or backward \((0,3,2,1)\). These are often called *gracefully* inscribed configurations.

The other **four** components represent non-graceful placements of square vertices along the Jordan parameter order. The standard unrestricted Toeplitz conjecture does not require gracefulness. Its definition allows these sectors; a reduction to one or two requires proof. This issue is not resolved by local source-arc degree or by the four-periodic-inscription theorems.

## 6. Collision strata: where a square can shrink

Allow nonnegative cyclic gaps summing to one, giving the naive compactification \(S^1\times\overline\Delta_3\). Every **proper boundary face** has at least one gap zero, so at least two adjacent curve-parameter points coincide. For an embedded curve, any zero of \(F_\sigma\) with two coinciding points must have *zero* side and therefore **all four** points coincide: a nondegenerate square has four distinct vertices.

The all-coincident parameter configurations occur precisely at the **four vertices** of the closed gap simplex, where exactly one gap is one and the other three are zero. Hence

\[
F_\sigma^{-1}(0)\cap\partial(S^1\times\overline\Delta_3)
\ \subseteq\ S^1\times\{\text{four gap vertices}\},
\]

and equality holds for every continuous embedding: each such vertex maps to four copies of \(\gamma(t)\), a collapsed square. In particular, **every** naively compactified cyclic sector has unavoidable zeroes on its boundary. One cannot simply use the ordinary unmodified boundary degree from v216 on the whole closed domain.

This compactification is a coarse *gap simplex*, not the Fulton–MacPherson/cyclohedral blow-up. It records no relative collision scales or tangent directions. More sophisticated compactifications (Vrećica–Živaljević; Cantarella–Denne–McCleary) treat those distinctions explicitly, but no boundary-free invariant for arbitrary rough Jordan curves is asserted here.

### Why smooth curves are different (standard local cone argument)

For a fixed regular \(C^1\) embedding, \(\gamma'(t)\ne0\) and its derivative is continuous. Choose a parameter neighborhood on which \(\|\gamma'(s)-\gamma'(t)\|<\|\gamma'(t)\|/4\). Every nontrivial chord of parameters in that neighborhood has a longitudinal component at least \(3|\gamma'(t)||\Delta t|/4\) and transverse component at most \(|\gamma'(t)||\Delta t|/4\) relative to \(\gamma'(t)\). All its *unoriented chord lines* therefore lie within \(\arctan(1/3)<\pi/4\) of the tangent line. Two such chord lines cannot be perpendicular. Thus no square has all four vertices sufficiently near that parameter value. By compactness and uniform continuity of the inverse embedding, there is a positive minimum side length of any square on that **fixed** \(C^1\) regular curve.

This excludes collapse on a fixed regular curve, **not uniformly along a sequence of rough approximations**. A sequence of smooth curves approximating an arbitrary Jordan curve in \(C^0\) can lose derivative control and make the permitted size threshold tend to zero. That is precisely the gap preventing a naive compactness proof for the general conjecture.

## 7. Independent exact polygon tests and scientific boundaries

A deliberately narrow exact polygon auditor enumerates four **distinct open polygon edges**, a square orientation, and all four square-vertex labels. Each choice reduces to a rational \(4\times4\) affine linear system. It retains only isolated nondegenerate solutions with all four parameters strictly inside their respective edges. For the rational convex quadrilateral

\[
(0,0),\ (2,0),\ (3,3),\ (0,2),
\]

it finds **one geometric square orbit** with **eight labeled descriptions**: four cyclic rotations in each of two geometric orientations. The signs in each fourfold cyclic orbit alternate and sum to zero. The labeled mod-2 count is zero, while the illustrative geometric-orbit count is odd.

**A negative or zero output from this finite auditor is not an absence theorem:** it omits squares with vertices on the polygon's vertices, multiple points on one edge, positive-dimensional families, and singular affine charts. The universal square-peg conjecture is *not* decided by this code. In particular, even a simple square polygon can evade the distinct-open-edge diagnostic because its obvious square lies at edge endpoints.

## 8. What v219 should attempt

1. Specify the correct **orientation local system / equivariant target bundle** for the free \(C_4\) quotient, rather than counting labeled zeroes with ordinary integer signs.
2. Identify the **collision link** in a compactification that includes relative scales. Check whether a proposed mod-2 orbit count can change through zero-sized squares under mere \(C^0\) approximation; the v211 tangent-escape and v210 coil families remain adversarial inputs.
3. Treat all six cyclic-order components and distinguish graceful versus non-graceful configurations; do not silently assert a universal graceful-square theorem.
4. Compare with existing cyclohedral/Fulton–MacPherson configuration-space approaches before claiming a new topological invariant. Stop if the new object is equivalent to a known smooth-curve parity argument without a new noncollapse theorem.

The **main lesson** is a no-go for the simplest global signed-degree sum, together with a sharp description of where degenerate boundary zeros arise. It narrows the search for a usable global invariant but does not provide one.

## 9. Source references

- Siniša Vrećica and Rade Živaljević, *Fulton–MacPherson compactification, cyclohedra, and the polygonal pegs problem*, arXiv:0810.1439, https://arxiv.org/abs/0810.1439. Existing configuration-space and cyclohedral framework.
- Jason Cantarella, Elizabeth Denne, and John McCleary, *Transversality for Configuration Spaces and the “Square-Peg” Theorem*, arXiv:1402.6174v2, https://arxiv.org/abs/1402.6174. Smooth-transversality and generic odd-square count background.
- Benjamin Matschke, *On the Square Peg Problem and Its Relatives*, arXiv:1001.0186v2, https://arxiv.org/abs/1001.0186. Generic and open-class conclusions; no universal-Jordan deduction here.
- Internal research: [v216 boundary-only degree](v216-boundary-degree-global-arc-selection.md), [v217 multiscale conservation](v217-multiscale-relative-degree-ledger.md), and archived v210–v211 counterexamples.

Reproduce: `python -m unittest discover -s tests -v` and `python code/configuration_symmetry_v218.py`. No formal proof-assistant certification or Floer homology implementation is claimed.