# v219 — A symmetry-aware quotient, a vanishing absolute obstruction, and the collision-link problem

**Research program:** Benjamin John Schulz  
**Date:** 18 November 2026  
**Scientific status:** Elementary quotient and bundle calculations, a conditional regular-curve collision-link estimate, an explicit collapsing-square polygon family, and exact-rational regression tests. **Not a proof of the unrestricted Toeplitz square-peg conjecture.** This is not a constructed relative Euler class or Floer homology.

## 1. Answer to the v218 question

v218 showed that the naive **integer signed count of *labeled*** squares vanishes in each complete cyclic orbit: cyclic relabeling reverses the source orientation and preserves the target orientation. That does *not* rule out the well-known **mod-2 orbit parity** available under smooth-curve transversality and compactification hypotheses (Cantarella–Denne–McCleary; Matschke). It also does *not* automatically provide such a parity for arbitrary Jordan curves.

Here we test a more basic possibility: might the cyclic quotient, with the natural square-residual target representation, force zeros by an *absolute* mod-2 Euler/Stiefel–Whitney obstruction? **No.** On **one** cyclic-order component it deformation-retracts equivariantly onto a circle, and we construct a globally nonzero equivariant section *explicitly*. Thus the top-degree absolute bundle obstruction vanishes. Any successful global proof must use *relative boundary data*, a more sophisticated class, or a geometric non-collapse estimate. This conclusion does not exclude a relative/equivariant parity theorem once correct collision data are specified.

## 2. The free C4 action and an equivariant center coordinate

Let one positively ordered configuration component be

\[
M=S^1\times \Delta_3^\circ,\quad
(t;g_0,g_1,g_2,g_3),\quad g_i>0,\quad\sum_i g_i=1.
\]

The point parameters are \(t,t+g_0,t+g_0+g_1,t+g_0+g_1+g_2\) modulo one. Cyclic relabeling is

\[
T(t;g_0,g_1,g_2,g_3)
=(t+g_0;g_1,g_2,g_3,g_0).
\]

This also acts on the naive **closed** simplex \(\bar M=S^1\times\Delta_3\) allowing zero gaps. Define

\[
\boxed{c=t+\frac{3g_0+2g_1+g_2}{4}\pmod1.}
\]

Direct substitution gives \(c(Tx)=c(x)+1/4\pmod1\). Keeping \(c\) fixed while interpolating gaps,

\[
g_i(s)=(1-s)g_i+s/4,\qquad 0\le s\le1,
\]

and reconstructing \(t(s)=c-(3g_0(s)+2g_1(s)+g_2(s))/4\), gives an **equivariant strong deformation retraction** of both \(M\) and \(\bar M\) onto the equal-gap circle.

The action on that circle is a quarter-turn \(c\mapsto c+1/4\); it is free. Therefore

\[
\boxed{M/C_4\simeq S^1,\qquad\bar M/C_4\simeq S^1,}
\]

where \(\simeq\) means homotopy equivalence. The action remains free on the closed simplex: a fixed configuration for \(T^k\), \(k=1,2,3\), would have periodic gaps and center shifted by \(k/4\), a contradiction. Thus the statement does not assume free action disappears upon allowing collisions.

**Scope:** This is only the *single positive cyclic-order component*. All labeled four-point configurations have six cyclic-order sectors (v218). Their interactions, especially the nongraceful ones, are not silently identified with this single quotient.

## 3. Residual representation and an explicit nowhere-zero equivariant section

Use v218's oriented-square map \(F_\sigma\) and its exact target action \(F_\sigma(Tx)=B_\sigma F_\sigma(x)\) for \(\sigma=\pm1\), where \(B_\sigma^4=I\), \(\det B_\sigma=+1\). We form the rank-four associated real bundle

\[
E_\sigma=(M\times\mathbb R^4)/\{(x,v)\sim(Tx,B_\sigma v)\}\longrightarrow M/C_4.
\]

A nonzero section can be written down without any knowledge of a Jordan curve. Set

\[
u=(1,0,1,0),\qquad v_\sigma=B_\sigma u=(0,-\sigma,0,-\sigma).
\]

The exact matrices give \(B_\sigma v_\sigma=-u\). Define on the covering space

\[
\boxed{s_\sigma(x)=\cos(2\pi c(x))u+\sin(2\pi c(x))v_\sigma.}
\]

Then

\[
s_\sigma(Tx)=s_\sigma(c+1/4)=B_\sigma s_\sigma(x),
\quad \lVert s_\sigma(x)\rVert_2^2=2.
\]

Hence \(s_\sigma\) descends to a **nowhere-zero section of \(E_\sigma\)**, even across the naive collision boundary. Since \(E_\sigma\) is oriented (\(\det B_\sigma=1\)) and its base retracts onto a circle, its absolute integral Euler class is zero and its mod-2 top Stiefel–Whitney class \(w_4(E_\sigma)\) is zero. Indeed \(H^4(M/C_4;\mathbb Z/2)=0\) purely on dimension/homotopy grounds.

**No-go conclusion:** A generic argument asserting that *equivariance or the absolute top Euler class alone* must force an inscribed square in this configuration component is false. This is a modest but rigorous obstruction to an oversimplified v219 strategy. It does **not** show that the actual square-residual section can be homotoped to this section while respecting collision boundary conditions. It cannot: the original residual vanishes at the all-coincident configurations on the naive boundary.

### Orientation local system

By v218, \(\det DT=-1\), while \(\det B_\sigma=+1\). The quotient configuration manifold is nonorientable even though the target bundle is oriented. Its fundamental local coefficient system changes sign under the generator. Working modulo two removes this sign, but still leaves the collision-boundary behavior unresolved. An ordinary global *labeled* integer degree is unsuitable; the presence of a mod-2 quotient count in generic smooth cases is not in dispute.

## 4. The normalized total-collision link for a regular C1 embedding

Consider the simplex vertex with \(g_3\to1\). Write

\[
(g_0,g_1,g_2,g_3)=(\varepsilon a,\varepsilon b,\varepsilon c,1-\varepsilon),
\quad a,b,c\ge0,\quad a+b+c=1.
\]

For a regular \(C^1\) embedding \(\gamma\) and \(v=\gamma'(t)\ne0\), the four points are, uniformly in \((a,b,c)\),

\[
\gamma(t),\quad\gamma(t)+\varepsilon av+o(\varepsilon),\quad
\gamma(t)+\varepsilon(a+b)v+o(\varepsilon),\quad
\gamma(t)+\varepsilon v+o(\varepsilon).
\]

Consequently

\[
\lim_{\varepsilon\downarrow0}\frac{F_\sigma}{\varepsilon}
=\left(bv-\sigma aJv,\;v-\sigma aJv\right).
\]

The **second** real two-vector alone satisfies

\[
\boxed{\lVert v-\sigma aJv\rVert_2^2=(1+a^2)\lVert v\rVert_2^2
\ge\lVert v\rVert_2^2>0.}
\]

This gives a nonvanishing normalized collision link. Other total-collision vertices follow by cyclic transformation of the residual with the invertible matrix \(B_\sigma\). The link is well behaved for every fixed regular \(C^1\) embedding, agreeing with known smooth-curve methods rather than improving on them.

### Quantitative regularity bound

For a fixed regular \(C^1\) embedding, let

\[
m=\min_t\lVert\gamma'(t)\rVert_2>0,
\qquad \omega(r)=\sup_{|s-t|\le r}\lVert\gamma'(s)-\gamma'(t)\rVert_2.
\]

The first-order error in the second residual block is at most \(2\varepsilon\omega(\varepsilon)\). If \(\omega(\varepsilon)\le m/4\), then

\[
\boxed{\lVert F_\sigma^{(2)}\rVert_2\ge (m/2)\varepsilon}
\]

uniformly in \(t,a,b,c\). This is the precise local non-collapse condition *for that fixed curve*. It is **not uniform over arbitrary C0-smooth approximations**, whose derivative lower bound and modulus need not remain controlled. A continuous Jordan curve has no canonical derivative from which to form the same link in general.

## 5. Exact C0 collapsing-square family: the microscopic notch

For each integer \(n\ge2\), put \(h=1/n\) and consider the simple closed polygon with vertices, in boundary order,

\[
(0,0),(1,0),(1,-h),(1+h,-h),(1+h,0),(2,0),(2,1),(0,1).
\]

This replaces a short segment on the lower side of the rectangle \([0,2]\times[0,1]\) with an outward rectangular detour. The four consecutive vertices

\[
(1,0),\quad(1,-h),\quad(1+h,-h),\quad(1+h,0)
\]

form an **exact nondegenerate Euclidean square of side \(h\)**. All polygons are embedded. Their boundary sets converge to the rectangle with Hausdorff distance at most \(h\): the detour is within \(h\) of the original bottom side; every point of the missing bottom subsegment is within \(h/2\) of a remaining vertical detour side. By assigning the detour an interval of parameter length \(O(h)\), one also obtains uniformly convergent embedding parameterizations.

Thus a sequence of Jordan polygons can converge in \(C^0\) to a fixed Jordan polygon while retaining a *chosen* square witness whose side length tends to zero:

\[
\boxed{s_n=1/n\longrightarrow0.}
\]

This **does not imply that the limit rectangle has no square**. It already has macroscopic inscribed squares. The example instead prevents using *mere existence of a square at each approximation level* to infer a nonzero-sized limiting square from that particular sequence. It exposes the missing selection/compactness requirement precisely.

All coordinates and side-square identities in this family are exact rationals; the checker validates polygon simplicity and square orientation for multiple values of \(n\). The geometric Hausdorff bound follows analytically, not from grid samples.

## 6. What is proved, what is conditional, and what is still open

| Claim | Status |
|---|---|
| Equivariant center and strong deformation retraction of one cyclic-order quotient to S1 | Proved algebraically |
| Free C4 action persists on naive closed gap simplex | Proved |
| Explicit nonvanishing equivariant target section and vanishing absolute top Euler/w4 obstruction | Proved for this associated bundle |
| Uniform nonzero normalized collision residual for each fixed regular C1 embedding | Proved conditionally on derivative regularity |
| Arbitrarily small selected square witnesses in simple polygons approaching a rectangle | Exact constructed family |
| A relative Euler class with a well-defined all-Jordan boundary value | **Not constructed** |
| Global mod-2 parity invariant under arbitrary rough Jordan regularization | **Not proved** |
| General Toeplitz square-peg conjecture | **Open** |

An absolute homotopy calculation is not a substitute for a *relative* obstruction. The associated target bundle has a nonzero section, but imposing physically meaningful boundary data on a collision blow-up may yield a nontrivial relative class. Establishing such data with insufficient regularity, while controlling the nongraceful sectors, is a distinct and hard problem.

### Next legitimate milestone (v220)

1. **Analyze the collision-link boundary class**, rather than the absolute bundle. On a regular \(C^1\) curve, calculate the normalized residual along each cyclohedral collision face and its equivariant clutching class. Compare with existing compactification proofs; avoid claiming a new parity theorem already known for smooth curves.
2. **Find a C0-robust boundary substitute** or a geometric scale floor derived from a *global* invariant, not from a single arbitrarily shrinking square. Test the v210 coil and v211 tangent-escape counterexamples plus this v219 notch.
3. **Track the five other cyclic-order components** and explain their relation to geometric square orbits. A proof for only one graceful sector is not automatically the universal Toeplitz result.
4. **If relative boundary continuity cannot be established for arbitrary Jordan curves, stop that lane honestly** and return to the Liouville-regular / Floer-action restrictions of the 2026 papers rather than silently assuming derivative control.

## 7. Code and sources

Run:

```bash
python -m unittest discover -s tests -v
python code/collision_quotient_v219.py
```

The exact-rational suite tests the group action, all four collision vertices, equivariant deformation retraction, the algebraic nowhere-zero section for both square orientations, the first-order residual for rational tangent data, and the shrinking-notch family. Two optional integration tests compare to the independently implemented v216 and v218 modules when checked out from the full repository. The Python program checks exact formulas; it is **not** Lean/Coq certification or a proof of a global invariant.

**Primary references, pre-existing work:**

- S. Vrećica and R. Živaljević, *Fulton–MacPherson compactification, cyclohedra, and the polygonal pegs problem* (2008), https://arxiv.org/abs/0810.1439 . Treats cyclohedra and configuration-space compactification.
- J. Cantarella, E. Denne, J. McCleary, *Transversality for Configuration Spaces and the “Square-Peg” Theorem* (2014/2021 revisions), https://arxiv.org/abs/1402.6174 . Proves a C1-dense smooth transversality/odd-square result under the stated hypotheses.
- B. Matschke, *On the Square Peg Problem and Its Relatives*, https://arxiv.org/abs/1001.0186 . Known parity and square-peg contexts; no newly claimed all-Jordan result.
- Internal starting point: [v218](v218-global-configuration-symmetry-and-collision-audit.md), [v217](v217-multiscale-relative-degree-ledger.md), archived v210–v211.
