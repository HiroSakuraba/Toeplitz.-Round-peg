# v248: A quantitative square-persistence theorem for rough boundary arcs

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact rational structural lemma, explicit examples and negative controls. **This does not prove the unrestricted Big Peg conjecture or Toeplitz.** It does establish an explicit positive square-size bound for an infinite-dimensional family of arbitrarily rough *continuous* Jordan curves.

## Why shift from optimizing decagons?

v244–v247 certified increasingly smaller *upper bounds* on any possible universal Big Peg constant by constructing individual polygons. Those records do not tell us whether a single positive lower bound works for **all** polygons. The missing ingredient is **noncollapse under increasingly irregular boundary perturbations**.

The v248 question is narrower and provable: given a particular nonsingular square of a polygon, can we certify a nonzero square-size bound for every sufficiently small continuous deformation of the four boundary arcs that support its corners, without assuming differentiability, root-finding convergence, or a bound on the oscillation frequency?

Yes, under an explicit graph and transversality hypothesis. The result is a *local* quantitative topological stability theorem, not a global square peg theorem.

## Theorem (normal-graph persistence)

Let P be a Jordan polygon with four **distinct ordered edges** e_k(t)=A_k+t E_k, 0<=t<=1 (k=0,...,3). Write a square by a first corner p=(x,y) and side vector S=(u,v), giving four corners

\[
C_0=p,\quad C_1=p+S,\quad C_2=p+S+JS,\quad C_3=p+JS,
\quad J(u,v)=(-v,u).
\]

The four line-incidence conditions

\[
\det(E_k,C_k-A_k)=0
\]

give a **4-by-4 rational linear system** M z=b in z=(x,y,u,v). Suppose (i) M is invertible, (ii) its unique square solution z_* has all edge coordinates t_k strictly in (0,1), and (iii) S_* is nonzero.

Now replace each selected edge by an **arbitrary continuous normal graph**

\[
\Gamma_k(t)=A_k+tE_k+g_k(t)\,\frac{JE_k}{\|E_k\|},
\qquad |g_k(t)|\le\varepsilon.
\]

The four graph functions may be nowhere differentiable. Other parts of the Jordan curve may be arbitrary, provided it contains all four specified arcs.

Write B_k=|E_{kx}|+|E_{ky}|, which is a rational upper bound for \|E_k\|. Let

\[
K=\max_j\sum_{k=0}^3|(M^{-1})_{jk}| B_k,\qquad
\rho=2K\varepsilon.
\]

Define L_k as the exact sum of absolute values of the four coefficients in the linear expression for the tangential coordinate

\[
t_k(z)=\frac{E_k\cdot(C_k(z)-A_k)}{\|E_k\|^2}.
\]

If, for every k,

\[
\rho L_k<\min(t_k(z_*),1-t_k(z_*))
\]

and

\[
s_{\min}^2=
\bigl(\max\{|u_*|-\rho,0\}\bigr)^2+
\bigl(\max\{|v_*|-\rho,0\}\bigr)^2>0,
\]

then **every** choice of the four continuous functions g_k within the given uniform amplitude bound has an inscribed square with side **at least s_min**.

### Proof

The signed normal line equation for the perturbed graph reads

\[
M z=b+h(z),
\quad h_k(z)=\|E_k\|g_k(t_k(z)).
\]

On the closed cube Q={z: \|z-z_*\|_\infty<=rho}, the tangential-margin inequalities ensure every t_k(z) lies strictly inside (0,1), so each continuous g_k is evaluated in its domain. The map

\[
T(z)=z_*+M^{-1}h(z)
\]

is continuous. Since \|E_k\|<=B_k and |g_k|<=epsilon, the exact bound on the inverse gives

\[
\|T(z)-z_*\|_\infty\le\varepsilon K=\rho/2<\rho.
\]

Thus T maps the entire compact convex cube into itself. **Brouwer's fixed-point theorem** produces a z in Q with T(z)=z. That equation is exactly the condition that all four square corners lie on the four perturbed graph arcs: the tangential and normal coordinates jointly reconstruct each corner. The side components differ from u_*,v_* by at most rho, so the exact strictly positive s_min lower bound proves the square is nondegenerate. QED.

**Trust base:** elementary finite-dimensional linear algebra, exact rational inequalities and Brouwer's fixed-point theorem. No numerical root solver, implicit differentiability or topological degree calculation is assumed.

## A quantitative family of rough curves

Take the following explicit rational octagon in cyclic order (coordinates divided by 20):

| vertex | X | Y |
|---:|---:|---:|
| 0 | 39 | 21 |
| 1 | 18 | 40 |
| 2 | -17 | 39 |
| 3 | -38 | 20 |
| 4 | -43 | -19 |
| 5 | -20 | -38 |
| 6 | 21 | -38 |
| 7 | 39 | -23 |

The four selected edges are indices (0,2,4,6). The exact root is

\[
z_*=(21561/15136,\;808131/529760,\;-44719/15136,\;-20205/105952),
\]

and the four exact edge parameters are

\[
(13257/26488,\;8577/13244,\;18673/26488,\;8317/13244).
\]

For **every** continuous normal perturbation on those four edges with Euclidean normal displacement <= **1/512**, v248 checks that the certificate cube has exact radius

\[
\rho=367709/3390464,
\]

and the resulting square has squared side at least

\[
\boxed{s_{\min}^2=
46593827703305/5747623067648> (27/10)^2}.
\]

Therefore **every** Jordan curve contained in the disk of radius 3 and containing these four perturbed arcs satisfies the **quantitative restricted-class Big Peg bound**

\[
\boxed{J=s_{\max}/r_{\mathrm{in}}>9/10.}
\]

Indeed, the whole bounded Jordan domain lies in the radius-3 disk, so r_in<=3, while the square side >27/10. This is an explicit positive bound for an infinite-dimensional family that includes continuous, nonsmooth, arbitrarily high-frequency or fractal normal graph perturbations. It is **not** a universal bound for arbitrary Jordan curves, because requiring these four arcs is a substantial structural hypothesis.

## Fragility warning from the sharp Big Peg counterexamples

The four-edge square of Pochuev's original decagon is so close to an edge endpoint that this conservative normal-graph certificate works at epsilon=1/10^6 but cannot automatically extrapolate to coarse perturbations.

For the optimized v247 polygon the maximizing square is **much closer to a vertex**: its smallest edge-parameter margin is less than 1/10^7. Our explicit transversality certificate works at epsilon=1/10^9 but fails at epsilon=1/10^8. A failed certificate **does not mean the perturbed curve has no square**: it only means this specific four-edge Brouwer box cannot guarantee one uniformly.

This illustrates the true difficulty in turning local square persistence into a universal Big Peg proof. As polygon geometry degenerates, the available stability margin may tend to zero. To solve Toeplitz through Big Peg, we need a theorem that replaces these polygon-dependent constants with **one positive inradius-scaled constant independent of edge count and near-vertex degeneracy**.

## Reproducibility

Pure Python standard library, exact fractions:

    python code/big_peg_persistence_v248.py
    python -m unittest -v tests.test_big_peg_persistence_v248

The tests check the exact M inverse/root, strict tangential margins, the universal perturbation bound, explicit constant normal-graph examples, singular matrices and bad edge assignments, rejection of inadmissible epsilons, and both near-vertex Pochuev/v247 cases. These are **universal graph-amplitude certificates**, not random-sample tests.

## Prior research (primary sources)

- B. Matschke, *A Survey on the Square Peg Problem* (Notices AMS, 2014): https://pure.mpg.de/rest/items/item_3120610/component/file_3120611/content . Explicitly identifies a positive universal quantitative square-size bound as a route from polygonal/smooth approximation to Toeplitz.
- R. E. Schwartz, Big Peg Conjecture 3.5 (Bull. AMS, 2022); see David Bernier's direct statement and 2026 sharp-constant challenge: https://dbernier.ca/2026/06/16/the-sharp-big-peg-conjecture/ .
- Grisha Pochuev's original exact ten-sided sharp-constant counterexample: https://github.com/Grisha-Pochuev/sharp-big-peg-counterexample .
- Greene and Lobb, *Square pegs between two graphs*, Comment. Math. Helv. (2026), DOI 10.4171/CMH/619. Their existence theorem for a more general graph class is distinct from the explicit local quantitative inequality proven here: https://doi.org/10.4171/cmh/619 .

## Recommended next milestone

Use the robust graph certificate as a **localization tool** in a covering argument. Identify shape invariants that give lower bounds on edge-parameter margins and inverse-incidence condition numbers, possibly under controlled turning or chord-arc hypotheses. Classify what happens when these invariants degenerate: near-vertex collisions, nearly singular incidence matrices, and oscillations at arbitrarily small scales.

Only a mechanism that controls all those degeneracies, with a constant tied to inradius but not polygon complexity, could close the general Big Peg conjecture and hence Toeplitz.
