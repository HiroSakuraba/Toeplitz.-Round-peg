# v260: A rigorous signed-degree gluing law through an actual positive-square pair death

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact, rationally certified, local four-dimensional Brouwer-degree cancellation for the **genuine noncollinear polygon-vertex square-pair death of v254**. The key new result is **uniform outer-boundary nonvanishing for the whole real interval containing the irrational collision**, including an explicit continuous perturbation budget. This makes the degree additivity and homotopy conservation rigorous, rather than inferring it from two sampled square counts.

**Most important conclusion:** The two disappearing squares have **opposite local indices**, \(-1\) and \(+1\). Their sum is zero, so the pair may disappear without violating degree. The protected local degree is **zero**, and therefore does **not** force an inscribed square to exist. This is a valuable obstruction to naive degree-based attacks on Toeplitz, not a general proof of a nonzero global square index.

## 1. A general, correct degree-gluing statement

Let \(U\subset\mathbf R^4\) be a bounded open parameter rectangle, and let \(F_t:\overline U\to\mathbf R^4\) be continuous in both its parameter \(s\) and homotopy time \(t\in[a,b]\).

Assume

\[
\boxed{0\notin F_t(\partial U)\quad \text{for ALL }t\in[a,b].}
\]

Then the integer Brouwer degree \(\deg(F_t,U,0)\) is well-defined and independent of \(t\) by homotopy invariance.

Now suppose \(U=U_-\cup\Sigma\cup U_+\), where \(\Sigma\) is a three-dimensional interior interface and \(U_\pm\) are the two open halfrectangles on either side. At a regular time for which \(F_t\) is nonzero on \(\Sigma\), degree additivity gives

\[
\boxed{\deg(F_t,U,0)
=\deg(F_t,U_-,0)+\deg(F_t,U_+,0).}
\]

If in each half there are isolated differentiable roots \(s_k\) with nonsingular Jacobian, their contributions are \(\operatorname{sgn}\det D_sF_t(s_k)\). The orientations here use **one shared ordered set of boundary parameters and the same ordering of the four square-defect equations**.

At a singular event where a zero lies **on the shared interior interface** \(\Sigma\), the separate half-domain degrees need **not** be defined: each half-domain has a boundary zero. But the degree of the **whole** \(U\) remains defined if the outer boundary is zero-free.

This is standard additivity and homotopy invariance of the Brouwer degree, not a novel general topological theorem. The new contribution of v260 is to verify *all* these hypotheses by exact rational polynomial signs for a real polygonal square-pair death.

## 2. The ACTUAL four-dimensional square-defect across the v254 corner

Reuse the strictly convex affine pentagon \(P_0(t),\dots,P_4(t)\) from [v254](v254-genuine-square-pair-death.md). Only \(P_0(t)\) moves. Its complete Jordan property is proved for all \(3/5\le t\le13/20\) by the existing exact v254 certificate.

Focus on the rational isolating interval for the **one algebraic collision time** \(\alpha\):

\[
a=\frac{64722}{100000},\qquad
b=\frac{64723}{100000},\qquad
a<\alpha<b.
\]

Let the 4D parameter \(s=(s_0,h,s_2,s_3)\) specify actual points on four ordered **polygon boundary arcs**, with the **second** arc passing continuously across the genuinely noncollinear vertex \(P_2(t)\):

\[
\begin{aligned}
p_0(t,s_0)&=P_0(t)+s_0(P_1(t)-P_0(t)),\\
p_1(t,h)&=
\begin{cases}
P_2(t)+h(P_2(t)-P_1(t)),&h\le0,\\
P_2(t)+h(P_3(t)-P_2(t)),&h\ge0,
\end{cases}\\
p_2(t,s_2)&=P_3(t)+s_2(P_4(t)-P_3(t)),\\
p_3(t,s_3)&=P_4(t)+s_3(P_0(t)-P_4(t)).
\end{aligned}
\]

At \(h=0\), both formulas for \(p_1\) equal the same polygon vertex, so this map is **continuous but piecewise affine** across the genuine corner.

Define the *same* direct square-defect map as v259, using \(J(x,y)=(-y,x)\):

\[
F_t(s)=
\begin{pmatrix}
p_2-p_1-J(p_1-p_0)\\
p_3-p_0-J(p_1-p_0)
\end{pmatrix}\in\mathbf R^4.
\]

A zero gives a counterclockwise square on the four actual boundary arcs. The v254 positive side² \(>1\) result proves that the collision itself is nondegenerate.

The negative half \(h<0\) is precisely v254 square assignment \(A=(0,1,3,4)\), and the positive half \(h>0\) is \(B=(0,2,3,4)\). Unlike a method that separately follows two Cramer square coordinates, **both are now embedded in one genuine source-parameter rectangle**, with its two halves sharing exactly the same \(h=0\) interface.

## 3. The exact outer box and ALL-TIME nonvanishing proof

The new [exact verifier](../code/square_relative_degree_v260.py) uses one fixed rational parameter rectangle

\[
\boxed{
\begin{aligned}
9/10\le s_0\le93/100,\\
-1/1000\le h\le1/1000,\\
43/100\le s_2\le46/100,\\
73/100\le s_3\le78/100.
\end{aligned}}
\]

Every value in this rectangle parameterizes actual interior portions of the intended polygon boundary edges, and both sides of its \(h\)-coordinate traverse the two adjacent genuine noncollinear edges. In particular, **the interface \(h=0\) is inside the box**, not on its outer boundary.

For a fixed real \(t\), restrict the direct four-dimensional map \(F_t\) to each halfbox. It is affine in \(s\) there:

\[
F_t(s)=M_\pm(t)s+c_\pm(t).
\]

The exact moving-polygon chart engine from v224 independently computes the determinant \(D_\pm(t)=\det M_\pm(t)\) and the unique extended Cramer solutions \(r_\pm(t)\) whenever \(D_\pm(t)\ne0\).

**Uniform exact polynomial checks** over the *entire real interval* \(t\in[a,b]\) prove all the following:

1. The determinants have strictly nonzero, opposite signs. In fact,
   \[
   D_-(t)<-1/10,\qquad D_+(t)>1/10
   \]
   for every real time. The determinant polynomial is independently compared to a separate 24-term exact permutation determinant evaluated at both time endpoints.
2. All 16 moving matrix entries satisfy \(|(M_\pm)_{ij}(t)|\le2\). Every entry is affine in \(t\), so it suffices to check both rational time endpoints.
3. Every coordinate of both **extended** Cramer roots \(r_\pm(t)\) lies **strictly inside** its corresponding outer-box interval by at least the positive exact gap
   \[
   \boxed{g=1/10000}.
   \]
   These are **eight strict rational-polynomial inequalities per halfbox**, verified on the whole real interval using the exact Bernstein coefficient procedure from v249. The extended roots are checked even when their \(h\)-coordinate has the *wrong sign* for their halfbox; in that case, they are not feasible geometric square roots, but the bound is still useful for excluding boundary zeros.

**Consequently, \(F_t\) never vanishes anywhere on the outer box boundary, for ANY real \(t\in[a,b]\).**

Proof: take a point \(s\) on the outer boundary. Depending on the sign of its \(h\)-coordinate, the map equals one of the affine functions \(M_\pm(s-r_\pm)\). If that function vanished, invertibility would force \(s=r_\pm(t)\), contradicting the exact strict distance of the branch root from all outer faces. At \(h=0\), the two affine formulas agree and the same argument applies.

This is the crucial proof obligation: it applies **even at the irrational collision** \(\alpha\), where the map is nonsmooth, rather than relying on regular-time sampling.

## 4. A quantitative positive exterior residual bound

The preceding argument is strengthened to a **strict explicit residual bound**.

With every \(4\times4\) matrix entry bounded by \(B=2\), every three-by-three cofactor has magnitude at most \(6B^3\) by the determinant's six-term Leibniz formula. Every determinant has magnitude greater than \(d=1/10\). The matrix inverse therefore satisfies

\[
\|M_\pm(t)^{-1}\|_\infty <
\frac{4\cdot6B^3}{d}
=\frac{24\cdot8}{1/10}
=\boxed{1920}.
\]

For any point on the *outer* boundary, some coordinate differs from either root \(r_\pm(t)\) by at least \(g=1/10000\). Thus

\[
\|F_t(s)\|_\infty
\ge\frac{\|s-r_\pm(t)\|_\infty}{\|M_\pm(t)^{-1}\|_\infty}
>\boxed{\frac1{19200000}}.
\]

No inverse condition number is assumed without proof; the uniform bound follows from exact determinant and entry certificates.

Now perturb each of the four actual boundary arcs by an **arbitrary continuous** function with each Cartesian coordinate altered by at most

\[
\boxed{\varepsilon=\frac1{10^9}}.
\]

Each square-defect component is the signed sum of at most four input Cartesian coordinates, so its change is at most \(4\varepsilon=4/10^9\), strictly less than \(1/19200000\).

Therefore **the perturbed four-arc square-defect map also has no zero on the OUTER boundary**, and its integer total degree remains **zero** through the entire time interval. No differentiability or fixed oscillation frequency is required.

**Important distinction:** Preserving a **zero** degree does **not** establish that an exact square survives such an arbitrary perturbation. The result is robust *degree bookkeeping*, not robust existence.

### Noncollapse in the entire source box

The exact vertex-extremum method also proves a positive side-size floor independent of square-branch identities. Across the whole time/parameter box, the first-to-second boundary vector \(u=p_1-p_0\) satisfies

\[
u_x<-74/100,\qquad u_y<-79/100.
\]

This follows by evaluating a multi-affine expression on the rational endpoint combinations of \(t,s_0,h\), split at \(h=0\). Hence any square root, even under the allowed coordinatewise perturbations, has

\[
\operatorname{side}^2\ge
\left(\frac{74}{100}-2\varepsilon\right)^2+
\left(\frac{79}{100}-2\varepsilon\right)^2>1.
\]

So the observed disappearance is **not** a collision of four square corners toward a point.

## 5. Exact matching at the algebraic inner interface

v254 proved that the collision time \(\alpha\) is the unique real zero in the interval of

\[
Q(t)=72468080000000t^2
-496738335230000t
+291143052416293.
\]

The v260 program re-expresses this polynomial in the normalized time \(\tau=(t-a)/(b-a)\), preserving its **exact rational coefficients**.

For the negative-side Cramer solution, its vertex-crossing coordinate is

\[
h_-(t)=u_{\text{edge 1}}(t)-1;
\]

for the positive side it is \(h_+(t)=u_{\text{edge 2}}(t)\).

The verifier checks that the **numerators** of both \(h_-,h_+\) are exactly divisible by \(Q(a+(b-a)\tau)\), with quotient having a strictly positive sign throughout the interval. Since the denominators have opposite signs, the roots obey

\[
\begin{array}{c|ccc}
 & t<\alpha & t=\alpha & t>\alpha\\
h_-&<0&0&>0\\
h_+&>0&0&<0.
\end{array}
\]

Thus both roots are in their **correct halfboxes before** the event, both lie on the **shared inner interface at** the event, and each is on the **wrong halfbox afterward**.

Moreover, for all four source parameter coordinates \(j\), the verifier checks the rational polynomial identity

\[
Q(a+(b-a)\tau)\ \bigm|\
N_{-,j}D_+-N_{+,j}D_-.
\]

As both denominators stay nonzero, the four parameters of the two roots **coincide at exactly the algebraic event**. The pre-existing v254 certificate independently establishes the same fact for all four geometric square coordinates and proves that the common second square vertex is exactly the genuine polygon vertex \(P_2(\alpha)\).

At \(t=\alpha\), the two sources **represent one geometric square**. The root is interior to the **whole** outer box but lies on the **boundary of each individual halfbox**.

## 6. A rigorous integer degree ledger

The actual signs of the four-dimensional Jacobian determinants, in the common parameter and residual orientation, are

\[
\boxed{\operatorname{sgn}\det M_-=-1,\qquad
\operatorname{sgn}\det M_+=+1.}
\]

**Before the collision \(t<\alpha\):** Both branches have a single regular root strictly inside their own halfbox, with one shared outer domain and a root-free internal interface. Additivity yields

\[
\deg(F_t,U,0)=
\underbrace{(-1)}_{h<0}+
\underbrace{(+1)}_{h>0}
=\boxed0.
\]

**At \(t=\alpha\):** The two roots coincide as one positive-size square on the interface \(h=0\). Individual halfbox degrees are undefined because of the interface root, but the **outer box boundary remains zero-free**; homotopy invariance gives

\[
\boxed{\deg(F_\alpha,U,0)=0.}
\]

**After the collision \(t>\alpha\):** The only candidate roots in the two affine extensions each violate the required sign of \(h\); no square-defect root exists in the local box. Hence

\[
\boxed{\deg(F_t,U,0)=0.}
\]

The certified local degree is therefore **zero for every real time**, including the irrational collision.

| Real time | Actual local geometric squares in the source box | Two-halfbox degree | Whole-box degree |
|---|---:|---|---:|
| \(t<\alpha\) | **2** | \(-1,+1\) | **0** |
| \(t=\alpha\) | **1** (shared kink square) | individual degrees undefined | **0** |
| \(t>\alpha\) | **0** | \(0,0\) | **0** |

This is a genuine **relative-degree gluing law**. The presence of a nonsmooth edge-chart transition does not break the topological invariant, but the algebraic signs cancel and the invariant **does not protect the disappearing pair**.

## 7. Independent geometric cross-check and the surviving square

The original independent v223 complete **unfiltered** geometric-square census is rerun at both rational interval endpoints:

- At \(t=64722/100000<\alpha\), the **whole pentagon** has exactly **3** positive geometric inscribed squares.
- At \(t=64723/100000>\alpha\), the **whole pentagon** has exactly **1**.

The local box accounts for the **two** disappearing squares. v254 independently proved that a **third** positive-size square remains on edge assignment \((0,1,2,4)\) throughout the entire deformation. The v260 degree-zero result does not claim to protect this persistent square; it lies outside the selected pair-death local box.

This is distinct from the v259 example, where the local degree was **+1** and therefore guaranteed a square even under controlled rough perturbations. The contrast is instructive:

\[
\boxed{\text{Nonzero local degree protects existence; zero local degree permits pair cancellation.}}
\]

The difficulty for a general Toeplitz proof is finding a **globally relevant, nonzero, size-controlled degree**, not merely constructing a local degree around a known square.

## 8. Reproducibility and exactness

The new standard-library-only program replays exact v254 Jordan and collision certificates; independently computes its oriented affine Jacobian determinants, certifies 16 outer-root coordinate bounds by rational Bernstein signs, and verifies all four parameter-matching identities via exact polynomial division. The inverse-norm and boundary-margin calculations are exact rational inequalities.

~~~sh
python code/square_relative_degree_v260.py
python -m unittest -v tests.test_square_relative_degree_v260
~~~

The regression suite covers the before/at/after algebraic event classification, the nested normal-root box, full-time nonvanishing, opposite oriented Jacobians, exact gluing identities, positive square-side bound, independent all-orbit censuses, arbitrarily continuous perturbation budget, and negative controls for altered box, fake geometry and overly large perturbations.

### Scope and next milestone

This theorem proves a **local degree on one four-arc parameter box**, not an invariant on the entire space of four ordered boundary points. In particular, it does not enumerate all positive squares at every irrational parameter value; the independently *complete whole-polygon counts* are proved at two rational snapshots, while the **local box classification** holds at *every real time*.

The degree remains **zero** through this pair death. This is entirely consistent with the classical generic odd-parity phenomenon for simple polygons. It does **not** show that the global positive square count can become zero.

The **v261** priority is to choose and orient a **finite atlas of disjoint isolating boxes** around *all* positive-size squares over a moving-polygon family, certify their external boundary nonvanishing, and determine whether the sum of local degrees produces a useful nonzero global class. Genuine rank drops, vertex gluing, repeated-edge configurations and collapse at tiny scales remain major potential obstructions.

**No universal positive Big Peg constant or unrestricted Toeplitz theorem is claimed.**

## References

- R. Schwartz, *An Integration Approach to the Toeplitz Square Peg Problem*, *Forum of Mathematics, Sigma* 6 (2018), for rigorous geometric-analytic methods on specified Jordan-curve classes.
- D. Salamon, *Lectures on Differential Topology*, chapter on Brouwer degree and oriented local root counts (available at https://people.math.ethz.ch/~salamon/PREPRINTS/difftop.pdf). The degree additivity and homotopy properties used here are standard.
- B. Matschke, *A Survey on the Square Peg Problem*, *Notices of the AMS* 61 (2014).
- [v254 exact positive square-pair death](v254-genuine-square-pair-death.md), [v259 exact local degree across a genuine polygon kink](v259-piecewise-arc-local-degree.md).
