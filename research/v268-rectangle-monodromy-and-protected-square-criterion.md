# v268: Rectangle-locus monodromy, a sharp Möbius no-go, and an exact odd-loop square certificate

**Research program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** A general **conditional topological forcing proposition**, a reusable **exact rational-function certificate** for continuous nondegenerate inscribed-rectangle loops, and two rigorously distinguished nontrivial loops on a **single fixed circle**. The first shows that two individual Möbius chord twists need **not** force a square; the second uses an **odd cyclic monodromy** to force a unique square at an irrational parameter. This is not an all-Jordan Toeplitz proof.

## 1. From secant spin to the *actual* rectangle locus

[v267](v267-secant-spin-rectangle-obstruction.md) proved that every continuous Jordan embedding has nonzero antipodal secants of odd winding and a nontrivial unoriented-chord Möbius line. It also showed that the square-perpendicularity equation lives in a **tensor sign line** of *two* diagonal chords and that the new diagonal residual's absolute rank-four obstruction still vanishes.

The missing hypothesis was concrete: on every Jordan curve, does the space of **actual inscribed nondegenerate rectangles** contain a closed component whose diagonal-perpendicularity sign line has nonzero first Stiefel–Whitney class?

v268 distinguishes two algebraically explicit closed rectangle loops with **different** monodromy. This proves the individual line twists are insufficient and isolates the *additional* topological property needed for a square-forcing route.

Throughout, a rectangle means **four distinct Euclidean corners** on one fixed Jordan curve. The corners are ordered in **geometric counterclockwise order**, not necessarily the same order as their parameters on an arbitrary, possibly nonconvex Jordan boundary. The geometric diagonal pairing is \((p_0,p_2)\mid(p_1,p_3)\); the other two pairings of curve-ordered labels remain relevant when trying to cover the *entire* rectangle locus of a nonconvex curve.

## 2. Conditional monodromy theorem for any Jordan curve

Let \(C:[0,1]\to\mathcal R(\gamma)\) be a continuous loop of **genuine nondegenerate geometric rectangles** inscribed in a fixed Jordan curve \(\gamma\). Suppose it admits a continuous choice of its four corners in geometric positive cyclic order

\[
p(t)=(p_0(t),p_1(t),p_2(t),p_3(t)).
\]

Because the geometric rectangle at \(t=1\) equals the one at \(t=0\), the endpoint differs only by one of the four cyclic relabelings

\[
\boxed{p(1)=T^k p(0),\quad T(p_0,p_1,p_2,p_3)=(p_1,p_2,p_3,p_0)}
\]

for some \(k\in\mathbb Z/4\). Reflection is excluded: the corners stay distinct and their **geometric** orientation does not change under a continuous path.

Use the two oriented diagonals and their dot product:

\[
d_0(t)=p_2(t)-p_0(t),\qquad
d_1(t)=p_3(t)-p_1(t),\qquad
f(t)=d_0(t)\cdot d_1(t).
\]

Under \(T\),

\[
(d_0,d_1)\longmapsto(d_1,-d_0),\qquad
f(Tp)=-f(p).
\]

Iterating yields

\[
\boxed{f(1)=(-1)^k f(0).}
\]

**Theorem 1 (odd-loop square forcing).** If \(k\) is odd, then the fixed Jordan curve \(\gamma\) contains a **positive-size inscribed square**.

**Proof.** If \(f(0)=0\), the first rectangle is already a square. Otherwise \(f(1)=-f(0)\), and the intermediate value theorem gives \(t_*\) with \(f(t_*)=0\). The configuration \(C(t_*)\) is by assumption a nondegenerate rectangle, so its two diagonals have the same midpoint and length. With \(d_0\cdot d_1=0\), the rectangle has perpendicular diagonals and is a positive-size square. QED.

**Corollary 2 (sharp non-square parity obstruction).** Any closed loop consisting *entirely of non-square rectangles* has even \(k\in\{0,2\}\). Conversely, an even-monodromy loop does not have to contain a square; an exact non-square \(k=2\) loop is constructed below.

For a generic transverse dot-product section on such a loop, an odd \(k\) forces an **odd number of simple square crossings modulo two**. The theorem does not require transversality or count multiplicities.

### Tensor-line interpretation

The dot product is a section of a sign line \(L_\perp\) associated with the character

\[
\chi:C_4\to\{\pm1\},\qquad \chi(T^k)=(-1)^k.
\]

The pulled-back line has

\[
\boxed{\langle w_1(L_\perp),[C]\rangle=k\bmod 2.}
\]

A nowhere-vanishing real section trivializes a real line bundle; thus odd \(k\) is an obstruction to avoiding perpendicular diagonals.

**Important correction to a common shortcut:** If a loop has odd \(k\), it **exchanges the diagonal labels**, so separate globally named diagonal lines \(L_0,L_1\) may **not exist individually** on the quotient loop. Their tensor/sign representation **still exists**, because exchanging the factors preserves the dot product. The formula \(w_1(L_0\otimes L_1)=w_1(L_0)+w_1(L_1)\) may be used only where the two diagonal line bundles are separately well-defined.

For even \(k=2\), both individual diagonal orientations can reverse after one circuit. Then each \(w_1(L_i)\) may be 1, while their tensor has \(w_1=1+1=0\). This is exhibited exactly in §3.

**The missing Toeplitz theorem:** Nothing here constructs an odd-\(k\) nondegenerate rectangle loop for **every** Jordan curve. The theorem is **conditional on such a loop actually existing**; the nontriviality of *arbitrary unordered chord space* cannot replace that hypothesis.

## 3. Sharp negative example: a non-square loop with two genuine Möbius diagonal twists

Use the **fixed** circle

\[
\Gamma=\{(x,y):x^2+y^2=25\}.
\]

Take the rectangle with semiwidth \(a=3\), semiheight \(b=4\), with cyclically ordered corners

\[
P_0=(3,4),\quad P_1=(-3,4),\quad
P_2=(-3,-4),\quad P_3=(3,-4).
\]

Rotate this rigid rectangle continuously through an angle \(\phi\in[0,\pi]\). Every corner remains on \(\Gamma\); every intermediate configuration is a nondegenerate rectangle of side lengths **6** and **8** and diagonals of length **10**.

The diagonal-dot product is invariant under rotation and equals

\[
d_0=(-6,-8),\quad d_1=(6,-8),
\qquad
\boxed{f=d_0\cdot d_1=28}
\]

for **every** real \(\phi\). Therefore **no square** occurs along this rectangle loop.

After a half-turn, the **unlabeled geometric rectangle returns to itself** but each corner advances two labels:

\[
P_j(\pi)=P_{j+2}(0),\qquad
\boxed{k=2\pmod4.}
\]

The two individually tracked diagonals both reverse:

\[
d_0(\pi)=-d_0(0),\quad
d_1(\pi)=-d_1(0).
\]

They therefore each carry a **nontrivial Möbius chord-orientation sign** over the loop, but

\[
\boxed{w_1(L_0)=w_1(L_1)=1,\qquad
w_1(L_\perp)=1+1=0.}
\]

This is a rigorous counterexample to the tempting inference that *individual chord-Möbius twists*, even both of them, force a square.

### Exact continuous algebraic parametrization

The code does not merely sample the rotation numerically. It uses the rational parametrization

\[
c(u)=\frac{1-u^2}{1+u^2},
\qquad
s(u)=\frac{2u}{1+u^2},\qquad 0\le u\le1,
\]

which rotates through \(0\le\phi\le\pi/2\). Two exact rational patches give the full half-turn:

\[
R_0(u)=\begin{pmatrix}c(u)&-s(u)\\s(u)&c(u)\end{pmatrix},
\qquad
R_1(u)=\begin{pmatrix}-s(u)&-c(u)\\c(u)&-s(u)\end{pmatrix}.
\]

At \(u=1\), \(R_0(1)=R_1(0)\), so the pieces glue **continuously**; their other endpoints are the identity and \(-I\). The checker proves the circle, midpoint, equal-diagonal, fixed-dot, endpoint-gluing, size-floor and \(T^2\) monodromy identities **exactly as rational function equalities for all real parameters**.

Thus the half-turn no-square loop is not an artifact of finite sampled frames.

## 4. Exact positive example: an odd-twisted loop with a unique irrational square

On the **same circle of radius 5**, choose a continuously varying rectangle whose semiwidth \(a\) rises from 3 to 4 while semiheight \(b\) falls from 4 to 3, preserving

\[
a^2+b^2=25.
\]

Use the exact rational stereographic parameter

\[
p(u)=\frac{3-u}{6},\quad
a(u)=5\frac{1-p(u)^2}{1+p(u)^2},\quad
b(u)=10\frac{p(u)}{1+p(u)^2},
\qquad 0\le u\le1.
\]

The code proves

\[
\boxed{a(u)^2+b(u)^2=25},\quad
\boxed{3\le a(u),b(u)\le4}
\]

for **all real \(u\in[0,1]\)** using exact rational polynomial and Bernstein sign certificates. Consequently every rectangle has strictly positive side lengths, with

\[
\boxed{\text{side}^2\ge36}
\]

and both diagonals of length **10** at every time.

Now rotate these rectangles through the **same continuous quarter-turn** \(R_0(u)\). At \(u=0\), the corners are those of the original \(6\times8\) rectangle. At \(u=1\), the unrotated rectangle is \(8\times6\), so the extra 90-degree rotation makes its **geometric corner set exactly the original one**, but with each corner advanced **one** label:

\[
\boxed{P_j(1)=P_{j+1}(0),\quad k=1.}
\]

Thus this closed geometric rectangle loop has the **nontrivial tensor monodromy** that the secant-spin proposal actually needs.

The signed dot product is

\[
f(u)=4\bigl(b(u)^2-a(u)^2\bigr).
\]

At the endpoints,

\[
\boxed{f(0)=+28,\qquad f(1)=-28.}
\]

It must vanish. The condition \(f(u)=0\) simplifies to

\[
p(u)^2+2p(u)-1=0,
\]

whose **unique root in** \(p\in[1/3,1/2]\) is \(p=\sqrt2-1\). Therefore

\[
\boxed{u_*=9-6\sqrt2}
\quad\text{and}\quad
\boxed{u_*^2-18u_*+9=0.}
\]

The rational-function numerator is proved divisible by the exact algebraic polynomial \(u^2-18u+9\); the Sturm routine independently verifies it has **exactly one irrational root in (0,1)**. Thus there is exactly one distinct square event in the loop. At this instant the two diagonals have length 10 and are perpendicular, giving

\[
\boxed{\text{square side}^2=10^2/2=50.}
\]

This is a protected **positive-size** square crossing, not a zero-size singularity.

**Critical distinction:** The example is *not evidence that arbitrary Jordan curves admit such an odd loop*. It constructs one on the already square-rich **circle** to validate the forcing criterion and demonstrate its monodromy behavior. The even loop in §3 lives on the *same* curve, showing that the existence of nontrivial chord twists or non-square rectangle loops is not equivalent to this odd condition.

## 5. A reusable exact certificate, rather than only two examples

The function

`certify_circle_rectangle_loop(patches, radius_squared, strict_side_floor)`

accepts any finite collection of **rational-function patches**, each describing four continuously varying corners of a candidate inscribed rectangle on a **fixed circle centered at the origin**.

It is an exact sufficient checker with the following fail-closed obligations:

1. **Every real parameter:** all four corners satisfy \(x^2+y^2=R^2\) **identically as rational functions**. Every patch satisfies the midpoint and equal-diagonal rectangle equations **identically**, without fitting sampled points.
2. **No poles:** every denominator is strictly nonzero on its entire closed parameter interval, certified by rational Bernstein inequalities.
3. **Noncollapse:** every consecutive-corner squared distance is strictly greater than the user-specified **positive rational side floor**, at every real parameter in every patch.
4. **Correct orientation:** the initial corners of each patch are geometrically CCW. Distinctness on every real parameter prevents a change of orientation.
5. **Continuous gluing:** the endpoint of each patch exactly equals the start of the next patch, in the same four-corner lift.
6. **Closed geometric loop:** the final corners exactly match a single \(T^k\) relabeling of the initial corners. Wrong identifications or nonclosed paths are rejected.
7. **Complete algebraic root coverage:** the signed perpendicularity is a rational function with certified nonzero denominator; exact Sturm isolation finds **all** distinct positive-size square event parameters. A patch where perpendicularity vanishes **identically** is marked **inconclusive as a finite event census** rather than assigned zero or one.
8. **Monodromy cross-check:** its endpoint dot-product values must transform by \((-1)^k\), and an odd-\(k\) loop must contain at least one certified polynomial zero.

The generic checker validates the even loop as \((k,\text{square-event count})=(2,0)\) and the odd one as \((1,1)\). Negative tests force rejection of an incorrect circle radius, a false side floor, a collapsed endpoint, a continuously square-valued interval, a loop that does not close, a pair of patches that do not glue, and a clockwise corner sequence in the CCW certificate.

Because the checker is specifically for the **known circle**, it **does not purport to construct** rectangle components on other Jordan curves. Extending it to exact rational-parametric polygons and computed rectangle components is a sensible next engineering step.

## 6. What the April 2026 Greene–Lobb paper adds

The most relevant recent primary source is Joshua Evan Greene and Andrew Lobb, [*Jordan curves inscribe a positive measure of rectangles*](https://arxiv.org/abs/2604.17116), submitted **18 April 2026**.

For **every** Jordan curve of diameter \(2R\) enclosing area \(A>0\), they prove there is a set \(S\subset(0,\pi)\) of **diagonal angles** supporting inscribed rectangles with

\[
\boxed{|S|\ge \frac{A}{R^2}.}
\]

Moreover, their **Corollary 2.2** yields a quantitative noncollapse statement: for any \(\varepsilon>0\), there is a set of such angles of measure at least

\[
\boxed{\frac{A-\pi\varepsilon}{R^2}}
\]

for which some inscribed rectangle has **squared circumradius at least \(\varepsilon\)**. This can be strictly positive whenever \(0<\varepsilon<A/\pi\).

Their paper also proves a **Genericity Theorem (3.1)**: for an open dense set of smooth Jordan embeddings, an appropriately **compactified** space of rectangles **and binormals** is a smooth closed one-dimensional manifold. The **binormal boundary states** at angles \(0,\pi\) matter and are not automatically loops of genuine rectangles at angles in \((0,\pi)\).

These facts are significant for our program. They show that (i) arbitrary continuous Jordan curves have far more inscribed rectangles than Vaughan's one-rectangle theorem detects, and (ii) a **size-selected** subset can be protected quantitatively from total collapse at a **positive-measure set of angles**, even for \(C^0\) curves. Their Floer spectral invariant has an action derivative connected to **rectangle radius squared**.

**But:** Positive-measure angle support, even with a radius floor, does **not** guarantee the specific square angle \(\theta=\pi/2\) belongs to that support. In elementary set theory, a set \(S\subset(0,\pi)\) can have measure **arbitrarily close to \(\pi\)** while omitting \(\pi/2\), for example \(S=(0,\pi/2-\delta)\cup(\pi/2+\delta,\pi)\) for small \(\delta>0\). This is not claimed to be a realizable rectangle-angle set of a Jordan curve; it only shows why **the measure lower bound alone** cannot logically force the square angle.

And size-selected rectangles at many **different angles** need not lie on **one closed connected component with odd monodromy**. Proving that connection is still missing.

For **rectifiable** Jordan curves, Greene–Lobb's Theorem B is stronger: the curve inscribes \(\theta\)-rectangles for **all** \(0<\theta<A/R^2\), which in particular includes \(\pi/2\) when \(A/R^2>\pi/2\). That restricted theorem is already covered by their paper. We must not claim it anew.

## 7. What would make this a genuine Toeplitz proof strategy?

The research now has a sharp **relative topology versus geometric size** separation.

**Topological obligation T.** From the Jordan embedding, produce an actual **connected** loop or relative cobordism of inscribed rectangles with **odd cyclic monodromy**. It must be a geometric rectangle family, not merely a loop in the large free unordered-chord Möbius strip. The §3 even-loop counterexample proves that individual chord-line \(w_1\) is insufficient.

**Geometric obligation G.** Ensure that the odd-monodromy component cannot pass through total-collision/zero-side strata. For a *fixed, compact, entirely nondegenerate loop*, its minimum side length is automatically positive by continuity. The difficulty is demonstrating such a loop or an equivalent boundary-relative class **uniformly under arbitrary \(C^0\) approximations**.

**Multiplicity/sector obligation S.** Capture nongraceful rectangles and squares, singular and positive-dimensional rectangle components, and pairing changes under the full geometric \(D_4\) action without silently treating a fourfold labeled orbit as a mod-2 obstruction.

**Literature-link obligation L.** Test whether Greene–Lobb's **quantitatively radius-protected spectral rectangles** can populate a single connected component with an odd cyclic monodromy, or whether the available compactified 1-manifolds can close only through binormal or collision states. No result found in the 2026 paper directly proves T.

Proving these would turn the elementary monodromy observation into a substantive all-Jordan argument. Until then, the odd-loop theorem is a valid **sufficient** criterion, not an existence theorem for arbitrary Jordan curves.

## 8. Reproduction and status

New source: [`code/rectangle_monodromy_v268.py`](../code/rectangle_monodromy_v268.py).  
New tests: [`tests/test_rectangle_monodromy_v268.py`](../tests/test_rectangle_monodromy_v268.py).  
Dedicated workflow: [`v268-rectangle-monodromy.yml`](../.github/workflows/v268-rectangle-monodromy.yml).

```sh
python code/rectangle_monodromy_v268.py
python -m unittest -v tests.test_rectangle_monodromy_v268
python -m unittest -v tests.test_moving_polygon_events_v224 tests.test_configuration_symmetry_v218
```

The code uses **exact rational arithmetic and rational polynomial identity checks**, with exact Bernstein and Sturm certificates for all closed real-parameter intervals. No floating-point comparison or numerical sampling is used as proof of an event or noncollapse. The topological odd-loop criterion is justified analytically by continuity and the intermediate value theorem, and the full-Jordan obligations are **not** machine-proved by finite tests.

### Related project results

- [v218](v218-global-configuration-symmetry-and-collision-audit.md): cyclic source-sign cancellation and the six curve-order sectors.
- [v219](v219-equivariant-quotient-collision-link.md): vanishing *absolute* rank-four obstruction, total-collapse strata, and microscopic-notch selected-square shrinkout.
- [v262](v262-cyclic-quotient-relative-parity.md): nonzero **relative** mod-2 quotient degree in a selected convex-pentagon atlas through square-pair annihilation.
- [v265](v265-generic-irrational-vertex-gluing.md): exact all-time algebraic vertex-pair gluing and the 3→2→1 pentagon root census.
- [v266](v266-two-corner-stability-and-convex-cutoff.md): independent hidden-corner/convex noncollapse arguments.
- [v267](v267-secant-spin-rectangle-obstruction.md): universal antipodal secant winding, the Möbius chord-line twist and the diagonal-tensor conditional theorem.

**Conclusion:** The odd-monodromy square-forcing test is mathematically exact, but the same circle also contains a non-square rectangle loop with **both individual chord orientation lines twisted** and **trivial perpendicularity tensor twist**. To prove the unrestricted Toeplitz conjecture using this framework, we must establish **a protected odd-monodromy rectangle family on every Jordan curve**, or a genuinely stronger relative/Floer boundary invariant. **The unrestricted Toeplitz conjecture remains open.**
