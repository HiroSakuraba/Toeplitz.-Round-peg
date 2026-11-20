# v257: A quantitative near-square on any approximately centrally symmetric Jordan curve

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Result:** An unconditional right-angle boundary-pair theorem for every Jordan curve, and a **quantitative approximate-square** theorem for curves that are near centrally symmetric. The theorem imposes **no smoothness, convexity, or local regularity hypothesis**. A fully rational verifier constructs explicit four-point witnesses for polygonal examples.

**Scientific status:** A theorem about four actual boundary points **near** the corners of a large square. It does **not** assert those four boundary points form an exact square. It therefore does not solve Toeplitz or the general weak Big Peg problem.

## 1. Two exactly orthogonal boundary points exist on every Jordan curve

Let \(\gamma\) be any simple closed continuous plane curve, \(c\) any point in its bounded interior, and \(R\) rotation by \(90^\circ\) about \(c\).

The **quarter-turn intersection lemma** from v256 gives

\[
\gamma\cap R\gamma\ne\varnothing.
\]

Take \(v\) in this intersection, and put \(w=R^{-1}v\). Both \(v,w\) belong to the original boundary. In position vectors centered at \(c\), \(w=-Rv\), so

\[
|v-c|=|w-c|,\qquad (v-c)\cdot(w-c)=0.
\]

If the bounded Jordan domain contains the closed disk \(D(c,r)\), then \(|v-c|\ge r\). Therefore **every Jordan curve** with that centered disk has two boundary points at the same distance at least \(r\) from \(c\), separated by a right angle.

This is a universal **two-point** statement, not an inscribed square. The remaining two square vertices need not lie on the curve.

## 2. The symmetry defect and the approximate square theorem

Let \(-\gamma\) denote the half-turn of \(\gamma\) about \(c\), and define its *Hausdorff central-symmetry defect*

\[
\delta=d_H(\gamma,-\gamma).
\]

Because \(\gamma\) is compact, all nearest points required in this argument exist, and \(\delta\) is finite.

Given \(v,w=-Rv\in\gamma\) from the quarter-turn intersection, the ideal centered square would have ordered corners

\[
S_0=v,\quad S_1=Rv=-w,\quad S_2=-v,\quad S_3=w.
\]

Two corners \(S_0,S_3\) already belong **exactly** to \(\gamma\).

By the definition of \(d_H(\gamma,-\gamma)\), there exist points \(A_1,A_2\in\gamma\) satisfying

\[
|A_1-(-w)|\le\delta,\qquad |A_2-(-v)|\le\delta.
\]

Set \(A_0=v\), \(A_3=w\). We obtain four **actual boundary points** with

\[
\boxed{|A_k-S_k|\le\delta\quad (k=0,1,2,3),}
\]

and in fact \(A_0=S_0\), \(A_3=S_3\).

The ideal side of this square is

\[
\boxed{s_{\rm ideal}=\sqrt2\,|v-c|\ge\sqrt2\,r}.
\]

For any pair of vertices, the reverse triangle inequality gives

\[
\big||A_i-A_j|-|S_i-S_j|\big|\le 2\delta.
\]

Thus the four boundary points form a **quantitatively near-square**: their side lengths differ from an ideal square's side by no more than \(2\delta\), and their positions are each within \(\delta\) of the ideal corners.

If \(\delta<r/2\), all four chosen points are **distinct**. Indeed every two ideal square corners are separated by at least \(\sqrt2\,r\), while each is moved by at most \(\delta\). Their final separation is at least \(\sqrt2\,r-2\delta>0\).

In particular, \(\delta=0\) recovers v256's exact square because the curve is centrally symmetric. A small positive \(\delta\), however, provides only *nearby* points, not four exact square corners.

### Why this isn't Toeplitz

A sequence of four boundary points forming arbitrarily accurate approximate squares **on a fixed Jordan curve**, with a uniform nonzero side-size floor, would yield an exact square by compactness.

But the theorem here gives an error \(\delta\) fixed by the given curve. For a generic asymmetric curve \(\delta>0\), and we have not proved we can make \(\delta\to0\) while keeping the same boundary. Reducing the symmetry defect generally changes the curve. Hence no compactness shortcut proves an exact square on the original nonsymmetric curve.

Similarly, without a control on \(\delta/r\), this result cannot provide an approximate square of arbitrarily good relative accuracy for every Jordan curve.

## 3. A fully rational polygonal certificate

The [combined v256/v257 verifier](../code/central_symmetry_square_v256.py) operates on rational Jordan polygons with **an even number \(n=2m\) of vertices**. This even-vertex requirement is a limitation only of the convenient **explicit indexed correspondence algorithm**, not of the general Jordan theorem.

Associate each edge \(i\) to edge \(i+m\). Define the **exact rational squared upper defect**

\[
\boxed{\Delta^2=\max_i |P_i+P_{i+m}|^2}
\]

with the origin at the intended symmetry point.

For any \(t\in[0,1]\), the corresponding boundary points

\[
p_i(t)=P_i+t(P_{i+1}-P_i),
\qquad
p_{i+m}(t)=P_{i+m}+t(P_{i+m+1}-P_{i+m})
\]

satisfy the exact identity

\[
p_i(t)+p_{i+m}(t)
=(1-t)(P_i+P_{i+m})+
t(P_{i+1}+P_{i+m+1}).
\]

The Euclidean norm of this convex combination is at most the maximum of its two endpoint norms, hence \(|p_i(t)+p_{i+m}(t)|\le\Delta\) for **every real interpolation parameter**.

Because the edge pairing is bijective on the entire polygon, this proves the **global** bound

\[
\boxed{d_H(\partial P,-\partial P)\le\Delta}.
\]

There is no claim that the index pairing always gives the *smallest possible* Hausdorff bound.

A rational intersection of \(\partial P\) and \(R(\partial P)\) supplies exact boundary points \(v,w\). Their opposite partners are constructed by interpolating the indexed opposite edges at the **same rational interpolation parameters**. Therefore the four boundary points are exactly rational, have certified corner error squared \(\le\Delta^2\), and are distinct whenever \(\Delta^2<r^2/4\).

The code checks each displayed corner lies on the actual polygonal boundary, independently checks the exact square geometry of the ideal corners, verifies the disk containment, and refuses an approximate-square certificate with insufficient margin.

## 4. An asymmetric, nonconvex example

Start with v256's centrally symmetric **12-sided polygon with inward notches**, which is neither convex nor star-shaped about the origin and contains a centered disk of radius 2.

Displace one vertex, initially \((-1,-2)\), by the rational vector

\[
\left(\frac1{20},\frac1{40}\right).
\]

The resulting polygon is a genuine asymmetric Jordan curve. The unit-by-unit antipodal vertex mismatch is zero for all but two partners, so

\[
\boxed{\Delta^2=(1/20)^2+(1/40)^2=\frac1{320}}.
\]

The closed radius-2 disk is still contained in the new polygon, independently certified by exact edge-to-center squared-distance tests.

The verifier constructs **four distinct rational boundary points**, two of which lie **exactly** on a centered square of squared side at least 8. The other two have exact squared error at most \(1/320\) from the corresponding ideal square vertices.

This is a direct, rational example of the general theorem in a **nonconvex, asymmetric** setting.

## 5. Research implications and next steps

- The v256 theorem gives a **sharp positive bound for exact squares** in a nonconvex, arbitrarily rough, *centrally symmetric* family.
- v257 proves that **two of the four large-square corners can always be found on any Jordan curve**, and gives a controlled **approximate** four-corner configuration under near symmetry.
- Exact central symmetry, not convexity, is sufficient for the full four-point result. A small violation of symmetry does **not** automatically preserve an exact centered square. The other two corners may fail to lie at precisely the right locations.
- A meaningful v258 would seek an additional *transversality or degree* assumption under which these near-square witnesses perturb to **exact boundary squares** without losing a size floor. Such a theorem must explain how the two missing incidence constraints are solved, not merely invoke compactness for fixed nonzero \(\delta\).

For an unrestricted Toeplitz proof, one ultimately needs a size-controlled existence result that does not rely on symmetry, convexity, or prescribed local chart conditioning.

### Reproduce

~~~sh
python code/central_symmetry_square_v256.py
python -m unittest -v tests.test_central_symmetry_square_v256
~~~

The exact tests check both the perfectly symmetric and perturbed examples, the opposite-edge Hausdorff upper bound, general rational quarter-turn edge intersections, segment overlap degeneracies, disk containment, four distinct corner certificates, and fail-closed invalid inputs.

**References:** The v256 [quarter-turn intersection proof](v256-centered-square-jordan-symmetry.md) is self-contained modulo the Jordan curve theorem. For broader square-peg context, see B. Matschke, *A Survey on the Square Peg Problem*, *Notices AMS* (2014), and the symmetry and regularity subclasses discussed in *On the Square Peg Problem and Its Relatives*, arXiv:1001.0186. No historical priority claim is made.
