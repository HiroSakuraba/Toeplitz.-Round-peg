# v271: Complete discrete sector/pairing ledger, not a global Jordan atlas

**Program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Exact universal symmetry ledger **plus a fully solved two-branch rectangle locus in one four-edge cell of an actual nonconvex Jordan polygon**. Not an all-cells global Jordan atlas.

## 1. Distinguishing three different quotients

Suppose four distinct points of a Jordan embedding are listed in source-parameter order `q0<q1<q2<q3`. Their *geometric* cyclic order need not match this listing when the curve is nonconvex. Treat their geometric order as a separate choice; never silently assume that all rectangles are gracefully ordered.

There are `4! = 24` orderings of four distinct geometric vertices. Rotation of a word by `T` gives a free four-element `C4` orbit, so **six** cyclic order sectors. Enlarging to the eight symmetries of a square (`D4`, including reversal) yields **three** eight-element combinatorial orbits. This quotient is useful bookkeeping, not a claim that curve-parameter sectors can be identified by a sign-preserving geometric homotopy.

For the four source positions there are exactly three candidate opposite-vertex pairings:

- `(0,1)|(2,3)`;
- `(0,2)|(1,3)`;
- `(0,3)|(1,2)`.

A supplied genuine nondegenerate rectangle has **exactly one pairing** satisfying coincident diagonal midpoints and equal diagonal lengths. Indeed its unique center is the common midpoint of its two actual diagonals; adjacent vertex pairs have distinct midpoints. The exact enumeration checks all `24` orderings of the axis-aligned `6x8` rectangle, with eight successful choices for each of the three pairings. The same uniqueness holds for an actual square.

## 2. Exact cyclic transformation of the square equation

Use v262's CCW square residual

`F(p0,p1,p2,p3) = (p2-p1-J(p1-p0), p3-p0-J(p1-p0)) in R4`.

Cyclic relabeling satisfies `F(Tp)=R F(p)` with

```
R = [[-1, 1, 1, 0],
     [-1,-1, 0, 1],
     [ 0, 1, 0, 0],
     [-1, 0, 0, 0]]
```

The code verifies this **universal linear identity on all eight Cartesian basis vectors**, not merely on a set of sampled rectangles. It also computes exactly `det R=+1`, `R^4=I`, and the source four-cycle's determinant `-1`.

At a transverse labeled zero, the cyclic relabeling reverses the local signed index, exactly as v218 and v262 observed. The sign change is not a new nonzero global Euler class: in the full fourfold labeled orbit, naïve signed indices cancel. Passing to a cyclic quotient requires an associated bundle and genuinely consistent local trivializations, and it still does not guarantee nonzero *absolute* rank-four obstruction (v219).

## 3. An actual nonconvex four-arc rectangle algebraic locus

Unlike a fixed-rectangle permutation audit, this example constructs and **exhausts a genuine one-dimensional rectangle locus within a specified four-edge assignment cell**.

Take the simple rational polygon whose eight CCW vertices are

```
(0,0), (1,0), (1,1), (17/20,1),
(17/20,11/10), (3/4,11/10), (3/4,1), (0,1).
```

The outward bump makes the polygon nonconvex; no self-intersections occur. Select four distinct boundary arcs: right edge, upper-left horizontal segment, left edge, bottom edge. On the compact parameter chart `2/5<=t,u<=3/5`, their corners in geometric CCW order are

`p0=(1,t), p1=(u,1), p2=(0,1-t), p3=(1-u,0)`.

Every point lies strictly in its designated straight boundary arc, including at `t,u=1/2`. The equal-midpoint identity is automatic:

`p0+p2 = p1+p3 = (1,1)`.

The remaining equal-diagonal equation is the **complete polynomial identity**

`|p2-p0|² - |p3-p1|² = 4(t-u)(t+u-1)`.

Therefore *every rectangle in this two-variable centered four-arc cell*, not merely sample paths, lies on exactly one of the two algebraic branches

`u=t` or `u=1-t`.

The diagonal dot is exactly

`(p2-p0)·(p3-p1)=2(t+u-1)`.

On `u=t`, it equals `4t-2` and crosses zero at `t=1/2`; on `u=1-t` it **vanishes identically**, so that entire branch consists of positive inscribed squares. The branches meet at `(t,u)=(1/2,1/2)` and cannot be treated as unrelated simple-zero orbits. For both, each side has Cartesian components of absolute value at least `2/5`, giving exact **squared side length >=8/25** throughout the cell. The geometry is uniformly noncollapsed.

This is a direct stress test for **positive-dimensional square loci and singular rectangle-branch intersections on a nonconvex Jordan polygon**. The analytic factorization proves all-time completeness *within this selected four-edge cell*; the exact implementation checks the polygon, each branch, and its side floor at adversarial endpoints. It does **not** enumerate other four-edge cells or every square on the entire polygon.

## 4. What this checks and does not check

**Verified exactly:** all permutation orbit sizes and diagonal pairings; the universal target representation; and all rectangle solutions in one nonconvex four-arc cell, including a positive-dimensional square branch and a singular crossing. This is a stronger test than only a permutation census.

**Not verified:** traversal of *every* polygon edge-assignment cell, all six global curve-order branch connections, transitions across vertex glue planes, whole-curve all-time root census, boundary binormal weights, or a nonzero all-Jordan relative class.

## 5. Next geometric atlas milestone

1. Start from one rational nonconvex Jordan polygon; enumerate four-edge assignment cells and all three potential diagonal pairings.
2. Solve the rectangle midpoint and equal-diagonal equations as algebraic sets in four edge parameters. Remove all-four-corner collision and classify flat endpoints using v272.
3. Glue crossing branches only after certifying their shared geometric corners and local orientation transition.
4. Add quotient-relative boundary weights, then independently compare against a complete square census before making any global index claim.

The current implementation covers the complete permutation ledger **and one actual nonconvex cell**, but the whole-polygon cell complex and global gluing are still missing. This is partial completion of the originally ambitious v271 global-atlas milestone, not a universal substitute.

Code and tests: [`code/relative_rectangle_v269_v273.py`](../code/relative_rectangle_v269_v273.py), [`tests/test_relative_rectangle_v269_v273.py`](../tests/test_relative_rectangle_v269_v273.py).