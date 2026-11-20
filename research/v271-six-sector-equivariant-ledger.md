# v271: Complete discrete sector/pairing ledger, not a global Jordan atlas

**Program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Exact universal linear symmetry check and complete `S4` permutation census of *one fixed actual rectangle*. **Not** an enumeration of rectangle branches of a general Jordan curve.

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

## 3. What this checks and does not check

**Verified exactly:** permutation orbit sizes, three pairings, the relevant universal target representation, and source orientation sign. This supplies a dependable sector-completeness gate for later polygonal rectangle-atlas algorithms.

**Not verified:** whole-curve existence of six 1-dimensional branches, traversal of every polygon edge-assignment cell, transitions across vertex glue planes, all-time root census, boundary binormal weights, or a nonzero all-Jordan relative class. Those tasks cannot be inferred from an `S4` census of a fixed rectangle.

## 4. Next geometric atlas milestone

1. Start from one rational nonconvex Jordan polygon; enumerate four-edge assignment cells and all three potential diagonal pairings.
2. Solve the rectangle midpoint and equal-diagonal equations as algebraic sets in four edge parameters. Remove all-four-corner collision and classify flat endpoints using v272.
3. Glue crossing branches only after certifying their shared geometric corners and local orientation transition.
4. Add quotient-relative boundary weights, then independently compare against a complete square census before making any global index claim.

The current implementation intentionally stops at the rigorous finite algebraic symmetry layer. This is *partial completion* of the originally ambitious v271 global-atlas milestone, not a substitute for its geometric branch computation.

Code and tests: [`code/relative_rectangle_v269_v273.py`](../code/relative_rectangle_v269_v273.py), [`tests/test_relative_rectangle_v269_v273.py`](../tests/test_relative_rectangle_v269_v273.py).