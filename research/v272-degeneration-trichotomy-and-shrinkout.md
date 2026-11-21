# v272: Exhaustive rectangle degeneration trichotomy and an explicit Jordan shrinkout family

**Program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Exact elementary geometric classification; rational Jordan polygon construction; not a universal noncollapse theorem.

## 1. Rectangles admit only three limit strata

Consider a convergent sequence of Euclidean rectangles represented by *cyclically geometrically ordered* vertices, and let the limiting adjacent side vectors be `u` and `v`. The rectangle equations imply `u·v=0` and, after translation, its vertices take the form `p, p+u, p+u+v, p+v`, also in the limit. Both squared side lengths are nonnegative. Exactly three cases are possible:

1. **Nondegenerate rectangle:** `|u|>0` and `|v|>0`. All corners remain distinct, the diagonal radius is positive.
2. **Two-pair flattening:** exactly one of `|u|, |v|` is zero. The four corners converge to two distinct points, each with multiplicity two, and the diagonal length remains positive. This is the ellipse's boundary type in v268.
3. **Total four-corner collision:** `|u|=|v|=0`. All four points converge to the same position and both diagonals vanish.

In particular, **there is no genuine one-pair-only or three-corner limit of actual Euclidean rectangles**. If one adjacent pair collides, the opposite adjacent pair must also collide. If two opposite corners collide, equal diagonal lengths force all four to collide. The classification follows directly from the rectangle algebra, without genericity or curve regularity.

**Square corollary:** since a square has `|u|=|v|`, a convergent sequence of squares can only have a positive square limit or a total four-corner collision. It can never converge to a nonzero two-pair flattening. Thus a relative sign crossing at a *flat end with nonzero diagonal dot* is safe as a crossing at that end, but it still gives no uniform guarantee against escape into a **different**, total-collision stratum under varying curves.

## 2. Exact selected-square shrinkout on genuine Jordan curves

Define `epsilon=1/n`, `n>=4`, `h=1/3`. Take the unit square boundary and replace the short top segment from `(h,1)` to `(h+epsilon,1)` with an *outward square bump*. The eight CCW polygon vertices are

```
(0,0), (1,0), (1,1), (h+eps,1),
(h+eps,1+eps), (h,1+eps), (h,1), (0,1).
```

These vertices describe a **simple Jordan polygon**, as can be proved by the two monotone top subarcs and confirmed by an exact nonadjacent-segment intersection test. Its four bump corners form a positive inscribed square of side `epsilon`; `F=0` exactly. Its squared circumradius is `epsilon²/2`.

The polygon boundary converges in Hausdorff distance to the unit square boundary with distance at most `epsilon`. Nevertheless the *selected bump square* has side squared `1/n² -> 0` and all four corners converge to `(h,1)`. Rescaling around that point by `1/epsilon` produces the fixed unit-square configuration

`(1,1), (0,1), (0,0), (1,0)`.

The construction does **not** show that the limiting unit square has no positive inscribed square (it obviously has one), nor does it show that *every* square on the approximants shrinks. It is a counterexample to a purported positive size floor on an **arbitrarily selected** branch of inscribed squares, even under convergence to a fixed nondegenerate Jordan boundary. The explicit blow-up demonstrates why a local microlimit does not by itself preserve the large-scale topological class.

## 3. Research implications

- The prospective compactified rectangle boundary should separate two-pair flats with nonzero diagonal from total-collision configurations. Merging them as one boundary type destroys the distinguishing relative sign mechanism of v269.
- Any continuation theorem based on a selected square needs a reason that its branch cannot approach the total-collision stratum. Topology or winding of an unrelated secant line is insufficient (cf. v267).
- A blow-up is only useful if one can control the *rescaled curves* and relate the local configuration to a protected global class. Arbitrary Jordan curves need not have a unique or locally compact tangent model.

**Unsolved:** a compactified geometric rectangle cobordism with quantitative noescape through total collision on all Jordan boundaries.

Code and tests: [`code/relative_rectangle_v269_v273.py`](../code/relative_rectangle_v269_v273.py), [`tests/test_relative_rectangle_v269_v273.py`](../tests/test_relative_rectangle_v269_v273.py).