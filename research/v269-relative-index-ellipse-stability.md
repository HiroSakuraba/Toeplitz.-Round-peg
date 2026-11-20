# v269: A relative rectangle sign class and local degree at an ellipse square

**Program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Proved restricted relative-interval and local-C0 persistence statements; exact rational determinant certificate. **No global compactified-rectangle index for arbitrary Jordan curves.**

## 1. Correct topological object

Let `M` be a compact connected 1-manifold with boundary, representing a *given* (not hypothetical) family of nondegenerate rectangles except at controlled endpoints. Let `L -> M` be the real sign line for the diagonal dot product `f=d_02·d_13`. Over a trivialized interval `M=[0,1]`, if `f(0) f(1) < 0` and `f` does not vanish on the boundary, then an interior square occurs. For a generic transverse section, the number of zeros mod 2 equals the mismatch of endpoint signs.

If `M=S^1` and `w1(L)[M]=1`, every continuous section has a zero; generic transverse zeros are odd mod 2. These two mechanisms have a common conceptual description as a relative first Stiefel-Whitney / zero-divisor class, *provided the line, compactification and boundary trivialization actually exist*. In particular, it is incorrect to evaluate a global `w1` on a made-up 1-cycle, or to add endpoint signs from unrelated local trivializations. For a cobordism, zeros may reach a boundary where `f` vanishes, invalidating the relative count.

The interval theorem is elementary: whenever all zeros are simple, crossing a zero reverses sign, so their count is odd iff the two end signs disagree. The nongeneric conclusion (at least one zero when signs disagree) follows from the intermediate value theorem. The zero-free-boundary hypothesis is essential.

The v268 noncircular ellipse has a genuine rectangle locus `M=(0,a)`, whose two ends naturally compactify to *two-pair-flat* configurations of nonzero diagonals. Its natural geometric cyclic labeling has `f=-4a²` and `f=+4b²` at the respective ends. This gives a relative parity `1` even though `H1(M;F2)=0` and no essential odd rectangle loop exists.

## 2. A second, independent persistence certificate: nonzero square Jacobian

Take the noncircular ellipse

`gamma(phi) = (a cos phi, b sin phi)`, `a,b>0`, `a != b`.

At its unique geometric inscribed square, define

`s = ab/sqrt(a²+b²)`.

In geometric CCW order, the square corners are

`p0=(s,s), p1=(-s,s), p2=(-s,-s), p3=(s,-s)`.

Define the exact square residual

`F1 = p2-p1-J(p1-p0)`, `F2 = p3-p0-J(p1-p0)`

with `J(x,y)=(-y,x)`. We view `F=(F1,F2)` as a map from four independent nearby ellipse parameters `(phi0,phi1,phi2,phi3)` into `R^4`. The normalized tangent columns at the four corners are `gamma'(phi_j)/s`:

- `(-a/b, b/a)`;
- `(-a/b,-b/a)`;
- `( a/b,-b/a)`;
- `( a/b, b/a)`.

The derivative columns of `F` from these tangent vectors `t_j` are respectively

- column 0: `(J t0, (J-I)t0)`;
- column 1: `((-I-J)t1, -J t1)`;
- column 2: `(t2, 0)`;
- column 3: `(0, t3)`.

Taking the determinant and restoring the four omitted factors of `s` gives

`det(D F) = -4 a² b² (a²-b²)/(a²+b²)`.

This is **nonzero precisely when the axes differ**. For `a=3,b=4`, it is exactly `4032/25`; the square has side squared `576/25`. The code computes the 4-by-4 rational determinant independently and compares it against the simplified closed formula for arbitrary positive rational unequal `a,b`. The circle (`a=b`) is excluded because rotated squares form a nonisolated family and this four-variable Jacobian degenerates.

### Local-C0 stability corollary

Choose four disjoint closed parameter arcs around the four square corners, small enough that `F` has just the one regular zero inside their product `U` and no zero on its boundary. The local Brouwer degree of `F` is `sign(det D F)=±1`.

Let `gamma_t` be a Jordan embedding uniformly close to the given ellipse *in this fixed parameterization*. Its corresponding square residual `F_t` is uniformly close to `F` on `U` because it depends linearly on the four corner positions. For sufficiently small perturbation, `F_t` remains nonzero on `∂U` and has the same degree. Hence `F_t` has a zero inside `U`. The separated arcs and the original positive intercorner distances imply a positive square side at that zero. **No differentiability of the perturbed embedding is required.**

This is a qualitative open-neighborhood result, not an explicit numerical perturbation budget and not a new proof of the known smooth Jordan theorem. It is also not an arbitrary-Jordan argument: an arbitrary curve need not be close to a noncircular ellipse in the required parameterized uniform sense.

## 3. Research outcome

**Achieved:** two correct realizations of square persistence on the same ellipse: a relative opposite-endpoint sign parity and an isolated local Brouwer degree. These do not rely on an odd closed rectangle loop. Exact tests distinguish zero-free relative boundaries, circle rank loss, and nonzero determinant.

**Unsolved obligation:** Construct, for arbitrary Jordan approximation families, a well-defined *global* relative class with canonical boundary trivialization/weights at binormal strata, and prove that contributions cannot disappear into all-four-corner collapse. Gluing selected local `±1` indices does not automatically give such a class; see v262 and v218–v219.

## Reproduction

`python -m unittest discover -s tests -p 'test_relative_rectangle_v269_v273.py' -v`

Code: [`code/relative_rectangle_v269_v273.py`](../code/relative_rectangle_v269_v273.py).