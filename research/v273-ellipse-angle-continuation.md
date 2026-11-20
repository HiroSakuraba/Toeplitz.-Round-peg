# v273: Exact full rectangle-angle continuation on ellipses and a sharp protected angle band

**Program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Complete analytic **noncircular ellipse** angle locus and exact size floors; sharp separation from the unrestricted Jordan problem.

## 1. A solvable angle-continuation model

By the v268 global classification, every actual inscribed rectangle in the noncircular ellipse `x²/a²+y²/b²=1` has axis-parallel vertices `(±x,±y)` with `x,y>0`. Choose the two directed diagonals

`d0=(-2x,-2y), d1=(2x,-2y)`.

Their lengths agree and the cosine of the directed diagonal angle is

`c = (d0·d1)/|d0|² = (y²-x²)/(x²+y²)`.

Every `c` strictly between `-1` and `+1` occurs at **exactly one** geometric rectangle. Solving the ellipse and angle equations gives

`D(c) = b²(1-c)+a²(1+c)`;

`x² = a² b² (1-c)/D(c)`;

`y² = a² b² (1+c)/D(c)`;

`rho²(c) = x²+y² = 2a² b²/D(c)`.

The expressions are strictly positive for all `-1<c<1`, and the ellipse incidence and cosine identity follow algebraically. The diagonal angle ranges through the **whole open interval `(0,pi)`** continuously and bijectively, even though the geometric rectangle-locus itself is an open interval with no essential closed loops.

At `c=0`, the rectangle is the unique square with side squared

`4a²b²/(a²+b²)`.

This confirms that an angle-continuation strategy can work without odd-monodromy loops. It does not prove that *arbitrary Jordan curves* have the same angle support or a size-protected branch.

## 2. Sharp positive-size protection near 90 degrees

Fix a rational `0<=eta<1` and limit to angles `|cos(theta)|<=eta`. The denominator obeys

`D(c) <= a²+b²+|a²-b²| eta`.

Therefore the squared circumradius has the *exact attained* lower bound

`rho² >= 2a² b² / [a²+b²+|a²-b²| eta] > 0`.

The squared lengths of the two sides are `4x²,4y²`; since `x²(c)` decreases and `y²(c)` increases strictly with `c`, the sharp minimum over the angle band is

`4 min(x²(eta), y²(-eta)) > 0`.

For `a=3`, `b=4`, and `eta=1/2`, the exact bounds are both `192/19` for squared circumradius and the smaller squared side (in this numerical case they coincide). These values are verified with rational arithmetic, including adversarial endpoints and rejection of degenerate angles `c=±1`.

## 3. A conditional square theorem useful beyond ellipses

Let `gamma_n` be Jordan curves converging uniformly in compatible parameterizations to a Jordan curve `gamma`. Suppose there are genuine inscribed rectangles `P_n` whose diagonal angles `theta_n -> pi/2` and whose squared circumradii satisfy `rho(P_n)² >= epsilon>0`, uniformly in `n`. If all curves are contained in a common bounded region, compactness of the four ordered corner tuples gives a convergent subsequence on the limit curve.

The radius floor prevents all-corner collapse. A rectangle with diagonal angle approaching `pi/2` has nonzero diagonal length; because its diagonals remain equal and become perpendicular, its limit is a **positive square**. The square assertion does not require the rectangle's individual side floors: any two-pair flat has diagonal angle `0` or `pi`, not `pi/2`.

Thus **size-protected rectangles with angles approaching 90 degrees are sufficient for Toeplitz**. This compactness theorem is elementary and conditional; it is not an existence theorem for the required rectangles.

## 4. Missing global ingredient

Greene–Lobb (2026) prove large *measure* of protected rectangle angles for arbitrary Jordan curves, but this need not force a protected angle sequence converging to `pi/2`. An entire interval on a particular curve, as in the ellipse, cannot be assumed universal. One must establish either

- a protected angle branch/component whose closure contains `pi/2`, or
- protection in every neighborhood of `pi/2` with a **single uniform positive size floor**.

The current ellipse formulas provide a reference calibration and exact unit tests. They do **not** supply a geometric angle-gap theorem for arbitrary Jordan curves.

## References

- Greene–Lobb (2026), [arXiv:2604.17116](https://arxiv.org/abs/2604.17116), Theorem A, Corollary 2.2.
- [v268 rectangle-locus classification](v268-rectangle-monodromy-and-protected-square-criterion.md).

Code and tests: [`code/relative_rectangle_v269_v273.py`](../code/relative_rectangle_v269_v273.py), [`tests/test_relative_rectangle_v269_v273.py`](../tests/test_relative_rectangle_v269_v273.py).