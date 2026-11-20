# v270: The standard spectral axioms do not themselves force radius at the square angle

**Program:** Benjamin John Schulz  
**Date:** 20 November 2026  
**Status:** Exact analytic *axiom-model no-go*, not a Floer spectral computation or a counterexample Jordan curve.

## 1. Source and precise question

Greene–Lobb, [*Jordan curves inscribe a positive measure of rectangles*](https://arxiv.org/html/2604.17116v1), Theorems A/B and §2, establish a monotone increasing spectral invariant `ell(gamma,theta)` with endpoint values `0` and enclosed area `A`, an `R²` Lipschitz bound (where `2R` is the curve diameter), and the triangle inequality `ell(x+y) <= ell(x)+ell(y)` for `x+y<=pi`. For generic smooth curves the derivative is, at almost every angle, the squared radius of some inscribed rectangle. The paper proves size-protected rectangle angles on a positive-measure set. These statements **do not** assert that the derivative is positive at `theta=pi/2`.

We ask whether the three *elementary analytic* properties (strict monotonicity, Lipschitz control, triangle inequality), plus the endpoint values, already force `ell'(pi/2)>0`.

## 2. Countermodel satisfying all these axioms

Fix `R>0` and `A` with `0<B=A/R²<=pi/2`. Define, for `0<=theta<=pi`,

`L(theta) = (A/pi) [theta + 1/2 sin(2 theta)]`.

Then:

1. `L(0)=0`, `L(pi)=A`. Its derivative is
   `L'(theta)=(2A/pi) cos²(theta) >= 0`. It is strictly increasing because this derivative is positive in every nonempty open interval apart from an isolated zero at `pi/2`, and hence is a continuous increasing bijection onto `[0,A]`.
2. `0 <= L'(theta) <= 2A/pi <= R²`, so `L` is `R²`-Lipschitz.
3. For `x,y>=0` with `x+y<=pi`, the **exact** trigonometric identity gives

   `L(x)+L(y)-L(x+y) = (2A/pi) sin(x+y) sin(x) sin(y) >= 0`.

   Therefore the spectral triangle inequality holds on its entire specified domain.
4. Yet `L'(pi/2)=0`. There is no positive angle-local radius bound at `pi/2` derivable from these axioms alone.

The proof is analytic for all real angles, not a consequence of finite floating-point tests. Code samples the formulas only as a regression and checks the exact parameter threshold `B<=pi/2` in dimensionless form.

**Logical scope:** This `L` is an *abstract admissible function*. We neither claim nor need that it is realized by a Jordan curve's Floer spectral invariant. The model refutes one inference from the listed axioms; it does **not** refute any additional spectral, geometric, or Floer constraint that could distinguish actual curves from abstract functions.

## 3. Filtered noncollapse and what is missing

Greene–Lobb's Corollary 2.2 yields, for any `epsilon>0`, a set of angles of measure at least `(A-pi epsilon)/R²` for which a rectangle can be selected with squared circumradius at least `epsilon`. Global measure or integral control does not locate a radius-protected rectangle at the single square angle.

A viable filtered continuation result would require an extra condition that links *the same nontrivial class* across angle `pi/2`, or a nonescape theorem for a particular relative rectangle component. Merely substituting `L'(pi/2)` into a generic derivative statement would be invalid if the derivative did not exist or if the a.e. spectrality result did not apply at that angle.

**v270 milestone outcome:** Proven no-go for an axiom-only spectral argument. The filtered-Floer conjecture remains open and should be formulated with an additional geometric/relative class, not with an unproved pointwise derivative lower bound.

## References

- Greene–Lobb (2026), [arXiv:2604.17116](https://arxiv.org/abs/2604.17116), §2 and Corollary 2.2.
- Greene–Lobb, [*Floer homology and square pegs*](https://maths.dur.ac.uk/users/andrew.lobb/JF.pdf), for the underlying Jordan Floer program.

Code and tests: [`code/relative_rectangle_v269_v273.py`](../code/relative_rectangle_v269_v273.py) and [`tests/test_relative_rectangle_v269_v273.py`](../tests/test_relative_rectangle_v269_v273.py).