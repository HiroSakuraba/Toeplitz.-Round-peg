# v243: Exact polygonal Big Peg square maxima and disk certificates

**Research program:** Benjamin John Schulz  
**Date:** 19 November 2026  
**Status:** Exact-arithmetic verifier for **rational simple polygons**, with a complete finite maximum-square census and an independently certified contained disk. The reported ten-sided counterexample has **not** yet been verified here because its original coordinates have not been obtained.

## Correction to v242 and source provenance

The v242 numerical search found no counterexample **among the curves it tested**. This must not be interpreted as evidence that the sharp sqrt(2) Big Peg conjecture is true.

Daniel Bernier's original June 16, 2026 prize challenge was updated **August 17, 2026**: he says Grisha Pochuev submitted a valid **10-sided polygonal counterexample** to the sharp conjecture, and he awarded the prize. The original source is https://dbernier.ca/2026/06/16/the-sharp-big-peg-conjecture/ . The organizer's exploratory floating-point checker lives at https://github.com/mariotrevi/strong-big-peg-check .

The public post says the polygon appears below, but the **original 10 rational/algebraic coordinates and its verification artifact were not available as machine-readable text** from the sources inspected for this PR. We do not manufacture, digitize, or attribute substitute coordinates, and cannot claim a verified numerical ratio for that polygon.

The **unrestricted Big Peg problem** is whether there exists *some* universal c>0 such that every Jordan polygon containing a disk of radius r has an inscribed square of side at least cr. The reported example refutes the sharp c=sqrt(2), not every positive c.

## Exact polygonal maximum-square algorithm

For each choice of four polygon edges (i,j,k,l) and an orientation sigma=+1 or -1, parameterize two *adjacent* square vertices by

    P=A_i+tE_i, Q=A_j+uE_j, 0<=t,u<=1.

The other two square corners are forced:

    R=Q+sigma*J(Q-P), S=P+sigma*J(Q-P),

where J(x,y)=(-y,x). Requiring R to lie on edge k and S on edge l imposes two affine equations

    cross(E_k,R-A_k)=0, cross(E_l,S-A_l)=0.

Membership in the closed *segments*, not merely their infinite lines, imposes four additional affine inequalities by projection:

    0 <= dot(E_k,R-A_k) <= dot(E_k,E_k)
    0 <= dot(E_l,S-A_l) <= dot(E_l,E_l).

The two variables t,u are also bounded by [0,1]. Consequently, the feasible locus is a **compact convex rational polytope in two dimensions**, potentially empty, a point, a line segment or a polygon.

The squared side length

    |Q-P|^2=|A_j-A_i+uE_j-tE_i|^2

is a **convex quadratic function** of (t,u). Its maximum on this compact polytope occurs at an extreme point. Therefore, it suffices to enumerate exactly the candidate vertices, including singular cases:

- **Rank 2:** intersect the two affine constraints; accept only if all inequalities hold.
- **Rank 1:** compute the complete line interval from the inequalities and test both endpoints. If inconsistent, reject.
- **Rank 0:** intersect every pair of boundary lines and retain every feasible polygon vertex.

All coordinates and comparisons are rational Fractions. The algorithm checks **2*n^4** ordered charts. Polygon vertices may belong to two edges, edges may be reused, and rank-deficient charts may contain a continuum of squares. These cases are **not silently discarded**. An attained witness square, edge assignment and rational parameters are returned.

This complements rather than supplants the existing v223 complete orbit census. The v243 target is the **largest** square, including the maximum over continuous square families, not parity or counts of isolated square orbits.

## A rigorous inradius lower bound

Given a rational centre C and positive rational radius squared r2, the verifier:

1. Checks Jordan polygon simplicity exactly using v223.
2. Certifies C strictly inside the polygon by rational ray-crossing.
3. Computes exact squared distances from C to every closed polygon edge, including its endpoints.
4. Checks all distances squared >= r2.

The **closed disk** with centre C and radius sqrt(r2) then lies in the polygon closure. This proves r_in^2 >= r2 without claiming to compute the optimal inradius.

If the exhaustive maximum squared square side M2 satisfies

    M2 < 2*r2,

the polygon is a **rigorously verified counterexample to the sharp polygonal Big Peg constant**. Note the directions: a certified *upper* bound on **all** square sides and a certified *lower* bound on inradius. If the inequality does not hold, the program does **not** infer that the conjecture holds for the candidate.

## Tests and trust base

Pure Python standard library, including fractions.Fraction. The test suite exercises full-rank, rank-one, rank-zero, and inconsistent affine systems; exact square/diamond maxima; a triangle cross-check with v223; reversed orientations; exact scaling; invalid Jordan polygons; disk tangency, invalid disks, and a negative-control polygon that cannot be called a counterexample. CI runs all tests with no optional numerical dependencies.

An exact square with corners (+/-1,+/-1) contains the unit disk and has maximum squared side 4, so the sharp counterexample claim must be **false** for this control.

The proof of completeness rests on convex quadratic maximization over each feasible polytope, **not** successful convergence of a numerical search.

## Reproduction

    python code/big_peg_exact_v243.py
    python -m unittest -v tests.test_big_peg_exact_v243

For any independent candidate polygon, create a JSON file with **strings** for every rational coordinate:

~~~json
{
  "name": "illustrative square; NOT the Pochuev decagon",
  "vertices": [["-1","-1"],["1","-1"],["1","1"],["-1","1"]],
  "disk": {"center": ["0","0"], "radius_squared": "1"}
}
~~~

Then run:

    python code/big_peg_exact_v243.py --polygon-json path/to/polygon.json

Decimal strings are interpreted as **exact rationals**. JSON floating-point numbers are refused to avoid concealing approximate coordinates as exact. If no disk is provided, only the complete maximum squared square-side is computed.

## Immediate next steps

**First:** obtain the published Pochuev polygon's ten vertices and original certification materials, preserving source provenance and rational/algebraic uncertainty, then run the exact checker to certify the sharp violation. An image trace is **not** an exact witness.

**Second:** add fast interval or bounding-box pruning and algebraic-coordinate support. The complete O(n^4) edge-chart census is practical for small polygons but is not claimed performant for every 64-gon.

**Third:** use exact verified ratios to study which smaller constants might remain universally valid. A positive universal c, if established in the necessary generality, would help prevent squares from collapsing in polygonal approximations and would be directly relevant to Toeplitz.
