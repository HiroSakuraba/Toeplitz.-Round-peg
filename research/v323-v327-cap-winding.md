# v323–v327. Local signed cap area, spiral obstruction, and spectral carrier localization

Benjamin John Schulz · November 2026

**Scope:** exact cap identity, a rational shrinking-square/nonvanishing-action counterexample, a conditional quantitative carrier localization theorem, and a two-obstruction reduction. No general noncollapse result or universal central-mass bound is claimed. In particular the example is about the action of an individual rectangle generator, NOT its selection as a Floer spectral carrier.

Inputs: Greene–Lobb, [Floer homology and square pegs](https://arxiv.org/html/2404.05179v2), §§2–5, especially Lemma 5.1; [Jordan curves inscribe a positive measure of rectangles](https://arxiv.org/html/2604.17116v1), §2, Theorem 5.1 and its proof. We use the conventions H(z,w)=|z-w|²/4, positive rotation by t of the relative coordinate, and action equal to integrated Hamiltonian minus cap area. Write A for the positive area of the bounded domain, not the signed area of the chosen registration. All action identities below are **modulo A** unless explicitly stated. For generic smooth curves use the actual real action and spectrality from those papers; a PL action residue is the elementary cap expression, not an assertion of PL Floer theory.

## v323. Exact two-arc reduction, with no length estimate

Consider a nondegenerate t-rectangle, 0<t<pi. Its diagonal endpoints p1,p2 rotate to p1',p2', about c with radius q. Choose two paths E_i on the Jordan boundary returning p_i' to p_i. For a smooth Jordan loop, use the cap construction of Greene–Lobb Lemma 5.1: paths in the deleted product may be chosen by compatible lifts, and adjusting their winding about the diagonal changes the symplectic area by a multiple of A (the core contributes 2A). Arbitrary choices of the individual projected return paths differ from these by integer traversals of the boundary, also invisible modulo A. Thus we may evaluate using the short-in-parameter arcs without requiring their simultaneous lift to avoid the diagonal.

Let lambda_c=((x-c_x)dy-(y-c_y)dx)/2. The integral over each forward circular trajectory is tq²/2. The two together cancel the Hamiltonian integral tq² **exactly**. Stokes on the two projected cap boundaries yields

    action = - integral_E1 lambda_c - integral_E2 lambda_c   (mod A).       (323.1)

The sum is unchanged by translating the primitive's origin, because the total endpoint displacement is zero. Equivalently, in terms of a boundary primitive F(s)=integral from a fixed basepoint to gamma(s) of lambda_0 on a parametrization lift,

    action = F(s1')+F(s2')-F(s1)-F(s2)   (mod A).

F(s+1)-F(s) is the signed area ±A. So no branch of the lift is needed for the residue. This formula is not a new construction of the diagonal-avoiding cap: it is an exact reduction using the cited construction.

Close E_i by the straight segment from p_i back to p_i'. Let W_i be the winding function of this closed loop on its complement (values on the trace are irrelevant for smooth/PL paths). Its signed area S_i=integral W_i dArea equals integral lambda_c around that loop, by the winding-number form of Green's theorem. Each closing chord contributes q² sin(t)/2. Therefore

    action = q² sin(t) - S_1 - S_2   (mod A).                            (323.2)

If both return paths and closing chords lie in a disk B(c,rho), their winding functions vanish outside the disk. We obtain the bound

    dist(action,A Z) <= q² + integral_B |W_1+W_2| dArea.                 (323.3)

In particular if |W_1+W_2|<=m almost everywhere, the right side is <=q²+pi rho² m. This bounds *joint signed winding*, not length, total enclosed area, nor the sum of unsigned path lengths. A length-free, sufficient finite test: if each chord meets its return arc only at its endpoints, the closed loop is Jordan and |W_i|<=1 off the trace; then dist(action,A Z)<=q²+2pi rho². More generally, if there are finitely many distinct transverse interior chord intersections, split the arc at these hits. Each piece with its own straight closing segment bounds a simple loop, unless the piece is a straight chord (zero area). Signed areas add, because all closure segments lie in the same line and telescope as oriented chains. With at most m_i interior hits on each arc, |W_i|<=m_i+1 almost everywhere. Collinear overlaps and infinitely many intersections are not covered by this finite test.

For registered uniform approximations gamma_n→Gamma, define the inverse modulus

    nu(d)=sup{circular parameter distance(s,t): |Gamma(s)-Gamma(t)|<=d}.

It tends to zero as d→0 by compactness and embeddedness. If the approximant registration error is e, rectangle contact parameters all lie in a parameter arc of length <=2 nu(2q+2e) provided nu<1/4 (center the arc at one contact). Let omega(s) be Gamma's forward modulus. The short return paths are contained in the disk centered at c of radius rho=q+2e+omega(2 nu(2q+2e)). The closing chords are there by convexity. Thus rho→0 as q,e→0. It follows from (323.3) that a sequence of collapsing rectangles with action residues at distance >=eta>0 from A_n Z must satisfy

    esssup |W_1+W_2| >= (eta-q²)/(pi rho²),                             (323.4)

and integral |W_1+W_2| >= eta-q². In particular the essential multiplicity tends to infinity. This is the precise winding payment complementary to the v321 length payment. It does not follow from Jordan simplicity that the two return paths cancel.

## v324. Exact interleaved diamond spiral: joint cancellation is false

We give shrinking elegant squares on simple rational polygons converging in registration to the boundary of [0,4]². Their two-return-path cap action stays in the interior modulo the enclosed area. Thus even area convergence, a persistent interior disk, and cyclic contact order do not imply uniform cap cancellation for all rectangle generators.

### Unscaled local gadget

Fix an integer N>=1, m=4N, d=1/(16N). Let e_j be the periodic sequence (1,0),(0,1),(-1,0),(0,-1) indexed modulo four. Define

    r_j=3/4-j d
    u_j=r_j e_j                       (0<=j<=m)
    v_j=(r_j-2d)e_j                   (-1<=j<=m-1).

The ordered arc is

    (1,0), u_0,...,u_m, (0,0), v_(m-1),...,v_(-1),
    (0,-1),(-1,-1),(-1,1),(0,1),(1,1).

Its four square contacts, in curve order, are (1,0),(0,0),(0,1),(1,1). Its endpoints are on x=1 and its remaining points are strictly to the left of x=1 except the first portion of the last horizontal segment's endpoint. The four square vertices are in cyclic geometric order. The two spirals do not intersect: inside each open quadrant their straight segments have strictly ordered endpoints on the two bounding rays (successive turns on a ray differ by 4d, and the second spiral is offset by 2d). Chords with that ordering cannot cross; rays shared between quadrants are met only at the listed distinct vertices. The radial entrance on the positive x axis is outside all turns; the inner positive x segment ends at radius 1/2, below every outbound positive x contact. The inner negative y segment reaches radius 1/2-d, below all inbound negative y contacts. The outer negative y segment starts at radius 3/4-d, above all inbound negative y contacts. The remaining outside path stays outside the radius-3/4 diamond. These checks also preclude self-intersection of each spiral. All segments have nonzero length.

Choose the starting diagonal ((0,0),(1,1)); positive quarter rotation gives ((1,0),(0,1)). One return path follows the full inbound spiral and its radial ends, the other is the upper horizontal segment. The primitive lambda_0 on a straight segment a→b integrates to det(a,b)/2. Consequently the joint return integral is

    J_N = (1/2) sum_(j=0)^(m-1) r_j r_(j+1) - 1/2
        = 19N/24 - 1/(384N) - 1/2.                                  (324.1)

Derivation: r_0=3/4, r_m=1/2. The mean integral of the squared linear profile is 19/48. The trapezoid sum of r² exceeds m(19/48) by m d²/6; using r_j r_(j+1)=(r_j²+r_(j+1)²-d²)/2 gives the stated sum and formula.

The reverse winding mostly cancels the inbound winding in the **full Jordan boundary integral**, but is not part of either selected return arc. This is the failure of joint cap cancellation. The signed primitive has large oscillation between the square contacts.

### Scaling and registration

For integer k>=2 choose N=k² and scale by eps=1/k; translate by (-eps,0). The gadget replaces the left-hand side segment (0,0)→(0,eps) of the clockwise boundary of [0,4]², with all other gadget points at x<0 except its endpoints. Close from (0,eps) through (0,4),(4,4),(4,0) to (0,0). The result P_k is a simple polygon. The modification lies in [-2eps,0]×[-eps,eps]. Register the whole gadget against the replaced segment, with endpoints matched by a piecewise-linear parametrization; the uniform error is <=sqrt(8)eps (the coordinate sup bound is <=2eps). The interiors coincide outside the modification box, by winding; A_k→16, with |A_k-16|<=4eps². The base interior disk of radius 1 persists. The marked square has side eps and q²=eps²/2.

By (323.1) and (324.1) its action residue is -eps² J_(k²) modulo A_k, with

    eps² J_(k²)=19/24 - 1/(2k²) - 1/(384k⁴)  → 19/24.

For k>=2 this is >1/2 and <1, whereas A_k>=15, so the distance to A_k Z is the displayed positive number. A complete Jordan boundary has winding 0 or ±1 away from it, even as the joint return winding grows as order k². The marked square is *not* asserted to carry either spectral invariant.

### Smooth and fixed-target versions

For each finite polygon one can regularize the marked four contacts using the marked-polygon lemma in general-noncollapse-goal.md, with disks arbitrarily small. For a fixed k the signed primitive along the two tagged paths changes arbitrarily little: the old and replacement arcs in each finite disk have length tending to zero (choose the connectors of bounded complexity for that contact), so their closed difference has signed integral tending to zero. Smooth the finitely many other corners within disjoint disks; preserve the four straight contacts and their nonsingular square matrix. Real analytic embeddings are C¹-dense in smooth embeddings by periodic Fourier approximation; the regular marked square persists by the implicit function theorem, with contact tags converging, so its action residue is continuous under these approximations. For each k choose errors in registration, area and action less than 1/k², or smaller if needed for positive scale and action gap. Genericity at the square angle can be obtained by a further sufficiently small perturbation; the regular square persists. This supplies the same obstruction for smooth curves, not just PL residues. Scale each smoothed loop to area 16 if exact area is required; scale factors tend to one.

The phenomenon can also be placed on one fixed nonrectifiable Jordan curve, by putting independent gadgets of size eps_j=2^(-j-3) and N_j=eps_j^(-2) along the left side at heights 2^(-j), j>=1. Their modification boxes are disjoint and converge to (0,0); concatenate the unchanged edge parts and the gadgets in vertical order. The limiting parametrization is continuous and injective on the circle because each gadget has distinct contacts, stays within its own box, and diameters tend to zero at the only accumulation point. Finite truncations are registered approximations; their interiors differ from the limit only inside a union of boxes with total area tending to zero. On each truncation use its smallest gadget's square, then regularize and smooth as above with sufficiently small errors. The action gaps again stay >1/2 for all sufficiently large j. This version does not rule out some other better approximation sequence; nor is it a square-free curve (it has the marked squares and many more). An obstruction to an arbitrary-generator bound is not an obstruction to a spectral-carrier-specific theorem.

## v325. Positive-scale spectral-carrier localization without chain-level overclaim

For a given generic smooth approximant gamma_n and an angular interval I, define the action-exclusion modulus

    E_n(s,I) = sup dist(action(P), A_n Z)

over actual inscribed t-rectangles of half-diagonal <=s and t in I; assign 0 to the empty supremum. This is a small-radius geometric quantity, expressible by (323.1) or (323.2), not the spectral function. Bound E_n by s² plus a uniform joint winding-area bound in (323.3) when such a bound is available.

Suppose for some fixed s,eta>0 and a nontrivial closed I inside (0,pi), for arbitrarily fine selected approximants,

    E_n(s,I)<eta,   eta<=ell_n(t)<=A_n-eta for all t in I.              (325.1)

Spectrality implies every action carrier of ell_n(t) has radius >s on I. At the full-measure regular differentiability angles, the proof of Greene–Lobb 2026 Theorem 5.1 identifies a *single* rectangle that realizes ell_n(t) as action and ell_n'(t) as squared radius (not just two unrelated rectangles). Hence ell_n'(t)>s² almost everywhere on I. For any a<b in I,

    ell_n(b)-ell_n(a)>=s²(b-a).

In a uniform subsequential limit the same inequality holds. In particular if I contains pi/2 in its interior this gives v322's central mass >=s², and a square of side >=sqrt(2)s on Gamma, with odd local degree via v321A. Even if I is only a singleton {pi/2}, spectrality and compactness give the same square-size conclusion (not a derivative conclusion). In this way action localization can solve the actual square problem without first proving central mass; a window upgrades it to the stronger slope statement.

There is also a direct necessary property of any representing Floer cycle: when the spectral value lies in the stated action interval and a cycle attains it, at least one generator of maximal action in that cycle is outside the small-radius set. A filtered chain complex exists in the external smooth theory. Nothing here proves the span of small-radius generators is a subcomplex, or asserts that the full Floer class survives in a quotient by that span. Radius does not itself define an action filtration; a differential may connect generators of different radii. A positive-radius continuation theorem requires additional strip control.

## v326. Two independent obstructions and a quantitative regularity criterion

The duality between the degree-1 and degree-2 invariants relates the two values at the square angle by sum A (Greene–Lobb 2024 §4; proof in [Square pegs between two graphs](https://arxiv.org/abs/2407.07798)). This does not supply a uniform interior gap for general curves.

Indeed in the normalized scalar setting of v319 with area ratio b<=1/2, the upper model f(x)=min(x,b) and lower model g(x)=b-f(1-x)=max(0,x-(1-b)) obey upper subadditivity, monotonicity, slope cap, endpoints, lower<=upper and duality. At the square angle upper=b, lower=0. Strict approximating models f_lambda=(1-lambda)f+lambda b x and their duals still obey these inequalities. Their gap goes to zero at x=1/2. Thus even using both scalar values cannot establish the interior gap. The earlier symmetric plateau example models the different failure: a positive interior gap *with* zero limiting central derivative. No claim is made that these scalar functions arise from curves.

A useful positive-scale criterion on other approximants follows directly from (323.1): if the centered primitive oscillation on subarcs with physical diameter <=rho tends uniformly to zero, all collapsing action residues vanish. Equivalently it suffices to bound the chord-closed area modulus uniformly by a function w(rho)→0; use (323.3) or the triangle inequality in (323.2). Combined with an interior spectral gap at pi/2 this precludes shrinkout. On an angle interval with a uniform gap, choose s and at sufficiently fine stages e so that s²+2w(rho(s,e))<eta. Then (325.1) gives a quantitative derivative and side bound. This is not restricted to finite length.

Check against current literature: Asano–Ike [Liouville regular Jordan curves](https://arxiv.org/abs/2609.11598) (2026 preprint) and Li–Pan [Rectangular Pegs on Jordan Curves of Finite p-Variation](https://arxiv.org/abs/2609.04918) (2026 preprint) obtain stronger aspect-ratio conclusions for classes with suitable controlled primitives. Do not claim novelty or general rectangle existence for those regularity classes from our criterion. The precise new purpose here is the winding-payment test and its relationship to the angular derivative carrier, and the explicit counterexample to an unconditional all-generator estimate.

## v327. Universal noncollapse attempt: a dichotomy usable for future tests

Use area-convergent generic approximants of any Jordan Gamma as in v318, with a uniform spectral subsequential limit L. Write M(t) for the maximal half-diagonal defined in v318. Suppose Gamma has no square. Then M(pi/2)=0; all rectangles on fine approximants with angles tending to pi/2 have radii tending to zero. By v318 the central slope of L is zero. By spectrality choose square action carriers P_n, or at regular angles theta_n→pi/2 choose the simultaneous action/slope carriers of Theorem 5.1. Along the selected subsequence, their real actions tend to L(pi/2), by the uniform Lipschitz cap. Subadditivity gives L(pi/2)>=A/2. Two possibilities:

* Endpoint saturation: L(pi/2)=A (dual lower spectral value tends to zero). The action residue can tend to zero; even perfect local winding cancellation would not contradict this possibility. The capped scalar model in v326 admits it.
* Interior action: A/2<=L(pi/2)<A. Put eta=min(L(pi/2),A-L(pi/2))/2>0. For all sufficiently fine selected stages the carriers have action distance >=eta from A_n Z. Their joint winding areas must satisfy (323.4) in shrinking contact disks. In particular their essential winding multiplicities must diverge at least on the scale rho_n^(-2). The v324 spiral demonstrates why this is not a contradiction using Jordan simplicity alone.

For A>pi R²/2, the spectral Lipschitz cap excludes the first possibility, leaving only the second. For arbitrary A, both must be addressed. The dichotomy is conditional on absence of a square and makes no assertion that one can realize the first or second behavior by an actual spectral carrier. Excluding both for at least one carefully selected area-convergent approximation sequence would prove Toeplitz; obtaining a bound E_n(s,I) and an action gap independently of a known square would meet the quantitative v318–322 goal.

**Research decision:** do not try to prove universal joint cap cancellation for arbitrary rectangles: v324 refutes that assertion with elegant squares. The next meaningful issue is a *spectral selection* constraint distinguishing homology-carrying rectangles from local high-winding generators, or an approximation rule preventing the relevant residues without requiring Liouville regularity. For the endpoint alternative, a separate geometric bound on the upper spectral value below A is needed. A chain-level proposal must exhibit an actual filtered continuation or prove localization of the relevant differential, not merely remove small generators from the complex.

### Exact replay

code/cap_winding_v323_v327.py computes signed area sums with Fraction, checks the entire polygon for self-intersection for moderate N, evaluates the closed form at larger k without generating the polygon, verifies the winding multiplicity at a rational point, and audits the two scalar endpoint models. The finite verifier is not a Floer solver, nor a proof of a continuum universal estimate. The analytic simplicity/registration proof above handles arbitrarily fine stages.
