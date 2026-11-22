# v318–v322. Support-free angular spectral noncollapse

Benjamin John Schulz · November 2026

**Status:** proved compactness and spectral-envelope reductions, a sharp obstruction to the scalar-budget approach, a quantitative angular budget, and two conditional routes to the existing odd-local-degree target. **No universal positive square-size estimate or proof of Toeplitz is claimed.** This replaces the *search for four affine support lines* with an intrinsic analytic question; it does not solve the new question.

## Inputs and conventions

Let Gamma be a Jordan curve, D its bounded complementary domain, A=area(D)>0, R=diam(Gamma)/2, K=R² and r=inradius(D). The angle t lies in (0,pi); a rectangle has vertices c±u,c±v with |u|=|v|=q, v=Rot_t u. At t=pi/2 its square side is sqrt(2)q. We do not impose cyclic ordering of the vertices on Gamma. The area means the area of the open interior, not the area of the closure (the boundary itself can have positive area).

External input: Greene and Lobb, [Jordan curves inscribe a positive measure of rectangles](https://arxiv.org/abs/2604.17116), v1, 18 April 2026, Section 2 and Theorem 5.1. For generic smooth Jordan curves the upper-grading spectral function ell is nondecreasing, K-Lipschitz, subadditive for sums <=pi, has endpoints 0 and A, and at almost every t its derivative equals q² for some inscribed rectangle of angle t. We use the derivative property, not an assertion that its carrier is the largest rectangle or that a single branch persists. See also [Floer homology and square pegs](https://arxiv.org/abs/2404.05179), v2, Lemma 5.1, for the distinct action-based shrinkout argument. Both references are cited as preprints here.

Use registered generic smooth approximations gamma_n with enclosed areas A_n→A. Existence even for a positive-area Jordan boundary: a Riemann map from the disk to D extends homeomorphically to the closed disk by Carathéodory's theorem. Images of circles of radii tending to 1 from below converge uniformly in the induced registration, enclose nested domains exhausting D, and are smooth embedded curves. Approximate each by a generic smooth embedding in the smooth topology with registration and area errors tending to zero (area is continuous in C¹ on smooth loops). Reparameterize by the boundary correspondence if a different registration is given. One can also rescale each curve by sqrt(A/A_n) to have exact area A if needed, since scaling preserves genericity of the rectangle and binormal configuration spaces. None of this assumes rectifiability of Gamma. Arbitrary C⁰ approximations without area control are not substituted for this construction.

## Relation to earlier versions and to the cited papers

Several statements below restate or sharpen earlier results; they are not new in substance.

- **Central slope gives a square (v318, v322).** This is the transfer theorem of [v277](v277-localized-spectral-slope.md). Upper semicontinuity of the envelope sharpens v277's side bound from sqrt(c) to sqrt(2m).
- **Scalar no-go (v319).** [v270](v270-spectral-axiom-no-go.md) already showed that monotonicity, the Lipschitz cap, subadditivity and the endpoint values do not force a positive derivative at pi/2. Its countermodel is strictly increasing with an isolated zero derivative there. v319 adds a model with a flat central interval and uses it to show that the v320 window bound is sharp.
- **Angular budgets (v318, v320).** The bound integral_0^pi M_Gamma² >= A is a size-weighted form of the Greene–Lobb proof of Theorem A (arXiv:2604.17116, §2). The window bound repeats that argument on a central interval. The prefix bound L(t) >= A/ceil(pi/t) is the subadditivity step in their proof of Theorem B.
- **Action gap and length payment (v321B).** The action gap for A > pi K/2 is the case theta = pi/2 of the Greene–Lobb argument for Theorem B. The length payment quantifies the cap estimate in their Lemma 5.1 (arXiv:2404.05179). They note in arXiv:2407.07798, §1, that their earlier work already turned an action gap into a rectangle-size bound depending on length.

What is new here is the upper-semicontinuous envelope M_Gamma with the domination inequality of Theorem 318.2, the plateau model and its sharpness, and the explicit form liminf rho_n T_n >= eta of the length payment.

## v318. Largest-radius envelope and domination theorem

For t in (0,pi), define M_Gamma(t) as the largest q for a t-rectangle on Gamma, allowing q=0 so the compact feasible set is never empty. The optimizing center and vectors live in a common compact box; all four vertices must lie on Gamma. Thus the maximum exists and lies in [0,R].

**Theorem 318.1.** If compact curves Gamma_n converge in Hausdorff distance to Gamma and t_n→t in (0,pi), then

    limsup M_Gamma_n(t_n) <= M_Gamma(t).

Proof: choose maximizers along a subsequence realizing the limsup. The center lies in the convex hull of the bounded union, and radii are bounded by half the diameters. Extract converging centers and diagonal vectors. The equal-length and rotation equations are closed, and each limiting vertex is on Gamma. A positive limit is a nondegenerate rectangle because t is strictly interior. The zero case is immediate. This also proves upper semicontinuity and Borel measurability of M_Gamma.

Let ell_n be the Greene–Lobb spectral functions of the registered generic approximations above. Uniform boundedness and equicontinuity give subsequential uniform limits L on [0,pi]. Every such limit is absolutely continuous, nondecreasing, subadditive, K-Lipschitz, with endpoints 0,A. We assert neither uniqueness nor continuity of L with respect to arbitrary changes of Gamma.

**Theorem 318.2 (envelope domination).** For every interior interval [a,b],

    L(b)-L(a) <= integral_a^b M_Gamma(t)^2 dt.

Proof: for almost every fixed t the derivative assertion holds simultaneously for the countably many ell_n (intersect their full-measure sets). Then ell_n'(t)<=M_Gamma_n(t)^2. Theorem 318.1 gives the pointwise limsup bound. Since all derivatives are bounded by one common constant, reverse Fatou gives the stated inequality by passing to the subsequence in ell_n(b)-ell_n(a). Consequently 0<=L'<=M_Gamma² almost everywhere, by differentiation of interval integrals. Taking limits at the endpoints gives the full angular budget integral_0^pi M_Gamma² >= A. The latter is a size-weighted version of the positive-measure input.

For any interior t0, upper semicontinuity gives

    limsup_(h→0+) [L(t0+h)-L(t0-h)]/(2h) <= M_Gamma(t0)^2.

Indeed for any eta>0, all M_Gamma(t)^2 in a sufficiently small neighborhood are at most M_Gamma(t0)^2+eta. If Gamma has no square, M_Gamma(pi/2)=0, so L is differentiable at pi/2 with derivative zero (both one-sided quotients are nonnegative and bounded by twice the central difference). More strongly, the angle-local Lipschitz constants of L tend to zero there. By the same compactness argument, in a sufficiently small angular neighborhood all rectangles on all sufficiently fine approximants have q<epsilon, for any chosen epsilon>0. These are necessary conditions for square-freeness, not contradictions yet.

A positive central limsup m for *any* L supplies a square with side >=sqrt(2m). It also supplies an approximant selection and certificate as in v321. No hypothesis about support lines, monotone arcs, coarse comparisons, square parity, or length enters 318.1–2.

## v319. Sharp scalar-axiom obstruction

We test whether the scalar spectral properties alone force a positive central derivative. Work in normalized coordinates x=t/pi, f(x)=L(pi x)/(pi K), b=A/(pi K) in (0,1]. Thus f(0)=0, f(1)=b, and the slope cap is one. Let a=b/2 and define

    f_b(x) = min(x,a) + max(0,x-(1-a)).

For b<1 this function has a flat central interval [a,1-a] of positive length, slopes 1,0,1, and symmetry f_b(x)+f_b(1-x)=b. It obeys every scalar property in v318, and even this extra reflection identity (we do not assume that a single actual upper-grading invariant satisfies the identity).

**Subadditivity proof.** Put g=x-f_b(x)=min(max(x-a,0),1-2a). We prove g(x+y)>=g(x)+g(y) for nonnegative x,y with x+y<=1. If either g vanishes, monotonicity suffices. Otherwise x,y>a. If either x or y is >=1-a, the other is <=a, contradiction. Both are therefore on the rising segment, and g(x)+g(y)=x+y-2a. This is at most x+y-a and at most 1-2a, hence at most g(x+y). Subtract from x+y to obtain subadditivity.

For 0<lambda<=1 define f_(b,lambda)=(1-lambda) f_b +lambda bx. All the same scalar constraints hold, and these functions are strictly increasing, with central slope lambda b→0. Also f_(b,lambda) has values strictly between 0 and b in the interior. Hence imposing strict monotonicity at each approximating stage, an interior action value at the square angle, or reflection symmetry still does not give a uniform positive central slope from scalar data alone. These are abstract functions, not asserted to be Floer spectral invariants of curves. The example does not disprove Toeplitz, a geometric inradius inequality, or a constraint involving actual rectangle carriers.

## v320. Quantitative angular budget and interior disk audit

For a measurable angular interval I=[a,b] of length w>0, endpoints and the global slope cap give

    integral_I M_Gamma² >= L(b)-L(a) >= max(0, A-K(pi-w)).

In particular on I_h=[pi/2-h,pi/2+h], 0<h<pi/2,

    max_(t in I_h) M_Gamma(t)^2 >= max(0, K-(pi K-A)/(2h)).

The maximum exists by upper semicontinuity. If the right side is positive, this is an explicit lower radius bound for a near-square, not necessarily an exact square. The threshold h>(pi-A/K)/2 cannot be improved from the scalar constraints alone: the f_b model has precisely the saturated missing central interval. For fixed t in (0,pi], subadditivity and monotonicity give L(t)>=A/ceil(pi/t), hence integral_0^t M²>=A/ceil(pi/t); this recovers a quantitative prefix budget, not a central one. The corresponding limits at endpoints are understood via interior intervals.

For an interior disk of radius r, A>=pi r². Thus the same lower bound follows with A replaced by pi r² (and R fixed). An interior disk persists under registered approximations, as shown in the general noncollapse goal. But the resulting bound is only positive for a fixed nonzero angular window whenever r<R. It cannot be sent to h=0. The v319 models survive even with b arbitrarily close to one. At equality A=pi R², the isodiametric equality case is a circle anyway. No local lower bound has been derived from disk inclusion alone.

**Necessary deficit inequality.** If M_Gamma(t)<=q0 on a central window I_h, then

    pi K-A >= 2h(K-q0²).

This is a falsifiable restriction on a candidate collapsing spectrum. There is no contradiction for an arbitrarily narrow window when the area deficit is positive.

## v321. Near-angle repair and action-payment connection

### 321A. Transfer positive angular slope to the actual odd local degree target

Suppose a uniform subsequential limit L satisfies a central limsup >=m>0. Fix any 0<m0<m. Choose h_j→0 with central increments per 2h_j >m0, allowing a strict margin by taking m0 smaller than m. Uniform convergence lets us pick increasing n_j with registered errors tending to zero and the same strict increment inequality for ell_nj. Absolute continuity and the derivative property then give a carrier t_j in I_hj of radius q_j with q_j²>m0. No interchange of the h and n limits is allowed: the stage index is chosen *after* the window and with spectral uniform error o(h_j).

Rotate the second diagonal by pi/2-t_j while fixing the first, center and their equal lengths. The resulting square has side sqrt(2)q_j, and the maximum vertex displacement is

    2q_j sin(|t_j-pi/2|/2) <= q_j |t_j-pi/2|.

It tends to zero. Thus all four corners are in a shrinking tube around Gamma (not just one vertex). Take registered polygonal approximations to Gamma of error tending to zero; choose the four contacts with sup-distance error e_j→0. Eventually sqrt(2)q_j>16e_j, so the explicit v315 PL homeomorphism places all four exact square corners on the polygon. Apply the marked polygon lemma in general-noncollapse-goal.md to regularize this square, with further displacement tending to zero. This supplies odd *local* degree above any fixed cutoff <sqrt(2m0) on arbitrarily fine approximations. No parity of the entire square census and no support box on Gamma are required. The local isolating boxes are built after finding the carriers. A spectral estimate is sufficient even if a global degree above the cutoff equals zero.

This uses the tube-witness repair part of v315, not a claim that angular spectral variation is equal to the mod-two degree or that an arbitrary homotopy's wall zeros are near the target. v310's fixed-cutoff instability result remains an independent sufficient route.

### 321B. A different necessary payment for action-carrying shrinkout

The spectral derivative and spectral action are distinct: the derivative controls radius; the action may remain finite even as the radius vanishes. We quantify the obstruction in Greene–Lobb 2024, Lemma 5.1, for a sequence of shrinking action-carrying rectangles of a fixed interior angle t. Use its cap and projection construction, with areas A_n and action considered modulo A_n.

Assume registered uniform approximations to an embedding. If q_n→0, the four parameter contacts cluster in a shrinking connected arc of the parameter circle (uniform continuity of the inverse of Gamma on its image). Let C_n denote the corresponding short-in-parameter boundary subarc spanning all four contacts. It is spatially contained in a disk B_n of radius rho_n→0; enlarge B_n by a vanishing amount to include the two circular orbit arcs. Let T_n be the arclength of C_n. As in the cited cap argument, each of the two projected paths can be reduced modulo a full traversal of the Jordan curve to a subarc of C_n, followed by the relevant circular orbit arc. The cap correction around the diagonal changes its area by an integer multiple of A_n, which is invisible modulo A_n. Each subarc has length <=T_n, and each circular orbit arc length <=pi q_n. Use the primitive ( (x-p_x)dy-(y-p_y)dx )/2 with p the disk center; this differs from x dy by an exact one-form and is bounded in norm by rho_n/2 on B_n. The Hamiltonian term is t q_n². Hence

    dist(action_n, A_n Z) <= t q_n² + rho_n (T_n + pi q_n).

If the action remains a distance >=eta>0 from A_n Z, it follows that

    liminf rho_n T_n >= eta.

This is a quantitative local length payment for collapsing rectangles (conditional on the cap/action setup), not just an observation that lengths become unbounded. It implies an action-based noncollapse criterion if one can construct approximants with rho_n T_n→0 for every potentially collapsing carrier sequence and with a uniform action gap at t=pi/2. Length may diverge globally under that criterion. It does not follow from uniform approximation: a simple arc can execute arbitrarily many tight turns in a small disk. A bound on length times disk radius or a sharper *signed* projection-area cancellation is still required.

The action gap itself is separate. At t=pi/2, subadditivity gives ell_n>=A_n/2, while the slope cap gives ell_n<=pi R_n²/2. Hence a fixed interior action gap follows when A>pi K/2. For general A, the upper gap does not follow from these estimates. This local payment argument extends the form of the known bounded-length estimate; it is not an independent derivation of the Floer cap or an unconditional shrinkout theorem.

## v322. Closure audit and next target

We now have a support-free sufficient condition, with explicit quantifiers:

    For each Gamma, construct generic area-convergent registered gamma_n
    and a subsequential spectral limit L with
    limsup_(h→0+) (L(pi/2+h)-L(pi/2-h))/(2h) >= m(Gamma)>0.

Then v318 gives a square with side >=sqrt(2m(Gamma)), and v321A constructs arbitrarily fine registered simple polygonal approximations with odd local degree above any smaller fixed cutoff. If m(Gamma)>=c r² with absolute c>0, the quantitative target is met (with side coefficient sqrt(2c)). This is not shown by v318–321. An alternative sufficient input is a sequence of positive-radius near-square rectangles from *any* construction, with angles→pi/2 and all contacts in the shrinking target tube; it need not be a spectral carrier.

A uniformly positive derivative is stronger than square existence and is not known to be necessary even if a curve has large squares. A flattening spectral limit therefore does not prove a curve is square-free. Do not interpret v319 as a curve counterexample, or v320's near-square estimate at fixed h as a vanishing-error tube witness.

### Specific next proof obligations

1. Derive an actual geometric carrier constraint unavailable to the scalar plateau model. The envelope inequality gives a pointwise upper bound, so the new input must give central mass from the geometry or topology of the rectangle space, not merely another scalar endpoint estimate. A chain-level continuation method localized away from the diagonal would need to certify that a positive-radius portion of the class survives at the square angle.
2. For the action alternative, use winding of the two projected caps jointly, rather than bounding them independently by arclength; prove a uniform cancellation or identify an actual Jordan obstruction. Control the upper action gap separately. Do not equate action with q².
3. If borrowing v315's cutoff ledger, show that all unpaired crossings have all four vertices in shrinking target tubes, rather than assuming a certified reference box. The angular repair theorem removes this requirement only if positive-scale near-square carriers are independently produced.
4. Compare any conjectured additional spectral constraint with v319 before investing in numerical exploration. Test circles, thin rectangles, asymmetric lobes, and high-winding finite coils. A finite search cannot certify the missing continuum estimate.

**Conclusion.** The coarse support selection is no longer a prerequisite for this route. The universal positive central spectral mass (or another positive-scale near-square selector) is the unresolved analytic core. We have not proved the original general noncollapse target. The unrestricted square-peg problem remains open.

### Reproduction

The companion exact-rational checker uses normalized angular coordinates to verify scalar plateau formulas, subadditivity at all vertices of their finite piecewise-affine subdivision of x>=0,y>=0,x+y<=1, and the angular budget inequalities. It is an algebraic negative control, not a Floer solver or a numerical proof of the geometric estimate.
