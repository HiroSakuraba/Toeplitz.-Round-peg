# v277 — A local spectral-slope criterion and a geometric source of the bound

**Status.** The transfer theorem is conditional. The uniform secant-cone estimate below supplies its premise for a restricted class of approximating families. No identification with a geometric relative index or sheaf filtration is asserted.

## Spectral normalization and hypotheses

Use the grading-two homogeneous class in Greene–Lobb's Jordan Floer homology, and let \(\ell(g,\theta)\) be its spectral invariant. In [Greene–Lobb, §2 and Theorem 5.1](https://arxiv.org/html/2604.17116v1#S5) (2026 preprint), for an embedding in their open dense smooth class \(\mathrm{Emb}^{\infty}_{\mathrm{reg}}\) from Theorem 3.1, a full-measure set of angles has

\[
\ell'(g,\theta)=R(P)^2=\frac14|z-w|^2
\]

for an associated inscribed rectangle \(P\) with opposite vertices \(z,w\). This \(R(P)\) is the circumradius. It is distinct from \(R(g)\), the half-diameter of the entire Jordan curve which gives the Lipschitz bound for \(\ell(g,\cdot)\). Theorem 5.1 selects a branch by spectrality and the action derivative; it imposes no derivative requirement at exactly \(\pi/2\). We use the entire full-measure set on nearby angles.

## Local spectral slope transfers to a positive square

Let \(g_n:S^1\to\mathbb C\) be embeddings in the class above, converging uniformly as maps to a Jordan embedding \(g\). Write \(\theta_*=\pi/2\). Suppose there are \(c>0\) and \(0<\delta_n\to0\) such that

\[
\frac{\ell(g_n,\theta_*+\delta_n)-\ell(g_n,\theta_*-\delta_n)}{2\delta_n}\ge c
\]

for all sufficiently large \(n\), or along an infinite subsequence. Then \(g\) inscribes a square of side at least \(\sqrt c\).

**Proof.** Each Lipschitz spectral function is absolutely continuous. Its derivative has mean at least \(c\) in the displayed interval, so the set where the derivative is greater than \(c/2\) has positive measure. It intersects the full-measure set of Theorem 5.1. Choose an angle there and the corresponding rectangle, with radius squared at least \(c/2\). Label its corners by a common midpoint and half diagonals. Compactness of the four circle parameters and uniform convergence of the embeddings yield a subsequence of the four vertices converging to vertices on \(g\). The limiting diagonals have the same midpoint and nonzero common half-length, and their angle tends to \(\pi/2\). Hence their endpoints give four distinct points forming a square of side \(\sqrt{2}R\ge\sqrt c\).

One usable variant is a limiting spectrum \(\ell_\infty\) with strictly positive lower central slope at \(\theta_*\), and \(\ell(g_n,\cdot)\to\ell_\infty\) uniformly on a fixed angle neighborhood. For a central interval of radius \(\delta\), replacing the limit function with the \(n\)-th function changes the quotient by at most \(\|\ell(g_n,\cdot)-\ell_\infty\|_\infty/\delta\). Choosing \(\delta_n\to0\) with that ratio tending to zero proves the premise. This checks the rate at the shrinking window and still requires control of the limiting central slope itself.

## A testable geometric source covering all spectral carriers

Let \(d\) be the circular distance on \(\mathbb R/\mathbb Z\). Assume that for some \(0<\eta<1/4\) and \(0\le\alpha<\pi/4\), and for all sufficiently large \(n\), the unoriented directions of every secant between points of \(g_n\) whose parameters lie in any single ball \(B_d(t,\eta)\) belong to a cone of angular half-width \(\alpha\) about one line. The line may depend on the ball and \(n\); the same line must work for all secants in that ball. Let

\[
\Delta=\min_{d(s,t)\ge\eta}|g(s)-g(t)|>0.
\]

Choose \(0<\delta_0<\pi/2-2\alpha\). Every inscribed rectangle with diagonal angle in \([\theta_*-\delta_0,\theta_*+\delta_0]\) on \(g_n\), for large \(n\), has circumradius at least \(\Delta/4\). Indeed, if the radius were smaller, every pair of vertices would be less than \(\Delta/2\) apart. With \(\|g_n-g\|_\infty<\Delta/4\), all their parameters then lie in the \(\eta\)-ball of the first vertex. Their two diagonals are secants in that ball, so the acute angle between their lines is at most \(2\alpha\), whereas it is at least \(\pi/2-\delta_0>2\alpha\) by the angle band.

Consequently the grading-two spectral derivative is at least \(\Delta^2/16\) almost everywhere on the fixed band, by Theorem 5.1. Integration proves the central slope premise with \(c=\Delta^2/16\) for every smaller window. This legitimately applies the geometric estimate to the spectral carrier since the estimate covers **all** geometric rectangles on the band. It is unnecessary to recognize a separately computed index class.

For example, every family converging in \(C^1\) to a regular \(C^1\) Jordan embedding satisfies the cone hypothesis. On sufficiently short parameter intervals, all derivatives of the limiting curve make an arbitrarily small oriented angle with the tangent at the interval center and have speed bounded away from zero. Uniform convergence of the derivatives preserves that cone. It is convex as an oriented cone, so integrating the derivatives in increasing parameter order gives the required secant statement. This yields a quantitative spectral check in a class with existing affirmative theorems, rather than an extension to all rough curves. Merely uniform convergence of the embeddings does not supply the cone.

## Countermodels and the missing comparison

The v270 model is \(L_0(\theta)=\tfrac A\pi(\theta+\tfrac12\sin(2\theta))\). For \(0<e_n\to0\), set \(L_n=(1-e_n)L_0+e_n A\theta/\pi\). These preserve the endpoints, strict increase, a common Lipschitz bound \(2A/\pi\), and the triangle inequality (subadditivity whenever the angle sum is at most \(\pi\)). For the last property use

\[
\sin(2x+2y)-\sin(2x)-\sin(2y)=-4\sin x\sin y\sin(x+y)\le0.
\]

Every \(L_n'(\theta_*)=e_nA/\pi>0\), and \(L_n\to L_0\) uniformly. Nevertheless the centered quotient at radius \(\delta_n\to0\) is

\[
\frac A\pi\left(1-(1-e_n)\frac{\sin(2\delta_n)}{2\delta_n}\right)\longrightarrow0.
\]

Thus individual positive square-angle slopes, even with a convergent inhomogeneous family of these analytic axioms, do not give a uniform lower bound. For \(L_0\) itself the sets where its derivative exceeds any fixed positive threshold avoid a neighborhood of \(\theta_*\), even though small thresholds yield sets of positive measure. These are abstract functions, with no claim of spectral realization.

- Total angle support or total enclosed area does not locate uniformly large spectral carriers near the square angle.
- Uniform convergence of curves and period integrals does not control the paired local cap defects. The v276 diamond spirals give an exact obstruction to the stronger global defect modulus.
- A geometric parity class and a spectral action class require a comparison that preserves the identity of the relevant carrier. A large rectangle elsewhere in the atlas cannot supply that comparison.

The remaining premise for an arbitrary Jordan limit is a uniform size or localized slope bound on an appropriate carrier, or a verified off-diagonal action condition tied to that carrier. [v276's cap criterion](v276-liouville-sheaf-rough-frontier.md) is one conditional formulation of the latter. These two routes are alternatives and no equivalence between their classes has been proved here.
