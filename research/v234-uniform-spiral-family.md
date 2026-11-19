# v234 — A uniform square-size bound on a spiral family

For \(0\le\lambda\le1\), set \(\theta_\lambda(r)=r^{-2}+\lambda r^{-1}\). The arm \(r e^{i\theta_\lambda(r)}\), its negative, and the unit semicircle between their outer endpoints form a Jordan curve. Each arm meets each radius once, the two arms meet only at the origin, and the semicircle meets the arms only at its endpoints.

**Result.** For every \(\lambda\in[0,1]\) the first arm contains a square of side strictly greater than \(1/2\). The same bound holds if that arm is continuously perturbed by at most \(10^{-9}\) in each Cartesian coordinate on the four source intervals of a certificate for \(\lambda\). In particular, an arbitrary continuous additive phase perturbation bounded there by \(10^{-9}\) radians satisfies the coordinate bound and preserves the Jordan construction.

This covers a continuum of curves outside Liouville regularity. On the first arm the pullback of \((x\,dy-y\,dx)/2\) is \((-r^{-1}-\lambda/2)\,dr\). Its primitive diverges logarithmically at the origin for every \(\lambda\). The same obstruction still holds for perturbations supported away from the origin. See [Asano–Ike, Lemma 6.2](https://arxiv.org/html/2609.11598v3#S6) for the agreement of a Liouville primitive with integration on smooth subarcs.

## Verification

The certificate table partitions \([0,1]\) into 512 closed rational intervals. On each interval write \(m\) for its midpoint, \(h\) for its half-width, and use \(u=\sqrt{r}\) on the spiral. Four proposed source parameters have the form \(p(\lambda)=c+v(\lambda-m)\), with four independent deviations \(z_i\in[-w_i,w_i]\). The rows contain the interval endpoints, the square orientation, and the dyadic numerators of \(c,v,w\) over \(2^{40}\). The verifier checks that the four source ranges in \(u\) stay disjoint and inside \((0,1)\) throughout the interval.

Let \(F\) be the four-coordinate square residual from v231, \(A\) the source Jacobian at the midpoint, and \(B\) a rational approximation to \(A^{-1}\). The verifier checks \(B\) is invertible over the rationals. For each fixed \(\lambda\), use

\[
T_\lambda(z)=z-B F(p(\lambda)+z,\lambda).
\]

At \(z=0\), enclose \(B F(p(\lambda),\lambda)\) by Taylor's theorem in \(\lambda\), with a midpoint value, midpoint derivative, and a second derivative bound along the moving center. The identity underlying the latter bound is

\[
\frac{d^2}{d\lambda^2}(u^2 e^{i\theta})=
e^{i\theta}\!\left[2v^2-u^2\dot\theta^2+
i(4uv\dot\theta+u^2\ddot\theta)\right],
\]

where \(u=c+v(\lambda-m)\), \(\dot\theta=u^{-2}-v(4u^{-5}+2\lambda u^{-3})\), and
\(\ddot\theta=v^2(20u^{-6}+6\lambda u^{-4})-4vu^{-3}\).

In the \(z\) directions, the mean value bound uses interval entries of \(I-BD_u F\) over the whole moving box. If \(R\) is the linear map from the eight vertex coordinates to the four residuals, coordinate perturbations of size \(\varepsilon\) add at most \(\varepsilon\sum_k |(BR)_{ik}|\) to row \(i\). This last bound uses continuity and the coordinate error alone. SageMath checks that the sum of the three bounds for each row is strictly less than \(w_i\). Thus the perturbed \(T_\lambda\) maps the box into its interior. Brouwer's fixed-point theorem gives a zero, and interval bounds on the first two vertices, with the perturbation included, certify a side strictly greater than \(1/2\).

Numerical continuation only supplies the proposed boxes. Replay checks every interval and the complete cover. It does not need to follow a particular square through a shared interval boundary.

## Reproduce the result

With the SageMath Conda environment active, run `python code/spiral_family_v234.py`. The optional `--generate --manifest PATH` mode reruns candidate generation and verifies a fresh table before writing it.

This theorem controls one specified family and perturbations on its certified source arcs. A useful next challenge is to find comparable boxes on rough curves whose geometry differs by more than the certified tolerance on the witness arcs.
