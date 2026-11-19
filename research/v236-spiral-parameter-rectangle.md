# v236 — A uniform square-size bound with variable roughness amplitude

Take the v235 base-32 Takagi series \(T\) and its cutoff \(\chi\), and use \(u=\sqrt r\). The first spiral arm has phase

\[
\theta_{\lambda,a}(u)=u^{-4}+\lambda u^{-2}+a\,\chi(u)T(u),
\qquad (\lambda,a)\in[0,1]\times[-1/100,1/100].
\]

For **every pair** in this closed rectangle, the first arm contains an inscribed square of side strictly greater than \(1/2\). The certificate also allows a further continuous Cartesian perturbation of \(10^{-9}\) per coordinate on the four witness arcs. The cutoff leaves the original smooth spiral near the origin, so the divergent Liouville primitive still applies throughout the family. When \(a\ne0\), the phase is nowhere differentiable on the witness arcs. The value \(a=0\) gives the smooth witness arcs of v234.

## How to verify the rectangle

The compressed JSON proof table has 30,995 leaves and encodes a binary partition of the parameter rectangle. An internal node `[axis,left,right]` splits \(\lambda\) (axis 0) or \(a\) (axis 1) at its exact rational midpoint. Both children are mandatory and cover the split boundary in both closed subrectangles. Each leaf gives an orientation and rational data for four moving source centers and a common box width. Replay starts from the whole rectangle and verifies every leaf. It can reject a missing tile without comparing floating-point endpoints.

On a tile, use the source centers \(p_j(\lambda,a)=c_j+v_j(\lambda-\lambda_0)+w_j(a-a_0)\). For each source interval, v235's chord and tail bounds give an affine reference \(L_j(u)\) for the unit Takagi series and an exact upper bound on its phase error. The verifier scales that bound by the largest absolute amplitude on the tile and by the maximum physical radius on the source interval. It adds the allowed Cartesian error separately. The reference phase \(u^{-4}+\lambda u^{-2}+a L_j(u)\) is smooth on that interval.

To bound the preconditioned square residual at the moving centers, apply Taylor's theorem along the straight ray from the tile midpoint to the chosen parameter pair. With ray direction \((d_\lambda,d_a)\), let \(d_u=v_jd_\lambda+w_jd_a\). On a reference arc:

\[
\dot\theta=\theta_u d_u+u^{-2}d_\lambda+L_j(u)d_a,\qquad
\ddot\theta=(20u^{-6}+6\lambda u^{-4})d_u^2+
2(-2u^{-3}d_\lambda+L_j'(u)d_a)d_u.
\]

Differentiate \(u^2e^{i\theta}\) as in v234 using these ray derivatives. The interval calculation combines the midpoint linear terms before multiplying by their parameter displacements, then bounds the second derivative for the Taylor remainder. The source Jacobian bounds and the Cartesian error enter the same fixed-point inclusion used in v234 and v235. A zero for every parameter pair follows from Brouwer's theorem. The squared side is bounded on the whole corresponding source box with the curve error included.

Candidate generation starts from v234's local squares. It uses floating-point refinement only to propose the dyadic center, two speeds, and width. The proof stage checks each proposed tile using exact rationals and outward-rounded SageMath intervals.

## Reproduce

With SageMath's environment active, run `python code/spiral_parameter_rectangle_v236.py`. Use `--generate --manifest PATH` to regenerate the partition, with `--jobs N` to search independent coefficient slabs in parallel. The result concerns this parameter rectangle and this specified phase law.
