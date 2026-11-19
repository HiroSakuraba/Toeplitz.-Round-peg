# v233 — Square certificates on mixed inverse-power spirals

Each of the following double spirals contains a square on its first arm. The physical radius is \(r\). The bounds are strict.

| Phase | Square side length | Radius of every vertex |
| --- | ---: | ---: |
| \(r^{-2}+r^{-1}\) | \(1/2\) | \(1/4\) |
| \(r^{-5/2}+\tfrac12r^{-2}\) | \(1/4\) | \(9/25\) |
| \(r^{-3}+\tfrac32r^{-3/2}\) | \(2/5\) | \(9/20\) |

Each phase differs without bound near zero from every single inverse power plus a constant. Even if the leading power and coefficient match, the second term still diverges. The witness intervals were discovered on the mixed phases themselves.

## Curves and irregularity

Write \(r=u^2\) and \(\theta(u)=\sum_n c_n u^{-n}\) for positive rational \(c_n\) and distinct positive integers \(n\). Set \(g(u)=u^2 e^{i\theta(u)}\) for \(0<u\leq1\), and \(g(0)=0\). Each radius occurs once per arm. The arms \(g\) and \(-g\) meet only at the origin, and either unit semicircle between their outer endpoints closes a Jordan curve.

On compact smooth subarcs, a continuous Liouville lift must agree with integration of the Liouville form up to a constant. In the convention \(\lambda=(x\,dy-y\,dx)/2\), its pullback is \(u^4\theta'(u)\,du/2=-\tfrac12\sum_n n c_n u^{3-n}\,du\). The primitive diverges at zero when any \(n\ge4\) occurs, as in all three examples. The terms have the same sign. Changing the primitive of the area form by an exact smooth differential adds only a bounded endpoint term. See [Asano–Ike, Lemma 6.2](https://arxiv.org/html/2609.11598v3#S6).

## Exact certificate

On a proposed rational interval \([a,b]\subset(0,1]\), each endpoint phase is rational, so the v231 Machin and Taylor calculation gives rational sine and cosine bounds. Grid rounding adds an explicit coordinate error.

\[
g''(u)=e^{i\theta(u)}[2-u^2\theta'(u)^2+i(4u\theta'(u)+u^2\theta''(u))].
\]

Here \(-u\theta'(u)=\sum_n n c_n u^{-n}\) is decreasing, and the imaginary term is \(\sum_n n(n-3)c_n u^{-n}\). Hence on \([a,b]\)

\[
\sup |g''|\le2+[a\theta'(a)]^2+\sum_n |n(n-3)|c_n a^{-n}.
\]

Multiplying by \((b-a)^2/8\) bounds the coordinatewise interpolation error. Add the larger endpoint error and take the maximum \(E\) over all four intervals. Let \(A\) be the v231 affine chord residual matrix and \(m\) its interior root's source margin. The strict test \(4E<m/\|A^{-1}\|_\infty\) keeps the straight-line residual homotopy away from zero on the boundary. Nonzero Brouwer degree gives a zero on the true arcs.

At that zero, the parameter of each vertex can differ from its value at the affine zero. The checker includes this displacement by bounding the distance from the affine vertex to **every point on its chord**, then adds the curve-to-chord error. For the first two arcs, subtract the sum of these two allowances from the affine side length to get a uniform side bound. Each vertex radius is at least its interval's \(a^2\). The positive side bound makes the four vertices distinct.

## Reproduction

Run the exact replay with `python code/mixed_power_spiral_v233.py`. Add `--rediscover` to run angle-spaced numerical sampling and bounded candidate search before exact verification. The stored source endpoints are integer numerators over \(2^{30}\) in \(u\). Run the focused regressions with `python -m unittest tests.test_mixed_power_spiral_v233 -v`.

With the SageMath Conda environment active, run `python code/verify_mixed_power_spiral_sage_v233.py` for a second calculation. It uses 160-bit interval trig and curvature bounds and solves the affine system over SageMath's rational field.

The proved instances are the three specified phase laws. The degree criterion for other laws remains conditional on finding a box that passes its inequalities. The unrestricted square-peg conjecture remains open.
