# v235 — Squares on a spiral family with a nowhere differentiable phase

Use \(u=\sqrt r\) and define
\[
T(u)=\sum_{k=0}^\infty 32^{-k}\operatorname{dist}(32^k u,\mathbb Z),
\qquad
\theta_\lambda(u)=u^{-4}+\lambda u^{-2}+\frac{\chi(u)}{100}T(u),
\qquad 0\le\lambda\le1,
\]
where \(\chi(u)=\min(1,\max(0,12(u-1/4)))\). The first arm is \(u^2e^{i\theta_\lambda(u)}\) and extends continuously to the origin. Its negative and the outer unit semicircle close it to a Jordan curve by the same strict radial separation as in v234.

The 938-box certificate gives an inscribed square of side **strictly greater than \(1/2\)** for every \(\lambda\) in the interval. The bound survives an additional continuous Cartesian error of \(10^{-9}\) per coordinate on the certified arcs.

All four arcs of every box lie in \(1/3\le u\le5/6\), where \(\chi=1\). There \(T(u)\ge \operatorname{dist}(u,\mathbb Z)\ge 1/6\). The phase displacement is at most \(32/(6200)<1\) radians, so the Euclidean displacement from the v234 arm is at least \(u^2 T(u)/200\) by \(2\sin(x/2)\ge x/2\) for \(x\in[0,1]\). At least one Cartesian component changes by at least \(1/21600\) at each source point. That exceeds the v234 perturbation budget by more than 46,000 times.

The series is uniformly convergent and nowhere differentiable. On each base-32 cell of level \(n\) the first \(n\) tent terms are affine, and the terms at levels \(n\) and higher vanish at its endpoints. The cell secant slopes at successive levels therefore differ by \(+1\) or \(-1\). If a finite derivative existed at a point, the secant slopes of the nested cells containing it would converge to that derivative. At grid endpoints one can take the nested cells to the right. A smooth local argument of \(u^{-2}g(u)\) for the arm \(g\) would recover the phase up to an additive constant, so the arm has no finite derivative in the radial parameter on the interiors of the certified source intervals. See [Allaart on the Takagi–van der Waerden family](https://arxiv.org/abs/1312.2119).

Near zero \(\chi=0\). The smooth original spiral there retains the logarithmically divergent action primitive from v234, which precludes a continuous Liouville lift. [Asano–Ike, Lemma 6.2](https://arxiv.org/html/2609.11598v3#S6)

## Proof certificate

For each source range \([a,b]\), form the affine chord of the first six terms of \(T/100\). The maximum deviation of each retained tent term from its chord occurs at an endpoint or a tent knot. At the knots the tent term has two alternating values. Within each parity the signed difference from the chord is affine in the knot index, so checking the first and last knot of each parity suffices even when a source interval contains thousands of knots.

Sum the six exact rational chord errors and add the geometric tail bound \(32^{-6}/[200(1-1/32)]\). Multiplication by \(b^2\) bounds the Cartesian deviation between the true rough arm and the arm with the affine reference phase. This uses \(|e^{ix}-e^{iy}|\le |x-y|\).

The v234 fixed-point test runs on the smooth reference for each source box. It includes the four independent Cartesian error bounds when enclosing the residual and the square side. The derivative and Taylor calculations apply only to the reference phases. Continuity of the infinite series is enough to apply Brouwer's theorem to the true residual. Replay verifies every rational parameter interval and rejects gaps in the cover.

## Reproduce

With SageMath's environment active, run `python code/rough_spiral_family_v235.py`. The optional `--generate --manifest PATH` starts from the v234 locator data, refines candidates for the new phase, and writes a table only after it passes replay.
