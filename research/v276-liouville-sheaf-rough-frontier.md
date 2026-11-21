# v276 — Local action defects and exact counterwinding controls

**Status.** Restricted primitive estimates and exact positive and negative controls, with a conditional Floer cap criterion. The filtered comparison with the geometric relative class is open.

## Coverage and imported results

Here zero area means zero planar measure of the *image*. The bounded complementary domain still has positive area. Use \(\alpha=(x\,dy-y\,dx)/2\). Since \(y\,dx=d(xy)/2-\alpha\), changing between the two conventions preserves locally uniform convergence of normalized primitives under uniform convergence of parametrizations.

| Hypothesis on the Jordan curve | Result used | Primary reference |
|---|---|---|
| Smooth embedding | Rectangles in every prescribed similarity class | [Greene–Lobb, Annals 194 (2021)](https://doi.org/10.4007/annals.2021.194.2.4) |
| Continuous Legendrian lift, via smoothly approximated primitives converging locally uniformly on the parameter lift | Prescribed-angle rectangles | [Asano–Ike v4, Definition 1.1 and Theorem 1.3](https://arxiv.org/html/2412.21057v4#S1), preprint |
| Rectifiable or locally monotone | Continuous Legendrian lift | [Asano–Ike v4, Propositions 5.8 and 5.11](https://arxiv.org/html/2412.21057v4#S5), preprint |
| Finite \(p\)-variation for \(p<2\) | Smooth approximation in \(q\)-variation for \(q>p\), with Young primitive convergence when \(q<2\) | [Li–Pan, Theorems 1.1–1.2 and proof in §4](https://arxiv.org/html/2609.04918#S4), preprint |
| Asymmetric variation and coordinate bounded variation on finitely many pieces, with their stated area-preserving bi-Lipschitz variants | Liouville regularity | [Asano–Ike v3, Theorem 1.6](https://arxiv.org/html/2609.11598v3#S1), preprint |
| Uniformly vanishing conformal radial-tail lengths, including the stated Hölder, bounded Jordan John and quasicircle cases | Liouville regularity | [Asano–Ike v3, Theorem 1.7 and §4](https://arxiv.org/html/2609.11598v3#S4), preprint |
| Positive area of the image | Prescribed-angle rectangles via density | [Asano–Ike v4, Remark 5.5](https://arxiv.org/html/2412.21057v4#S5), preprint |

## Rational diamond spirals

On the \(k\)-th block, \(k\ge0\), the \(\ell^1\) radius decreases linearly between successive vertices from \(R_k=2^{-k}\) to \(R_k/2\), in \(N_k=4\cdot4^k\) steps. Every vertex is on one of the four coordinate rays. Its successor is on the next ray in the chosen quarter-turn direction. There are two choices for this direction:

- Monotone winding uses a positive quarter turn at every step.
- Counterwinding reverses sign every \(B_k=4\cdot2^k\) steps within the block, beginning with positive turns.

Both block and burst lengths are multiples of four. Every block begins and ends on the positive \(x\)-axis. Each open segment lies in a quadrant where its \(\ell^1\) radius is affine and strictly decreasing. Parametrize the resulting arm by that radius and extend it continuously to zero. The opposite arm is its central reflection. Equality of points on the two arms at positive radius would give a point equal to its negative, because each arm meets a given radius once. The two arms therefore meet only at zero. On the outer diamond of radius one, append the two lower edges between the outer endpoints. This gives a Jordan curve of planar measure zero, as a countable union of segments and a point.

Truncation after finitely many complete blocks replaces the two remaining tails by the radial diameter through their last endpoints. The approximating polygons are Jordan by the same reasoning and converge uniformly in compatible parametrizations, since the changed part of each curve lies in a shrinking diamond. Traverse inward on one arm and outward on its reflection. The two arm integrals cancel in every truncation by central reflection invariance of \(\alpha\), the connecting diameter contributes zero, and the two lower outer edges contribute exactly \(1\). Thus the enclosed area is \(1\) for every truncation. The limiting enclosed area is also \(1\): away from the zero-measure limit curve the winding numbers eventually agree by uniform convergence, and all bounded domains lie in a fixed compact set.

The absolute integral on one quarter-turn segment with radii \(r_j,r_{j+1}\) is \(r_j r_{j+1}/2\). It lies between \(R_k^2/8\) and \(R_k^2/2\). In particular:

1. **Obstruction.** The monotone winding increases the primitive by at least \(N_kR_k^2/8=1/2\) in every block, although the whole block lies in a diamond of radius \(R_k\). Any continuous primitive would agree with the segment integrals on every open segment by [Asano–Ike v3, Lemma 6.2](https://arxiv.org/html/2609.11598v3#S6), and at the limiting endpoints of each segment by continuity. Summing the first \(k\) blocks contradicts continuity at zero. This is a closed, zero-measure, non-Liouville-regular Jordan example with exactly area-preserving PL approximations. The same failure of primitive control holds for an area-preserving smooth approximation sequence: for each finite truncation round corners with uniform and primitive error smaller than \(1/n\), then rescale by a factor tending to \(1\). Smoothing can be taken arbitrarily fine on each fixed polygon, so the primitive error after rescaling can also be made smaller than \(1/n\) despite unbounded lengths over the sequence.
2. **Positive control.** Each one-sign burst contributes at most \(B_kR_k^2/2=2^{1-k}\), and the unsigned burst sums decrease because the radii decrease. Alternating sums of this decreasing sequence stay in \([0,2^{1-k}]\). Primitive values within a burst lie between the values at its ends, so the same bound holds throughout the block. The sum over all remaining blocks tends to zero. Therefore the primitives of the closed truncations converge uniformly on one period. Rounding each fixed polygon in disjoint corner disks preserves the embedding and changes the primitive by arbitrarily little: each replacement is made by trimming both incident segments and inserting a short tangent-matching arc in their triangle, and the sum of the affected lengths tends to zero for the fixed polygon. Choose that smoothing error to vanish along the sequence. The one-period convergence and convergence of period integrals give locally uniform convergence on the lifted parameter line. Thus this counterwinding curve admits a continuous Legendrian lift, and Asano–Ike's existing theorem supplies the rectangles.

Both curves have infinite \(p\)-variation for every \(1\le p\le2\): the sum of \(p\)-th powers of the edge lengths in block \(k\) is at least \(N_k(R_k/2)^p\ge1\). Both wind through all four coordinate directions in every block, so every nonzero projection takes both positive and negative values arbitrarily close to the endpoint of the arm with projection zero. This precludes local monotonicity at zero. They also fail the bounded-turning condition: two consecutive returns to the positive ray in a single burst have separation \(4(R_k/2)/N_k\), but the diameter of the intervening arc is at least \(R_k/2\), while the complementary arc contains the outer join. The ratio tends to infinity, ruling out a quasicircle by the [arc condition in Ahlfors, Theorem 1](https://archive.ymsc.tsinghua.edu.cn/pacm_download/117/5956-11511_2006_Article_BF02391816.pdf). For the positive example we do not use a conformal criterion or classify its domain as John or Hölder. For the negative example the conformal sufficient hypotheses and other sufficient lift hypotheses in the table are excluded by nonregularity.

Run the finite rational audit with:

~~~text
python code/liouville_primitive_experiments_v276.py
python -m unittest discover -s tests -p test_liouville_primitive_experiments_v276.py
~~~

## A conditional cap criterion that only tests the selected rectangles

For a smooth positively oriented embedding \(g_n:\mathbb R/\mathbb Z\to\mathbb C\), and an oriented parameter arc lifted to \([s,t]\), define the chord-closed signed defect

\[
D_n(s,t)=\int_s^t g_n^*\alpha-\tfrac12 g_n(s)\wedge g_n(t).
\]

The wedge is the planar determinant. This is invariant under translation of the curve: the translated line integral and chord term differ by the same endpoint contribution. A short graph in a thin strip can make this signed defect small even when its length is large.

Suppose \(g_n\to g\) uniformly in compatible parametrizations, with \(g\) a Jordan embedding, and that prescribed-angle nondegenerate rectangles on \(g_n\) carry actual normalized Greene–Lobb actions \(\mathcal A_n\) with \(\operatorname{dist}(\mathcal A_n,A_n\mathbb Z)\ge\varepsilon>0\), where \(A_n\) is the enclosed area of \(g_n\). On each subsequence with total collapse, require that the sum of the two chord-closed defects of the return connectors tends to zero. Here the connectors run from the terminal diagonal of the Hamiltonian trajectory back to the initial diagonal, with both parameters taken in the shrinking common lifted arc containing all four vertices. Both terms are signed. A sufficient, stronger check is a uniform modulus \(D_n(s,t)\to0\) on all short arcs in the possible collapse region.

This prevents shrinkout. For a collapsing subsequence the four parameters cluster at a single parameter after extraction since the limit is an embedding. Their lifted common arc lies in a disk of radius \(h_n\to0\) by uniform convergence and continuity of \(g\). In [Greene–Lobb, §5, proof of Lemma 5.1](https://arxiv.org/html/2404.05179#S5), the cap is represented by the two rotating coordinate arcs and a return path on \(g_n\times g_n\) disjoint from the diagonal, with a core-winding correction. Modulo \(A_n\) the winding correction has zero contribution. The two projected return paths have the same integrals as the chosen short connectors modulo \(A_n\), because their parameter windings differ by full periods. Center the primitive at the disk center. Each connector integral differs from its defect by at most \(h_n^2/2\). Both rotation arcs lie in the disk of radius \(2h_n\), each has length at most \(\pi h_n\), and the Hamiltonian term is at most \(\pi h_n^2\). Hence

\[
\operatorname{dist}(\mathcal A_n,A_n\mathbb Z)
\le |D_n(s_1,t_1)+D_n(s_2,t_2)|+(3\pi+1)h_n^2.
\]

This contradicts the gap. The argument also permits angles in \((0,\pi)\) converging to a fixed interior angle: the estimates are uniform, and any positive-radius limit has that angle. It is a criterion on the very same selected action carriers, with no implication that all rectangles stay large.

There is a sharp scope distinction. A uniform vanishing modulus for \(D_n\) on *all* short parameter arcs gives equicontinuity of the normalized primitives, since the endpoint term has a uniform modulus by compatible uniform convergence. A finite subdivision bounds the primitives uniformly on one period, so Arzelà–Ascoli gives a convergent subsequence and period convergence, hence a continuous Legendrian lift. Conversely, uniform primitive convergence implies such a global defect modulus. Thus the global condition is a repackaging of the existing regularity mechanism, not an extension. Only the localized, paired check can potentially be weaker. The monotone diamond spiral shows why geometric closeness and fixed area alone do not give the stronger modulus: at each small block, the endpoints are collinear with the origin and its defect is at least \(1/2\).

## Sheaf gate

[Asano–Ike v4, Theorem 4.1 and Remark 4.2](https://arxiv.org/html/2412.21057v4#S4) locate the change at a persistence value \(a_0\) through the corresponding microlocal Hom. At a fixed angle they give a conditional route via vanishing of its derived sections over the diagonal preimage, at the *same* critical value. The numerical cap criterion above does not prove that vanishing. In particular it neither identifies a geometric parity class with their generator nor supplies a persistence class or action gap on a new rough family. A comparison would need to identify the precise nontrivial class, its selected off-diagonal carriers and the action normalization, and then verify the paired defect for those carriers. This obligation remains open.
