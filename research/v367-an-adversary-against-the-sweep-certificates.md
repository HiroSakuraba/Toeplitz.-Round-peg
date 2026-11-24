# v367. An adversary against the sweep certificates

Can successively placed folds stop every [v366](v366-sweeps-read-from-the-curve.md) certificate from reaching two stretches a fixed gap apart along the original loop? This note runs an adversary on the [v247 decagon](v247-true-big-peg-ratio-improvement.md) that bumps the boundary exactly where a certificate is weakest. Inserted vertices inherit registration times, so the stretches keep their meaning on the original loop. The runs are finite, so they can expose a failure mechanism but cannot prove that none exists.

**Results.**
1. **Attacks move parity sideways.** A bump that carries the boundary across the fourth corner at a locus end kills that certificate. It makes the neighbouring triple across the attacked face odd instead. The certificate on the union of the two neighbours keeps degree one and its clearance, because its outer faces were never touched.
2. **Attacking the best certificate** for seven stages alternates the parity between two neighbouring triples. The shared face's clearance falls by a factor of about 0.45 per stage, from 0.095 to 0.00076. The union certificate keeps clearance 0.0772 throughout.
3. **Attacking the most robust certificate** for eleven stages never empties the set. The number of certificates rises from 20 to as many as 78. The best side floor rises from 1.40 to 2.57, because the bumps complete near-squares of the decagon. The best margin falls from 0.095 to 0.017.
4. **Greedy pushing** moves the parity through five triples at a total bump size of 0.245. It then stops at its parity check: the remaining separated move, across the union box's outer face, does not transfer parity. That one bump, about 0.42 in registered distance, destroys the union certificate, which until then had kept its clearance of 0.211. 64 other certificates survive it.

At every stage of every run, some certificate had its stretches separated by at least one empty bin, and its side floor stayed at or above 1.4006. What the adversary can drive to zero is the clearance of the faces it attacks.

## Why attacks only move parity

Take two stretch triples that differ only in one stretch, where the two versions of that stretch are adjacent bins. Their boxes share a face, and together they form the box of the union triple.

**Lemma 367.1.** If the boundaries of both boxes and of their union are zero-free, the union's degree is the sum of the two degrees modulo two. A redraw that changes $F$ only near the shared face, keeping the union's outer faces zero-free, cannot change the union's degree. So it can only move odd parity from one half to the other.

The mod-two degree is additive over boxes glued along a zero-free face. The union's degree depends only on $F$ near its own boundary, by homotopy. ∎

An attacker can move parity only across a face where the locus actually ends, and each move needs a bump at least the clearance of the fourth corner there. To collapse a certificate, the attacker must keep pushing until the union of all visited triples brings two stretches together. Until then, that union is a certificate whose outer clearance has not changed.

## The runs

The bump is a triangle on the edge nearest the attacked fourth corner $D$. Its apex lies beyond $D$, at up to twice the clearance. The checker tries widths from three times the clearance downward and keeps the first simple polygon that puts $D$ on the other side. All search settings are those of v366: 16 bins, spans one and two, one empty bin between stretches, and both orientations.

**Attack the best.** At stage 0 the best certificate is $(7,3,10)$ with clockwise completion. Its weakest end is on the face between Q-bins 2 and 3, with clearance 0.0954.

| stage | odd triple | clearance at the shared face | union $(7,\{2,3\},10)$ |
|---|---|---|---|
| 0 | $(7,3,10)$ | $9.54\times10^{-2}$ | odd, $7.72\times10^{-2}$ |
| 1 | $(7,2,10)$ | $4.27\times10^{-2}$ | odd, $7.72\times10^{-2}$ |
| 2 | $(7,3,10)$ | $1.91\times10^{-2}$ | odd, $7.72\times10^{-2}$ |
| 4 | $(7,3,10)$ | $3.82\times10^{-3}$ | odd, $7.72\times10^{-2}$ |
| 6 | $(7,3,10)$ | $7.63\times10^{-4}$ | odd, $7.72\times10^{-2}$ |

Each bump adds two squares inside the union box, which has 13 by stage 6. The bump sizes shrink geometrically, so the polygons converge. In the limit the attacked fourth corner lands on the curve. That would put a zero on the shared face and block both single-bin triples, but not the union, whose outer faces keep their clearance.

**Attack the most robust.** This policy chooses the certificate whose smallest end clearance is largest, trying to drive the best margin over all certificates to zero.

| stage | vertices | certificates | best side floor | best margin |
|---|---|---|---|---|
| 0 | 10 | 20 | 1.4006 | $9.5\times10^{-2}$ |
| 3 | 19 | 70 | 2.5606 | $3.8\times10^{-2}$ |
| 6 | 28 | 74 | 2.5690 | $1.8\times10^{-2}$ |
| 10 | 40 | 66 | 2.5690 | $1.7\times10^{-2}$ |

The large floors are real. After four bumps the polygon has an exact square of side about 2.599, with all four corners on the boundary. It comes from a near-square of the original decagon, whose fourth corner was about 0.04 off the boundary, and the adversary's own bump completed it. The margin decreases by roughly 0.83 per stage on average, but not monotonically, because new certificates appear with fresh clearance.

**Greedy push.** Starting from $(7,3,10)$, the greedy adversary attacks the cheapest locus end whose face leads to an unvisited single-bin triple that keeps the gaps. The path is
$$
(7,3,10)\to(7,2,10)\to(6,2,10)\to(6,2,11)\to(6,1,11)\to(6,1,12),
$$
with costs $0.095$, $0.077$, $0.025$, $0.0025$ and $0.045$, totalling $0.245$. The union triple, bins $\{6,7\}$, $\{1,2,3\}$, $\{10,11,12\}$, has degree one before and after. Its clearance is exactly the same rational number, about 0.211, and the number of squares inside it grew from 1 to 11.

The run stops because of parity, not separation. From $(6,1,12)$ the locus still has an end on the face toward $(5,1,12)$, and that triple keeps the required gaps. That end lies on the union box's outer face, with clearance 0.211. The bump across it is simple, and its registered distance is about 0.42. After it, $(5,1,12)$, $(6,1,12)$, the union triple, and the union widened to bin 5 all have degree zero. 64 other certificates remain, with best side floor about 2.29. An attack on an outer face of the union box does not move its parity; it destroys it. [v368](v368-certificate-margin-and-one-stage-criterion.md) shows that any redraw that changes the union's degree must exceed its margin, about 0.025 here.

## What this says about the noncollapse question

In this family, folds placed against the certificates do not make them retreat into smaller stretches. Attacks on inner faces only push parity across a fixed bin grid, and parity is conserved on unions of stretches. The danger is different. Each attacked face's clearance can be driven geometrically to zero, and a bump at an outer face of a union box, comparable to its clearance there, can destroy that family's parity outright. In every run, other separated certificates survived.

A proof along these lines would need a lower bound, valid at every depth, on the clearance at the outer faces of some union box whose stretches stay separated. The circle shows that such a bound can fail for a particular target. On near-circles the clearance is proportional to the distance from quarter-turn symmetry (v366). The natural next test is an adversary that keeps pushing across outer faces with a fixed total bump budget. Bounding how many bins it can cross would show whether separation survives.

[Code](../code/refinement_stress_v367.py) · [Tests](../tests/test_refinement_stress_v367.py) · [Results](../results/v367/)
