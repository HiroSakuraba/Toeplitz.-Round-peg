# v371. How much symmetry the half-turn mechanism needs

[v370](v370-half-turn-wild-points-cannot-absorb-parity.md) proved that an exactly half-turn-symmetric wild point cannot absorb the odd square parity. This note asks how much asymmetry the mechanism tolerates. The answer: the large graceful degree is stable under every deformation that does not create a square of side exactly $s$. So the tolerance is set at the scale of the cutoff, not at the wild point. Exact counts on twisted double spirals show this directly.

**Results.**
1. **Stability (Proposition 371.1).** An explicit uniform-distance bound keeps the large aligned and reversed degrees fixed. It depends only on how far the residual is from zero on rows whose side is near $s$. Collapse and tiny squares do not enter.
2. **Path principle (Corollary 371.2).** Along a continuous family of closed curves with no square of side exactly $s$, the large degrees do not change. Intermediate curves need not be simple.
3. **Beyond symmetry (Corollary 371.3).** Suppose a Jordan curve can be joined to an exactly half-turn-symmetric wild curve without passing a square of side $s$. Then it has a square of side greater than $s$. This includes every curve within an explicit uniform distance of a symmetric one, however wild it is anywhere.
4. **Twisted spirals.** Turn arm two by $\pi+\phi$ instead of $\pi$:
   - at $\phi=0.001$, every count is unchanged;
   - at $\phi=0.01$, the small aligned and reversed counts differ by up to 6, but their sum stays even at every depth and both cutoffs, so the large graceful count stays odd;
   - at $\phi=0.1$, cutoff $1/20$, one aligned square is lost at depths 3 through 6. That is a single cutoff event at the outer scale, not at the wild point. At cutoff $1/40$ nothing changes.

## Stability of the large degree

Use the rows and residual $G$ of [v369](v369-relabeling-dichotomy-and-collapse-localization.md) for any registered closed curve $C$, simple or not. Write $\mathrm{side}_C(x)=|p_1-p_0|$. The *large aligned degree* at cutoff $s$ is the mod-two degree of $G_C$ over the rows with $\mathrm{side}_C>s$ in the aligned quotient. The *large reversed degree* uses the reversed quotient.

Both are defined when $C$ has no square of side exactly $s$. Zeros on $\partial X$ have two equal consecutive points, hence all four points equal, hence side $0$. So they never lie in the large region, and this holds for non-simple curves too.

**Proposition 371.1.** Let $C$, $C'$ be registered closed curves with $\sup_t|C(t)-C'(t)|<\delta$. Suppose every row $x$ with $|\mathrm{side}_C(x)-s|\le2\delta$ has $|G_C(x)|>4\sqrt2\,\delta$. Then $C$ and $C'$ have the same large aligned degree and the same large reversed degree at cutoff $s$.

*Proof.* Take the straight homotopy $C_\lambda$. Each point moves by less than $\delta$. So the side changes by less than $2\delta$, and each residual vector by less than $4\delta$, which is $4\sqrt2\,\delta$ in $\mathbb R^4$. A zero of $G_\lambda$ with $\mathrm{side}_\lambda=s$ would need $|\mathrm{side}_C-s|<2\delta$ and $|G_C|<4\sqrt2\,\delta$, which is excluded. So the moving region $\{(x,\lambda):\mathrm{side}_\lambda(x)>s\}$ of $X\times[0,1]$ has no zeros on its boundary. The zero sets at $\lambda=0$ and $\lambda=1$ are cobordant, and the degrees agree. ∎

**Corollary 371.2 (path principle).** Let $C_\lambda$, $\lambda\in[0,1]$, be a uniformly continuous family of registered closed curves, none of which has a square of side exactly $s$. Then all $C_\lambda$ have the same large degrees.

*Proof.* The set of pairs $(\lambda,x)$ with $\mathrm{side}_\lambda(x)=s$ is compact, and $|G|>0$ on it. So there is a uniform gap, and Proposition 371.1 applies to short steps. ∎

**Corollary 371.3.** Let $\gamma_0$ satisfy the hypotheses of Theorem 370.2, and let $s\le s_0(\gamma_0)$ be such that $\gamma_0$ has no square of side exactly $s$. Then:
- the large graceful degree (aligned plus reversed) of $\gamma_0$ is odd;
- every Jordan curve joined to $\gamma_0$ as in Corollary 371.2 inscribes a square of side greater than $s$;
- in particular, so does every Jordan curve within the uniform distance $\delta$ of Proposition 371.1 for $C=\gamma_0$.

*Proof.* The symmetric drawings $D_n$ of v370 converge uniformly to $\gamma_0$, and their large graceful degrees are odd. For large $n$, Proposition 371.1 applies with the positive gap of $\gamma_0$. The rest is Corollary 371.2. ∎

The nearby curve needs no symmetry, no tameness and no special point. It may be wild everywhere within $\delta$. What it must avoid is a square of side exactly $s$ along the way. That is a condition at the scale $s$, away from any collapse.

This is the classical deformation strategy. Two things differ here. The base curves may have wild points, and the large region never sees collapse, so no control of tiny squares is needed along the path. The difficulty of the general problem reappears exactly as it should. A curve with no square would have to be separated from every such base curve, for every small $s$, by curves carrying squares of side exactly $s$.

## Twisted double spirals

Arm one is $r=\theta^{-1/2}$ from $\theta=3$ inward, with 24 rational vertices per turn. Arm two is arm one turned by $\pi+\phi$, by an exact rational rotation. A chord joins the inner ends. For every $\phi\in(-\pi,\pi)$ the two arms stay disjoint and interleaved. All drawings are checked to be Jordan, and no degenerate candidates occurred. By Lemma 370.1, the large graceful parity is $1+(\text{small aligned})+(\text{small reversed})$.

| depth | $\phi=0.001$, side $<1/20$ | $\phi=0.01$ | $\phi=0.1$ | $\phi=0.5$ |
|---|---|---|---|---|
| 3 | 5, 5 | 5, 1 | 4, 9 | 0, 16 |
| 4 | 29, 29 | 29, 23 | 28, 33 | 21, 40 |
| 5 | 53, 53 | 53, 47 | 52, 57 | 45, 64 |
| 6 | 78, 78 | 78, 74 | 77, 82 | 70, 88 |
| 7 | 104, 104 | 104, 100 | 102, 108 | 95, 114 |

Each entry is small aligned, small reversed.

The large graceful count is odd exactly when the two entries have the same parity:
- at $\phi=0.001$ and $\phi=0.01$, at every depth, for both cutoffs $1/20$ and $1/40$;
- at $\phi=0.1$ for cutoff $1/40$, but not for cutoff $1/20$ at depths 3 to 6.

The $\phi=0.1$ defect is a single aligned square present from depth 3 on. It has its corners on the outer turns, where squares of side about $1/20$ live. It is not near the wild point. At $\phi=0.5$ several such events occur. The parity that symmetry protected changes only through squares of side exactly $s$, as Corollary 371.2 predicts.

## What remains

Every result in this strand reduces to whether some path from a protected curve avoids squares of side exactly $s$. For the twisted spirals that can be checked at each depth, and small twists pass. For a general wild curve, the cutoff gap at scale $s$ along any natural path is the quantity to control. That is the same quantity as the v368 margin, now at the cutoff instead of at a box face.

[Code](../code/twist_stability_v371.py) · [Tests](../tests/test_twist_stability_v371.py) · [Scan](../results/v371/twist_scan.json)
