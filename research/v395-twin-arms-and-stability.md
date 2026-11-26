# v395. Twin arms, and stability under slow angular deviations

This is the first step of the reduction from model wild points to real ones. Two things are shown:
- simplicity of the curve alone forces a radially monotone wild point to be a two-arm spiral with twist below contact;
- the existence theorem of v393–v394 survives every slowly varying angular deviation from the power-law model, however large the deviation gets.

**Results.** Items 2–4 are computer-assisted, in interval arithmetic with outward rounding, for the power-law profile with $a=1/2$ and level-one, comparable-offset labels as before.

1. **Twin-arm lemma (Lemma 395.1).** Let $p$ lie on a Jordan curve, with arcs $\alpha$ and $\beta$ ending at $p$. Suppose each arc is *radially monotone* near $p$: its distance to $p$ is strictly monotone, so it meets each small circle around $p$ once. Then:
   - their angle functions satisfy $\theta_\beta(r)-\theta_\alpha(r)\in(2\pi k,2\pi(k+1))$ for one fixed $k$ and all small $r$;
   - so if either arm winds infinitely often around $p$, both do, in the same direction and at the same rate up to less than a full turn;
   - choosing the lift, $\theta_\beta=\theta_\alpha+\pi+\varphi(r)$ with $|\varphi(r)|<\pi$.

   **Proof.** The arcs are disjoint away from $p$, so the difference is never a multiple of $2\pi$. It is continuous on the connected interval $(0,r_0)$, so it stays in one gap. ∎

   A radially monotone wild point is therefore *automatically* a two-arm spiral with twist below contact (v394's Lemma 394.1 was the special case of power-law arms). The only freedom left is the common winding profile $\theta_\alpha(r)$ and the regularity of the twist.
2. **Common deviations only count through their slope (Theorem 395.2).** Any radially monotone two-arm curve can be written as
   $$\theta_1(r)=r^{-1/a}+e_1(r),\qquad\theta_2(r)=r^{-1/a}-\pi+e_1(r)+\varphi(r),$$
   with $e_1$ the common deviation from the power law and $\varphi$ the relative twist.
   - v390's anchored label map sees angles only through differences between corners. So $e_1$ enters only through its variation across a square.
   - If $|re_1'|\le L_1$, that variation adds $L_1|\log(|X|/|A|)|$ to each corner's lift deviation, and $2L_1\log\frac{1+u}{1-u}$ to the twist sum in v386's identity.
   - **The size of $e_1$ never matters.** The two arms may rotate together without bound, provided they do it slowly.
   - All of v392's estimates take the pair $(L,L_1)$. v393/v394's degree step involves only the plain spiral and is unchanged.
3. **Thresholds with a common deviation (Theorem 395.3).** Every one of the 60 tabulated triples $(\Phi,L,L_1)$ gets a threshold $c$ such that every curve as in item 2, with $|\varphi|\le\Phi$, $|r\varphi'|\le L$ and $|re_1'|\le L_1$ on the label's annulus, inscribes a square of every label with $m\ge c$ (and $m\ge11$), $M\le2m$, $|j|=1$. The triples run over gaps $\delta=\pi-\Phi$ from 1 to 0.03, $L\le20$ and $L_1\le10$. Infinitely many squares accumulate at the point. Near contact, $c\approx(L+2L_1+\kappa)/\delta$ with $\kappa$ between about 1 and 4.
4. **Example: unbounded common rotation.** Take $e_1=5\log\log(1/r)+\sin(3\log1/r)$, which grows without bound but has slope at most $5/\log(1/r)+3$, and twist $\varphi=2\sin(4\log1/r)$.
   - Every label with $m\ge45$ whose annulus lies inside $r\le0.08$ is covered: infinitely many squares.
   - Plain squares at four labels in three sign patterns, up to offset 310, were followed through the full deformation. All four stay strictly inside their certified regions, keep their degree, and end with their corners on the deformed curve to $10^{-10}$, so they really are inscribed squares.
   - The common rotation on their annuli is 8–10 radians, and still growing at smaller scales.

## Why simplicity does so much

The model of v380–v394 looked like a strong assumption: two arms, the same spiral, offset by half a turn plus a twist. Lemma 395.1 shows most of that is forced. Once each arm meets small circles one point at a time:
- the arms can't wind at different rates, in opposite directions, or one without the other, since any of these makes the relative angle grow until the arms cross;
- what's left is one shared profile plus a twist trapped in a window of width $2\pi$.

The numerics illustrate each case:
- arms winding in opposite directions meet (first at $r\approx0.38$ in the test);
- a winding arm and a straight arm meet;
- arms following $r^{-2}$ and $r^{-1/0.45}$ meet;
- the twin pair from item 4 never meets on $[10^{-4},0.3]$, and its relative twist stays in $(-2,2)$.

## Why the size of a common deviation doesn't matter

A constant common rotation $e_1=\alpha$ just rotates the whole picture, and v390's anchored map is built from differences, $|X|^{-1/a}-\Theta_A-\arg(X/A)-\pi o_X+(E_X-E_A)$. So it is exactly unchanged; a test checks this with $\alpha=123.4$. A non-constant $e_1$ changes by at most $L_1\cdot|\log(|X|/|A|)|$ between two corners. Every estimate already carried such terms for the twist, so $e_1$ just adds to them:
- **Lift deviation.** Each corner's deviation gains $L_1(-\log(1-q))$ for $B$ and $D$, and $L_1(-\log(1-\sqrt2q))$ for $C$.
- **Twist sum.** The signed sum $e_1(A)+e_1(C)-e_1(B)-e_1(D)$ is two differences, so the identity's slack gains $2L_1\log\frac{1+u}{1-u}$.
- **Side bound and floor.** They use the same lift deviations, so they pick up the same terms.

That's all. The certificates are otherwise unchanged, which is why the threshold grows by about $2L_1/\delta$.

## The table

$c(\Phi,L,L_1)$, with gap $\delta=\pi-\Phi$:

| $\delta$, $L$ \ $L_1$ | 0 | 1 | 2 | 5 | 10 |
|---|---|---|---|---|---|
| 1.0, 0 | 20 | 20 | 20 | 23 | 34 |
| 1.0, 5 | 20 | 21 | 22 | 29 | 41 |
| 1.0, 20 | 34 | 37 | 40 | 45 | 55 |
| 0.3, 0 | 30 | 31 | 32 | 51 | 85 |
| 0.3, 5 | 32 | 37 | 46 | 65 | 105 |
| 0.3, 20 | 85 | 93 | 101 | 126 | 148 |
| 0.1, 0 | 59 | 59 | 59 | 114 | 224 |
| 0.1, 5 | 66 | 91 | 103 | 169 | 265 |
| 0.1, 20 | 225 | 234 | 255 | 317 | 421 |
| 0.03, 0 | 143 | 144 | 152 | 354 | 704 |
| 0.03, 5 | 179 | 249 | 319 | 529 | 879 |
| 0.03, 20 | 705 | 775 | 845 | 1055 | 1405 |

The column $L_1=0$ reproduces v394's table. All 60 entries satisfy the degree conditions ($s_{\rm lo}\ge\tfrac14$, $B^*$ inside the region, $m\ge11$).

## What this leaves of the reduction

The reduction from v394's discussion needs a wild point to look like the model at infinitely many scales. After v395, the remaining gaps are:

1. **Radial monotonicity.** Arms that swing back out in radius aren't covered, and the twin lemma doesn't apply to them.
2. **The winding rate.** $\theta_1(r)=r^{-2}+e_1(r)$ with $|re_1'|$ bounded means the arms wind like $r^{-2}$ up to relative accuracy $O(r^2)$ in the rate. A spiral with another exponent, $r^{-1/a'}$, has $re_1'$ unbounded, so it isn't covered. That's step 3: certificates with the exponent $a$ as an interval parameter.
3. **Twist with unbounded slope**, and arms approaching contact faster than about $r^{2/3}$ (v394's heuristic).
4. **Non-winding and rough-everywhere behaviour.** These are outside this line entirely.

**In short.** Simplicity forces radially monotone wild points to be twin-arm spirals with twist below contact. Any such spiral whose shared winding follows $r^{-2}$ up to a slowly varying angle, and whose twist stays a fixed angle short of contact at a bounded rate, inscribes infinitely many squares accumulating at the point.

[Code](../code/twin_arms_v395.py) · [Tests](../tests/test_twin_arms_v395.py) · [Data](../results/v395/twin_check.json)
