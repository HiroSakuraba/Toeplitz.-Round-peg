# v364. When three patches force the fourth corner onto the target

Three shrinking edit patches can all meet a Jordan target while an odd three-way change in square parity keeps its fourth corner a fixed distance away. The target can even miss the limiting cutoff entirely. There is a positive result if the patches can be moved across the target. Odd parity can be carried from completions inside to completions outside without choosing a continuous branch of squares. Below, an exact polygon replay forces this transfer through a pair creation and a pair loss.

## The sweep lemma at one positive side

Fix a target Jordan boundary $G$, its inside-negative signed distance $\rho$, and a side $s>0$. Use the geometric cutoff wall from [v363](v363-cutoff-parity-under-disjoint-edits.md): a point $Q$ specifies the center and orientation of a physical square of side $s$, with the angle taken modulo a quarter turn. Take common large outer center bounds with no wall zeros on the outer faces throughout.

There are two parameters. The edit time is $t\in[0,1]$. The sweep position is $l\in[0,1]$. At each position take three nonempty disjoint closed patches $K_i(l)$ varying continuously in Hausdorff distance. Their maximum diameter over the sweep is strictly smaller than $s$. For each of the eight subsets $S$ the field along the edit has the form
$$
 f_S(l,t,x)=f_0(l,t,x)+\sum_{i\in S}h_i(l,t,x),
 \qquad \operatorname{supp}_x h_i(l,t,\cdot)\subset K_i(l).
$$
All fields are continuous jointly. They are common at $t=0$ (take all increments zero there). During the edits their zero sets need not be Jordan curves. A wall zero means that all four values of the field at the vertices of $Q$ vanish. Call a zero *marked* here if one vertex lies in each of the three patches. These three vertices are distinct. The patch diameters ensure that the remaining vertex is outside all three. Call it $q_*$.

The following conditions suffice.

1. At the four parameter corners all endpoint fields are cutoff-free. At $l=0$ the XOR of the eight endpoint degrees above $s$ at $t=1$ is one. Every marked wall zero at $l=0$ has $q_*$ strictly inside $G$.
2. At $l=1$ every marked wall zero has $q_*$ strictly outside $G$.
3. Across the sweep, the finished drawings ($t=1$) have no marked wall zeros. Cutoff clearance for all finished drawings is a stronger, often simpler check. The common starting face cancels without a clearance assumption in its interior.

**Conclusion.** At some sweep position, time and mask there is a side-$s$ square with one corner in each patch and its fourth corner **on $G$**. At finite patch width the first three corners are only guaranteed to lie in the stated patches. In particular, if such sweeps exist at $s_n\to s>0$ and throughout each sweep all three patches lie within $\eta_n\to0$ of the same target, then $G$ contains a square of side **exactly $s$**. The rest of the comparison boundary can remain far away. Compactness of the four physical corners proves this last assertion.

### Why crossings cannot all disappear in the sweep

At a marked configuration the distance of $q_*$ from the union of the patches is at least $d=s-\max_{i,l}\operatorname{diam}K_i(l)>0$: in each patch compare with its marked vertex. To distinguish an inside completion without labeling the square, use the continuous scalar
$$
 H(Q,l)=\sum_{q\text{ corner of }Q}
 \min\big(1,\operatorname{dist}(q,\cup_i K_i(l))/d\big)\,\rho(q).
$$
At a marked configuration this equals $\rho(q_*)$. At other configurations its sign is immaterial, provided it is the same for all masks.

The cutoff wall including $t,l$ has dimension five. The four corner evaluations form a section of rank four. After regularization its zeros are compact one-dimensional manifolds. Cut them at $H=0$ and sum the boundary counts of their portions with $H<0$, modulo two over the eight masks. A boundary point whose square misses one patch occurs with and without that edit, at the same parameters. Flip the first missed patch. This cancels that point, including on the cut $H=0$. At $t=0$ all eight copies cancel. Condition 3 removes marked endpoints at $t=1$. At $l=0$ the inside contribution is odd: the ordinary total wall count is the XOR of changes in endpoint degrees, the common starting degrees cancel, and every marked crossing is inside. At $l=1$ there is no marked inside contribution. A compact one-dimensional manifold has even boundary count, so an odd contribution is required on $H=0$. Its surviving square must see all three patches.

For continuous and singular data use the same regularization for every mask. Smoothly approximate the additive fields with each increment supported in a shrinking space-time neighborhood of its patch graph, choosing the same approximants for masks sharing an increment. Keep the common start, or take a collar with vanishing increments there after a common time reparametrization. Smooth the scalar $H$ and use a common small perturbation of the section and scalar, transverse on the bulk and on parameter faces. Finite local bundle frames and bump functions span the evaluation fibers; a simultaneous generic choice works for the finite family, also for the cut section of rank five. The perturbed sections then coincide when toggling an unseen edit. Rank and dimension exclude corner intersections after perturbation; original endpoint clearance preserves the degrees. Compactness and strict endpoint signs exclude wrongly signed marked crossings on the corresponding end faces for small perturbations. The absence of marked finished zeros similarly gives a margin there. Sending the perturbations and support enlargements to zero gives the claimed marked zero on the cut. A boundary hit along one particular branch is never assumed. If the finished face is left uncontrolled, a marked inside crossing there is an additional possible outcome of the same boundary count.

## A second test when the comparison can be corrected at the fourth corner

There is also a static consequence of v363. Keep the three marked edits $I$ fixed. Suppose the background near all feasible fourth-corner positions can be corrected to the target by a finite list $J$ of auxiliary edits. Assume all their compact supports and the three marked supports are pairwise disjoint. Assume the cutoff degrees used below are defined, with all hybrid endpoints cutoff-free for the Boolean expansion. Write $c(B)=D_Ip(B)$ for the cubic difference on background $B\subset J$. Then
$$
 c(J)=c(\varnothing)\ \oplus\ \bigoplus_{j\in J}
    \left(c(\{j\})\oplus c(\varnothing)\right).
$$
Indeed $p$ has degree at most four in all edit bits, so after three finite differences $c$ has degree at most one in the remaining bits. Each correction term is a fourth difference. By v363 it vanishes if a side-$s$ square cannot touch the marked patches and that correction patch at once. Intersecting the computed fourth-corner enclosures with the correction patches gives a quick sufficient exclusion.

If the reconstructed $c(J)$ is odd, apply v363 on the corrected background while varying only the three marked edits. The missing corner lies outside the marked patches, so it is a zero of the **fixed corrected background**. Agreement of that background's zero set with the target on the feasible fourth-corner region puts the corner on the target. A one-sided inclusion in a shrinking target tube suffices for the limit. This test needs no sweep and can use one background correction at a time. The uncorrected cubic difference alone does not determine its outcome, as the example below illustrates.

## An exact failure of the static three-patch implication

Use the triangle from v363 with middle square side $L=12/7$ and corners, starting at lower left,
$$
 (4/7,0),\quad(16/7,0),\quad(16/7,L),\quad(4/7,L).
$$
Edit only its first three contacts, in the v363 tangent directions. At depth $n$ take displacement $e=2^{-18-n}$ and taper halfwidth $w=2^{-n}/64$. Each support fits in a box of coordinate radius $3w$ about its old contact. Fix the stage cutoff $s_n=(12+13e/2)/7$. If $a,b,c$ are the edit bits and $z\in[0,1]$ is the amount installed, the central square side vector is
$$
 u=(12+ze(a+3b+3c))/7,\qquad v=ze(a-b).
$$
If a bit is absent, $u\le(12+6e)/7$, $|v|\le e$, and $u^2+v^2<s_n^2$. With all bits active the unique crossing is at $z=13/14$, and the fourth corner there is **exactly $(4/7,L)$ for every depth**. The regular four-line system and the nonzero derivative of the side give local crossing parity one.

The complete strip check is the [v363 chart cover](../code/cutoff_interactions_v363.py) with the new cutoff. Its mixed matrices are fixed. At depth zero each rejected chart's interval excludes the entire squared-size segment from $L^2$ to $s_0^2$ (verified for the closed segment in the tests). All these error intervals shrink around their respective nominal values as $n$ increases. On surviving charts the bounds on corner displacement, the connector offsets and the taper widths scale by the same $2^{-n}$. Thus the four cyclic lifts of the one central chart are the only possibilities at all depths. The starting and final faces are clear and the cubic difference is one.

Take as target the fixed v363 triangle with **only the fourth edit** at its original unscaled width and displacement $E=2^{-18}$. It still contains the first three displayed points. Its altered fourth support is $3x-y=-3E$. The old fourth vertex is at Euclidean distance $3E/\sqrt{10}$ from the shifted straight piece (the perpendicular foot lies inside the plateau). All other segments lie farther away. The same complete cutoff cover at $L$, followed by the fourth-edit square formula $u=(12+3E)/7$, $v=0$, excludes every square of side $L$ on this target. Thus even a persistent odd cubic difference, shrinking marked patches and a positive limiting cutoff do not supply the last corner from one static comparison.

Correcting the background at this missing corner changes the answer. Use the target with the fourth edit installed, and superimpose any subset of the three shrinking edits. All these boundary motions converge uniformly as sets to the fixed target. Since the target has no side-$L$ square, compactness excludes side $s_n$ throughout all eight motions for sufficiently large $n$. The corrected cubic difference is therefore zero, whereas its uncorrected value is one. The corresponding fourth difference detects the remote correction, whose support has fixed width. Its forced incidence need not approach the target.

## Transport on a polygon with a sloping roof

Translate the comparison triangle and its three patches horizontally by $l/2$. The target now has a flat bottom from $(0,0)$ to $(4,0)$ and a roof at height $L$ over $[2,4]$. Over $[0,3/2]$ its roof is $y=L+23/28-x$. Join the two roof pieces linearly and use vertical end walls. As $l$ varies, the first two patch centers stay on the bottom and the third stays on the flat roof. All three boxes lie within $3\sqrt2 w$ of the target.

The fourth completion starts below the sloped roof and ends above it. The checker independently certifies strict signs from the boxes alone for all $n\ge1$: for each possible right-angle corner $B$ among the marks $A,B,C$, exclude impossible assignments by the side and diagonal distance intervals, then enclose $A+C-B$ by interval arithmetic. For the surviving assignment the completion box around the translated fourth point has radius $9w$. At depth one and finer its possible vertical and horizontal displacements total less than the endpoint roof margin $1/4$; the other assignments remain excluded as the boxes shrink (their pair-distance bounds separate from both possible cutoffs already at depth one). Intersection with a closed polygon edge makes the sign test inconclusive. The depth-zero enclosure really does fail this test.

To test persistence through branch changes, grow the edits with a piecewise linear schedule in $t$ at knots $0,1/3,2/3,1$. The prescribed values of $z$ are
$$
 0,\quad 9/10+l/10,\quad 17/20+l/10,\quad 1.
$$
At $l=0$ there is one crossing on the final time segment. At $l=2/7$ an additional pair is born. At $l=11/14$ a pair dies. At $l=1$ the surviving crossing lies on the first time segment. In the intervening open range there are three crossings. There is no continuous choice $t(l)$ across the entire sweep: the original last-segment root is isolated from the first root up to the event where the last two meet and cease to exist. All finished masks stay clear by translation and the chart check. The lemma applies. As an independent geometric check, the limit square runs from $x=23/28$ to $x=23/28+12/7$ between heights zero and $12/7$. All four corners lie on the stated boundary.

The open selection problem is to obtain sweep endpoints with different inside contributions for general Jordan data and control marked crossings on the finished face. The static counterexample identifies why an odd cubic difference on its own cannot do this.

[Completion test and exact fixtures](../code/three_patch_sweep_v364.py) · [Tests](../tests/test_three_patch_sweep_v364.py)
