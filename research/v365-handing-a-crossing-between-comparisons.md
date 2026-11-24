# v365. Handing a crossing from one comparison to the next

[v364](v364-three-patch-fourth-corner-sweep.md) has one uncontrolled exit. A marked crossing can reach the finished face $t=1$ of the edit comparison before the sweep forces its fourth corner onto the target. This note chains comparisons in edit time: the next comparison starts from the finished drawings of the previous one, mask by mask. There are three findings.
- **Same patches.** The finished face of a non-final stage is an interior face of the chain. A crossing that leaves one stage there enters the next, and the two counts cancel. An exact replay on the v364 fixture has a crossing move between stages while its fourth corner is still inside. Neither stage satisfies v364 alone, but the chain does.
- **Changed third stretch.** If the next stage edits a different stretch for the same corner, the chain is a v364 comparison with a disconnected patch. The parity of the joint switch is the XOR of the two single-stretch parities plus a term from squares meeting all four patches. When no such square exists, the chain proves nothing that one of the two stretches does not already prove alone.
- **Removing an edit.** Trading one edit for another by blending forces a marked crossing on the finished face whenever the two single-stretch parities differ.

So handoffs make intermediate escapes harmless, but they cannot move parity between separated squares. The finished face of the last stage remains the open exit.

## 1. Gluing stages with the same patches

Keep the setting of v364. There are three patches $K_i(l)$, a side $s>0$, a target $G$ and the inside-sensitive scalar $H$. Stage one has fields
$$
f_S(l,t)=f_0(l,t)+\sum_{i\in S}h_i(l,t),\qquad t\in[0,1],
$$
with increments vanishing at $t=0$. Stage two starts where stage one ends, mask by mask:
$$
g_S(l,t)=f_S(l,1)+\sum_{i\in S}h'_i(l,t),\qquad t\in[0,1],
$$
with $h'_i(l,0)=0$ and $\operatorname{supp}_x h'_i\subset K_i(l)$.

**Lemma 365.1 (gluing).** The concatenation is a v364 comparison with increments $h_i$ followed by $h_i+h'_i$, each supported in $K_i$. Its boundary count uses the common start of stage one, the finished face of stage two, and the sweep ends $l=0,1$. The finished face of stage one is the start face of stage two with the same field for every mask. It appears in both counts and cancels modulo two, including its marked crossings and including those with $H<0$.

v364 allows any continuous joint field whose increments vanish at $t=0$ and stay inside the patches. Concatenation preserves both properties. ∎

Taken alone, each stage can fail. Stage one may start with even parity and may have marked zeros on its own finished face. Stage two has no common start, so its first face does not cancel.

## 2. An exact escape and re-entry

Use the v364 fixture unchanged: the triangle, its three shrinking edits, the cutoff $s_n$, the sloping-roof target, and the translation by $l/2$. Only the installation schedule changes. With knots $t=0,1/3,2/3,1$, the full-mask installed amount $z$ takes the values
$$
0,\qquad \tfrac{13}{14}-\tfrac{3}{56}+\tfrac{l}{8},\qquad \tfrac{13}{14}-\tfrac{1}{56}+\tfrac{l}{14},\qquad 1 .
$$
The full-mask square has side $s_n$ exactly at $z=13/14$, at every depth (v364). The two interior levels pass $13/14$ at $l=3/7$ and $l=1/4$. No pair is born or lost: there is exactly one crossing for every $l$.

| $l$ | crossing time $t$ | stage |
|---|---|---|
| $0$ | $11/15$ | two |
| $1/4$ | $2/3$ | the seam |
| $3/7$ | $1/3$ | one |
| $1$ | $13/42$ | one |

The fourth-corner enclosure of v364 certifies the side from the patch boxes alone. At depth $n\ge2$ it is inside at $l=0$ and $l=1/4$, and outside at $l=3/4$ and $l=1$. At depth 3 it is also inside at $l=3/7$. In the limit the fourth corner meets the roof at $l=1/2$.

Split at $t=2/3$ and count crossings on each face with $H<0$:
- **Stage one alone** ($t\le2/3$) has no crossing at $l=0$, so its start parity is even and v364 does not apply. Its finished face carries one marked crossing, at $l=1/4$ and inside. Its end at $l=1$ carries one outside. If that finished-face crossing is counted, the count does force a cut crossing. But that count is exactly the information unavailable for a general curve.
- **Stage two alone** ($t\ge2/3$) has no common start. Its start face carries the same crossing at $l=1/4$, inside. Its $l=0$ face carries one inside crossing, and its end is empty. The count forces nothing. That is correct: stage two never meets the target, since its crossing exists only for $l<1/4$.
- **The chain** is the v364 sweep with this schedule. Its finished face is v364's at $z=1$, and its chart cover is v364's because the cover bounds every partial installation. The start at $l=0$ is odd and inside, the end at $l=1$ is outside, and the seam cancels. So a square with its fourth corner on the target is forced. It occurs in stage one near $l=1/2$, and its limit is the v364 square.

The crossing starts in the later stage. At $l=1/4$ it leaves that stage backward through the seam, still inside, and meets the target in the earlier stage. Nothing about the seam had to be evaluated.

**For refinement sequences.** The comparison for stage $n$ may contain every earlier stage, and only the last finished face needs control. An escape at an earlier stage costs nothing if later edits in the same patches re-absorb it. The last finished face is still an exit, and gluing does not remove it.

## 3. Changing the third stretch

Now let stage two edit a different stretch $K_3'$ for the third corner, while keeping $h_3$ in place. Bit 3 then switches $h_3+h_3'$, supported in $K_3\cup K_3'$. This patch need not be connected, and its diameter need not be small.

**Theorem 365.2 (disconnected patches).** v364 holds with a third patch $K_3^*=K_3\cup K_3'$ if the diameter bound is replaced by this condition: every side-$s$ square with one corner in each of $K_1,K_2,K_3^*$ has its remaining corner at distance at least $d>0$ from $K_1\cup K_2\cup K_3^*$.

In v364 the diameter bound was used only for two things: to make the three marked corners distinct, and to give $H=\rho(q_*)$ at marked zeros. The stated condition gives both. Cancelling by switching an unseen edit uses only disjointness of the patches. ∎ The checker's union separation test bounds this distance with the v364 completion enclosures.

**Parity of the joint switch.** For any parity function $p$ of the four edit bits,
$$
p(S\cup\{3,3'\})\oplus p(S\cup\{3\})\oplus p(S\cup\{3'\})\oplus p(S)=D_{\{3,3'\}}p(S).
$$
Summing over $S\subseteq\{1,2\}$ gives
$$
D_{12,3^*}p=D_{123}p\ \oplus\ D_{123'}p\ \oplus\ D_{1233'}p .
$$
Under the condition of Theorem 365.2, no square meets all four patches, so the last term vanishes by [v363](v363-cutoff-parity-under-disjoint-edits.md). **The joint switch carries odd parity exactly when one of the two stretches does.** Two stretches that each carry odd parity cancel.

**Theorem 365.3 (reduction).** Assume Theorem 365.2's condition, an odd start with every marked crossing inside, an end with every marked crossing outside, and a last finished face without marked zeros. Then exactly one of the two systems carries odd parity: either $(K_1,K_2,K_3)$ with edits $h_1,h_2,h_3$, or $(K_1,K_2,K_3')$ with $h_1,h_2,h_3'$. That system satisfies the same hypotheses by itself, with the separation condition in place of the diameter bound.

*Proof.* A marked zero of the chain has corners in $K_1$, $K_2$ and exactly one of $K_3,K_3'$. Corners in both would put the last corner inside a patch. If its third corner is in $K_3$, its last corner is at least $d$ from $K_3'$. So it never sees $h_3'$, and the four field values agree with those of the system without $h_3'$, at every parameter. So the chain's face counts split into the two systems' face counts. An odd start makes exactly one system odd, and the end and finished conditions pass to both. ∎

**When the separation fails.** Put $K_3$ and $K_3'$ at the two upper contacts of the v363 triangle. Both are corners of the same square, so a side-$s$ square meets all four patches, and the union separation test returns no bound. v363 found odd degree on the fifteen proper masks and even degree on the full mask. Then $D_{123}p=D_{123'}p=0$ while $D_{12,3^*}p=1$. All of the joint parity comes from the four-patch term. It is carried by a square whose "fourth" corner lies in a marked patch, so it says nothing about the target.

## 4. Trading edits forces a finished crossing

Suppose instead that the handoff removes $h_3$ while installing $h_3'$. On the finished face, bit 3 switches $(1-\mu)h_3+\mu h_3'$ for $\mu\in[0,1]$.

**Proposition 365.4.** Assume the separation condition. The number of marked zeros on the finished face along the blend is congruent to $D_{123}p\oplus D_{123'}p$ (both evaluated on the seam's finished drawings). In particular, if the parity sits on $K_3$ and not on $K_3'$, some finished blended drawing has a marked square of side exactly $s$.

At $\mu=0$, squares through $K_3'$ do not see bit 3, so they cancel in pairs. The marked parity of the finished drawings is $D_{123}p$. At $\mu=1$, squares through $K_3$ do not see bit 3, and the parity is $D_{123'}p$. The parity changes only at marked zeros on the finished face. ∎

The square that needs $h_3$ to exceed the cutoff must pass back through the cutoff as $h_3$ is removed. Parity cannot be moved between separated squares by trading edits.

## What remains

Gluing makes every intermediate escape free, and the last finished face is the only exit left. A useful handoff would have to connect two squares that share all three marked stretches, which is v364 itself. The next two notes test the other two ideas: a search that reads sweeps off the drawing ([v366](v366-sweeps-read-from-the-curve.md)), and an adversary that keeps breaking the best one ([v367](v367-an-adversary-against-the-sweep-certificates.md)).

[Code](../code/escape_handoff_v365.py) · [Tests](../tests/test_escape_handoff_v365.py)
