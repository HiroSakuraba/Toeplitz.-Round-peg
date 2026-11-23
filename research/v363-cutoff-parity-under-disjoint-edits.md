# v363. How many separate edits can change square parity together?

Fix a physical side length $s>0$. Draw several edits in disjoint patches and make a new loop by selecting any subset of them. If the cutoff is clear, the parity above $s$ is a Boolean polynomial of degree at most **four** in those choices. A mixed term can be nonzero only if a single side-$s$ square touches every patch in that term during a comparison. An explicit triangle requires a fourth-order term. This gives a finite restriction on the distant crossings in [v315](v315-cutoff-payment.md).

## The locality theorem

Take a continuous defining function $f$ of the starting loop, negative inside and positive outside. Let $h_i$ be continuous changes with pairwise disjoint compact supports $K_i$. Write $f_S=f+\sum_{i\in S}h_i$. When speaking of loops, assume each referenced endpoint field defines the corresponding Jordan boundary with the correct signs. Only its zero section is needed during a homotopy. This model includes replacements made in disjoint disks which retain their exterior boundary and inside/outside data: take defining functions of the edited loops and blend each into the old function in a collar away from the changed arc where their signs and zero sets already agree. The literal signed-distance differences need not have compact support.

Let $p(S)$ denote the mod-two degree of physical squares above $s$ on the endpoint. Use centers, sizes and orientations modulo a quarter turn as in [v310](v310-moving-degree-and-scale-loss.md). All planar zero sets throughout the comparisons are bounded. Choose common outer faces large enough to be zero-free. Assume first that no endpoint realizes $s$. For a fixed background $B$ and nonempty set $I$ disjoint from $B$, form the mixed difference
$$
 D_Ip(B)=\bigoplus_{S\subset I}p(B\cup S).
$$
If this is one, at some time of the comparisons $f_B+t\sum_{i\in S}h_i$ there is an exact side-$s$ square whose four vertices are zeros of the intermediate field and which has a vertex in **each** $K_i$ for $i\in I$. In particular $D_Ip(B)=0$ for $|I|>4$. Already for two patches the difference is zero if their pair-distance range contains neither $s$ nor $\sqrt2 s$. For four patches with an odd difference the four contacts must belong to four different patches.

**Proof.** The cutoff wall consists of center, angle modulo a quarter turn, and time. At each point $z=(Q,t)$ its section evaluates the four field values at the vertices of the physical square $Q$. For a regular wall section, the endpoint degree change is the number of wall zeros modulo two. The sum of starting degrees over the cube cancels.

For continuous inputs, smoothly approximate the base field and the increments on the common bounded domain, the latter supported in arbitrarily small disjoint open enlargements of the $K_i$. Preserve all endpoint and outer zero-free margins. Add **one common small smooth section perturbation** on the compact wall, making all finitely many sections transverse at once. This can be done using finitely many local bundle frames and bump functions spanning each fiber, and choosing a regular parameter for all the resulting families. Keep boundary collars free of zeros. This preserves each original wall degree by homotopy. The value of this common perturbation is independent of $S$.

At a fixed regularized wall zero, if an edit in $I$ is unseen by all four vertices, toggling the first such edit pairs two zeros at the identical $z$. For this counting argument membership in an enlarged support depends only on $Q$, so the pairing is an involution. An odd total leaves a zero seeing all the enlarged patches. Shrink the enlargements and both perturbation errors to zero. Compactness and the finite number of masks give a limiting cutoff square with each required incidence and all four field values zero. The endpoints are clear. Disjointness of the closed patches then gives the order bound. This counts geometric squares once, including when the original unperturbed wall had a continuum.

## Recovering the full degree

Every Boolean function on masks has the expansion $p(S)=\bigoplus_{I\subseteq S}D_Ip(\varnothing)$ with $D_\varnothing p(\varnothing)=p(\varnothing)$. Thus only subdrawings with at most four active edits enter the full-degree reconstruction. The weight on the value of a subdrawing with $j\le4$ active edits out of $m$ is
$$
 \sum_{a=0}^{\min(m,4)-j}{m-j\choose a}\pmod 2.
$$
For this reconstruction alone it suffices that the full, empty, and contributing small endpoints are cutoff-free. To see this without assigning a degree to a singular omitted endpoint, use the *same* small perturbation for all wall comparisons. At each of the finitely many transverse wall zeros in their union, the binary indicator of a hit, as a function of the mask, depends only on the edits seen by its four corners. Apply the Boolean identity directly to the indicator there and sum. The wall counts at the participating regular endpoints equal their degree changes. For omitted endpoint fields if necessary use a common small generic perturbation also at the time faces: each face has dimension three and fiber rank four, so the perturbed face has no zeros. Zero-free margins at the participating ends protect their original degrees.

If no side-$s$ square can visit any pair of patches, this specializes to
$$
 p(\mathrm{all})=p(\varnothing)\ \oplus\ \bigoplus_i\big(p(\{i\})\oplus p(\varnothing)\big).
$$
Independent odd changes cancel in even numbers. The patches can have different shapes, growth rates and locations. Equal dynamics and a spatial symmetry are unnecessary. This is an identity at a prescribed cutoff; computing a small subdrawing's degree can still require its entire boundary.

**Four contacts on the target.** Suppose on a sequence of comparisons an odd fourth difference is certified on four disjoint patches, the cutoffs tend to $s>0$, and every one of those patches is contained in the $\eta_n$-neighborhood of one fixed Jordan target with $\eta_n\to0$. The locality theorem produces a square with all four corners in that neighborhood, even if the rest of the comparison stays far away. A convergent subsequence gives a square of side $s$ on the target. Patches of diameter tending to zero that meet the target suffice. An odd lower-order term certifies only the corresponding number of patch incidences. That distinction matters for the unrestricted problem.

## Choosing a cutoff from the patch geometry

For two, three or four closed axis-aligned boxes the checker computes **all** squared side lengths for which a square can visit every box at a distinct corner. For each assignment of corners, write its vertices using a free translation and a side vector $(u,v)$. The offsets from the first vertex are $(0,0),(u,v),(u-v,u+v),(-v,u)$. In each coordinate, a common translation exists exactly when the intervals permitted by all the boxes intersect. Comparing each lower bound with every upper bound eliminates the translation and leaves closed linear inequalities on $(u,v)$. Two assigned distinct corners already bound this polygon. Empty, point and segment cases are included.

On a nonempty feasible polygon the attainable squared sides form a closed interval by connectedness. The maximum of $u^2+v^2$ occurs at a vertex. The minimum occurs at the origin, a vertex, or the foot of the perpendicular from the origin to a supporting line. These tests and the union over corner assignments are rational. In a proposed positive band the code chooses the midpoint of the widest gap outside the forbidden size intervals of the specified groups. The output is a squared cutoff. It proves that those groups cannot contribute mixed terms there. Endpoint clearance and the remaining parity contributions still require verification.

For example put boxes of halfwidth $1/100$ around $(0,0),(1,0),(1/2,7/8)$. Every pair can be hit by a unit square. The exact joint test excludes a simultaneous hit on all three at every side length and selects squared cutoff one in the band $[99/100,101/100]$. This rules out a three-way interaction that no pairwise distance test at that cutoff excludes.

## Four distant crossings in an asymmetric folded family

Start from $[-12,12]\times[-1,1]$. Read the height $H=2$ from the end walls of the target and take $s=H/5=2/5$. The unedited rectangle has zero-free cutoff and odd degree by [v313](v313-half-turn-degree.md). On the top edge install the v313 overhanging profile at $x=-9,-4,1,6$, with final scales $19/20,24/25,49/50,1$, respectively. Grow the four from zero proportionally to their own final scales, then translate all sites two units to the right. At zero scale the ears are absent. This gives continuous signed defining fields, for example from signed distances of the growing Jordan domains. Add a shallow backtracking fold to the lower baseline near the origin.

Each top event is confined to one bay and its baseline. Indeed the opposite row and end walls are out of reach for a cutoff square meeting the protruding part. The old sites are five units apart. The horizontal and vertical spans of a square are equal, so its horizontal reach from the attachment is at most $7/16+17/32<1$ at scale one, and scales proportionally. Thus the finite baseline segments in the v313 profile census suffice. The only two positive normalized squared sides are $1386334197/9194892100$ and $533/2890$. Only the larger can reach $s$. For final scale $a$, it crosses once at time $t^2=s^2/(a^2(533/2890))$. Both the interior contact system and the scale crossing are regular. Translation produces no events. A square entirely on the original rectangle would have side two.

The bottom edit may be any Jordan replacement supported in $[-1,1]\times[-1-1/4,-1+1/16]$ agreeing with the old boundary and sides outside the patch. A side-$s$ square meeting that edit sees only the bottom boundary in this neighborhood, and all such contacts are within a horizontal band of height $5/16<s$. Its vertical span would be at least $s$. The same exclusion holds along a defining-field interpolation, since any new zeros off the baseline are confined to the patch. Thus arbitrarily many folds there are allowed.

Every one of the four transverse top crossings has a vertex of height at least $1+s$. Its horizontal distance from the final top protrusions is at least $2-(7/16+17/32)-7/16=19/32$. The end walls are farther away. Consequently each crossing has a corner at sup distance at least $2/5$ from the final boundary. These genuinely distant events have even total degree despite four distinct crossing times and four unequal profiles. The final cutoff is clear and the degree above it is odd, for all the allowed bottom replacements. The result here is the global parity accounting. This wide rectangle also contains simple directly visible boundary squares.

## Four edits really can interact at a single cutoff

Use the triangle $(0,0),(4,0),(1,3)$ from v310. Its horizontal-base square has side $12/7$. Start at its lower left corner and label the four corners counterclockwise. In small disjoint neighborhoods install parallel straight pieces, with short sloping connectors back to the original lines. Bits $a,b,c,d$ activate the four edits. Take $e=2^{-18}$. The bottom two pieces shift down by $e$, the upper right piece shifts up by $e$, and the upper left one shifts left by $e$. The explicit polygon and all sixteen endpoint checks are in the replay. At fractional edit time $t$, the first side vector is
$$
 u=(12+te(a+3b+3c+3d))/7,\qquad v=te(a-b).
$$
This follows by solving the two bottom height equations and the line equations $x+y=4+tce$, $3x-y=-3tde$. Set $s=(12+19e/2)/7$. With any bit absent, $u\le(12+9e)/7$, $|v|\le e$, giving $u^2+v^2<s^2$ throughout. With all bits present, $v=0$ and the square crosses transversely once at $t=19/20$.

**Why this is a global parity test.** Every swept vertex is within line-residual error $3e$ of one of the three original supports $y=0$, $x+y=4$, $3x-y=0$. The checker covers all $3^4=81$ assignments of these strips to square corners. One strip alone has insufficient projection width. For every mixed assignment the four-by-four line matrix is invertible. Its inverse bounds all coordinate errors, and hence the interval of squared side lengths. Only the four cyclic choices of the horizontal-base square can reach the cutoff. There each corner's horizontal deviation stays below $1/64-e$, excluding the sloping connectors and placing the corners on their respective straight pieces. The original census has exactly one square above this cutoff, the square of squared side $360/121$. The global degree stays odd on all fifteen proper subsets and becomes even on the full set. The fourth difference is therefore one. A rule using only additive, pairwise and triple interactions would fail on actual Jordan polygons, even for these very small edits.

## Remaining proof obligation

For arbitrary input we still need a positive cutoff and a decomposition whose distant contributions are even or whose uncancelled contributions put all four corners near the target. The order bound removes interactions of order five and higher for disjoint edits at a fixed cutoff. It supplies neither such a decomposition for overlapping refinements nor parity control of the remaining one-, two- and three-patch contributions.

[Exact checker](../code/cutoff_interactions_v363.py) · [Tests](../tests/test_cutoff_interactions_v363.py) · [Previous cutoff transfer](v315-cutoff-payment.md)
