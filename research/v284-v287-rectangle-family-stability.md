# v284–v287 — Stability of the rectangle family

The v283 diamond witness extends to a two-parameter family with moving supports. An explicit smoothing at its four flat contacts retains the entire framed path. A separate degree argument protects a connected band near the square under continuous perturbations of all four contact arcs. The band also applies to a concave polygon with a finite-coil singularity.

## v284 — A connected near-square band under continuous perturbations

Give four reference segments their affine parameters in $[0,1]$. Write $p=(p_0,p_1,p_2,p_3)$ for their contacts. For $x\in(0,1)$ define the four-coordinate residual

$$
 F_x(p)=(p_0+p_2-p_1-p_3,\;x(p_3-p_0)-(1-x)J(p_1-p_0)).
$$

Here $J$ is the counterclockwise quarter turn. A zero with $p_2\ne p_0$ is a genuine CCW rectangle. Its two side lengths $w=|p_1-p_0|$ and $h=|p_3-p_0|$ satisfy $x=w/(w+h)$. The dot product of the ordered diagonals has sign $\operatorname{sgn}(1/2-x)$ because it equals $h^2-w^2$.

Suppose the affine chart $F_{1/2}$ has invertible parameter matrix $A$ and a root $u_*$ in the open cube. Around $u_*$ take a closed coordinate box $D$ of halfwidth $\eta$ contained in $(0,1)^4$. For each $i$ let $\widetilde p_i(u_i)$ be a continuous replacement of its contact arc, with $\|\widetilde p_i(u_i)-p_i(u_i)\|_\infty\le\varepsilon$ on its tagged window. The replacements need only be continuous. When they are portions of a Jordan embedding, the resulting zeros give inscribed rectangles on that embedding.

Let $C=\|A^{-1}\|_\infty$, and let $M$ bound $\|p_3-p_0+J(p_1-p_0)\|_\infty$ on $D$. The difference between the preconditioned perturbed residual and the translated identity satisfies

$$
 \bigl\|A^{-1}F_x(\widetilde p(u))-(u-u_*)\bigr\|_\infty
 \le C\bigl(4\varepsilon+M|x-1/2|\bigr).
$$

The midpoint error is at most $4\varepsilon$. In the rotated-side equation the error is at most $2\varepsilon$, since $x$ and $1-x$ are positive weights summing to one. If $C(4\varepsilon+M\delta)<\eta$, the straight homotopy to $u-u_*$ avoids zero on the side of the cylinder $[1/2-\delta,1/2+\delta]\times D$. The slice degree is $+1$ in the preconditioned coordinates at every $x$. The lower bound for $|p_2-p_0|$ follows by subtracting both contact window displacements and $2\varepsilon$ from the reference diagonal's infinity norm.

For completeness, the degree bounds give a single compact connected set meeting both ends of the cylinder. Otherwise the compact zero set would split into two disjoint compact subsets containing the respective end fibers, by the compact-space equality of components and quasicomponents (Whyburn separation lemma). Isolate the first subset by an open neighborhood in the cylinder with zero-free relative boundary. The slice degree on this neighborhood is locally constant in $x$: near any given slice the isolated zeros fit in a fixed open subdomain, and compactness excludes other zeros from the neighborhood there, so excision and homotopy apply. At the first end its degree is $+1$, by excision from $D$, and at the other end its degree is zero. This contradiction gives the connected set. Its projection contains the whole $x$ interval, and the two end sections have opposite diagonal-dot signs in the chosen labeling. These end frames are on interior angle sections of the rectangle locus. The conclusion about the flat frame in v286 uses the additional end construction below.

The exact bounds are computed by `code/protected_rectangle_band_v284.py`. The reference square in the diamond chart is the v283 square. The second chart uses the interior square of the v274 concave quadrilateral in forward edge order.

| Four-arc chart | Coordinate sup tolerance $\varepsilon$ | Long-to-short ratios covered | Strict circumradius lower bound |
|---|---:|---:|---:|
| Diamond spiral outer arms and cap | $1/30000$ | $[1,10001/9999]$ | $9/10$ |
| Concave quadrilateral | $1/2000$ | $[1,501/499]$ | $7/10$ |

The certificates require compatible source tags on the four windows. They impose no smoothness there and no restrictions on the rest of the curve beyond the Jordan hypothesis.

## v285 — Independent motion of the cap and the first arm

Put the cap bottom at $V=(0,-h)$, with $1\le h\le101/100$. The left support joins $A=(-1,0)$ to $W=(0,-k)$, with $7/8\le k\le177/200$. The first arm and its reflection have successively decreasing $\ell^1$ radii $1,k,3/4,5/8,1/2$, on the same rays as in v283. Join their inner ends by any Jordan arc in the closed diamond of radius $1/2$. Each outer arm meets each of its radii once, and the two arms at equal radius are negatives of each other. They lie strictly inside the unit diamond except at the outer endpoints. The cap lies outside its interior. These facts verify the Jordan construction, including an arbitrary inner replacement.

Let $q=(a,b)$ traverse the complementary arc from $(1,0)$ to $W$ as in v283. Set

$$
 \begin{gathered}
 B=h+b,\quad L=k(1+a)+b,\quad D=h-k,\quad M=h+k,\\
 p_1=(-\ell,h(\ell-1)),\qquad p_2=(r,h(r-1)),\\
 p_0=(a-\ell-r,b+h(\ell-r)),\qquad \ell=(Mr-L)/D.
 \end{gathered}
$$

The last identity puts $p_0$ on the left support line. Perpendicularity of the adjacent sides is the equation $G(r)=0$, where

$$
 G(r)=-2h(1-kh)r^2+\bigl(2ha-2khB+(1-h^2)L\bigr)r+L(hB-a).
$$

Throughout the open arc $L>0$ and $hB-a>0$. For the latter inequality use $hB-a=(B-ha)+(h-1)(B+a)$, with $B-ha>0$ and $B+a\ge0$. The leading coefficient is negative. Thus the larger root $r$ is the unique positive root in the open arc. At $(1,0)$ it is $1$, with $G_r=-2D$. At $W$ it is $0$, with $G_r=-2khD$. Thus both endpoint roots are simple.

The proof of actual segment membership uses the following five comparisons. With $v=b+k$, their factorizations are

$$
 \begin{aligned}
 G(L/M)&=\frac{DL}{M^2}\bigl((1-kh)(ha-b)+h^3-k\bigr)\ge0,\\
 G((L+D)/M)&=\frac{D(B-ha)}{M^2}
       \bigl(-(1-kh)(ka+b)-hk^2-2h-k\bigr)\le0,\\
 G((ha+b+k)/(2h))&=\frac{D}{2h}
       \bigl(h(a-v)^2+v(2hD-(h-1)^2 a)\bigr)\ge0,\\
 G((ha+b+h)/(2h))&=\frac{D}{2h}(B-ha)(hb-a-1)\le0,\\
 G(B/h)&=-\frac{(B-ha)(2h-k+b-ka)}{h}<0
          \qquad\text{on the open arc}.
 \end{aligned}
$$

All five arguments are nonnegative. The checks use $0\le L\le2k$, $b\ge-k$, $ha+b+k\ge0$, $ha-b>-1$, $hb-a<1$, and $|a|+|b|\le1$. At a fixed $h,k$, it suffices to verify the affine bounds at the outer arc vertices and the four vertices of the inner diamond. At those vertices the expressions are separately affine in $h,k$, so the four parameter-box corners suffice. The only allowable zeros of $L$ and $B-ha$ are respectively the terminal and initial vertices. The remaining scalar bounds follow from $kh<1$, $2h(h-k)>(h-1)^2$, and $2h-k>1$. For example, the first bracket exceeds $h^3-k-(1-kh)\ge0$. The second bracket is negative by $|ka+b|\le1$, and the last bracket's second factor is positive by $|b-ka|\le1$.

The signs give $0\le\ell\le1$, $0\le r\le1$, and $-1\le p_{0,x}\le0$. The open arc has $r>0$ and $B-hr>0$, so $p_2-p_1$ has positive horizontal component and $q-p_2$ has positive vertical component. Their perpendicularity gives a genuine CCW rectangle. The endpoints pair $(A,(1,0))$ and $(W,V)$, with diagonal dot products $4$ and $-(h-k)^2$. Each has positive radius. The bounded-radius labeled path therefore has relative index $1$ against its two flat ends for every allowed $h,k$ and inner arc. It covers every aspect ratio, since the normalized side-length coordinate runs from $0$ to $1$.

The diagonal between $p_0$ and $p_2$ crosses from the half-plane $hx-y\le k$ to the line $hx-y=h$. Thus the squared circumradius is at least $(h-k)^2/(4(h^2+1))\ge529/323216>1/625$, uniformly along the path and throughout the parameter box.

## v286 — Smooth both flat collars without losing the path

There is a local surgery on a rectangle path when each of its two pairs of collapsing corners approaches its flat point in a common inward longitudinal direction. Write the two pair midpoints as $A(s),B(s)$, with $s$ the short half-side length. The original four corners near a flat consist of

$$
 A(s)\pm s n(s),\qquad B(s)\pm s n(s),
$$

where $n(s)$ is a smooth unit normal to $B(s)-A(s)$. Suppose the short-side vector has nonzero derivative in the original path parameter pointing away from the flat. Then $s$ is a smooth local parameter on the one-sided path. For each pair choose a signed longitudinal axis parallel to the flat chord so the midpoint projection has positive derivative at $0$.

Choose $\phi(s)=s\chi(s/\zeta)$, for a smooth nondecreasing cutoff $\chi$ equal to zero near zero and to one for arguments at least $2/3$. Replace the four traces in a small collar by

$$
 A(\phi(s))\pm s n(\phi(s)),\qquad B(\phi(s))\pm s n(\phi(s)).
$$

The two short sides still match, and they are perpendicular to the same long vector as at the old parameter $\phi(s)$. The long length stays fixed at that reference value and the short length grows from $2\phi(s)$ to $2s$, so the squared circumradius does not decrease compared with the reference at $\phi(s)$. The trace joins the old path exactly for $s\ge 2\zeta/3$. At zero, each pair traces one smooth straight subarc through its old flat point, with the two halves occupying opposite sides.

Here is the embedding check for the replacement. Use the signed longitudinal coordinate $\xi$ in either pair and a transverse coordinate $\eta$. The reference midpoint is a smooth graph $\eta=H(\xi)$ near the flat. The longitudinal derivative of a replacement half is bounded below by $c\phi'(s)-K(\phi(s)+s\phi'(s))$ for some $c>0$ and finite $K$: the normal is transverse at zero and has bounded derivatives. Where $\phi>0$, the inequality $\phi\le s\phi'$ makes this positive for small $\zeta$. The two traced points lie on opposite sides of the midpoint graph since their signed deviations are

$$
 \pm s\bigl(n_\eta(\phi(s))-H'(\xi_*)\,n_\xi(\phi(s))\bigr)
$$

for intermediate coordinates $\xi_*$; the bracket has the sign of $n_\eta(0)$ for small $\zeta$. The two halves are therefore disjoint except at the flat point. This check applies where they meet the old traces as well. Each surgery can be contained in arbitrarily small neighborhoods of the two flat points. On a Jordan curve with the local contacts on the two incident segments, choose these neighborhoods away from the remaining compact subarcs and perform the replacement between the unchanged joins. The resulting curve is Jordan, with a smooth arc across each flat contact.

The exact endpoint derivatives give the hypotheses throughout the v285 parameter box. They are measured in the affine parameter on the segment followed by $q$, increasing into the rectangle path. Differentiate $G(r)=0$ and use $\ell=(Mr-L)/D$. At the first flat, $q'=(-1,k)$ gives $r'=-1$ and $\ell'=-M/D$. At the second, $q'=(3/4,k)$ gives $r'=\ell'=7/8$.

| Flat | Two inward longitudinal midpoint derivatives | Short-side derivative |
|---|---|---|
| $(-1,0),(1,0)$ | $M/D,\ 1$ | $(0,-M)$ |
| $(0,-k),(0,-h)$ | $k,\ 7h/8$ | $(7/4,0)$ |

The two pairs at each flat lie on four incident segments and trace initial subsegments with nonzero derivatives. The four flat points are distinct and are outside the inner diamond. Consequently the two surgeries can be performed independently for every member of the moving family, keeping the interior of the path. Choose the collars small enough that the retained interior still has both diagonal-dot signs and includes a square. The endpoint frames are unchanged. The rounded paths have the same uniform radius lower bound $1/25$ as in v285. The reference diamond keeps the v283 bound $1/512$ on squared radius and its exact square. This is an existence construction of smooth flat collars. Connecting the v284 continuum to flat ends under arbitrary perturbations of all contact arcs remains open.

## v287 — A finite-coil example with a different outer geometry

Use the concave quadrilateral $((0,0),(4,0),(1,1),(0,3))$ of v274. Replace the bottom segment from $(2,0)$ to $(3,0)$ by the translated finite-coil arc of [Asano–Ike, Proposition 6.6](https://arxiv.org/html/2609.11598v3#S6.SS2). Its coils lie in disks centered at $(2+2^{-n},0)$ with radii $2^{-3n}$, $n\ge1$. The arc agrees with the bottom line outside the disks. The disks lie in $2<x\le21/8$, $|y|\le1/8$. The next polygon edge has height at least $11/24$ there, and the remaining two edges miss the box. The first tagged contact window ends at $4(9/44+1/100)<2$. All four windows therefore retain the exact affine references used in v284. The polygon corners may be rounded in disjoint small neighborhoods away from the coil enclosure and the tagged windows.

The finite-coil arc is smooth away from its accumulation point and its primitive is unbounded there, by Proposition 6.6. Translation adds a bounded endpoint term to that primitive on the compact arc. Any continuous Liouville primitive of the resulting Jordan curve would restrict to the forbidden primitive on the arc, so the example lies outside Liouville regularity. Finite coil truncations give compatible uniform approximations, since all discarded modifications lie in shrinking disks near the accumulation point. On the limit and all the truncations the same four-arc box gives a connected band with circumradius greater than $7/10$ up to aspect ratio $501/499$. The bounds persist under the v284 tolerance on the four supported windows. This example certifies an interior band and leaves the full flat-to-flat continuation on its geometry open.

## Reproduction

Run `python code/protected_rectangle_band_v284.py --check` and `python code/moving_spiral_support_v285_v286.py --check`. The rational records cover the box estimates, coil enclosure, moving supports, and the hypotheses of the smoothing proof. The SageMath tests independently verify the polynomial comparison identities and replay true rectangles at moved-support samples.

The current results provide perturbation-stable near-square families in two rough geometries and a complete flat continuation for a controlled moving and rounded spiral family. A general criterion that selects such a family on an arbitrary Jordan curve remains open.
