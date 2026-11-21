# v309 — Separate the square degrees by their centers

Benjamin John Schulz · November 2026

The perturbed two-run polygon in [v308](v308-deep-dents-and-square-parity.md) has two actual large squares and total degree zero modulo two. Their centers have horizontal coordinates approximately $-4.000090$ and $3.999944$. Each square carries its own persistent odd degree. For **every** Jordan embedding within registered Euclidean sup distance $1/4$ of the original two-run polygon there is a square centered at $x<-2$ and another centered at $x>2$, each with side at least $13\sqrt2/8$. The existence proof does not need the square locations or transversality on the target. The reference coordinate frame and length units stay fixed.

## One-cap center-wall lemma

Put all trial corners in a horizontal slab of height $H$. Fix a vertical line $x=m$. In the part with $m+e<x\le m+H/2$, suppose every allowed point lies in the top cap of thickness $d$. Any square with allowed vertices and center on the line has half-diagonal
$$
\rho\le\sqrt2(e+d).
$$
Either axis can be reflected, so this applies when the cap is at the bottom or the restricted region is to the left. All conditions are on possible vertex locations.

**Proof.** The horizontal offsets of square vertices from their center are $\pm u,\pm v$ with $u\ge v\ge0$. Their vertical offsets are the same multiset and $\rho^2=u^2+v^2$. In particular $u\le H/2$. For $v>0$ the two vertices on the right have offsets $(u,\pm v)$ and $(v,\mp u)$, so their vertical separation is $u+v$, including when $u=v$. If $v>e$, both are in the top cap, giving $u+v\le d$. If $v\le e<u$, the rightmost vertex is in the cap and its vertical offset has magnitude $v$. The center is at depth at most $d+e$ from the top wall. The highest vertex has offset $u$, hence $u\le d+e$. The remaining case $u\le e$ satisfies the same bound. This proves the lemma.

## From one parity for the whole curve to one per run

Use the horizontal caps and ordered vertical lanes of v308. Write $w$ for maximum lane width, $g$ for minimum intervening gap, $d$ for cap thickness and
$$
a=\max(w,\sqrt2d),\qquad b=\tfrac12\min(g,H-2d).
$$
Assume $a<b$. The v308 geometric exclusion applies to every square whose corners belong to the allowed cover: none has half-diagonal in $(a,b)$. Choose any cutoff $q$ in that open interval. Thus every square counted above $q$ already has $\rho\ge b$.

Require occupied outer cells and one empty cell between consecutive occupied runs (a run may contain several cells). Across the open horizontal gap of each empty cell require that the curve use just one cap. Each gap may choose its own top or bottom cap. Let $[l_i,r_i]$ and $[l_{i+1},r_{i+1}]$ be the lanes bordering that empty cell, with midlines $m_i,m_{i+1}$. Choose dividing lines $x=L\le m_i$ and $x=U\ge m_{i+1}$. Require the midline separation $m_{i+1}-m_i>H/2$. The following additional strict inequalities certify the two walls:
$$
L+H/2<l_{i+1},\qquad U-H/2>r_i,\qquad
\sqrt2\big(d+\max(r_i-L,U-l_{i+1})\big)<b .
$$
For each run take the full-height reference rectangle between its bounding lane midlines. If its width is $W$, require each adjacent dividing line to be inset into it by less than $\min(H,W)/(2\sqrt2)$. Zero inset is permitted. The inequalities keep the center lines ordered and leave room for pure reference squares.

**Localization theorem.** The center region belonging to each occupied run has degree one modulo two above $q$. The intervening regions have degree zero. Consequently there are at least as many distinct squares of side at least $\sqrt2 b$ as occupied runs, including for an even number of runs. The center regions are open vertical slabs between successive dividing lines, unbounded at the two outer ends.

1. Take signed distances (negative inside) for the target and for the union of reference rectangle boundaries, with the sign for the latter taken from the union of the closed rectangles. They agree in sign on the open cells and outside the outer box. In an empty gap, the unused cap can also be removed from the cover: its points outside the retained cover connect to the exterior through that horizontal wall without touching either boundary. Thus the point zeros in the linear interpolation are confined to the restricted cover. This statement concerns the full homotopy; endpoint information by itself would not give local degree invariance.
2. On $x=L$ the one-cap lemma applies with $e=r_i-L$: even a square of maximum horizontal reach cannot reach the far lane on the right. Reflection gives the same argument on $x=U$. The wall bound is strictly below $b$, so there is no simultaneous four-corner zero on a dividing line above the cutoff. The size gap excludes the cutoff face. Center and radius are bounded since the swept point zeros lie in a compact set. Hence the localized geometric-square degree on each slab is invariant. As in [v307 §2](v307-convex-shell-square-gap.md), angles are taken modulo a quarter turn, with the residuals cyclically permuted across the seam.
3. At the reference, a square belonging to one rectangle has $\rho\ge\min(H,W)/2$ by the convex inball estimate. Its horizontal half-span is at least $\rho/\sqrt2$. This puts its center strictly on the run side of its dividing lines. No mixed square (vertices on different reference rectangles) can have its center in a run slab: to hit another rectangle would require horizontal reach greater than $H/2$ across a neighboring gap. The reference degree there is therefore the degree of that one convex rectangle, namely one. Excision is valid even for a continuum of reference squares: near its compact pure zero set all corner signed distances agree with those of the single rectangle.
4. For completeness, mixed reference squares can touch only a consecutive pair of rectangles. The horizontal distance between each successive pair of reference rectangles exceeds $H/2$ by the midline separation test. Two such gaps preclude a nonconsecutive pair because square horizontal span is at most $H$. A mixed center for a consecutive pair is in the intervening slab by the same reach test, using the near wall of each run. For this pair alone, moving the rectangles horizontally apart preserves the global degree above $q$: mixed diameter is at least the gap between rectangles, greater than or equal to $g>2q$, and pure squares stay above $q$. Once the gap exceeds $H$, only pure zeros remain, of degree one plus one. Excision before translation gives degree zero for the mixed set. This yields zero in each intervening slab. No finite-orbit assumption enters these degree computations.

## Two squares despite total parity zero

On the two-run fixture the levels are $y=\pm2$, the cap thickness is $1/32$, and the inner lanes have width $1/4$ centered at $x=\pm2$. Only the top cap crosses the central gap. Here $b=15/8$. The lines at $x=\pm3$ pass: the wall bound squared is $1369/512 < (15/8)^2$. Therefore two separate center regions force squares with side at least $15\sqrt2/8$, one strictly to the left of $-3$ and one strictly to the right of $3$.

The obvious midline $x=0$ would give a false boundary exclusion on the original polygon: $(\pm2,\pm2)$ are four actual boundary points forming a large square centered exactly there. The two dividing lines leave that square in a separate middle slab of degree zero. The unrestricted v308 strip cover also allows an axis-aligned side-four square centered on $x=-3$ if the missing bottom cap is restored. The one-cap assumption is needed for this wall proof.

**Registered perturbations.** Expand every lane by $\eta$ at each horizontal end, the slab by $\eta$ at its top and bottom, and the cap thickness to $d+2\eta$. Do the same on the unused-cap exclusions by shortening each open empty gap to the new lane boundaries. Euclidean perturbations of size at most $\eta$ remain inside that enlarged cover. The cell midpoints are unchanged and have original clearance at least $b>\eta$ whenever the enlarged gap test passes. Their winding numbers survive straight interpolation of the parametrizations. Recompute all wall tests in these enlarged dimensions; the new upper size-gap endpoint is $b-\eta$. Only the target must be a Jordan embedding.

For $\eta=1/4$ take walls $x=\pm2$. The new gap endpoint is $13/8$; the squared upper radius on a wall is $841/512<(13/8)^2$. The registered perturbation used for the exact v308 census has norm below $1/1024$. Relative to this *perturbed census polygon*, any further registered Jordan perturbation of at most $255/1024$ remains within the quarter-unit tube by the triangle inequality. Each side still has odd degree. Additional pairs may form in either region.

**Tracking under motion.** For a continuous path of parametrized Jordan embeddings staying in the certified tube, the signed-distance sections vary continuously: away from the boundary the winding is locally constant, and near it the distance tends uniformly to zero. On each outer slab the same compact continuation argument as [v302 §4](v302-connected-rough-rectangle-crossing.md) transports odd degree on the geometric-square bundle. If the compact zero set had no component joining the two time faces, separating those fibers into disjoint compact sets and excising a neighborhood of the first would contradict degree invariance. Starting at the census polygon there is just one large square per outer slab, so each belongs to a distinct spanning continuum. This proves persistence through possible folds or branching; uniqueness on intermediate curves would need separate hypotheses.

## Independent counts and replay

On the v308 perturbed two-run polygon the exact edge census above half-diagonal 1 yields slab counts **1, 0, 1**, with centers approximately $(-4.000090,0.000030)$ and $(3.999944,0.000031)$. Their squared side lengths are approximately 15.87720 and 15.99829. On the independently perturbed three-run example the slab counts are **1, 0, 1, 0, 3**. All five large squares from v308 are accounted for. This distinguishes degree one from uniqueness and improves the previous guaranteed multiplicity from one to three on the certified three-run class. The small $1/1024$ tube certifies the census center walls with unit inset. Although its wall bound can exceed 1, the open size gap keeps every square above the census cutoff on the large side.

Run `python code/center_walls_v309.py` to regenerate the [exact centers, side lengths and degree certificates](../results/v309/center-walls.json), and `python -m unittest discover -s tests -p 'test_center_walls_v309.py' -v`. The CLI accepts the v308 polygon format via `--input`, with optional `displacement` and `inset`. Reflection of a bridge and different caps on different gaps are permitted. The one-cap and reach hypotheses restrict this theorem to a class of deeply indented curves. Extending localization through a center wall containing actual zeros remains open.
