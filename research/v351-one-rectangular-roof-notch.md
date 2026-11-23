# v351. One rectangular roof notch

Cut a rectangle of width $W$ and height $H$ from the top down to a floor at height $h$, leaving towers of widths $a,d>0$ on the left and right. Write $n=W-a-d>0$ for the horizontal width of the cut and assume $0<h<H$. The boundary has eight edges and two inward corners. This is a first model where the [v350 disk translation](v350-inward-tilts-and-monotone-frontiers.md) can be blocked on one side.

**Theorem.** If $n>h$, the robust largest inscribed square side is at least $\sqrt2 r$, where $r$ is the inradius of the notched region. This holds for arbitrary real dimensions. The proof also gives sufficient conditions for narrower notches. Robustness means that each strict side cutoff below the bound survives sufficiently close registered Jordan redraws.

## An exact reduction for disks, at every notch width

Fill the bottom strip $[0,W]\times[0,h]$ and retain just one tower. This gives two comparison L shapes, both contained in the notched region. For any disk inside the region, its portion strictly above $h$ is convex, hence connected when nonempty. It must lie in one tower since the open notch separates the two towers above $h$. Its remaining portion lies in the bottom strip. Thus every disk lies entirely in at least one comparison L. The two L shapes become nested after reflecting one of them, so **the larger tower gives exactly the inradius of the full region**.

For instance, if $A=\max(a,d)$, a closed disk of radius $s$ fits precisely when $2s\le\min(W,H)$ and $(A-s)_+^2+(h-s)_+^2\ge s^2$, by [v345](v344-v345-convex-and-L-disk-bounds.md). This is a feasibility result about the actual disk; a continuous motion across the notch is unnecessary. In the $14$ by $10$ rectangle with tower widths $2,7$ and floor height $3$, the exact inradius is $7/2$. The disk of that radius grounded at the left corner would enter the notch, but it fits at the right corner.

## Transfer the square only after checking its four contacts

Reflect if necessary to put the wider tower on the left. The notched boundary and this L boundary agree locally everywhere along the L with $x<A+n$. An isolated labeled square component whose contact traces lie strictly in this range therefore has the same local residual, degree and minimum side on both boundaries. Parametrize the matching boundary arcs identically and use an isolating neighborhood inside their product. This permits arbitrary small redraws of the full notched boundary.

If $H>W$, the two end walls have height greater than $W$. The [v340 tall-wall theorem](v339-v340-end-wall-selection.md) gives robust side $W$ (its proof uses real lengths and allows the roof to dip).

Now take $H\le W$. If $A>H$, the comparison L has a horizontal sliding family of side $H$ with contacts at $x\le A$. If $A=H$, the rotated v340 plateau family from v345 has all contacts in $[0,H]^2$: exchange axes and reflect the width coordinate in the normalized plateau construction. Both families therefore transfer, independent of $n$. Their robust floors reach the comparison L's disk target.

If $A<H$, the v345 marked family has span $t=\min(W,H,A+h)$ and every contact satisfies $x\le t$. Thus **$t<A+n$ suffices**, with the explicit minimum side from v345. When $n>h$ it holds automatically. Since the disk target agrees exactly for both regions, this completes the theorem.

## A different family can span the notch

There is another sufficient test even when the wider-tower square touches the cut. Suppose $W>H$ and $W-H<\min(a,d)$. In the uncut rectangle, side-$H$ squares slide horizontally with their left contacts in $[0,W-H]$ and right contacts in $[H,W]$. These bands lie strictly clear of the cut at the roof. The v340 sliding component, rotated by a right angle, therefore transfers intact. Its robust side is exactly $H$ by the height constraint on the reference. Only contacts matter; part of the interior of a square can extend into the notch.

**A boundary case that loses the first certificate.** Set $W=8,H=7,a=d=3,n=h=2$. The left L selects $(2,0),(5,2),(3,5),(0,3)$, of side $\sqrt{13}$. Its second contact is exactly the foot of the right-hand notch wall. Near this square the other three contacts force the second point to run on $(u+3,u)$: the left anchor has ordinate $3$ because the opposite contact is on $x=3$. The incident horizontal floor and right vertical wall both lie on the same side of the line $x-y=3$. Move the right wall slightly to the left, keeping all other supporting lines fixed. The line meets neither incident segment locally: the floor ends before $x=5$ and the intersection with the new wall is below height $2$. Thus this isolated square can disappear. The full-height square above survives the same event with side approaching $7$ under redraws. The [exact tests](../tests/test_rectangular_notch_v351.py) recover zero signed degree for the first component and absolute degree one for the second.

The [selector](../code/rectangular_notch_v351.py) reports a squared floor or an unresolved result. The tests compare exact disk targets on 540 wide-notch parameter rows and independently trace the components on actual notched polygons. When the cut is narrow, the comparison square meets the restored tower, and the two outer contact bands meet the cut, these tests can leave a gap. For instance $W=14,H=8,a=7,d=5,n=2,h=3$ still passes the local L transfer since $t=8<9$. At $W=13,H=7,a=d=6,n=1,h=3$, both strict comparisons hit equality. A separate v341 trace from the axis square then finds 17 incident edge charts, minimum squared side $49/2$ and absolute degree one. This still clears the disk target (squared side $18$). The open question is how to select such a component throughout the remaining narrow-notch parameter range.
