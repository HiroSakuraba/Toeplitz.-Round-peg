# v354. A square survives when the second cut swallows the turning endpoint

In the [v353 example](v353-two-cuts-and-degree-transfer.md), widen the second roof cut $[3,q]$ past $q=7$. The detached turning interval collapses and disappears. At floor height 4, the different square $(5,0),(8,5),(3,8),(0,3)$ remains fixed before, at, and after this event. Its side is $\sqrt{34}$, above the disk target of at most $\sqrt{32}$. The replacement also works, with a change of contacts, at every possible second floor height after the endpoint is covered.

## Parameter family

Use box height $B$ and width $B+d$, first cut $[A,B]$ at floor height $h$, second cut $[p,q]$ at floor height $g$. Assume
$$
0<d<p,\quad h>0,\quad p+h<A<B,\quad 2h\le A,\quad p<B-h\le q<A,\quad 0<g<B .
$$
The old turning endpoint has top abscissa $B-h$. Every shape in this parameter family satisfies the robust disk bound: each strict side cutoff below $\sqrt2 r$ survives sufficiently small registered Jordan redraws. In particular $A=8,B=10,d=1,h=p=3$ permits every $7\le q<8$ and every $0<g<10$.

**What happens to the old family.** Just before the event, the detached interval of v353 has top abscissa $v\in[q,B-h]$. At equality this leaves one square. Near it the bottom and left anchors are interior, and completing them gives the other contacts $(t,u),(v,t)$ with $t=u+v$. At the first cut's right wall or adjacent floor, $t\le B$ and $u\ge h$, hence $v\le B-h$. At the second cut's right wall or adjacent top, $v\ge q$. Widening the cut therefore removes every square in this marked neighborhood. This also proves the singleton at equality has zero degree. All other parameters are held fixed, including the strict floor clearance $g<B$.

**Disk size.** Filling the second cut can only increase the inradius. By the [one-notch disk reduction](v351-one-rectangular-roof-notch.md) the filled shape has the inradius of the comparison L with tower width $A$ (the other tower has width $d<A$). That radius is exactly $A/2$ by the [orthogonal L formula](v344-v345-convex-and-L-disk-bounds.md): $A\ge2h$ and both outer spans exceed $A$. Thus it suffices to protect a square spanning at least $A$ horizontally, since its side is at least its horizontal span divided by $\sqrt2$.

## The replacement depends on which boundary pieces meet

Use the same bottom and left anchors $(u,0),(0,v)$, but now mark the **left** wall or adjacent floor of each cut. The two completed contacts are $(t,u),(v,t)$ with $t=u+v$.

| Second floor | Span and anchors | Two notch contacts |
|---|---|---|
| $g<A$ | $t=A$, $v=p$, $u=A-p$ | One on each vertical wall |
| $g>A$ | $t=g$, $u=h$, $v=g-h$ | One on each horizontal floor |
| $g=A$ | $t=A$, $h\le u\le A-p$, $v=A-u$ | A sliding interval on the first wall and second floor |

In the first row, $h<A-p<B$ and $g<A<B$ make the wall contacts interior. This row works even before $q=B-h$, since it only uses the second cut's left wall. In the second row $A<g<B$ and $p<g-h<B-h\le q$, so both floor contacts are interior. The anchor equations $(u+v=A,v=p)$ or $(u=h,u+v=g)$ in the respective rows have independent normals. Completing the square determines the other two edge parameters at nonzero speed. Thus the four-parameter residual has invertible derivative and absolute degree one.

At equality the selected interval is compact. Its outer contacts avoid the far ends of both cuts: $A<B$ and $A-h<B-h\le q$. Near it the first marked contact requires $t\ge A$, $u\ge h$, with equality in at least one; the second requires $t\ge g$, $v\ge p$, with equality in at least one. If $t>A=g$ both imply $t=h+p<A$, impossible. This isolates exactly the stated interval. Lower $g$ slightly. In the same neighborhood $t\ge A>g$ requires $v=p$; $t>A$ would again force $u=h$. Only the first row's regular square remains. Registering boundaries edge by edge and applying degree invariance proves odd degree for the whole interval. Uniform registered redraws preserve the degree on a small isolator. Every member has span at least $A$, which proves the claimed robust lower bound.

For the example at the start, the squared side is 34 for $g<8$, the minimum squared side of the sliding interval is 32 at $g=8$, and the squared side is $9+(g-3)^2$ for $8<g<10$. These are floors from the selected component; other squares can be larger.

[Exact selector](../code/swallowed_endpoint_v354.py) · [Independent component checks](../tests/test_swallowed_endpoint_v354.py). Thirty parameter rows check the component degree and minimum squared side, including both equalities and fractional dimensions. Three more traces follow the fixed wall square through the original endpoint event. The argument uses $2h\le A$ to compare sizes; it does not establish the disk conjecture for arbitrary two-notch rectangles.
