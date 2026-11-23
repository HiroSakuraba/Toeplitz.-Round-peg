# v353. A second roof cut can move the degree to a detached family

In [v352](v352-all-rectangular-roof-notches.md) a square sliding across a rectangle could turn at a notch wall. A second cut can split that turning interval. We give a family in which the separated piece carries the disk bound, and classify the rotating components for any finite number of cuts in a square box. All cuts are disjoint rectangular bites in the top edge, with positive floors, widths and intervening roof lengths. Robustness refers to arbitrary sufficiently small registered Jordan redraws and every strict cutoff below the stated size.

## Cutting through the turning interval

Write the first cut as $[A,B]$ with floor $h$, box width $B+d$ and variable height $H$. The second cut is $[p,q]$ with floor $g$. Assume
$$
0< B-A\le h<B,\qquad 0<d<p<q<B-h\le A,qquad 0<g<B.
$$
At $H=B$ the left end of the old sliding interval turns along the first notch wall $x=B$. Its moving top contact used to traverse $[0,B-h]$. The second notch removes its stretch $(p,q)$. Focus on the detached interval $q\le v\le B-h$, with $u=B-v$ and corners
$$
(u,0),\quad(B,u),\quad(v,B),\quad(0,v).
$$
It has odd local degree and minimum squared side $v_*^2+(B-v_*)^2$, where $v_*$ is the point of $[q,B-h]$ closest to $B/2$. In particular its minimum side is at least $B/\sqrt2$, which reaches the disk target since the height is $B$.

**Isolation.** All bottom and left contacts are interior. If these anchors are $(u,0),(0,v)$, the other contacts are $(t,u),(v,t)$ with $t=u+v$. In a neighborhood of the interval the first of these contacts is on the right wall of the first cut or its floor: $t\le B$, $u\ge h$, and $t<B$ implies $u=h$. The second is on the intervening tower's top or its incident walls: $q\le v\le A$, and $t<H$ implies $v=q$ or $A$. For $H=B$ and $t<B$ this would give either $t=h+q$ or $h+A$. The latter is at least $B$. The former is strictly less than $B$ by a fixed amount, so exclude it by taking a small neighborhood of the family. This leaves precisely the displayed closed interval at $t=B$.

**Counting the detached piece.** Increase $H$ slightly above $B$, raising $h$ slightly as well only if $h=B-A$. Use the same marked neighborhood and edge-by-edge registration. Now $t\le B<H$ forces $v=q$ or $A$. The latter is impossible since $u\ge h'$ and $A+h'>B$. At $v=q$, the floor alternative remains outside the neighborhood, so $t=B$, $u=B-q$. Both nonanchor vertices are strictly inside vertical edges: $g<B<H$ and $h'<B-q<H$. The reduced equations $u+v=B$, $v=q$ have independent normals. The edge parametrizations then recover the other two marked parameters at nonzero speed. This unique regular square gives degree of absolute value one; homotopy invariance returns it to the original interval. On a small isolator its boundary residual stays nonzero under uniform redraws. Shrinking it gives the minimum-side guarantee. This counts the detached piece independently of all squares still attached to the original slide.

**Both sides of the height transition.** For $\max(A,h+q,g)<H<B$, choose $u=h$, $v=H-h$. Then $q<v<A$ and $A<H<B$. The two nonanchor contacts are interior to the first floor and the intervening roof, fixing $u=h$, $u+v=H$ transversely. This square spans the entire height.

For any $H>B$, the regular wall square $u=B-q,v=q$ from the degree calculation is still present. To compare its size with interior disks even for tall boxes, put $G=\max(g,h)$ and $T=\max(p,A-q,d)$, the largest tower width. A disk whose center has height at most $G$ has radius at most $G$ by the bottom boundary. Above $G$ its horizontal diameter lies in a single tower, so $2r\le\max(2G,T)$ for every disk. Hence $B\ge\max(2G,T)$ suffices for the square to reach the disk target. This is a sufficient test on the tall side of the transition.

**Example.** For width 11, take cuts $[8,10]$ at height 3 and $[3,4]$ at height 4. At outer height 10 the detached interval has minimum squared side 50 and absolute degree one. The independent polygon trace from the axis square has minimum squared side 58 but degree zero: selecting the component with the larger minimum would lose this certificate. Below the transition the new square has side $\sqrt{9+(H-3)^2}$ for $8<H<10$. Above it, the square $(6,0),(10,6),(4,10),(0,4)$ has side $\sqrt{52}$ for every $H>10$. It reaches the disk target because the disk diameter is at most 8. Zero degree alone makes no claim that a component must disappear under a particular redraw.

## Any finite number of roof cuts in a square box

**Theorem.** Cut one or more disjoint rectangular notches from the top of a side-$S$ square, preserving positive roof lengths at both outer corners and between the cuts. The robust side is at least $S/\sqrt2$, hence at least $\sqrt2 r$, for arbitrary cut depths strictly between zero and $S$.

Start with the turning family $(u,0),(S,u),(v,S),(0,v)$, $u+v=S$. It meets the remaining roof precisely for $v$ outside the open cut intervals. The two outer roof intervals join through the axis square at $v=0,S$ after one cyclic change of labeling, as in v352. Take this one continuous labeling across that join. The cut gaps separate the other cyclic lifts. Each interior roof interval yields a separate labeled component. Each selected square spans $S$.

At a nonaxis member the bottom, right and left contacts are interior and force $t=u+v=S$. The top contact therefore cannot descend a notch wall even at a cut endpoint. Near the axis square all corners have the unchanged local geometry of the uncut box. These facts isolate each listed component on the four-parameter square residual. They do not assert that all other inscribed squares have been listed.

Increase only the outer width to $S+e$ for small positive $e$, keeping every cut fixed. Near any nonaxis member the three unaffected contacts would force the fourth point to height $S+e$, which is unavailable. Near the axis square the local rectangle has just the horizontal slide of side-$S$ squares (v352); for sufficiently small $e$ all contacts avoid all cuts. This sliding interval has odd degree by v340. By compactness every nearby zero in each isolator must tend to that isolator's reference component. Thus the joined component has odd degree and each interior roof component has degree zero. The minimum of the joined component is at least $S/\sqrt2$ and the same strict-cutoff transfer applies.

[Exact component tests](../tests/test_two_notches_v353.py) check the detached interval, single wall or floor square on either side, ties and rational heights, and the predicted degree distribution in square boxes with two and three cuts. General two-cut rectangles of unequal width and height remain open, in particular when a cut swallows the end of the turning interval rather than splitting its interior.
