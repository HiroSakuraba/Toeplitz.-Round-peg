# v360. Three traversals and a free fourth contact

[v359](v359-fixed-arcs-through-refinement.md) prescribed four narrow strips. Here only three boundary windows are marked. The fourth corner may land anywhere on the loop. Two inside/outside probes replace its strip. This specializes the [v312 winding approach](v312-two-arc-winding.md) to an explicit test that tolerates folds in all three marked paths. The colors refer to fixed parameter arcs across a sequence, even when the probes or the fourth contact change.

## Certificate

Choose orthonormal coordinates and $L>2e\ge0$. A Jordan curve has three disjoint marked traversals in the strips $|y|\le e$ (bottom), $|x-L|\le e$ (red) and $|x|\le e$ (blue). Write $X,Y,W$ for their tangential coordinates respectively. Order each traversal from low endpoint to high endpoint, independently of boundary orientation. Choose $a<b$ inside the range of the blue endpoint ordinates. Require the bottom endpoint coordinates to strictly bracket $[L-b-4e,L-a+4e]$ and the red endpoint ordinates to strictly bracket $[L-b-8e,L-a+8e]$. Tangential backtracking is allowed.

If the probes $(a,L)$ and $(b,L)$ are on opposite sides of the Jordan curve and each has Euclidean distance greater than $8e$ from it, there is a square whose opposite corners visit red and blue. Its side is at least $(L-2e)/\sqrt2$. No parametrization or strip condition is imposed on the fourth contact. The inside/outside and clearance hypotheses are global conditions on the boundary.

**Proof.** Trim the blue traversal to its last visit to ordinate $a$ before its first subsequent visit to $b$ (take the first visit to $b$ after the low endpoint). Throughout the trimmed arc $a\le W\le b$. For one point $P$ on bottom, $C$ on blue, and $T$ on red, the required red and fourth points in square order $P,T,D,C$ are $P-J(C-P)$ and $D=C-J(C-P)$, where $J$ is counterclockwise rotation by a quarter turn. On the three-parameter box use two coordinates of $T-[P-J(C-P)]$ and the signed distance of $D$ to the **whole** Jordan curve. A zero is a genuine square since red and blue are separated.

Homotope the three traversals parametrically to linear traversals of their projected supporting lines, preserving tangential endpoint coordinates and the transverse strip bounds. Keep the signed-distance function of the original curve fixed throughout this calculation. Intermediate traversals need not be embedded. If the first two residuals vanish, projection gives
$$
 |X+W-L|\le4e,\qquad |Y-X|\le4e.
$$
Indeed each projected point differs from its trial point by at most $e$, and quarter turns preserve distances. Thus neither the bottom nor red parameter can be at an endpoint. On either blue end, the fourth completion is at most $3e$ from $(W,X+W)$ by projection, hence at most $7e$ from its probe. It cannot meet the boundary. This excludes every boundary zero throughout the homotopy.

On the straight traversals the two residual equations reduce to $X=L-W$, $Y=X$. In the third residual we can replace the completion $(W,X+W)$ by $(W,L)$ using straight interpolation, since they already agree whenever the other two residuals vanish. The signed distance along the probe segment has opposite endpoint signs. Homotope it to the affine interpolation of those endpoint values: the probe-end signs remain fixed, and the first two coordinates still exclude the other faces. The final system has one nonsingular interior zero in linear traversal coordinates. Its mod-two degree is one. The original system therefore has a zero regardless of how many times the free fourth contact could occur, even along a continuum. The red-blue horizontal separation yields the claimed side bound.

A successful reference also protects every strict cutoff below its floor for sufficiently small registered Jordan redraws. On the fixed parameter box there is a positive boundary residual margin. Under a registered perturbation of size $\delta$ the distance magnitudes change by at most $\delta$ at fixed points. Outside the $\delta$-tube the interior sign is unchanged by winding invariance; inside the tube a sign change costs at most $3\delta$ in signed distance. Completions move uniformly as well. Thus the degree and any strict size margin persist. Across larger edits one can instead verify the three traversals and the two probe signs anew, keeping the colors fixed.

## An actual contact exchange

Keep the bottom and red windows centered at $(6,0),(11,6)$ and the blue window at $(0,5)$, each extending half a unit in both tangential directions as in v359. Take $L=11$ and $e=1/256$. In the v357/v359 polygon, reroute the upper boundary after $(8,12)$ along
$$
 (5.1,12),(5.1,10),(5,10),(5,k),(4.9,k),(4.9,9),(2,9).
$$
Then resume its original route. Keep the original parameter intervals for the three traversals. Use the same probes $(4.8,11),(5.2,11)$ for every $10<k\le12$. They are outside and inside respectively; both stay at least $0.1$ from the boundary. The entire family therefore keeps opposite corners on the same colored windows, with side at least $(11-2/256)/\sqrt2>7.77$. This remains true for any Jordan edits preserving the three traversal hypotheses and the probes, including arbitrarily many folds.

With the three marked traversals straight, the fourth point runs along $(w,11)$. Above $k=11$ it can hit $w=4.9,5,5.1$. At equality the first two hits merge into a horizontal contact interval, which has zero degree. Below equality that entire interval is gone. The $w=5.1$ contact remains odd and sits more than $1/32$ away from the supporting line $x=5$ of the former fourth window. Thus at the last stage no square from these three straight windows can use that old fourth strip. The three-variable degree counts all eligible fourth contacts together. It survives the loss without choosing a branch.

[Exact probe scanner](../code/free_contact_v360.py) · [Tests](../tests/test_free_contact_v360.py). Given rational polygon windows, the scanner checks strip bounds and scans midpoint probes between level-crossing events, using exact point-to-edge distances and ray parity. A miss is inconclusive. Tests run the scanner on the three stages with straight and independently folded windows and separately check the disappearing and surviving components on the straight cases with v341. The general problem is now to select three suitable traversals and find opposite probe signs, or to replace these conditions when every feasible probe has the same sign. No universal arc selector is established.
