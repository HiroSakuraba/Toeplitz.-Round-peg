# v296–v299 — A rectangle witness with four rough flat contacts

A scale-dependent positional estimate can replace derivative control all the way to a positive-radius flat. The resulting collar supports arcs that reverse their transverse projection at arbitrarily small scales. Applied at both ends of the rounded concave witness, this places finite-coil accumulation points at all four flat contacts. A first-order counterexample identifies the role of the quadratic information.

## v296 — A four-contact rescaling using only positions

Normalize by a prescribed flat chord of length $L>0$ and use axes $e,n$ with $Jn=e$ for the counterclockwise quarter turn $J$. In the normalized coordinates the two distinct flat points are $A=(0,0)$ and $B=(1,0)$. Suppose the two embedding arcs obey

$$
 \begin{aligned}
 P_A(t)&=(a t^2,t+b t^2)+E_A(t),\\
 P_B(t)&=(1+c t^2,t+d t^2)+E_B(t),\\
 \|E_i(t)\|_\infty&\le\varepsilon t^2\qquad (|t|\le t_0).
 \end{aligned}
$$

The signed $t$ is the compatible source parameter. It gives the common leading transverse direction at the flat, while the arc positions away from zero need only be continuous. The two arcs may reverse their transverse projection. In particular, the proof retains four independent source variables.

Prescribe short-to-long ratio $2r>0$, and take consecutive geometric corners $P_A(rz_0),P_A(rz_1),P_B(rz_2),P_B(rz_3)$ in a box about $z_*=(-1,1,1,-1)$. Use the midpoint difference $M=p_0+p_2-p_1-p_3$ and similarity difference $R=2r(p_3-p_0)-J(p_1-p_0)$. Scale their components in the order

$$
 H_r(z)=\bigl(M_e/r^2,\ M_n/r,\ R_e/r,\ R_n/(2r^2)\bigr).
$$

For the exact quadratic models this has the limit

$$
 G(z)=\left(
 a(z_0^2-z_1^2)+c(z_2^2-z_3^2),\
 z_0+z_2-z_1-z_3,\
 2+z_0-z_1,\
 z_3-z_0+\tfrac a2(z_1^2-z_0^2)
 \right).
$$

Its root at $z_*$ has Jacobian

$$
 J_0=\begin{pmatrix}
 -2a&-2a&2c&2c\\
 1&-1&1&-1\\
 1&-1&0&0\\
 a-1&a&0&1
 \end{pmatrix},\qquad
 \det J_0=4(c-a-2ac).
$$

This is twice the determinant of the chord-length Hessian from v292 in the normalized quadratic model. Assume $c-a-2ac\ne0$. Both rounded concave ends satisfy it. Indeed, if the v292 rounding parameter is $\zeta=L/K$ at the respective end, the jets in transverse source coordinates have $(a,c)=(24K/5,-12K/5)$ and $(K,2K/5)$, with $b=d=0$. Reparametrizing the smooth references by their transverse coordinate preserves the longitudinal quadratic terms, since the longitudinal first derivative vanishes at the flat. A common sufficiently large $K$ gives roundings small enough for v292. The determinant is nonzero at both ends for every $K\ge1$.

### Exact lateral degree and size bounds

Let $C=\|J_0^{-1}\|_\infty$, $B_2=\max(1,2(|a|+|c|),|a|)$, and $M_2=8\max(1,|a|+|b|+|c|+|d|)$. With source halfwidth $h\le1/4$ and $0<r\le1$, all $|z_i|$ are at most two. Expanding the quadratic model and bounding the position errors gives

$$
 \bigl\|J_0^{-1}H_r(z)-(z-z_*)\bigr\|_\infty
 \le C(B_2h^2+M_2r+24\varepsilon).
$$

Here the Taylor remainder of $G$ contributes at most $B_2h^2$. The exact model's three lower-order corrections are bounded by $M_2r$. In the most sensitive error row, $M_e/r^2$, four contact errors contribute at most $16\varepsilon$. The row $R_e/r$ contributes at most $16\varepsilon r^2+8\varepsilon r$, and the other rows have smaller bounds. Thus the bound is uniform down to arbitrarily small $r$.

The rational certificate uses

$$
 h=\min(1/4,1/(8CB_2)),\qquad
 \varepsilon\le h/(192C),\qquad
 r_*\le\min\!\left(1,\frac{h}{8CM_2},\frac{t_0}{1+h},
                    \frac{1}{16(1+|a|+|c|+2\varepsilon)}\right).
$$

Then the error is at most $3h/8$. The preconditioned residual is homotopic to the translated identity through maps without lateral zeros, and each positive ratio slice has degree $+1$. The longitudinal separation from the first to the fourth contact is at least $1-4r^2(|a|+|c|+2\varepsilon)\ge3/4$. The zero is a genuine rectangle of the prescribed ratio, with physical circumradius at least $3L/8$ along the whole box of solutions.

Apply the connected-continuation argument of v284 for $r$ in $[r_j,r_*]$ with positive $r_j$ tending to zero. The selected compact connected sets have a subsequential Hausdorff limit in the compact parameter box. Their limit remains connected, meets both $r=0$ and $r=r_*$, and contains only limits of genuine rectangles. Even if the scaled residual has no limit at zero, the geometric contacts converge to $(A,A,B,B)$ because $z$ stays bounded. Hence the zero-ratio endpoint is the stipulated geometric flat, reached with the same radius bound.

## v297 — Finite-coil arcs at the actual flats

Use the transverse graph description $(f(t),t)$ of each smooth v292 reference arc. Choose a small source length $\ell$ on its positive side. Insert the finite-coil arc of [Asano–Ike, Proposition 6.6](https://arxiv.org/html/2609.11598v3#S6.SS2) after restoring its first finitely many disks to diameters. Its $n$-th disk is at $a_n=2^{-n}$ with radius $R_n=2^{-3n}$, and the compatible parametrization differs from the straight one by at most $2R_n$ there. After scaling by $\ell$, transfer the arc through the local shear $(x,y)\mapsto(f(x)+y,x)$ in the $(e,n)$ coordinates. The shear has sup-norm Lipschitz constant at most two on sufficiently small source intervals; the reference projection derivative vanishes at the flat. Each coil lives in its own deformed disk and joins the unchanged graph near its disk boundary. The insertion remains an embedded arc and equals the reference outside the insertion interval.

For source parameters in the $n$-th disk, the positional error from the reference divided by squared source distance is bounded by

$$
 \frac{4R_n}{\ell(a_n-R_n)^2}
 =\frac{4}{\ell}\,\frac{2^{-n}}{(1-2^{-2n})^2}.
$$

This decreases to zero. Given the tolerance from v296 and any sufficiently small $\ell$, retain only late enough coils to spend at most half the positional error budget. Shrinking the reference collar allocates the other half to its smooth Taylor remainder. The bounds apply to the limit and every finite truncation, and the anchor itself stays fixed. All four original flat contacts can receive independent insertions in disjoint neighborhoods. This construction uses the full source-based lemma, since winding inside the disks prevents representing the perturbed arcs as graphs in the fixed transverse axis.

The primitive obstruction of the imported arc survives. Removing finitely many disks leaves its unbounded excursions on the remaining disks. The shear preserves the physical area form in the stated coordinate orientation, and the similarity scales it by a positive constant. Pulling back the primitive changes this scaled primitive by a smooth exact form on the chart. Its bounded endpoint contribution on the insertion cannot remove the excursions. Thus each of the four flat contacts can be a non-Liouville-regular accumulation point.

## v298 — Matching the rough collars to the same interior

Choose a small positive ratio gate at each end before performing the insertions. On the smooth reference, its solution expressed in the scaled $z$ coordinates tends to $z_*$ as the ratio tends to zero. This follows from the v292 parametrization by short half-length: the long side tends to the flat chord and the pair midpoint displacement is quadratic in the short half-length. The rescaled reference residual converges in $C^1$ on the fixed box to $G$ by the smooth Taylor expansion. By making the gate small, it is injective on that box after preconditioning, since $\|J_0^{-1}DG-I\|_\infty\le 2CB_2h\le1/4$. Thus the gate has exactly the reference rectangle there.

Every admissible source value at a gate has magnitude at least $r_*(1-h)$. The four insertions can have smaller support than this and can also avoid the central reference contact windows. The latter windows stay a positive source distance from all four flat points on the closed angle interval between the gates. With the original gate arc windows retained, every connected collar set meets the same exact gate rectangle as the interior. On the central interval use the uniform degree tube from v294. It allows small compatible $C^0$ changes of the interior when the gate windows are kept. If the gate windows move with sufficiently small $C^1$ changes, uniqueness within the two gate charts and continuity of the tag changes give the same matching as in v294, provided the collars still obey their quadratic envelopes.

The two rough collars and the interior give a single compact connected witness of positive-radius rectangles and their two flat limits. The collars have circumradius at least $3L/8$ with $L=5,\sqrt2$ respectively, and the interior retains the v294 bound $1/3$. In the original cyclic labels the flat diagonal dots have opposite signs. The joined witness therefore carries the nonzero relative sign obstruction and gives all aspect ratios. For the four-coil fixture with the central geometry unchanged, the reference square of squared side $145/121$ is on the same witness for the limit and the truncations.

## v299 — Positions and tangents alone can lose the pinned flat

In the same normalized coordinates consider two graph arcs $P_A(t)=(0,t)$ and $P_B(t)=(1+t^3,t)$. The flat chord and both tangent directions agree with the quadratic model. If the two short sides of a rectangle lie in these graph arcs, their nonzero common transverse span may be written $2r>0$. The left shared-secant slope is zero, while on the right it is

$$
 \frac{(v+r)^3-(v-r)^3}{2r}=3v^2+r^2>0.
$$

Thus no such genuine rectangle exists on the tagged arcs. In rectangles approaching the stipulated flat, the two corners near either endpoint must be adjacent. If they were diagonally opposite, the whole diameter would tend to zero, despite the positive distance to the other two corners. Thus a sequence approaching the pinned chord would need short sides on the two given arcs, which the secant inequality rules out. The replacement can be inserted in a shrinking neighborhood of a smooth quadratic flat by a graph cutoff, matching the old arcs outside. The error is $O(t^2)$ on the changed windows. For every $0\le\beta<1$ its weighted sup norm with weight $|t|^{1+\beta}$ tends to zero with the support size. This control concerns the stipulated flat. Other flat chords may move to nearby locations.

## Replay and scope

Run `python code/rough_flat_contacts_v296_v299.py --check` for exact scaled bounds, finite-coil tolerances, and the first-order control. The SageMath tests derive the rescaled polynomial limit and its determinant independently. The artifact gives exact sample constants for dimensionless curvature $K=1$. The global proof chooses sufficiently small roundings, computes the bounds from their jets, and takes the gates inside the corresponding source domains.

The result gives full continuation at rough flat contacts under quadratic positional envelopes and compatible joins. A theorem selecting suitable endpoint arcs on arbitrary Jordan curves remains open.
