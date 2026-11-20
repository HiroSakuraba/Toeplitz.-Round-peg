# v292–v295 — Moving flat endpoints on a rough Jordan curve

A curved rounding of the concave reference path gives two isolated flat chords and a regular fixed-angle path between them. There is a neighborhood of this rounded Jordan curve in which every curve has a connected rectangle witness spanning both flats, with circumradius at least $1/3$. The neighborhood requires $C^1$ control on four small endpoint arcs and compatible $C^0$ control on the rest of the curve. With $C^2$ control at the endpoint arcs, the two local branches are unique graphs over the aspect ratio. The neighborhood includes curves with two finite-coil singularities.

## v292 — A rounding with nondegenerate flat chords and regular transitions

Use the persistent concave path of v288 on the quadrilateral $((0,0),(4,0),(1,1),(0,3))$. Near a flat write its four corners as

$$
 A(s)-s n(s),\quad A(s)+s n(s),\quad B(s)+s n(s),\quad B(s)-s n(s).
$$

Here $s$ is the short half-side length, $A,B$ are the two pair midpoints, and $n$ is the unit short-side direction. The long vector is $B-A=L e$, with $L>0$ and $J n=e$ for the positive quarter turn $J$. A cyclic relabeling near the second flat places its shrinking side first. The exact short derivative of v288 is nonzero at both flats, so $s$ is a smooth one-sided parameter.

Take a fixed smooth nondecreasing cutoff $\chi$ equal to zero on $[0,1/3]$ and one from $2/3$ onward. For $0\le t\le1$ put $\psi(t)=t+(1-t)\chi(t)$ and extend it by one for $t\ge1$. On each sufficiently small collar use

$$
 \phi(s)=s\psi(s/\zeta).
$$

Replace the pair midpoints and their common normal by $A(\phi(s)),B(\phi(s)),n(\phi(s))$, and retain the short half-length $s$. Thus $\phi=s^2/\zeta$ near zero and $\phi=s$ at the join. In $0<s\le\zeta$ we have $\phi\le s$, $k=\phi'>0$, $k$ bounded independently of $\zeta$, and $s/k\le\zeta$. The two halves at each flat form a smooth common arc, since the signed transverse parametrization is $A(\phi(|s|))+s n(\phi(|s|))$ and $\phi(|s|)$ is even near zero. The v286 embedding check applies: $\phi\le sk$ keeps each half monotone in the inward longitudinal projection for small $\zeta$, and the two halves lie on opposite sides of the local midpoint graph. The roundings fit in disjoint neighborhoods of the original vertices and coincide with the old traces past their joins. The circumradius cannot decrease from the reference at $\phi(s)$, since the long side is the same and $s\ge\phi(s)$.

Write $a=A'(0)\cdot e(0)$ and $c=B'(0)\cdot e(0)$, where the derivatives are in the old short half-length. The exact endpoint data from v288 give

| Flat angle in the reference labels | $L$ | $a$ | $c$ |
|---|---:|---:|---:|
| $0$ | $5$ | $24/5$ | $-12/5$ |
| $1$ | $\sqrt2$ | $1$ | $2/5$ |

For example at the first flat the inward midpoint projections are $60,30$ on opposite chords and the old short-half-length derivative is $5/2$. The analogous projections at the second flat are $1,2/5$ on the same chord and the old short-half-length derivative is $\sqrt2/2$.

Write $\kappa=1/\zeta$. At each new flat the paired curves have common unit tangent $n(0)$ and longitudinal second derivatives $2\kappa a,2\kappa c$ in the signed transverse coordinate. The Hessian of half the squared chord length is

$$
 \begin{pmatrix}1-2\kappa La&-1\\ -1&1+2\kappa Lc\end{pmatrix}.
$$

Its determinants at the two flats are respectively $72\kappa(16\kappa-1)$ and $-(6\sqrt2/5)\kappa-(16/5)\kappa^2$. Both are nonzero for every $\kappa\ge1$. In particular the chords are isolated binormals. A binormal is a chord perpendicular to the tangent at each endpoint, equivalently an off-diagonal critical point of squared chord length. See [Greene–Lobb, Appendix A](https://arxiv.org/html/2604.17116v1#A1).

### Regularity at every strictly positive angle of the rounded path

The rounding also preserves regularity of the *selected fixed-angle solve*, including the transition to the original polygon edges. At a rounded rectangle with short half-length $s>0$, use the axes $e,n$ at the old parameter $\phi(s)$. Write $A'=(a,b)$, $B'=(c,d)$, and $n'=\omega e$ there. The contact velocities along each of the four traces, with the local nonnegative short half-length as parameter, are $kA'\mp T$ and $kB'\pm T$, where $T=n+sk\omega e$.

Let $m_A,d_A,m_B,d_B$ denote the means and half-differences of source variations within the pairs. Use the two vector equations saying the short half-vectors match and $s\,\Delta(B-A)=LJ\,\Delta U$, with $\Delta U$ the variation of the first short half-vector. These equations are equivalent to the differential of the fixed-angle rectangle residual after a row operation using short-vector matching. Divide the longitudinal short-matching row and the transverse similarity row by $k$. The resulting matrix is

$$
 \begin{pmatrix}
 s\omega&a&-s\omega&-c\\
 1&kb&-1&-kd\\
 -ska-L&-s^2k\omega-Lkb&skc&s^2k\omega\\
 -sb+Ls\omega&-s/k+La&sd&s/k
 \end{pmatrix}.
$$

As $\zeta$ tends to zero the entries converge uniformly on its shrinking collar to those obtained by setting $s=s/k=0$ and using the flat values for the geometry, with $k$ left in a fixed bounded interval. The determinant of this limiting matrix is $L^2ac$, independent of $k$, with values $-288$ and $4/5$. Thus every sufficiently small rounding has nonsingular fixed-angle equations throughout each punctured collar. Outside the collars the old affine determinant is $12(1+3x+x^2)>0$. Moreover $2s/(L(\phi(s))+2s)$ is strictly increasing into the path on both sufficiently small collars, by boundedness of $k$ and smoothness of $L$. The entire selected path therefore has a regular parametrization by the original angle coordinate on $0<x<1$. The orientation change at the second end sends the local small coordinate to $1-x$.

## v293 — Local continuation at a moving binormal

In the fixed long/transverse axes of a reference flat, represent the two near-end arcs by graphs $(f(u),u)$ and $(g(v),v)$. Small $C^1$ perturbations of the tagged arcs remain graphs on a common slightly smaller interval because the transverse projections have nonzero derivative. For projected short half-length $r$ define the chord midpoint and divided-difference operators

$$
 M_r f(u)=\tfrac12(f(u+r)+f(u-r)),\qquad
 D_r f(u)=\tfrac12\int_{-1}^1 f'(u+rt)\,dt.
$$

Let $D=M_r g(v)-M_r f(u)$, which stays positive on the collar, and let $\lambda$ be the short-to-long side ratio. The following three equations use only these graphs:

$$
 D_r f(u)-D_r g(v)=0,\qquad
 D\,D_r f(u)+v-u=0,\qquad
 2r-\lambda D=0.
$$

The first equation gives the matching short vectors. The second gives perpendicularity of the common short vector and the long side from one pair midpoint to the other. In these coordinates the short length is $2r\sqrt{1+(D_r f)^2}$ and the long length is $D\sqrt{1+(D_r f)^2}$, so the third imposes their ratio. At $\lambda=0$ we have $r=0$; the first two equations are exactly the two common-normal conditions. For $\lambda>0$, $r>0$ and the solution is a genuine rectangle, with the same ordered short-side convention throughout.

At the reference flat the $u,v$ Jacobian of the first two equations, with $r=0$, is

$$
 J_0=\begin{pmatrix}f''&-g''\\ Lf''-1&1\end{pmatrix}.
$$

The determinant of the chord-length Hessian above is $-L\det J_0$. The third row at $\lambda=0$ gives derivative $2$ in $r$, with zero center derivatives. Thus the full angular collar equation is nonsingular in $(u,v,r)$ at the reference flat. With $C^2$ control the integral formula for $D_r$ controls the first derivatives of the collar equations as well. The implicit function theorem, or preconditioned contraction on a sufficiently short tube, gives a unique local branch over $\lambda$ on each perturbed pair of arcs. Its endpoint is a moving binormal. Both perturbed ends preserve their reference flat type in the fixed labels.

For $C^1$ control the equations are uniformly close as continuous maps, including at $r=0$. At the smooth reference curve, enclose the short branch in a fixed-width box centered along its regular solution at each small $\lambda$. Precondition by the reference Jacobian and take the box small enough that the straight homotopy to the translated identity is zero-free on its sides. The uniform $C^1$ control on the arcs preserves this boundary estimate and the slice degree $+1$. Apply the connected-continuation argument of v284 on intervals $[\lambda_j,\lambda_{\rm gate}]$ with $\lambda_j>0$ tending to zero. The selected compact connected sets have a subsequential Hausdorff limit in the compact collar cylinder. That limit is connected, meets the gate and $\lambda=0$, and lies in the closure of genuine rectangles because every set in the sequence uses positive ratios. Its zero-ratio points are geometric flats and satisfy the binormal equations. With this weaker control there may be several solutions at a given collar parameter. The matching uses uniqueness at the positive gate from the next step.

## v294 — A mixed neighborhood with moving joins and a uniform radius

Choose two fixed interior angles $x_-,x_+$ deep inside the smooth reference collars, one near each end. Cut the reference path there. On the compact central interval $[x_-,x_+]$ the reference path is regular at every fixed angle by v292. Use compatible source parameters and a small box around its parameter vector $u(x)$. Let $A_x$ be the contact Jacobian at the reference zero, in these source coordinates.

By compactness, the box can have positive uniform halfwidth $h$. On it the derivative of the preconditioned reference residual differs from the identity by less than $1/4$ in the infinity norm. The perturbed central residual differs uniformly from the reference residual under small $C^0$ displacement of the tagged parametrization. Choose this displacement so the preconditioned error is less than $h/4$. The straight homotopy again has a zero-free lateral boundary and degree $+1$. The same compact continuation argument gives a connected set meeting both gate fibers.

At the two fixed gate angles all four contact windows can be chosen inside the endpoint arcs where $C^1$ control is assumed. The derivative of the preconditioned perturbed fixed-angle residual differs there from the identity by less than $1/2$ when that control is small enough. Integration along segments in each convex contact box gives injectivity, so each gate fiber has a unique zero. In the local collar construction of v293 take the boxes and arc perturbations small enough that the entire local gate fiber maps into this same fixed-angle box in ordered source coordinates. This is possible because the reference gate rectangle is common to both charts and the source changes to transverse graph coordinates are continuous under $C^1$ control. At each gate the local and central connected sets therefore share the same actual rectangle of the perturbed curve. Joining the three sets gives one compact connected witness.

The reference radius is at least $1/2$ along the whole path. Shrink the boxes and perturbation tolerances so the central opposite-corner distance is at least $2/3$. Near each flat the long side remains above $2/3$ because the reference chords have lengths $5,\sqrt2$. Thus the resulting witness has circumradius at least $1/3$ throughout, uniformly on a mixed-open neighborhood. The perturbed flats have opposite diagonal-dot signs in the continued labels. The witness sees every normalized angle from zero to one and consequently every prescribed aspect ratio. Comparing the endpoint signs also gives a nonzero restriction of the relative sign obstruction. All joins in this proof can move with the perturbation.

The neighborhood requires compatible $C^0$ control around the whole reference parametrization and $C^1$ control on fixed small source neighborhoods of the four flat contacts. One may impose the endpoint condition directly on the transverse graph functions by choosing compatible source coordinates at the gates. The proof supplies a positive tolerance by compactness and local invertibility. It leaves its numerical size unspecified. Requiring $C^2$ on the endpoint arcs gives the unique collar branches in v293 with the same central argument.

To obtain rough examples inside this neighborhood, choose the parametrization to keep the four square windows affine. Insert two arbitrarily small scaled finite-coil arcs there as in v290, leaving the endpoint neighborhoods alone. The tagged error divided by distance from the square parameter is bounded by $8/255$, independently of the insertion scale, as in v290. The affine preconditioner at the square still has norm at most two. The inequality $\|u-u_*\|_\infty\le (64/255)\|u-u_*\|_\infty$ gives the original square as the unique square in the central tube for these anchored copies. The sup error tends to zero with the size of the two insertion windows. This gives two-rough-contact curves inside the mixed-open neighborhood and covers further small $C^0$ changes of their interior, including moving contact positions at the square. Their local non-Liouville obstruction comes from the unchanged tail of the Asano–Ike construction. [Asano–Ike, Proposition 6.6](https://arxiv.org/html/2609.11598v3#S6.SS2)

## v295 — The straight tilt and the positive-radius time fold

In the earlier straight rounding take graph pieces $f(t)=0$, $g(t)=L$, with $L>0$. Their local binormals form an interval and the chord Hessian is singular. Tilt the second piece to $g(t)=L+\delta t$ with arbitrarily small nonzero $\delta$. The first divided-difference residual is now the constant $-\delta$, both at zero and at positive short length. Hence the entire selected pair of straight arc windows contains no flat and no rectangle with both short sides supported on those windows. A smooth cutoff outside the chosen windows can realize this as a small perturbation of the rest of the Jordan curve. This control concerns the tagged endpoints.

The v281 fold has the separate transverse model $|s|-t=0$. At positive time its two local degrees are $-1,+1$; the zero-time root has degree zero. Its positive radius and its nonzero relative angle classes coexist with zero degree in the time direction. This identifies the missing continuation hypothesis in that control.

## Replay and scope

Run `python code/moving_flat_endpoints_v292_v295.py --check` for the exact endpoint and fold records. The SageMath tests derive the normalized contact matrix directly from four corner variations, check its determinant at the limiting transition, and verify the chord Hessian by symbolic differentiation. The standard-library tests replay the records and check displaced parabolic collars and the tilted straight model.

The mixed-open result is obtained from the explicit rounded concave witness. A method for finding the required regular path and endpoint arcs on arbitrary Jordan curves remains open.
