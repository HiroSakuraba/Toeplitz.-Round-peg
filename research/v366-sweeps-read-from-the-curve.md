# v366. Sweeps read from the curve

The input is only a polygon's vertices and their registration. Slide a point $P$ along one stretch and a point $Q$ along another, and complete the segment $PQ$ to a square $P,Q,R,D$. Keep the configurations in which $R$ lies on a third stretch. These configurations form a piecewise linear locus. Along it, the fourth corner $D$ moves affinely on each piece and changes sides of the boundary exactly at inscribed squares. A count at the ends of the locus certifies a square that survives small registered redraws. This is a direct form of the [v360](v360-three-windows-free-fourth-contact.md) certificate: no strips, frames or probes are chosen by hand.

**Results.**
1. **Certificate.** An odd number of locus ends with $D$ inside, and no locus end with $D$ on the boundary, gives a square with $P,Q,R$ in the three stretches on every sufficiently close registered redraw. All tests are exact.
2. **The v364 roof without marks.** Searching all stretch triples on the bare sloping-roof hexagon finds the v364 square. With 32 bins it finds v364's own labeling: two corners on the bottom, one on the flat roof, and the free corner on the slope.
3. **The v247 decagon.** The search returns 20 certificates. The best puts its corners on bins 7, 3 and 10, has side floor above 1.400, and certifies the decagon's square of side about 1.410.
4. **Circle control.** Exactly quarter-turn-symmetric polygons give no certificate at all, for every size tested. That is the polygonal form of [v361](v361-blind-three-window-selection.md)'s obstruction for the round circle. One bump of relative size $b$ on the symmetric octagon creates certificates whose clearance is proportional to $b$. From 12 vertices up, one bump creates none.

## The certificate

Fix closed parameter arcs $A,B,E$ with pairwise disjoint closures and an orientation $\sigma=\pm1$. Write $J$ for the counterclockwise quarter turn and $\gamma$ for the registered boundary. On the box $\mathcal B=A\times B\times E$ put
$$
F(a,b,e)=\Big(\gamma(b)+\sigma J\big(\gamma(b)-\gamma(a)\big)-\gamma(e),\ \ \rho\big(\gamma(a)+\sigma J(\gamma(b)-\gamma(a))\big)\Big),
$$
where $\rho$ is the signed distance to the **whole** boundary, negative inside. A zero is an inscribed square with $P=\gamma(a)$, $Q=\gamma(b)$, $R=\gamma(e)$. The fourth corner may be anywhere on the boundary.

Let $L$ be the zero set of the first two coordinates. For a polygon it is a union of segments, one for each triple of edge pieces with a nonsingular $2\times3$ system. Its endpoints on the outer faces of $\mathcal B$ are where one of $P,Q,R$ reaches the end of its stretch.

**Theorem 366.1.** Suppose $D$ is off the boundary at every outer endpoint of $L$, no piece of $L$ runs inside an outer face, and no chart meeting $\mathcal B$ is rank-deficient. Then $F$ has no zero on $\partial\mathcal B$. Its mod-two degree equals the parity of the outer endpoints of $L$ with $D$ inside, counted once per piece end. If this parity is odd, every sufficiently close registered redraw has a square with $P,Q,R$ in $A,B,E$. Its side is at least the largest of $d(A,B)$, $d(B,E)$ and $d(A,E)/\sqrt2$, where $d$ is the Euclidean distance between stretch images.

*Proof.* A zero on $\partial\mathcal B$ would be a point of $L$ on an outer face with $D$ on the boundary, and the hypotheses exclude it. Perturb so that $L$ is a compact one-manifold meeting $\partial\mathcal B$ transversally at its ends. The zeros of $F$ are the points of $L$ where $\rho(D)$ changes sign. Their number on an arc of $L$ is odd exactly when its two ends have opposite signs. Closed components contribute an even number. Summing over components gives the parity of inside ends.

Where $L$ touches an outer face and returns, the two piece ends there cancel. A chart met at a single point contributes twice. Uniformly small registered redraws change $F$ uniformly, so the degree persists while $\partial\mathcal B$ stays zero-free. The side bound holds because $P,Q$ and $Q,R$ are adjacent and $P,R$ are opposite. ∎

The checker also counts independently the boundary crossings of $D$ along every piece of $L$. Each one is an exact square on the polygon. It rejects any triple where this count's parity disagrees with the end count. Collinear overlaps, rank-deficient charts and pieces lying in an outer face all make the triple inconclusive. The distance from $D$ to the boundary at the outer ends is reported as the *clearance*. It sets the scale of redraws the certificate can tolerate, but it is not a proved tolerance. [v368](v368-certificate-margin-and-one-stage-criterion.md) defines the margin of $F$ over all of $\partial\mathcal B$, proves the tolerance and computes it exactly.

## Search

Stretches are whole bins of the registered parameter: 16 equal bins with origin 0, spans of one or two bins, at least one empty bin between consecutive stretches, and both orientations. Every triple is evaluated, and results are sorted by the side floor.

**Sloping roof.** The hexagon is v364's target with equal times on its six edges. The best certificate has stretches 1, 7 and 11 (bottom, flat roof, slope) and side floor $(12/7)^2$. Its unique square is $[23/28,23/28+12/7]\times[0,12/7]$, the square that v364 forced in the limit. In this labeling the free corner is on the bottom.

With 32 bins, the triple of bins 1, 3 and 14 (two on the bottom, one on the flat roof) has degree one. Its square is the same, with the free corner $(23/28,12/7)$ on the slope, which is exactly v364's marking. The search reads the sweep off the drawing: as $P$ slides along the bottom, $D$ runs along $y=12/7$ and crosses the slope once.

**Decagon.** With equal times on its ten edges, the search returns 20 certificates. The best uses bins $(7,3,10)$ with span one and clockwise completion. Its side floor exceeds $1.4$, and its single square has side about $1.4102$. The tests also move every vertex by $\pm1/4096$ and keep the registration. The same triple still has degree one.

## The circle control

The polygon below has $4m$ rational points on the unit circle, invariant under the exact quarter turn $(x,y)\mapsto(-y,x)$. Every boundary point $p$ then gives the centred square $p,Jp,-p,-Jp$, so the centred squares form a continuum. Any box whose interior contains one of them is crossed by the whole family, which forces zeros on $\partial\mathcal B$. The search finds no certificate at all for $m=2,3,4,6,8$.

Push one vertex of the octagon ($m=2$) radially outward by a factor $1+b$. Eight certificates appear, all with side floor about $1.207$. Their best clearance is about $0.46\,b$ for $b=1/64,1/256,1/1024$, so it is linear in $b$. The tests check that quadrupling $b$ quadruples the clearance to within five percent. For $m\ge3$ a single bump leaves enough of the centred family to block every triple.

So on near-circles the certificate's own margin measures the distance to symmetry. On a refinement sequence converging to the circle, any certificates must have clearance tending to zero. The method loses control there, as v361 predicted, and it loses it in a measurable way.

## What this does and does not give

The search turns an unmarked drawing into finitely many exact, robust certificates with three named stretches. That is the input v364 and [v359](v359-fixed-arcs-through-refinement.md) need. It does not decide whether one triple of stretches works at every stage of a refinement. [v367](v367-an-adversary-against-the-sweep-certificates.md) tests that against an adversary.

[Code](../code/sweep_search_v366.py) · [Tests](../tests/test_sweep_search_v366.py)
