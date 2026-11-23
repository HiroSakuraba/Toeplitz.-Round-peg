# v359. Two fixed arcs through refinement

The goal is a square on an arbitrary Jordan curve. If two disjoint closed parameter arcs must each receive a square corner along arbitrarily fine registered approximations to one embedding, their images on that embedding have positive distance. For image distance D and uniform error e, any such square has side at least (D−2e)/√2. Compactness gives four distinct limiting square corners on the embedding. This needs no universal inradius coefficient. The [v247 decagon](v247-true-big-peg-ratio-improvement.md) already rules out √2 times the inradius as a universal bound.

## Exact arc test for whole traced families

Give each edge a parameter interval using rational knots from 0 to 1, with 0 and 1 identified. A polygon runs linearly on each interval. Inserted edges inherit subintervals of the old times. Never repartition the circle by the new edge count. An arc starts at offset + j/N modulo one and covers L consecutive bins. We require its **open interior** to receive a corner from every square of the selected family. This stricter test lets us shrink an isolator and transfer the visits to a close redraw. We choose pairs whose closures have a gap in each direction.

On a labeled v341 cell that is a point or segment, each of the four traversal times varies affinely along the cell. Split the segment wherever a time reaches a bin boundary. At every event and each open interval midpoint, record all arcs hit by at least one of the four corners. Intersect these records over the whole component; the corner supplying a hit may change. Take all separated pairs from this intersection if the signed local degree is odd. Testing endpoints alone could accept an arc hit by one corner at the first endpoint and a different corner at the last, with a gap in between. Higher-dimensional cells need a different arrangement test.

The odd degree gives existence in any sufficiently small isolator. The condition of hitting both open arcs is open in the four parameters, so compactness of the selected component gives such an isolator entirely within the hitting condition. Across stages one may take the union of certified pairs over available components at each stage, then intersect those sets across stages. One successful pair proves a finite replay. An infinite sequence needs a theorem.

## Replays

Use 16 equal traversal bins with origin 1/128 and baseline vertices at equal traversal times (the last time is 1). Numbers are zero-based arc starts. At span two each arc covers two consecutive bins.

| Reference | Span | Fixed pair | Refinement or obstruction |
|---|---:|---|---|
| v358 exact target tie, 10 by 9, cuts [1,2] at 4 and [5,9] at 8 | 2 | 5, 9 | Trace from (2,7),(3,0),(10,1),(9,8). No separated single-bin pair works on this component at this offset. |
| v247 decagon | 1 | 3, 7 | Exact seed from the recorded maximizing edge assignment. Here the proof needs no comparison with the inradius. |
| v357 simultaneous collision | 2 | 0, 9 | At floor heights 8.9, 9 and 9.1, use a disappearing tilted square at the first stage and the fifth square at the other two. |
| v357 surviving square | 1 | 4, 15 | Holds on the reference and with two and four folded tiles on each of the four contact arcs. The disappearing component at the collision has zero degree. |

On the tie and decagon, tests replace the midpoint of an edge outside the entire contact trace by a rational kink at two amplitudes. These are controls for registration and local edits. The four-contact test changes the marked arcs themselves, including the bottom and upper walls. Each tile backtracks tangentially between two parallel steps. The next argument covers every refinement in a continuous class, even when the original square moves off its old vertices.

## Four folded contact arcs, with arbitrary connectors

Use the 15 by 12 reference of [v357](v357-simultaneous-loss-of-two-square-pairs.md). Around its successive contacts (6,0),(11,6),(5,11),(0,5), mark the half-unit of straight edge on either side of each point. Retain those four parameter windows in the same order. Replace each by any continuous path inside the strip of normal half-width e=1/32 about its supporting line, whose tangential coordinates at its two endpoints still equal the original endpoint coordinates. The rest of the boundary can change arbitrarily provided the result is a Jordan curve with this registration. The paths may reverse tangential direction arbitrarily often. This is a sufficient geometric hypothesis on the four arcs. No entire bottom edge is required. It specializes the traversal idea of [v317](v317-neck-lens.md).

**Claim.** Every such curve has a square visiting arcs 4 and 15 with side at least (11−2e)/√2. This remains true at every stage of any sequence preserving these four windows. The stated bound is greater than 7.73 in the reference length units.

Parametrize the four trial points in the product of their windows. Project to the supporting lines, obtaining (x,0),(11,y),(5,z),(0,w). Interpolate the four marked paths to the reference straight traversals at matching parameters. The interpolants need not be embeddings; their normal error remains at most e, and their tangential endpoint coordinates stay fixed. In the residual for completing the adjacent first and fourth vertices, the four projected entries are 11−x−w, y−x, 5−w, z−x−w. A zero of the actual residual therefore bounds each of these by 4e in absolute value. Substitution gives errors in w,x,y,z of at most 4e,8e,12e,8e from 5,6,6,11. None reaches 1/2. On a parameter-box face at least one tangential coordinate is exactly 1/2 from its reference value, so the boundary is zero-free throughout. The straight residual has a single nonsingular zero, hence the terminal residual has odd degree. The windows on x=11 and x=0 lie strictly within bins 4 and 15 respectively. The two actual contacts there are opposite square corners and differ in abscissa by at least 11−2e. This proves the claim.

## A cutoff-crossing alternative for the general proof

Fix disjoint closed parameter arcs R and B for an embedding. On a compact four-parameter box whose first and third parameters range in R and B, use the ordinary square-completion residual. If its degree changes between two boundary-zero-free stages of a continuous registered deformation, there must be a square on the box boundary at some intermediate stage. Even at that exit the two corners are still in the **closed** colored arcs. If the deformation stays within e of a fixed target, the boundary square also has side at least (D−2e)/√2. Thus, on an infinite sequence of refinements uniformly approaching the target, infinitely many such exits already suffice by compactness. If there are only finitely many exits, one needs an odd box after the last one. A last exit can remove the only odd box; that is the remaining obstruction. This is the fixed-arc counterpart of [v310](v310-moving-degree-and-scale-loss.md) and [v315](v315-cutoff-payment.md).

A general proof still needs to **select** separated arcs from an arbitrary unmarked loop and prevent a final loss with no replacement. The strip theorem assumes its windows. A next test is to erase one of the four strip constraints while preserving just the two colors, then trace all families reaching those colors through the first exchange. A failed pair is useful data; changing the colors at every level does not prove noncollapse.

[Arc tracer](../code/separated_arcs_v359.py) · [Exact replays](../tests/test_separated_arcs_v359.py)
