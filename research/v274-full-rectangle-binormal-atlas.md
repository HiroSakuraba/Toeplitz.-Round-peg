# v274 — Complete rectangle atlas for a concave quadrilateral

**Proved restricted.** On \(P=((0,0),(4,0),(1,1),(0,3))\), the atlas exhausts all positive-radius rectangles, with their flat limits. There is no positive-radius cutoff. The graph describes the compactification by actual corner locations with labeled flat limits. It uses the parameterized polygon boundary, with edge \(i\) traced affinely on \([i,i+1]\).

## Coverage and chart equations

The rational segment-intersection audit verifies simplicity. The four edge directions are pairwise nonparallel, and the turn at \((1,1)\) is reflex. Give a rectangle counterclockwise labels \(p_0,p_1,p_2,p_3\) and let \(\theta\) be the counterclockwise angle from \(p_0-p_2\) to \(p_1-p_3\). Set \(x=\tan(\theta/2)/(1+\tan(\theta/2))\). The true rectangles have \(0<x<1\); flats lie at \(x=0,1\). The square angle is \(x=1/2\).

Every genuine rectangle uses at least three of this polygon's edge lines. Three distinct corners are never collinear. If two edge lines each carry a side, the sides are opposite and would be parallel. If instead each line carries a diagonal, both segments would contain their common midpoint in their interiors, which contradicts simplicity. This leaves the \(24+144=168\) ordered edge assignments using four distinct edges or one doubled edge. All counterclockwise labels occur, including all possible sectors of curve order. Starting with ordered assignments automatically covers whichever diagonal pairing is geometrically valid.

The chart equations are equal midpoints and

\[
((1-x)^2+x^2)(p_1-p_3)
=(1-2x)(p_0-p_2)+2x(1-x)J(p_0-p_2),
\]

where \(J\) is rotation by \(90^\circ\). The rotation clears a strictly positive denominator. With affine edge points this gives a \(4\times4\) system over rational polynomials in \(x\). On open intervals off the determinant roots its solution consists of rational functions. Their \(0,1\) wall roots, determinant roots, squared-radius zeros, the square angle and the flat angles cut out a finite common partition. Each sign is constant between consecutive distinct roots, so evaluating at a rational point in that chamber decides feasibility there. Half-open edge ownership prevents duplicate open branches when a corner remains on a polygon vertex.

At regular event values the exact endpoint solution is checked. At an interior determinant root the verifier separately enumerates the vertices of the affine slice of \([0,1]^4\) by activating independent cube walls. This catches configurations that would be missed by limits of the rational functions. All 120 singular chart slices in this fixture are empty. The implementation refuses unresolved persistent or positive-dimensional singular charts. The flat nodes are limits of genuine charts. Arbitrary pairs of boundary points at \(\theta=0,\pi\) need not be such limits.

## Gluing and independent checks

Each finite branch has two exactly evaluated endpoints. A boundary parameter \(i+1\) is the next edge's parameter \(0\), including the wrap from \(4\) to \(0\). Endpoints are identified by their exact angle and ordered normalized boundary parameters. The checker also requires their exact corner coordinates to match. The true positive branches are injectively represented within each open angle chamber, because the angle and the unique owned edge tuple determine the solution.

Cyclic relabeling sends \(x\) to \(1-x\). The checker matches every labeled endpoint and interval under this operation and checks the interior transitions by rational-function identities. Every nondegenerate rectangle has a counterclockwise labeling, so cyclic quotienting accounts for the dihedral relabelings after fixing that orientation.

The labeled graph has 680 vertices and 672 branches, including 16 labeled flat ends and eight labeled square vertices. It has no total-collision ends. There are 85 event angle values. The exact minimum of each branch's squared circumradius is computed using endpoints and derivative roots. The independent v223 square solver agrees on the two geometric squares, with squared sides \(1\) and \(145/121\). A separate adjacent-side ratio solver agrees on whole rectangle slices at \(x=1/3,2/3\) and at a genuine irrational-angle event. That solver uses the given adjacent-side ratio \((1-x)/x\) and never uses the atlas's determinant polynomials or gluing.

The artifact stores real algebraic values by minimal polynomial and a dyadic isolating interval, checked with the independent rational Sturm counter. The parallel-edge bump continua from v271 belong to a separate control category. This atlas makes no completeness claim on parallel-edge polygons or in a smooth blow-up of the flat locus.

Run `python code/rectangle_atlas_v274.py` in the SageMath environment to replay the stored artifact, or add `--write` to regenerate it. Run `python -m unittest discover -s tests -p test_rectangle_atlas_v274.py -v` there for the independent slice and failure checks.
