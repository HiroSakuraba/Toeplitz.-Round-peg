"""Real-valued actions of square generators from the preferred capping (v331, third step).

Greene-Lobb (arXiv:2404.05179, Lemma 2.5) give each generator x = R_theta(p) a preferred capping
that avoids the diagonal Delta(C). Under (z, w) -> (midpoint, w - z) the complement of Delta(C)
is C x C*, so a loop bounds a disk avoiding Delta(C) exactly when w - z has winding zero along it.
The Hamiltonian trajectory from p to x turns w - z by +theta, so the preferred return path is the
unique pair-path in gamma x gamma minus its diagonal along which arg(w - z) changes by -theta, and

    A(x) = integral over the pair-path p -> x of (lambda_z + lambda_w),
           with arg(w - z) changing by exactly +theta along it,

where lambda = (x dy - y dx)/2 (the sum is origin-independent). The Hamiltonian term cancels the
trajectory's own lambda-integral as in (323.1). On an ellipse this gives 2ab arctan((a/b) tan(theta/2))
on the major branch, whose theta-derivative is q1^2 (the derivative property) and which runs from 0 to
the area; the replay checks the polygon version against convex examples and against duality.

Polygons are parametrized by t in [0, n): t = m + u is the point P[m] + u (P[m+1] - P[m]). Along a
straight path in the (t_z, t_w) parameter plane, w - z is piecewise linear, so its total change of
argument is a finite sum of signed angles.
"""
from fractions import Fraction as Q
import math

from angular_windows_v328_v330 import area, primitive


class PolygonParam:
    def __init__(self, P):
        self.P = list(P)
        self.n = len(P)
        self.F = primitive(self.P+[self.P[0]])           # primitive at vertices, origin (0, 0)
        self.total = self.F[-1]                          # signed loop integral (= +-A)

    def point(self, t):
        n = self.n
        k = math.floor(t)
        m, u = k % n, t-k
        a, b = self.P[m], self.P[(m+1) % n]
        return (a[0]+u*(b[0]-a[0]), a[1]+u*(b[1]-a[1]))

    def F_ext(self, t):
        """Primitive along the lifted parameter, F(t + n) = F(t) + total."""
        n = self.n
        k = math.floor(t)
        m = k % n
        v = self.point(t)
        a = self.P[m]
        return self.F[m]+(a[0]*v[1]-a[1]*v[0])/2+(k//n)*self.total

    def param_of(self, v):
        """Parameter of a point on the polygon (first match; vertices give integer parameters)."""
        n = self.n
        for m in range(n):
            a, b = self.P[m], self.P[(m+1) % n]
            if v == a:
                return Q(m)
            d = (b[0]-a[0], b[1]-a[1])
            if d[0]*(v[1]-a[1])-d[1]*(v[0]-a[0]) == 0:
                s = ((v[0]-a[0])*d[0]+(v[1]-a[1])*d[1])
                dd = d[0]**2+d[1]**2
                if 0 < s < dd:
                    return m+Q(s, 1)/dd
        raise ValueError("point not on polygon")


def _arg_change(par, a0, b0, a1, b1):
    """Total change of arg(gamma(b) - gamma(a)) along the straight parameter path."""
    cuts = {Q(0), Q(1)}
    for lo, hi in ((a0, a1), (b0, b1)):
        if lo != hi:
            for k in range(math.floor(min(lo, hi))+1, math.ceil(max(lo, hi))):
                cuts.add((k-lo)/(hi-lo))
    total = 0.0
    prev = None
    for s in sorted(cuts):
        z, w = par.point(a0+s*(a1-a0)), par.point(b0+s*(b1-b0))
        v = (float(w[0]-z[0]), float(w[1]-z[1]))
        if v == (0.0, 0.0):
            raise ArithmeticError("pair path meets the diagonal")
        if prev is not None:
            total += math.atan2(prev[0]*v[1]-prev[1]*v[0], prev[0]*v[0]+prev[1]*v[1])
        prev = v
    return total


def real_action(par, p, x, theta=math.pi/2):
    """Real action of the generator x = R_theta(p); p and x are pairs of points on the polygon."""
    n = par.n
    a0, b0 = par.param_of(p[0]), par.param_of(p[1])
    if b0 <= a0:
        b0 += n
    a1, b1 = par.param_of(x[0]), par.param_of(x[1])
    # lift x into the same strip a < b < a + n of the parameter plane (the strip avoids the diagonal)
    while b1 <= a1:
        b1 += n
    while b1-a1 >= n:
        b1 -= n
    d0 = _arg_change(par, a0, b0, a1, b1)
    d1 = _arg_change(par, a0, b0, a1+n, b1+n)
    step = d1-d0
    if abs(abs(step)-2*math.pi) > 1e-6:
        raise ArithmeticError("core loop does not wind once around the diagonal", step)
    k = (theta-d0)/step
    if abs(k-round(k)) > 1e-6:
        raise ArithmeticError("no lift turns w - z by exactly theta", d0, step)
    k = round(k)
    a1, b1 = a1+k*n, b1+k*n
    return (par.F_ext(a1)-par.F_ext(a0))+(par.F_ext(b1)-par.F_ext(b0))


def square_actions(par, square):
    """Real actions of the generators x = (c1, c3) = R(c0, c2) and x' = (c2, c0) = R(c1, c3)
    of a counterclockwise square (c0, c1, c2, c3) at theta = pi/2."""
    c0, c1, c2, c3 = square
    return real_action(par, (c0, c2), (c1, c3)), real_action(par, (c1, c3), (c2, c0))


def ellipse_major_action(a, b, theta):
    """Closed form on the ellipse: 2ab arctan((a/b) tan(theta/2)) (major-axis branch)."""
    return 2*a*b*math.atan((a/b)*math.tan(theta/2))


# ---------------------------------------------------------------- fast path for large censuses

class FastParam(PolygonParam):
    """Same results as PolygonParam, with vectorized point location and argument sums.

    The argument change only selects an integer (it is checked to be within 1e-6 of one), so it
    is computed in floating point; the action itself stays exact."""
    def __init__(self, P):
        import numpy as np
        super().__init__(P)
        self.np = np
        self.V = np.array([[float(x), float(y)] for x, y in self.P])
        self.D = np.roll(self.V, -1, axis=0)-self.V
        self.vertex_index = {v: m for m, v in enumerate(self.P)}
        self.step = None

    def param_of(self, v):
        if v in self.vertex_index:
            return Q(self.vertex_index[v])
        np = self.np
        x, y = float(v[0]), float(v[1])
        cr = self.D[:, 0]*(y-self.V[:, 1])-self.D[:, 1]*(x-self.V[:, 0])
        dot = (x-self.V[:, 0])*self.D[:, 0]+(y-self.V[:, 1])*self.D[:, 1]
        dd = np.sum(self.D*self.D, axis=1)
        for m in np.nonzero((np.abs(cr) < 1e-9*(1+np.sqrt(dd))) & (dot > -1e-12) & (dot < dd+1e-12))[0]:
            m = int(m)
            a, b = self.P[m], self.P[(m+1) % self.n]
            d = (b[0]-a[0], b[1]-a[1])
            if d[0]*(v[1]-a[1])-d[1]*(v[0]-a[0]) == 0:
                s = (v[0]-a[0])*d[0]+(v[1]-a[1])*d[1]
                if 0 < s < d[0]**2+d[1]**2:
                    return m+Q(s, 1)/(d[0]**2+d[1]**2)
        raise ValueError("point not on polygon")

    def points(self, t):
        np = self.np
        k = np.floor(t).astype(int)
        u = t-k
        m = k % self.n
        return self.V[m]+u[:, None]*self.D[m]

    def arg_change(self, a0, b0, a1, b1):
        np = self.np
        a0f, b0f, a1f, b1f = map(float, (a0, b0, a1, b1))
        cuts = [np.array([0.0, 1.0])]
        for lo, hi in ((a0f, a1f), (b0f, b1f)):
            if lo != hi:
                ks = np.arange(math.floor(min(lo, hi))+1, math.ceil(max(lo, hi)), dtype=float)
                cuts.append((ks-lo)/(hi-lo))
        s = np.unique(np.concatenate(cuts))
        v = self.points(b0f+s*(b1f-b0f))-self.points(a0f+s*(a1f-a0f))
        if np.any(np.all(v == 0, axis=1)):
            raise ArithmeticError("pair path meets the diagonal")
        return float(np.sum(np.arctan2(v[:-1, 0]*v[1:, 1]-v[:-1, 1]*v[1:, 0], np.sum(v[:-1]*v[1:], axis=1))))

    def real_action(self, p, x, theta=math.pi/2):
        n = self.n
        a0, b0 = self.param_of(p[0]), self.param_of(p[1])
        if b0 <= a0:
            b0 += n
        a1, b1 = self.param_of(x[0]), self.param_of(x[1])
        while b1 <= a1:
            b1 += n
        while b1-a1 >= n:
            b1 -= n
        d0 = self.arg_change(a0, b0, a1, b1)
        if self.step is None:
            self.step = self.arg_change(a0, b0, a1+n, b1+n)-d0
            if abs(abs(self.step)-2*math.pi) > 1e-6:
                raise ArithmeticError("core loop does not wind once around the diagonal", self.step)
        k = (theta-d0)/self.step
        if abs(k-round(k)) > 1e-6:
            raise ArithmeticError("no lift turns w - z by exactly theta", d0)
        k = round(k)
        return (self.F_ext(a1+k*n)-self.F_ext(a0))+(self.F_ext(b1+k*n)-self.F_ext(b0))

    def square_actions(self, square):
        c0, c1, c2, c3 = square
        return self.real_action((c0, c2), (c1, c3)), self.real_action((c1, c3), (c2, c0))


def carrier_filter(N, workers=4, fast=True):
    """Residue survivors -> generators with real action in [A/2, pi R^2/2] -> even parity."""
    from residue_census_v331 import census, midside_fixture
    from floer_parity_v331 import calibration, square_generators
    from angular_windows_v328_v330 import pi_enclosure, radius_squared
    even = calibration()
    row = census(N, workers=workers)
    P = midside_fixture(N)
    par = FastParam(P) if fast else PolygonParam(P)
    A, F = area(P), primitive(P+[P[0]])
    top = pi_enclosure(30)[1]*radius_squared(P)/2
    tally = {"survivors": len(row["survivors"]), "no_real_action_in_window": 0, "odd_only": 0,
             "vertex_contact": 0, "candidates": 0, "duality_failures": 0}
    shifts = {}
    for _, sq in row["survivors"]:
        gens = square_generators(P, sq, F=F, A=A)
        acts = par.square_actions(sq) if fast else square_actions(par, sq)
        if sum(acts) != A:
            tally["duality_failures"] += 1
        for (r, _), act in zip(gens, acts):
            k = int((act-r)/A)
            shifts[k] = shifts.get(k, 0)+1
        window = [s for (r, s), act in zip(gens, acts) if A/2 <= act <= top]
        if not window:
            tally["no_real_action_in_window"] += 1
        elif even in window:
            tally["candidates"] += 1
        elif all(s is None for s in window):
            tally["vertex_contact"] += 1
        else:
            tally["odd_only"] += 1
    return {"N": N, **tally, "max_capping_turns": max((abs(k) for k in shifts), default=0)}


if __name__ == "__main__":
    import json
    import sys
    turns = (30,) if sys.argv[1:] == ["--large"] else (20, 24)
    print(json.dumps({"carrier_filter": [carrier_filter(N) for N in turns]}, indent=2))
