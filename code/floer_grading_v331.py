"""Integer Floer grading of square generators from the preferred capping (v331, fourth step).

The grading of a generator x = R_theta(p) is minus the Maslov index of the loop of Lagrangian
planes around its preferred capping (Greene-Lobb, arXiv:2404.05179, Lemma 2.5). Under the phase map
det^2 on the Lagrangian Grassmannian of C^2 the loop splits into three explicit pieces:

* along the Hamiltonian flow from p to x the plane is R_t(T_p L0); R_t is unitary with det = e^{it},
  so the phase rises by 2 theta;
* at x the canonical path from T_x L1 to T_x L0 rotates back through the Kahler angles
  theta_k in (0, pi) of the transverse pair, contributing -2 sum theta_k;
* along the preferred return path x -> p in L0 the plane is the product of the tangent lines at z
  and w, so the phase changes by twice their total turning (signed exterior angles on a polygon).

The total is checked to be an exact multiple of 2 pi. Calibration needs no constant: on every
ellipse tested the major-axis branch gets grading 2 and the minor-axis branch grading 1. Checks:
the parity reproduces the calibrated v331 signs, and the two generators of each square have gradings
summing to 3 (duality at theta = pi/2). A contact at a polygon vertex has no tangent line; under
smoothing its tangent is the incoming edge rotated by lam * (exterior angle), lam in (0, 1), and
possible_gradings collects the gradings over a grid of lam.
"""
import cmath
import itertools
import math

import numpy as np

from floer_actions_v331 import FastParam


def kahler_phase(A, B):
    """-2 sum theta_k for the canonical path from the Lagrangian B(R^2) to A(R^2)."""
    W = np.linalg.solve(A, B)
    angles = np.mod(np.angle(np.linalg.eigvals(W @ W.T)), 2*math.pi)/2
    if np.any(angles < 1e-9) or np.any(angles > math.pi-1e-9):
        raise ValueError("planes are not transverse")
    return -2*float(np.sum(angles))


def rotation(t):
    e = cmath.exp(1j*t)
    return np.array([[(1+e)/2, (1-e)/2], [(1-e)/2, (1+e)/2]])


def ellipse_grading(a, b, theta, branch):
    """Grading of the major- or minor-axis rectangle generator of the ellipse at angle theta."""
    base = 0.0 if branch == "major" else math.pi/2
    polar_z, polar_w = base-theta/2, base-theta/2+math.pi

    def param(psi):
        ph = math.atan2(a*math.sin(psi), b*math.cos(psi))
        while ph < psi-math.pi/2:
            ph += 2*math.pi
        while ph > psi+math.pi/2:
            ph -= 2*math.pi
        return ph

    def line(ph):
        v = complex(-a*math.sin(ph), b*math.cos(ph))
        return v/abs(v)

    def turning(ph0, ph1, steps=2000):
        total, prev = 0.0, cmath.phase(line(ph0))
        for k in range(1, steps+1):
            cur = cmath.phase(line(ph0+(ph1-ph0)*k/steps))
            total += (cur-prev+math.pi/2) % math.pi-math.pi/2
            prev = cur
        return total
    pz, pw = param(polar_z), param(polar_w)
    xz, xw = param(polar_z+theta), param(polar_w+theta)
    phase = (2*theta+kahler_phase(np.diag([line(xz), line(xw)]), rotation(theta) @ np.diag([line(pz), line(pw)]))
             + 2*(turning(xz, pz)+turning(xw, pw)))
    mu = phase/(2*math.pi)
    if abs(mu-round(mu)) > 1e-6:
        raise ArithmeticError("non-integer loop phase", mu)
    return -round(mu)


class Grader:
    """Gradings of square generators (theta = pi/2) on a rational polygon."""
    def __init__(self, P):
        self.par = FastParam(P)
        self.n = n = len(P)
        d = [complex(float(P[(m+1) % n][0]-P[m][0]), float(P[(m+1) % n][1]-P[m][1])) for m in range(n)]
        self.dir = [z/abs(z) for z in d]
        self.turn = [cmath.phase(self.dir[v]/self.dir[v-1]) for v in range(n)]   # exterior angle at vertex v

    def lifted(self, p, x, theta=math.pi/2):
        """Endpoint parameters (a0, b0, a1, b1) of the preferred pair-path p -> x."""
        par, n = self.par, self.n
        a0, b0 = par.param_of(p[0]), par.param_of(p[1])
        if b0 <= a0:
            b0 += n
        a1, b1 = par.param_of(x[0]), par.param_of(x[1])
        while b1 <= a1:
            b1 += n
        while b1-a1 >= n:
            b1 -= n
        d0 = par.arg_change(a0, b0, a1, b1)
        if par.step is None:
            par.step = par.arg_change(a0, b0, a1+n, b1+n)-d0
        k = (theta-d0)/par.step
        if abs(k-round(k)) > 1e-6:
            raise ArithmeticError("no lift turns w - z by exactly theta")
        k = round(k)
        return a0, b0, a1+k*n, b1+k*n

    def _line(self, t, lam):
        v = math.floor(t)
        if t != v:
            return self.dir[v % self.n]
        return self.dir[(v-1) % self.n]*cmath.exp(1j*lam*self.turn[v % self.n])

    def _turning(self, s, t, lam_s, lam_t):
        lo, hi, sign = (s, t, 1) if t >= s else (t, s, -1)
        total = sign*sum(self.turn[v % self.n] for v in range(math.floor(lo)+1, math.ceil(hi)))
        if s == math.floor(s):
            d = self.turn[math.floor(s) % self.n]
            total += (1-lam_s)*d if t >= s else -lam_s*d
        if t == math.floor(t):
            d = self.turn[math.floor(t) % self.n]
            total += lam_t*d if t >= s else -(1-lam_t)*d
        return total

    def grading(self, p, x, lams=None, theta=math.pi/2):
        """Grading of x = R_theta(p); None at a polygon-vertex contact unless smoothing parameters
        lams (one per endpoint a1, b1, a0, b0) are given, or at a non-transverse smoothing."""
        a1_0 = self.lifted(p, x, theta)
        a0, b0, a1, b1 = (float(t) for t in a1_0)
        ts = (a1, b1, a0, b0)
        if lams is None:
            if any(t == math.floor(t) for t in ts):
                return None
            lams = (0.5, 0.5, 0.5, 0.5)
        lines = [self._line(t, lam) for t, lam in zip(ts, lams)]
        ret = 2*(self._turning(a1, a0, lams[0], lams[2])+self._turning(b1, b0, lams[1], lams[3]))
        try:
            conn = kahler_phase(np.diag(lines[:2]), rotation(theta) @ np.diag(lines[2:]))
        except ValueError:
            return None
        mu = (2*theta+conn+ret)/(2*math.pi)
        if abs(mu-round(mu)) > 1e-6:
            raise ArithmeticError("non-integer loop phase", mu)
        return -round(mu)

    def possible_gradings(self, p, x, samples=(0.02, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.98)):
        """Gradings over smoothings of vertex contacts (a single value for generic contacts)."""
        ts = [float(t) for t in self.lifted(p, x)]
        ts = (ts[2], ts[3], ts[0], ts[1])
        grids = [samples if t == math.floor(t) else (0.5,) for t in ts]
        out = set()
        for lams in itertools.product(*grids):
            g = self.grading(p, x, lams)
            if g is not None:
                out.add(g)
        return out

    def square_generators(self, square):
        c0, c1, c2, c3 = square
        return ((c0, c2), (c1, c3)), ((c1, c3), (c2, c0))


def grading_filter(N, workers=4):
    """Window generators of the v331 residue survivors and their gradings."""
    from residue_census_v331 import census, midside_fixture
    from angular_windows_v328_v330 import area, pi_enclosure, radius_squared
    row = census(N, workers=workers)
    P = midside_fixture(N)
    G = Grader(P)
    A, top = area(P), pi_enclosure(30)[1]*radius_squared(P)/2
    generic, vertex, sums, outside = {}, [], {}, 0
    for _, sq in row["survivors"]:
        gens = G.square_generators(sq)
        acts = G.par.square_actions(sq)
        window = [g for g, act in zip(gens, acts) if A/2 <= act <= top]
        if not window:
            outside += 1
            continue
        g = G.grading(*window[0])
        if g is None:
            vertex.append(sorted(G.possible_gradings(*window[0])))
            continue
        generic[g] = generic.get(g, 0)+1
        partner = G.grading(*[q for q in gens if q != window[0]][0])
        if partner is not None:
            sums[g+partner] = sums.get(g+partner, 0)+1
    return {"N": N, "survivors": len(row["survivors"]), "real_action_outside_window": outside,
            "generic_window_gradings": {str(k): v for k, v in sorted(generic.items())},
            "generic_grading_two": generic.get(2, 0), "vertex_contact_cases": len(vertex),
            "vertex_possible_gradings": sorted(set().union(*map(set, vertex))) if vertex else [],
            "vertex_cases_admitting_two": sum(2 in v for v in vertex),
            "duality_grading_sums": {str(k): v for k, v in sorted(sums.items())}}


def report(large=False):
    calibration = {f"{a}:{b}": [[ellipse_grading(a, b, t, "major"), ellipse_grading(a, b, t, "minor")]
                                for t in (0.3, 0.8, 1.2, 1.5)] for a, b in ((2.0, 1.0), (5.0, 1.0), (1.3, 1.0))}
    return {"ellipse_major_minor": calibration,
            "spiral": [grading_filter(N) for N in ((30,) if large else (24,))]}


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps(report(large=sys.argv[1:] == ["--large"]), indent=2))
