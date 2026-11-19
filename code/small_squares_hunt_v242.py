"""v242: numerical hunt for curves whose inscribed squares are all small.

Question (Matschke's sufficient condition): does every smooth Jordan curve
containing a disk of radius r inscribe a square of side at least sqrt(2)*r?
If yes (or even with any constant c>0 in place of sqrt 2), the standard
approximation argument proves Toeplitz.  We study

    J(curve) = (side of the largest inscribed square) / inradius,

which is sqrt(2) for a circle, and try hard to make it smaller.

NUMERICAL EXPERIMENT, NOT A PROOF.  Requires numpy and scipy; the optional
search driver also needs `cma`.

    python code/small_squares_hunt_v242.py            # replay checks (fast)
    python code/small_squares_hunt_v242.py --reverify # re-verify every recorded run (slow)
    python code/small_squares_hunt_v242.py --search K SEED GENS  # new star-shaped CMA run
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

TWO_PI = 2 * np.pi
SQRT2 = np.sqrt(2.0)
DATA = Path(__file__).resolve().parents[1] / 'research' / 'v242-small-squares-hunt.json'


# ------------------------------------------------------------------ curves --
class Radial:
    """Star-shaped smooth curve r(θ) = exp(f(θ)), f = Σ_k a_k cos kθ + b_k sin kθ. Always Jordan."""

    def __init__(self, a, b):
        self.a, self.b = np.asarray(a, float), np.asarray(b, float)
        self.k = np.arange(1, len(self.a) + 1)

    @classmethod
    def from_vector(cls, x):
        x = np.asarray(x, float)
        return cls(x[:len(x) // 2], x[len(x) // 2:])

    def f(self, th):
        kt = np.multiply.outer(np.asarray(th, float), self.k)
        return np.cos(kt) @ self.a + np.sin(kt) @ self.b

    def gamma(self, th):
        th = np.asarray(th, float)
        return np.exp(self.f(th)) * np.exp(1j * th)

    def g(self, w):
        """Zero exactly on the curve; negative inside, positive outside."""
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.log(np.abs(w)) - self.f(np.angle(w))

    def log_radius_spread(self, n=20000):
        f = self.f(np.linspace(0, TWO_PI, n, endpoint=False))
        return float(f.max() - f.min())


class Polygon:
    """Any simple closed polygon (dense samples of a general curve). g = signed distance."""

    def __init__(self, pts):
        P = np.asarray(pts, complex)
        self.A, self.B = P, np.roll(P, -1)
        self.cum = np.r_[0, np.cumsum(np.abs(self.B - self.A))]
        self.L = self.cum[-1]

    def gamma(self, th):
        s = (np.asarray(th, float) % TWO_PI) / TWO_PI * self.L
        i = np.clip(np.searchsorted(self.cum, s, side='right') - 1, 0, len(self.A) - 1)
        t = (s - self.cum[i]) / np.maximum(self.cum[i + 1] - self.cum[i], 1e-300)
        return self.A[i] + t * (self.B[i] - self.A[i])

    def g(self, w):
        w = np.asarray(w, complex)
        flat = w.ravel()
        out = np.empty(flat.shape)
        AB = self.B - self.A
        nAB = np.maximum(np.abs(AB) ** 2, 1e-300)
        for k in range(0, flat.size, 2048):
            x = flat[k:k + 2048, None]
            t = np.clip(((x - self.A) * np.conj(AB)).real / nAB, 0, 1)
            d = np.min(np.abs(self.A + t * AB - x), axis=1)
            with np.errstate(divide='ignore', invalid='ignore'):
                wind = np.sum(np.angle((self.B - x) / (self.A - x)), axis=1) / TWO_PI
            out[k:k + 2048] = np.where(np.abs(wind) > 0.5, -d, d)
        return out.reshape(w.shape)


# ----------------------------------------------------------------- squares --
def _diag(curve, s, u):
    p, q = curve.gamma(s), curve.gamma(u)
    m, h = (p + q) / 2, (q - p) / 2
    return m + 1j * h, m - 1j * h


def _adj(sigma):
    def build(curve, s, u):
        p, q = curve.gamma(s), curve.gamma(u)
        d = q - p
        return q + sigma * 1j * d, p + sigma * 1j * d
    return build


def _newton(curve, s, u, build, iters=40, h=1e-7):
    for _ in range(iters):
        r1, r2 = build(curve, s, u)
        R = np.stack([curve.g(r1), curve.g(r2)], -1)
        J = np.empty(R.shape + (2,))
        for col, (ds, du) in enumerate(((h, 0), (0, h))):
            a1, a2 = build(curve, s + ds, u + du)
            J[..., col] = (np.stack([curve.g(a1), curve.g(a2)], -1) - R) / h
        det = J[..., 0, 0] * J[..., 1, 1] - J[..., 0, 1] * J[..., 1, 0]
        det = np.where(np.abs(det) > 1e-300, det, np.inf)
        ds_ = np.clip((J[..., 1, 1] * R[..., 0] - J[..., 0, 1] * R[..., 1]) / det, -.05, .05)
        du_ = np.clip((-J[..., 1, 0] * R[..., 0] + J[..., 0, 0] * R[..., 1]) / det, -.05, .05)
        s, u = s - ds_, u - du_
    r1, r2 = build(curve, s, u)
    return s, u, np.maximum(np.abs(curve.g(r1)), np.abs(curve.g(r2)))


def find_squares(curve, M=256, method='diag', min_side=1e-3):
    """Squares found from sign changes on an M×M grid of vertex pairs, refined by Newton.

    method='diag': the two grid points are opposite vertices; 'adj': adjacent
    vertices, both orientations.  Missing a square can only make the largest
    side look SMALLER, so low J values must be re-checked at high M.
    """
    th = np.linspace(0, TWO_PI, M, endpoint=False)
    dth = TWO_PI / M
    S, U = np.meshgrid(th, th, indexing='ij')
    circ = np.minimum(np.abs(S - U), TWO_PI - np.abs(S - U))
    mask = circ > 2.5 * dth
    builds = [_diag] if method == 'diag' else [_adj(1), _adj(-1)]
    if method == 'diag':
        mask &= S < U
    found = []
    for build in builds:
        r1, r2 = build(curve, S, U)
        G = [curve.g(r1), curve.g(r2)]
        cand = mask.copy()
        for Gi in G:
            C = np.stack([Gi, np.roll(Gi, -1, 0), np.roll(Gi, -1, 1), np.roll(np.roll(Gi, -1, 0), -1, 1)])
            cand &= np.all(np.isfinite(C), 0) & (C.min(0) < 0) & (C.max(0) > 0)
        idx = np.argwhere(cand)
        if not len(idx):
            continue
        s, u, res = _newton(curve, th[idx[:, 0]] + dth / 2, th[idx[:, 1]] + dth / 2, build)
        for si, ui in zip(s[res < 1e-10], u[res < 1e-10]):
            p, q = curve.gamma(si), curve.gamma(ui)
            a, b = build(curve, np.array(si), np.array(ui))
            V = np.array([p, complex(a), q, complex(b)]) if method == 'diag' else np.array([p, q, complex(a), complex(b)])
            if abs(V[1] - V[0]) >= min_side:
                found.append(V)
    uniq, seen = [], set()
    for V in found:
        key = tuple(sorted((round(z.real, 7), round(z.imag, 7)) for z in V))
        if key not in seen:
            seen.add(key)
            uniq.append(V)
    return uniq, max((abs(V[1] - V[0]) for V in uniq), default=0.0)


# ---------------------------------------------------------------- inradius --
def _seg_dist(x, A, B):
    AB = B - A
    t = np.clip(((x - A) * np.conj(AB)).real / np.maximum(np.abs(AB) ** 2, 1e-300), 0, 1)
    return np.min(np.abs(A + t * AB - x))


def incircle(curve, N=6000, G=140):
    """(radius, centre) of the largest inscribed disk: grid search then Nelder–Mead."""
    th = np.linspace(0, TWO_PI, N, endpoint=False)
    P = curve.gamma(th)
    A, B = P, np.roll(P, -1)
    tree = cKDTree(np.c_[P.real, P.imag])
    lo, hi = P.real.min() + 1j * P.imag.min(), P.real.max() + 1j * P.imag.max()
    X, Y = np.meshgrid(np.linspace(lo.real, hi.real, G), np.linspace(lo.imag, hi.imag, G))
    W = (X + 1j * Y).ravel()
    W = W[curve.g(W) < 0]
    d, _ = tree.query(np.c_[W.real, W.imag])
    best = (0.0, 0j)
    for i in np.argsort(-d)[:6]:
        def neg(z):
            w = complex(z[0], z[1])
            return 1e3 if not curve.g(np.array(w)) < 0 else -_seg_dist(w, A, B)
        r = minimize(neg, [W[i].real, W[i].imag], method='Nelder-Mead',
                     options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 500})
        if -r.fun > best[0]:
            best = (-r.fun, complex(r.x[0], r.x[1]))
    return best


def inradius(curve, N=6000, G=140):
    return incircle(curve, N, G)[0]


def verified_ratio(curve, high=True):
    """J with the largest side taken over several resolutions and both parametrisations."""
    if high:
        side = max(find_squares(curve, 768)[1], find_squares(curve, 1536)[1],
                   find_squares(curve, 1024, 'adj')[1])
        return side / inradius(curve, N=40000, G=400)
    side = max(find_squares(curve, 384)[1], find_squares(curve, 512, 'adj')[1])
    return side / inradius(curve, N=12000, G=250)


# ----------------------------------------------------------- shape families --
def regular_polygon(n, per_edge=40, phase=0.37):
    V = np.exp(2j * np.pi * np.arange(n) / n + 1j * phase)
    return Polygon([V[k] + (V[(k + 1) % n] - V[k]) * t
                    for k in range(n) for t in np.linspace(0, 1, per_edge, endpoint=False)])


def ellipse(a=2.0, b=1.0, n=600):
    th = np.linspace(0, TWO_PI, n, endpoint=False)
    return Polygon(a * np.cos(th) + 1j * b * np.sin(th))


def general_fourier(x, K, npts=240):
    """General (non-star) curve exp(iθ) + Σ_{k≠0,1, |k|<=K} c_k exp(ikθ), as used in the search."""
    th = np.linspace(0, TWO_PI, npts, endpoint=False)
    ks = np.r_[np.arange(-K, 0), np.arange(2, K + 1)]
    x = np.asarray(x, float)
    c = x[:len(ks)] + 1j * x[len(ks):]
    return np.exp(1j * th) + np.exp(1j * np.outer(th, ks)) @ c


# -------------------------------------------------------------- search (opt) --
def cma_search(K, seed, gens=150, sigma0=0.3, M=192, cap=5.0):
    """Star-shaped CMA-ES run minimising J, with the two guards learned the hard way:
    a cap on the radius spread (e^cap) and a high-resolution re-check of low values."""
    import cma
    grid = np.linspace(0, TWO_PI, 1500, endpoint=False)

    def obj(x):
        c = Radial.from_vector(x)
        f = c.f(grid)
        if f.max() - f.min() > cap:
            return 10.0 + f.max() - f.min()
        _, side = find_squares(c, M)
        if side <= 0:
            return 10.0
        J = side / inradius(c, N=3000, G=90)
        if J < 1.45:
            J = max(side, find_squares(c, 768)[1]) / inradius(c, N=12000, G=250)
        return J

    rng = np.random.default_rng(seed)
    k = np.r_[np.arange(1, K + 1), np.arange(1, K + 1)]
    es = cma.CMAEvolutionStrategy(rng.normal(0, sigma0, 2 * K) / k, sigma0,
                                  {'seed': seed, 'popsize': 12, 'verbose': -9, 'maxiter': gens})
    best = (np.inf, None)
    while not es.stop():
        X = es.ask()
        F = [obj(x) for x in X]
        es.tell(X, F)
        i = int(np.argmin(F))
        if F[i] < best[0]:
            best = (F[i], list(map(float, X[i])))
    return {'K': K, 'seed': seed, 'search_J': best[0], 'x': best[1]}


# ------------------------------------------------------------------- replay --
def load_data():
    return json.loads(DATA.read_text())


def replay(reverify=False):
    data = load_data()
    out = {}
    e = ellipse()
    s = find_squares(e, 256)[1]
    out['ellipse_2x1'] = {'side': s, 'exact_side': 4 / np.sqrt(5), 'inradius': inradius(e)}
    out['triangle_J'] = verified_ratio(regular_polygon(3), high=False)
    out['hexagon_J'] = verified_ratio(regular_polygon(6), high=False)
    floor = next(r for r in data['runs'] if r['role'] == 'floor shape')
    out['floor_shape_J'] = verified_ratio(Radial.from_vector(floor['x']), high=False)
    trap = next(r for r in data['runs'] if r['role'] == 'resolution trap')
    c = Radial.from_vector(trap['x'])
    ir = inradius(c, N=12000, G=250)
    out['resolution_trap'] = {'J_at_M192': find_squares(c, 192)[1] / ir,
                              'J_at_M768': find_squares(c, 768)[1] / ir}
    if reverify:
        for r in data['runs']:
            if r['verified_J'] is None:
                continue
            curve = (Radial.from_vector(r['x']) if r['class'].startswith('star')
                     else Polygon(general_fourier(r['x'], r['K'])))
            r['reverified_J'] = verified_ratio(curve) if r['class'].startswith('star') else \
                max(find_squares(curve, 512)[1], find_squares(curve, 768)[1]) / inradius(curve, N=8000, G=200)
        out['reverified'] = [(r['class'], r['K'], r['seed'], r.get('reverified_J')) for r in data['runs']]
    vj = [r['verified_J'] for r in data['runs'] if r['verified_J'] is not None]
    out['min_verified_J_recorded'] = min(vj)
    out['sqrt2'] = SQRT2
    return out


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--search':
        print(json.dumps(cma_search(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))))
    else:
        res = replay(reverify='--reverify' in sys.argv)
        print(json.dumps(res, indent=2, default=float))
        assert abs(res['ellipse_2x1']['side'] - res['ellipse_2x1']['exact_side']) < 1e-4
        assert abs(res['triangle_J'] - 6 * (2 - np.sqrt(3))) < 1e-3
        assert abs(res['hexagon_J'] - (2 * np.sqrt(3) - 2)) < 1e-3
        assert SQRT2 - 1e-5 <= res['floor_shape_J'] < SQRT2 + 1e-4
        assert res['resolution_trap']['J_at_M192'] < 1.2 < 1.6 < res['resolution_trap']['J_at_M768']
        assert res['min_verified_J_recorded'] >= SQRT2 - 1e-5
